"""
scrna_qc_local.py - Local fallback scRNA-seq QC tool.

Used when single-cell-rna-qc@life-sciences is unavailable or incompatible.
Computes core QC metrics, plots distributions, recommends filters with biological justification.
Does NOT automatically filter cells unless --apply-filters is explicitly set.
"""

import argparse
import csv
import gzip
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

# Soft imports
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

try:
    import scanpy as sc
    import anndata as ad
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

VERSION = "1.0.0"
TOOL_NAME = "scrna_qc_local.py"


def setup_logging(output_dir: Path, verbose: bool) -> logging.Logger:
    log_path = output_dir / "scrna_qc.log"
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stderr)],
    )
    return logging.getLogger(TOOL_NAME)


def mad_threshold(values, n_mads: float = 5.0):
    """MAD-based outlier thresholds (lower, upper)."""
    if not HAS_NUMPY:
        return None, None
    arr = np.array(values)
    median = np.median(arr)
    mad = np.median(np.abs(arr - median))
    lower = median - n_mads * mad
    upper = median + n_mads * mad
    return lower, upper


def get_mito_prefix(species: str, user_prefix: str) -> list:
    if user_prefix:
        return [user_prefix]
    species_lower = species.lower()
    if "human" in species_lower or "homo" in species_lower:
        return ["MT-"]
    elif "mouse" in species_lower or "mus" in species_lower:
        return ["mt-"]
    else:
        return ["MT-", "mt-", "Mt-"]


def load_data(input_path: Path, log: logging.Logger):
    """Load count matrix. Returns AnnData or None."""
    if not HAS_SCANPY:
        log.error("scanpy/anndata not available. Install: pip install scanpy anndata")
        return None

    path_str = str(input_path)
    log.info(f"Loading: {path_str}")

    if input_path.is_dir():
        log.info("Detected 10x directory format.")
        try:
            adata = sc.read_10x_mtx(path_str, var_names="gene_symbols", cache=False)
            log.info(f"Loaded 10x: {adata.shape}")
            return adata
        except Exception as e:
            log.error(f"Failed to read 10x directory: {e}")
            return None

    ext = input_path.suffix.lower()
    if ext == ".h5ad":
        try:
            adata = sc.read_h5ad(path_str)
            log.info(f"Loaded h5ad: {adata.shape}")
            return adata
        except Exception as e:
            log.error(f"Failed to read h5ad: {e}")
            return None

    if ext in (".csv", ".tsv", ".txt"):
        sep = "," if ext == ".csv" else "\t"
        try:
            df = pd.read_csv(path_str, index_col=0, sep=sep)
            log.info(f"Loaded CSV/TSV: {df.shape}")
            adata = ad.AnnData(X=df.values.T.astype(float))
            adata.obs_names = [str(i) for i in df.columns]
            adata.var_names = list(df.index)
            return adata
        except Exception as e:
            log.error(f"Failed to read CSV/TSV: {e}")
            return None

    log.error(f"Unsupported format: {ext}")
    return None


def compute_qc(adata, mito_prefixes: list, ribo_prefixes: list, log: logging.Logger) -> dict:
    """Compute per-cell QC metrics. Returns dict of metric arrays."""
    import numpy as np

    var_names = list(adata.var_names)
    var_names_lower = [v.upper() for v in var_names]

    # Mitochondrial genes
    mito_mask = [
        any(v.startswith(p.upper()) for p in [mp.upper() for mp in mito_prefixes])
        for v in var_names_lower
    ]
    adata.var["is_mito"] = mito_mask
    n_mito = sum(mito_mask)
    log.info(f"Mitochondrial genes found: {n_mito} (prefixes: {mito_prefixes})")

    # Ribosomal genes
    ribo_mask = [
        any(v.startswith(p.upper()) for p in [rp.upper() for rp in ribo_prefixes])
        for v in var_names_lower
    ]
    adata.var["is_ribo"] = ribo_mask
    n_ribo = sum(ribo_mask)
    log.info(f"Ribosomal genes found: {n_ribo} (prefixes: {ribo_prefixes})")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["is_mito", "is_ribo"] if any(ribo_mask) else ["is_mito"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    metrics = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_mito_genes": n_mito,
        "n_ribo_genes": n_ribo,
        "total_counts": list(adata.obs["total_counts"].astype(float)),
        "n_genes_by_counts": list(adata.obs["n_genes_by_counts"].astype(float)),
    }

    if "pct_counts_is_mito" in adata.obs.columns:
        metrics["pct_mito"] = list(adata.obs["pct_counts_is_mito"].astype(float))
    elif n_mito == 0:
        metrics["pct_mito"] = [0.0] * adata.n_obs
        log.warning("No mitochondrial genes found with provided prefixes. pct_mito = 0.")

    if "pct_counts_is_ribo" in adata.obs.columns:
        metrics["pct_ribo"] = list(adata.obs["pct_counts_is_ribo"].astype(float))

    # Top expressed genes
    top_genes = (
        adata.var["total_counts"].sort_values(ascending=False).head(20).index.tolist()
        if "total_counts" in adata.var.columns
        else []
    )
    metrics["top_expressed_genes"] = top_genes

    return metrics


def suggest_filters(metrics: dict, mad_n: float, log: logging.Logger) -> list:
    """Suggest filters using MAD-based outlier detection with filter decision framework."""
    suggestions = []

    filter_specs = [
        (
            "total_counts", "total_counts", "both",
            "Low counts may indicate empty droplets or dead cells. High counts may indicate doublets OR large/activated cells.",
            {
                "low": {
                    "technical_artifact": "Empty droplets, barcodes capturing ambient RNA, lysed cells",
                    "biological_signal": "Rare mature cell types with very low transcriptional activity (mature RBCs, platelets, quiescent stem cells)",
                    "evidence_for_filtering": "Counts near or below ambient RNA level; no detectable marker genes; concordant low n_genes; barcode not in cell-calling whitelist",
                    "evidence_against_filtering": "Detectable lineage markers (HBB, PPBP); cell count in tissue matches expected rare population size; consistent with known biology of sample",
                    "recommended_action_low": "Flag - inspect for known markers before filtering. If marker-negative and counts near 0, consider Filter.",
                },
                "high": {
                    "technical_artifact": "Multiplet / doublet barcodes capturing two or more cells",
                    "biological_signal": "Large cells (megakaryocytes, hepatocytes), highly activated immune cells (plasma cells, effector T cells), cycling tumor cells",
                    "evidence_for_filtering": "High doublet score from scrublet/DoubletFinder; co-expression of mutually exclusive lineage markers (e.g., CD3D + CD19); counts 2x the expected single-cell median",
                    "evidence_against_filtering": "No doublet score available; tissue known to contain large cells; high counts concordant with expected activated/cycling population",
                    "recommended_action_high": "Flag - run doublet detection (scrublet/DoubletFinder) before filtering. Do not remove without doublet score.",
                },
            }
        ),
        (
            "n_genes_by_counts", "n_genes_by_counts", "both",
            "Low genes may indicate poor-quality cells OR mature cell types (RBCs, platelets). High genes may indicate doublets OR actively transcribing cells.",
            {
                "low": {
                    "technical_artifact": "Poor lysis, low capture efficiency, barcodes in low-RNA regions",
                    "biological_signal": "Mature erythrocytes, platelets, or highly specialized cells with restricted transcriptomes",
                    "evidence_for_filtering": "Low counts and low genes together; no marker-gene expression; outside expected distribution for protocol/chemistry",
                    "evidence_against_filtering": "HBB/HBA1 or PPBP expression present; known erythroid/platelet biology in tissue; consistent with matched bulk RNA profile",
                    "recommended_action_low": "Flag - check for HBB/HBA1 expression (RBC markers) before filtering low-gene cells.",
                },
                "high": {
                    "technical_artifact": "Doublets (two cells merged in one droplet)",
                    "biological_signal": "Large transcriptionally active cells, tumor cells, activated immune cells in S/G2M phase",
                    "evidence_for_filtering": "High doublet score; co-expression of incompatible lineage markers; gene count 2x the modal distribution",
                    "evidence_against_filtering": "No doublet detection available; G2M marker signature present; consistent with tumor or activated immune context",
                    "recommended_action_high": "Flag - run doublet detection. Stratify by cell cycle if G2M markers present.",
                },
            }
        ),
        (
            "pct_mito", "pct_mito", "upper",
            "High mitochondrial % may indicate dying cells OR metabolic tissue, hypoxic/stressed biology, cardiac/muscle cells.",
            {
                "high": {
                    "technical_artifact": "Cell lysis (cytoplasmic RNA lost, mt-RNA retained), membrane damage during dissociation",
                    "biological_signal": "Highly metabolic tissue (heart, brown fat, liver), hypoxic tumor cells, activated immune cells under oxidative stress, ferroptotic cells, cardiac/muscle cells",
                    "evidence_for_filtering": "Concordant low counts and low genes alongside high mito%; no tissue-specific metabolic marker expression; pattern consistent across all batches/samples",
                    "evidence_against_filtering": "Tissue is cardiac, hepatic, brown fat, or skeletal muscle; high mito% cells express tissue-specific markers; high mito% cluster is unique to specific sample or condition",
                    "recommended_action_high": "Stratify - separate high-mito cells and check for tissue-specific markers. Filter only if marker-negative and count/gene metrics are also low.",
                },
            }
        ),
    ]

    for metric_name, values_key, direction, bio_note, actions in filter_specs:
        if values_key not in metrics:
            continue
        values = metrics[values_key]
        lower, upper = mad_threshold(values, mad_n)

        suggestion = {
            "metric": metric_name,
            "direction": direction,
            "mad_threshold": mad_n,
            "computed_lower": round(lower, 3) if lower is not None else None,
            "computed_upper": round(upper, 3) if upper is not None else None,
            "biological_note": bio_note,
            "filter_decision_framework": actions,
        }

        if direction in ("lower", "both") and lower is not None:
            n_flagged = sum(1 for v in values if v < lower)
            suggestion["n_flagged_low"] = n_flagged
            suggestion["pct_flagged_low"] = round(100 * n_flagged / len(values), 1)

        if direction in ("upper", "both") and upper is not None:
            n_flagged = sum(1 for v in values if v > upper)
            suggestion["n_flagged_high"] = n_flagged
            suggestion["pct_flagged_high"] = round(100 * n_flagged / len(values), 1)

        suggestion["validation_note"] = (
            "Plot distribution before applying. Consider tissue context. "
            "Run with different MAD thresholds. Inspect flagged cells for known marker expression. "
            "Perform sensitivity analysis: compare downstream results with and without filter."
        )
        suggestions.append(suggestion)

    return suggestions


def make_plots(adata, metrics: dict, output_dir: Path, log: logging.Logger):
    """Generate QC distribution plots."""
    if not HAS_MATPLOTLIB:
        log.warning("matplotlib not available. Skipping plots.")
        return []

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    generated = []

    def save_hist(values, title, xlabel, fname):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(values, bins=50, color="#2196F3", edgecolor="white", linewidth=0.3)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Number of cells")
        path = plots_dir / fname
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        generated.append(str(path))
        log.info(f"Plot saved: {path}")

    if "total_counts" in metrics:
        save_hist(metrics["total_counts"], "Total UMI Counts per Cell", "Total counts", "total_counts.png")

    if "n_genes_by_counts" in metrics:
        save_hist(metrics["n_genes_by_counts"], "Genes Detected per Cell", "Number of genes", "n_genes.png")

    if "pct_mito" in metrics:
        save_hist(metrics["pct_mito"], "Mitochondrial % per Cell", "% mitochondrial", "pct_mito.png")

    # Scatter: genes vs counts
    if "total_counts" in metrics and "n_genes_by_counts" in metrics:
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(metrics["total_counts"], metrics["n_genes_by_counts"],
                   s=1, alpha=0.3, c="#2196F3")
        ax.set_xlabel("Total counts")
        ax.set_ylabel("Genes detected")
        ax.set_title("Genes vs Counts")
        path = plots_dir / "genes_vs_counts.png"
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        generated.append(str(path))

    # Scatter: mito% vs counts
    if "pct_mito" in metrics and "total_counts" in metrics:
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(metrics["total_counts"], metrics["pct_mito"],
                   s=1, alpha=0.3, c="#E53935")
        ax.set_xlabel("Total counts")
        ax.set_ylabel("% mitochondrial")
        ax.set_title("Mitochondrial % vs Counts")
        path = plots_dir / "mito_vs_counts.png"
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        generated.append(str(path))

    return generated


def write_metrics_csv(metrics: dict, output_dir: Path) -> Path:
    path = output_dir / "qc_metrics.csv"
    n = metrics.get("n_cells", 0)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_index", "total_counts", "n_genes", "pct_mito", "pct_ribo"])
        for i in range(n):
            row = [
                i,
                round(metrics["total_counts"][i], 2) if "total_counts" in metrics else "",
                int(metrics["n_genes_by_counts"][i]) if "n_genes_by_counts" in metrics else "",
                round(metrics["pct_mito"][i], 3) if "pct_mito" in metrics else "",
                round(metrics["pct_ribo"][i], 3) if "pct_ribo" in metrics else "",
            ]
            writer.writerow(row)
    return path


def build_markdown_report(args, metrics: dict, suggestions: list, plots: list,
                           output_dir: Path) -> str:
    lines = [
        "# scRNA-seq QC Report (Local Fallback)",
        f"\nGenerated: {datetime.now().isoformat()}",
        f"\n**Tool**: `{TOOL_NAME}` v{VERSION}",
        f"\n**Note**: This report was generated by the local fallback tool. If `single-cell-rna-qc@life-sciences` is available and compatible, prefer that official skill.",
        "\n## Dataset Summary\n",
        f"- **Input**: `{args.input}`",
        f"- **Species**: {args.species}",
        f"- **Tissue**: {args.tissue}",
        f"- **Protocol**: {args.protocol}",
        f"- **Chemistry**: {args.chemistry}",
        f"- **Genome build**: {args.genome_build}",
        f"- **Cells**: {metrics.get('n_cells', '?')}",
        f"- **Genes**: {metrics.get('n_genes', '?')}",
        f"- **Mitochondrial genes detected**: {metrics.get('n_mito_genes', '?')}",
        f"- **Ribosomal genes detected**: {metrics.get('n_ribo_genes', '?')}",
        "\n## QC Metrics Summary\n",
    ]

    for key, label in [("total_counts", "Total counts"), ("n_genes_by_counts", "Genes per cell"),
                        ("pct_mito", "Mitochondrial %")]:
        if key in metrics and metrics[key]:
            vals = metrics[key]
            import statistics
            lines.append(f"**{label}**: median={round(statistics.median(vals), 1)}, "
                         f"mean={round(sum(vals)/len(vals), 1)}, "
                         f"min={round(min(vals), 1)}, max={round(max(vals), 1)}")

    if metrics.get("top_expressed_genes"):
        lines.append(f"\n**Top expressed genes**: {', '.join(metrics['top_expressed_genes'][:10])}")

    lines.append("\n## Recommended Filters (NOT applied unless --apply-filters is set)\n")
    for s in suggestions:
        lines.append(f"\n### {s['metric']}")
        lines.append(f"- MAD threshold: ±{s['mad_threshold']} MADs")
        if s.get("computed_lower") is not None:
            lines.append(f"- Suggested lower bound: {s['computed_lower']}")
        if s.get("computed_upper") is not None:
            lines.append(f"- Suggested upper bound: {s['computed_upper']}")
        if "n_flagged_low" in s:
            lines.append(f"- Cells flagged below lower bound: {s['n_flagged_low']} ({s['pct_flagged_low']}%)")
        if "n_flagged_high" in s:
            lines.append(f"- Cells flagged above upper bound: {s['n_flagged_high']} ({s['pct_flagged_high']}%)")
        lines.append(f"\n**Biological note**: {s['biological_note']}")

        fdf = s.get("filter_decision_framework", {})
        for bound, details in fdf.items():
            lines.append(f"\n**Filter decision ({bound} outliers):**")
            lines.append(f"- Observed metric: {s['metric']} ({bound} tail)")
            lines.append(f"- Possible technical artifact: {details.get('technical_artifact', 'unknown')}")
            lines.append(f"- Possible biological signal: {details.get('biological_signal', 'unknown')}")
            if "evidence_for_filtering" in details:
                lines.append(f"- Evidence FOR filtering: {details['evidence_for_filtering']}")
            if "evidence_against_filtering" in details:
                lines.append(f"- Evidence AGAINST filtering: {details['evidence_against_filtering']}")
            action_key = f"recommended_action_{bound}"
            if action_key in details:
                lines.append(f"- **Decision (Filter/Flag/Stratify/Preserve)**: {details[action_key]}")

        lines.append(f"\n**Validation**: {s['validation_note']}")

    lines.append("\n## Artifact vs. Plausible Biology\n")
    lines.append("| Observation | Possible Technical Explanation | Possible Biological Explanation | Confidence |")
    lines.append("|-------------|-------------------------------|--------------------------------|------------|")
    lines.append("| High mito% cells | Dying cells, membrane rupture | Metabolic tissue, hypoxia, tumor stress, cardiac/muscle | low - requires tissue context |")
    lines.append("| High UMI/gene cells | Doublets | Large cells, activated immune/tumor cells, plasma cells | low - requires doublet scoring |")
    lines.append("| Low gene count cells | Empty droplets, dead cells | Mature RBCs, platelets, sparse cell types | low - requires marker check |")

    lines.append("\n## Plots Generated\n")
    for p in plots:
        lines.append(f"- `{p}`")
    if not plots:
        lines.append("- No plots generated (matplotlib not available or --no-plots set).")

    lines.append("\n## Assumptions\n")
    lines.append(f"- Species: {args.species} (user-provided)")
    lines.append(f"- Tissue: {args.tissue} (user-provided)")
    lines.append(f"- Mitochondrial prefix(es): {get_mito_prefix(args.species, args.mito_prefix)}")
    lines.append(f"- MAD threshold: {args.mad_threshold}")
    lines.append("\n## Limitations\n")
    lines.append("- Doublet detection not performed (use scrublet or DoubletFinder).")
    lines.append("- Ambient RNA estimation not performed (use SoupX or DecontX).")
    lines.append("- Cell type annotation not performed.")
    lines.append("- Filters are recommendations only unless --apply-filters is used.")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Local scRNA-seq QC fallback tool.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Input: .h5ad, 10x dir, .csv, or .tsv")
    parser.add_argument("--output-dir", default="reports/scrna_qc", help="Output directory")
    parser.add_argument("--species", default="unknown", help="Species (human/mouse/other)")
    parser.add_argument("--tissue", default="unknown", help="Tissue type")
    parser.add_argument("--protocol", default="unknown", help="Protocol (10x_v3/dropseq/etc)")
    parser.add_argument("--chemistry", default="unknown", help="Chemistry version")
    parser.add_argument("--genome-build", default="unknown", help="Genome build (GRCh38/etc)")
    parser.add_argument("--sample-metadata", default=None, help="CSV with sample-level metadata")
    parser.add_argument("--mito-prefix", default="", help="Mitochondrial gene prefix (e.g., MT-)")
    parser.add_argument("--ribosomal-prefixes", default="RPS,RPL", help="Ribosomal gene prefixes (comma-separated)")
    parser.add_argument("--min-genes", type=int, default=None, help="Hard minimum genes (overrides MAD)")
    parser.add_argument("--max-mito-pct", type=float, default=None, help="Hard max mito pct (overrides MAD)")
    parser.add_argument("--mad-threshold", type=float, default=5.0, help="MAD multiplier for outlier detection")
    parser.add_argument("--recommend-only", action="store_true", help="Only recommend filters, do not apply")
    parser.add_argument("--apply-filters", action="store_true", help="Apply suggested filters (not recommended without review)")
    parser.add_argument("--json", action="store_true", help="Write JSON summary (default: always written)")
    parser.add_argument("--markdown", action="store_true", help="Write Markdown report (default: always written)")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logging(out_dir, args.verbose)

    log.info(f"{TOOL_NAME} v{VERSION}")
    log.info(f"Python: {sys.version}")
    log.info(f"scanpy available: {HAS_SCANPY}")
    log.info(f"Args: {vars(args)}")

    if not HAS_SCANPY:
        log.error("scanpy and anndata are required. Install: pip install scanpy anndata")
        summary = {
            "error": "scanpy not available",
            "tool": TOOL_NAME,
            "version": VERSION,
            "install": "pip install scanpy anndata",
        }
        with open(out_dir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2)
        return 1

    adata = load_data(Path(args.input), log)
    if adata is None:
        return 1

    log.info(f"Dataset shape: {adata.shape} (cells x genes)")

    mito_prefixes = get_mito_prefix(args.species, args.mito_prefix)
    ribo_prefixes = [p.strip() for p in args.ribosomal_prefixes.split(",") if p.strip()]

    metrics = compute_qc(adata, mito_prefixes, ribo_prefixes, log)
    suggestions = suggest_filters(metrics, args.mad_threshold, log)
    plots = make_plots(adata, metrics, out_dir, log)

    csv_path = write_metrics_csv(metrics, out_dir)
    log.info(f"Metrics CSV: {csv_path}")

    report_md = build_markdown_report(args, metrics, suggestions, plots, out_dir)
    report_path = out_dir / "qc_report.md"
    with open(report_path, "w") as f:
        f.write(report_md)
    log.info(f"Markdown report: {report_path}")

    summary = {
        "tool": TOOL_NAME,
        "version": VERSION,
        "python_version": sys.version,
        "scanpy_version": sc.__version__ if HAS_SCANPY else None,
        "generated": datetime.now().isoformat(),
        "primary_path": "single-cell-rna-qc@life-sciences",
        "fallback_note": (
            "This output was produced by the local fallback tool (tools/scrna_qc_local.py). "
            "The official Anthropic Life Sciences skill single-cell-rna-qc@life-sciences is the preferred "
            "and primary scRNA QC path. This fallback was used because the official skill was unavailable "
            "or the input format was incompatible. Re-run with the official skill when available."
        ),
        "input": args.input,
        "species": args.species,
        "tissue": args.tissue,
        "protocol": args.protocol,
        "chemistry": args.chemistry,
        "genome_build": args.genome_build,
        "n_cells": metrics.get("n_cells"),
        "n_genes": metrics.get("n_genes"),
        "n_mito_genes": metrics.get("n_mito_genes"),
        "mito_prefixes_used": mito_prefixes,
        "mad_threshold": args.mad_threshold,
        "filter_suggestions": suggestions,
        "plots": plots,
        "outputs": {
            "metrics_csv": str(csv_path),
            "report_md": str(report_path),
        },
        "notes": [
            "Filters are recommendations only. Review biologically before applying.",
            "Official single-cell-rna-qc@life-sciences should be preferred when available.",
            "Doublet detection and ambient RNA estimation not performed.",
        ],
    }
    json_path = out_dir / "summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log.info(f"JSON summary: {json_path}")

    print(f"\nQC complete. Outputs in: {out_dir}")
    print(f"  Cells: {metrics.get('n_cells')}, Genes: {metrics.get('n_genes')}")
    print(f"  Suggestions: {len(suggestions)} filter criteria recommended (not applied)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
