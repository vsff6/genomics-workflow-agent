"""
atac_qc_local.py - Local ATAC-seq / scATAC-seq QC summarization.

Computes metrics from fragments files, peaks BED, and optional BAM.
Uses bedtools when available for blacklist overlap and FRiP validation.
Detects chromosome naming mismatches between input files.
Never downloads references. Clearly reports skipped metrics and why.
"""

import argparse
import csv
import gzip
import json
import logging
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

VERSION = "2.0.0"
TOOL_NAME = "atac_qc_local.py"
MAX_FRAGMENTS_SCAN = 5_000_000


def setup_logging(output_dir: Path, verbose: bool) -> logging.Logger:
    log_path = output_dir / "atac_qc.log"
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stderr)],
    )
    return logging.getLogger(TOOL_NAME)


def tool_available(cmd: str) -> bool:
    try:
        result = subprocess.run([cmd, "--version"], capture_output=True, timeout=5)
        return result.returncode == 0
    except Exception:
        return False


def get_tool_version(cmd: str) -> str:
    try:
        result = subprocess.run([cmd, "--version"], capture_output=True, text=True, timeout=5)
        output = (result.stdout + result.stderr).strip()
        return output.split("\n")[0] if output else "unknown"
    except Exception:
        return "unknown"


# ──────────────────────────────────────────────────────────
# Chromosome naming helpers
# ──────────────────────────────────────────────────────────

def sample_chroms_from_bed(path: Path, n: int = 100) -> set:
    """Read up to n records from a BED/fragments/GTF file and return observed chromosome names."""
    chroms = set()
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if parts:
                    chroms.add(parts[0].strip())
                if len(chroms) >= n:
                    break
    except Exception:
        pass
    return chroms


def classify_chrom_style(chroms: set) -> str:
    """Return 'ucsc' (chr1), 'ensembl' (1), or 'unknown'."""
    for c in chroms:
        if c.startswith("chr"):
            return "ucsc"
        if re.match(r"^\d+$", c) or c in ("X", "Y", "MT", "M"):
            return "ensembl"
    return "unknown"


def check_chrom_mismatches(fragments_path, peaks_path, blacklist_path, gtf_path, log) -> list:
    """Detect chromosome naming mismatches between input files. Returns list of warning dicts."""
    warnings = []
    styles = {}

    sources = {}
    if fragments_path:
        sources["fragments"] = Path(fragments_path)
    if peaks_path:
        sources["peaks"] = Path(peaks_path)
    if blacklist_path:
        sources["blacklist"] = Path(blacklist_path)
    if gtf_path:
        sources["gtf"] = Path(gtf_path)

    for name, path in sources.items():
        if path.exists():
            chroms = sample_chroms_from_bed(path)
            style = classify_chrom_style(chroms)
            styles[name] = {"style": style, "example_chroms": list(chroms)[:5]}

    unique_styles = {v["style"] for v in styles.values() if v["style"] != "unknown"}
    if len(unique_styles) > 1:
        style_summary = {k: v["style"] for k, v in styles.items()}
        msg = (
            f"Chromosome naming mismatch detected across input files: {style_summary}. "
            "UCSC style (chr1) vs Ensembl style (1) will cause bedtools intersect and "
            "interval-based analyses to produce zero overlaps silently. "
            "Rename chromosomes to a consistent style before analysis."
        )
        log.warning(msg)
        warnings.append({
            "type": "chrom_naming_mismatch",
            "details": style_summary,
            "warning": msg,
        })
    return warnings


# ──────────────────────────────────────────────────────────
# bedtools integration
# ──────────────────────────────────────────────────────────

def run_bedtools_blacklist_overlap(fragments_path: Path, blacklist_path: Path, log) -> dict:
    """Compute fraction of fragments overlapping blacklist regions using bedtools intersect."""
    n_total = 0
    opener = gzip.open if str(fragments_path).endswith(".gz") else open
    try:
        with opener(fragments_path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                n_total += 1
    except Exception as e:
        return {"error": f"Could not count fragments: {e}"}

    if n_total == 0:
        return {"blacklist_fraction": None, "skipped_reason": "No fragments loaded"}

    cmd = ["bedtools", "intersect", "-u", "-a", str(fragments_path), "-b", str(blacklist_path)]
    log.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        out = {
            "command": " ".join(cmd),
            "return_code": result.returncode,
            "tool_version": get_tool_version("bedtools"),
            "stderr_snippet": result.stderr[:500] if result.stderr else "",
            "n_total_fragments": n_total,
        }
        if result.returncode == 0:
            n_blacklisted = sum(1 for line in result.stdout.splitlines() if line.strip())
            blacklist_frac = n_blacklisted / n_total if n_total > 0 else 0.0
            out["n_blacklisted_fragments"] = n_blacklisted
            out["blacklist_fraction"] = round(blacklist_frac, 4)
            out["note"] = (
                "Fragments overlapping blacklist regions are likely artifacts "
                "(satellite repeats, rDNA, centromeres, high-signal regions). "
                "These should be removed before accessibility analysis."
            )
        else:
            out["error"] = result.stderr[:500]
        return out
    except subprocess.TimeoutExpired:
        return {"command": " ".join(cmd), "error": "bedtools intersect timed out (300s)"}
    except Exception as e:
        return {"command": " ".join(cmd), "error": str(e)}


def run_bedtools_frip_validation(fragments_path: Path, peaks_path: Path, log) -> dict:
    """Validate FRiP using bedtools intersect (authoritative, handles large files)."""
    n_total = 0
    opener = gzip.open if str(fragments_path).endswith(".gz") else open
    try:
        with opener(fragments_path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                n_total += 1
    except Exception as e:
        return {"error": f"Could not count fragments: {e}"}

    cmd = ["bedtools", "intersect", "-u", "-a", str(fragments_path), "-b", str(peaks_path)]
    log.info(f"Running bedtools FRiP validation: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        out = {
            "command": " ".join(cmd),
            "return_code": result.returncode,
            "tool_version": get_tool_version("bedtools"),
            "stderr_snippet": result.stderr[:500] if result.stderr else "",
            "n_total_fragments": n_total,
        }
        if result.returncode == 0:
            n_in_peaks = sum(1 for line in result.stdout.splitlines() if line.strip())
            frip = n_in_peaks / n_total if n_total > 0 else 0.0
            out["n_fragments_in_peaks"] = n_in_peaks
            out["frip_bedtools"] = round(frip, 4)
            out["note"] = (
                "FRiP computed by bedtools intersect (authoritative). "
                "FRiP varies by tissue, protocol, and cell type - do not apply universal thresholds."
            )
        else:
            out["error"] = result.stderr[:500]
        return out
    except subprocess.TimeoutExpired:
        return {"command": " ".join(cmd), "error": "bedtools intersect timed out (300s)"}
    except Exception as e:
        return {"command": " ".join(cmd), "error": str(e)}


# ──────────────────────────────────────────────────────────
# Streaming fragment metrics (local fallback - no bedtools)
# ──────────────────────────────────────────────────────────

def read_fragments(fragments_path: Path, max_lines: int = MAX_FRAGMENTS_SCAN, log=None):
    """Read fragments file. Returns list of (chrom, start, end, barcode) tuples."""
    records = []
    opener = gzip.open if str(fragments_path).endswith(".gz") else open
    try:
        with opener(fragments_path, "rt", errors="replace") as f:
            for i, line in enumerate(f):
                if i >= max_lines:
                    if log:
                        log.warning(f"Fragment scan capped at {max_lines} lines.")
                    break
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 4:
                    continue
                try:
                    records.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
                except ValueError:
                    continue
    except Exception as e:
        if log:
            log.error(f"Error reading fragments: {e}")
    return records


def compute_insert_sizes(records: list) -> list:
    return [end - start for (_, start, end, _) in records if end > start]


def compute_frip(records: list, peaks_path: Path, log=None) -> dict:
    """Compute FRiP using simple Python interval overlap (local fallback)."""
    if not peaks_path or not peaks_path.exists():
        return {"frip": None, "skipped_reason": "peaks BED not provided"}

    peaks_by_chrom = defaultdict(list)
    try:
        opener = gzip.open if str(peaks_path).endswith(".gz") else open
        with opener(peaks_path, "rt") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                peaks_by_chrom[parts[0]].append((int(parts[1]), int(parts[2])))
    except Exception as e:
        return {"frip": None, "skipped_reason": f"Error reading peaks: {e}"}

    n_total = len(records)
    n_in_peaks = 0

    for (chrom, start, end, _) in records:
        for (ps, pe) in peaks_by_chrom.get(chrom, []):
            if start < pe and end > ps:
                n_in_peaks += 1
                break

    if n_total == 0:
        return {"frip": None, "skipped_reason": "No fragments loaded"}

    frip = n_in_peaks / n_total
    return {
        "frip": round(frip, 4),
        "n_total_fragments": n_total,
        "n_fragments_in_peaks": n_in_peaks,
        "note": (
            "FRiP computed by local Python interval overlap (streaming fallback). "
            "For large datasets, use bedtools intersect for authoritative FRiP. "
            "FRiP threshold varies by tissue and protocol. Do not apply universal cutoffs."
        ),
    }


def count_peaks(peaks_path: Path, log=None) -> dict:
    if not peaks_path or not peaks_path.exists():
        return {"n_peaks": None, "skipped_reason": "peaks BED not provided"}

    widths = []
    try:
        opener = gzip.open if str(peaks_path).endswith(".gz") else open
        with opener(peaks_path, "rt") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        widths.append(int(parts[2]) - int(parts[1]))
                    except ValueError:
                        pass
    except Exception as e:
        return {"n_peaks": None, "skipped_reason": str(e)}

    if not widths:
        return {"n_peaks": 0}

    return {
        "n_peaks": len(widths),
        "median_peak_width_bp": round(sorted(widths)[len(widths) // 2], 1),
        "mean_peak_width_bp": round(sum(widths) / len(widths), 1),
        "min_peak_width_bp": min(widths),
        "max_peak_width_bp": max(widths),
    }


def compute_barcode_stats(records: list) -> dict:
    bc_counts = defaultdict(int)
    for (_, _, _, bc) in records:
        bc_counts[bc] += 1

    if not bc_counts:
        return {"n_barcodes": 0}

    counts = list(bc_counts.values())
    import statistics
    return {
        "n_barcodes": len(counts),
        "median_frags_per_barcode": round(statistics.median(counts), 1),
        "mean_frags_per_barcode": round(sum(counts) / len(counts), 1),
        "min_frags_per_barcode": min(counts),
        "max_frags_per_barcode": max(counts),
    }


def make_insert_size_plot(insert_sizes: list, output_dir: Path, log=None) -> str:
    if not HAS_MATPLOTLIB or not insert_sizes:
        return None

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    sizes = [s for s in insert_sizes if 0 < s < 2000]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(sizes, bins=200, color="#4CAF50", edgecolor="none")
    for x, label in [(147, "mono"), (294, "di"), (441, "tri")]:
        ax.axvline(x, color="red", alpha=0.5, linestyle="--", linewidth=0.8)
        ax.text(x + 3, ax.get_ylim()[1] * 0.85, label, color="red", fontsize=7)
    ax.set_xlabel("Insert size (bp)")
    ax.set_ylabel("Count")
    ax.set_title("Insert Size Distribution (ATAC-seq)\nRed dashed: nucleosomal periodicity")
    path = plots_dir / "insert_size.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    if log:
        log.info(f"Insert size plot: {path}")
    return str(path)


def build_markdown(args, metrics: dict, skipped: list, plots: list,
                   blacklist_overlap: dict = None, frip_bedtools: dict = None,
                   chrom_warnings: list = None) -> str:
    lines = [
        "# ATAC-seq QC Report (Local)",
        f"\nGenerated: {datetime.now().isoformat()}",
        f"\n**Tool**: `{TOOL_NAME}` v{VERSION}",
        "\n## Inputs\n",
        f"- Fragments: `{args.fragments or 'not provided'}`",
        f"- Peaks BED: `{args.peaks or 'not provided'}`",
        f"- GTF/GFF: `{args.gtf or 'not provided'}`",
        f"- Blacklist BED: `{args.blacklist or 'not provided'}`",
        "\n## QC Metrics\n",
    ]

    if chrom_warnings:
        lines.append("\n## Chromosome Naming Warnings\n")
        for w in chrom_warnings:
            lines.append(f"- **{w['type']}**: {w['warning']}")

    if "fragments" in metrics:
        lines.append(f"- **Total fragments loaded**: {metrics['fragments'].get('n_fragments', '?')}")

    if "insert_sizes" in metrics:
        ins = metrics["insert_sizes"]
        lines.append(f"- **Median insert size**: {ins.get('median', '?')} bp")
        lines.append(f"- **Mean insert size**: {ins.get('mean', '?')} bp")

    if "frip" in metrics:
        frip_val = metrics["frip"].get("frip")
        if frip_val is not None:
            lines.append(f"- **FRiP (local)**: {frip_val:.4f}")
            lines.append(f"  - Fragments in peaks: {metrics['frip'].get('n_fragments_in_peaks')}")
            lines.append(f"  - Note: {metrics['frip'].get('note', '')}")
        else:
            lines.append(f"- **FRiP**: skipped - {metrics['frip'].get('skipped_reason')}")

    if frip_bedtools:
        if "frip_bedtools" in frip_bedtools:
            lines.append(f"- **FRiP (bedtools)**: {frip_bedtools['frip_bedtools']:.4f}")
            lines.append(f"  - Command: `{frip_bedtools.get('command', '')}`")
            lines.append(f"  - bedtools version: {frip_bedtools.get('tool_version', 'unknown')}")
            lines.append(f"  - Note: {frip_bedtools.get('note', '')}")
        elif "error" in frip_bedtools:
            lines.append(f"- **FRiP (bedtools)**: Error - {frip_bedtools['error']}")

    if "peaks" in metrics:
        lines.append(f"- **Number of peaks**: {metrics['peaks'].get('n_peaks', '?')}")
        lines.append(f"- **Median peak width**: {metrics['peaks'].get('median_peak_width_bp', '?')} bp")

    if "barcodes" in metrics:
        lines.append(f"- **Number of barcodes**: {metrics['barcodes'].get('n_barcodes', '?')}")
        lines.append(f"- **Median fragments per barcode**: {metrics['barcodes'].get('median_frags_per_barcode', '?')}")

    if blacklist_overlap:
        lines.append("\n## Blacklist Overlap (bedtools)\n")
        if "blacklist_fraction" in blacklist_overlap:
            lines.append(f"- **Blacklist fraction**: {blacklist_overlap['blacklist_fraction']:.4f}")
            lines.append(f"  - Blacklisted fragments: {blacklist_overlap.get('n_blacklisted_fragments')}")
            lines.append(f"  - Total fragments: {blacklist_overlap.get('n_total_fragments')}")
            lines.append(f"  - Command: `{blacklist_overlap.get('command', '')}`")
            lines.append(f"  - bedtools version: {blacklist_overlap.get('tool_version', 'unknown')}")
            lines.append(f"  - Note: {blacklist_overlap.get('note', '')}")
        elif "error" in blacklist_overlap:
            lines.append(f"- Error: {blacklist_overlap['error']}")

    lines.append("\n## Skipped Metrics\n")
    for item in skipped:
        lines.append(f"- **{item['metric']}**: {item['reason']}")
        if "missing_biological_conclusion" in item:
            lines.append(f"  - _Missing biological conclusion_: {item['missing_biological_conclusion']}")
        if "enable_with" in item:
            lines.append(f"  - Enable with: `{item['enable_with']}`")

    lines.append("\n## Plots\n")
    for p in plots:
        if p:
            lines.append(f"- `{p}`")
    if not any(plots):
        lines.append("- No plots generated.")

    lines.append("\n## Biological Interpretation Notes\n")
    lines.append("FRiP, TSS enrichment, and insert-size periodicity vary by:")
    lines.append("- Cell type composition (immune vs. epithelial vs. fibroblast)")
    lines.append("- Protocol (bulk ATAC vs. scATAC, nuclei vs. whole cell)")
    lines.append("- Tissue (open chromatin landscape differs dramatically)")
    lines.append("- Technical quality (transposition efficiency, lysis quality)")
    lines.append("\nDo not apply universal FRiP or TSS thresholds without context.")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Local ATAC-seq QC summarization tool.",
    )
    parser.add_argument("--fragments", help="Fragments file (.tsv or .tsv.gz)")
    parser.add_argument("--peaks", help="Peaks BED file (for FRiP)")
    parser.add_argument("--bam", help="BAM file (optional)")
    parser.add_argument("--gtf", help="GTF/GFF annotation (required for TSS enrichment)")
    parser.add_argument("--blacklist", help="Blacklist BED (required for blacklist fraction)")
    parser.add_argument("--genome-fasta", help="Genome FASTA (optional)")
    parser.add_argument("--chrom-sizes", help="Chromosome sizes file (optional)")
    parser.add_argument("--output-dir", default="reports/atac_qc")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--markdown", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logging(out_dir, args.verbose)

    has_bedtools = tool_available("bedtools")
    has_deeptools = tool_available("deeptools") or tool_available("bamCoverage")

    log.info(f"{TOOL_NAME} v{VERSION}")
    log.info(f"Python: {sys.version}")
    log.info(f"bedtools available: {has_bedtools}")
    log.info(f"deeptools available: {has_deeptools}")

    metrics = {}
    skipped = []
    plots = []
    blacklist_overlap = None
    frip_bedtools_result = None
    chrom_warnings = []

    # ── Chromosome mismatch check ─────────────────────────────────
    if any([args.fragments, args.peaks, args.blacklist, args.gtf]):
        chrom_warnings = check_chrom_mismatches(
            args.fragments, args.peaks, args.blacklist, args.gtf, log
        )

    # ── Fragments ─────────────────────────────────────────────────
    if args.fragments:
        frag_path = Path(args.fragments)
        if frag_path.exists():
            log.info("Reading fragments...")
            records = read_fragments(frag_path, log=log)
            metrics["fragments"] = {"n_fragments": len(records)}

            insert_sizes = compute_insert_sizes(records)
            if insert_sizes and HAS_NUMPY:
                metrics["insert_sizes"] = {
                    "n_sampled": len(insert_sizes),
                    "median": round(float(np.median(insert_sizes)), 1),
                    "mean": round(float(np.mean(insert_sizes)), 1),
                    "min": min(insert_sizes),
                    "max": max(insert_sizes),
                }
                plot_path = make_insert_size_plot(insert_sizes, out_dir, log)
                plots.append(plot_path)
            elif not HAS_NUMPY:
                skipped.append({
                    "metric": "Insert size statistics",
                    "reason": "numpy not available. Install: pip install numpy",
                    "missing_biological_conclusion": (
                        "Cannot assess nucleosomal banding pattern (mono/di/tri-nucleosomal periodicity). "
                        "Cannot detect over-digestion (short fragment enrichment) or under-digestion. "
                        "Cannot compare insert-size distributions across samples."
                    ),
                    "enable_with": "pip install numpy",
                })

            # Local FRiP (streaming Python fallback)
            metrics["frip"] = compute_frip(records, Path(args.peaks) if args.peaks else None, log)
            metrics["barcodes"] = compute_barcode_stats(records)

            # bedtools FRiP validation (authoritative, when available)
            if args.peaks and has_bedtools:
                peaks_path = Path(args.peaks)
                if peaks_path.exists():
                    log.info("Running bedtools FRiP validation...")
                    frip_bedtools_result = run_bedtools_frip_validation(frag_path, peaks_path, log)
            elif args.peaks and not has_bedtools:
                skipped.append({
                    "metric": "FRiP validation (bedtools)",
                    "reason": "bedtools not installed - local Python FRiP used only",
                    "missing_biological_conclusion": (
                        "FRiP computed by local Python interval overlap may be inaccurate for large datasets. "
                        "bedtools intersect provides authoritative, strand-aware overlap counting."
                    ),
                    "enable_with": "conda install -c bioconda bedtools",
                })

        else:
            log.warning(f"Fragments file not found: {frag_path}")
            skipped.append({
                "metric": "All fragment metrics",
                "reason": f"File not found: {frag_path}",
                "missing_biological_conclusion": (
                    "Cannot compute any ATAC QC metrics. Verify fragments file path and permissions."
                ),
            })
    else:
        skipped.append({
            "metric": "Fragment metrics",
            "reason": "No fragments file provided",
            "missing_biological_conclusion": (
                "Cannot compute fragment counts, insert sizes, FRiP, or barcode-level QC. "
                "Provide a fragments.tsv.gz file with --fragments."
            ),
        })

    # ── Peaks ─────────────────────────────────────────────────────
    if args.peaks:
        metrics["peaks"] = count_peaks(Path(args.peaks), log)
    else:
        skipped.append({
            "metric": "Peak count and width",
            "reason": "No peaks BED provided",
            "missing_biological_conclusion": (
                "Cannot summarize the accessible chromatin landscape. "
                "Cannot compute FRiP (signal-to-noise). "
                "Cannot assess peak count, peak width distribution, or peak-level QC. "
                "Provide a peaks BED file from MACS3 or nf-core/atacseq with --peaks."
            ),
            "enable_with": "--peaks /path/to/peaks.bed (from MACS3 or nf-core/atacseq output)",
        })

    # ── Blacklist fraction ────────────────────────────────────────
    if args.blacklist:
        blacklist_path = Path(args.blacklist)
        if args.fragments and blacklist_path.exists():
            frag_path = Path(args.fragments)
            if frag_path.exists():
                if has_bedtools:
                    log.info("Running bedtools blacklist overlap...")
                    blacklist_overlap = run_bedtools_blacklist_overlap(frag_path, blacklist_path, log)
                else:
                    skipped.append({
                        "metric": "Blacklist fraction (bedtools)",
                        "reason": "Blacklist BED provided but bedtools not installed.",
                        "missing_biological_conclusion": (
                            "Cannot remove artifactual peaks/fragments even though blacklist BED was provided. "
                            "Install bedtools to compute blacklist overlap fraction."
                        ),
                        "enable_with": "conda install -c bioconda bedtools",
                    })
            else:
                skipped.append({
                    "metric": "Blacklist fraction",
                    "reason": "Fragments file not found - cannot compute blacklist overlap.",
                    "missing_biological_conclusion": (
                        "Cannot identify artifactual fragments without a valid fragments file."
                    ),
                })
        else:
            skipped.append({
                "metric": "Blacklist fraction",
                "reason": "Blacklist BED provided but no fragments file - cannot compute overlap.",
                "missing_biological_conclusion": (
                    "Cannot identify artifactual fragments without a fragments file."
                ),
            })
    else:
        skipped.append({
            "metric": "Blacklist fraction",
            "reason": "Blacklist BED not provided",
            "missing_biological_conclusion": (
                "Cannot determine what fraction of reads or peaks overlap problematic genomic regions "
                "(satellite repeats, rDNA, centromeres, high-signal artifacts). "
                "Apparent accessibility in these regions is almost always artifactual, not biological."
            ),
            "enable_with": "--blacklist /path/to/blacklist.bed (ENCODE blacklist for hg38: https://github.com/Boyle-Lab/Blacklist)",
        })

    # ── TSS enrichment ────────────────────────────────────────────
    if not args.gtf:
        skipped.append({
            "metric": "TSS enrichment",
            "reason": "GTF/GFF annotation not provided",
            "missing_biological_conclusion": (
                "Cannot assess whether transposase insertion is enriched at gene regulatory elements. "
                "Cannot compare signal-to-noise across samples. "
                "Cannot identify globally low-quality libraries or globally accessible chromatin states."
            ),
            "enable_with": "--gtf /path/to/annotation.gtf (then use deeptools computeMatrix or nf-core/atacseq)",
        })
    elif args.bam and has_deeptools:
        # Document the command but do not auto-run on non-toy inputs
        skipped.append({
            "metric": "TSS enrichment",
            "reason": (
                "GTF and BAM provided and deeptools is available. "
                "TSS enrichment not auto-run to avoid unsafe computation on large files. "
                "Run manually: "
                "deeptools computeMatrix reference-point --referencePoint TSS "
                "-b 2000 -a 2000 -S <bigwig> -R <gtf> -o matrix.gz && "
                "deeptools plotHeatmap -m matrix.gz -o tss_enrichment.png"
            ),
            "missing_biological_conclusion": (
                "Cannot assess regulatory signal enrichment automatically. "
                "Run the deeptools command documented above on your BAM/bigwig."
            ),
            "enable_with": "deeptools computeMatrix reference-point --referencePoint TSS",
        })
    else:
        reason_parts = ["GTF provided but TSS enrichment calculation not implemented in local tool."]
        if not args.bam:
            reason_parts.append("BAM file not provided (required for deeptools TSS enrichment).")
        if not has_deeptools:
            reason_parts.append("deeptools not installed.")
        skipped.append({
            "metric": "TSS enrichment",
            "reason": " ".join(reason_parts),
            "missing_biological_conclusion": (
                "Cannot assess regulatory signal enrichment even though GTF was provided. "
                "Use nextflow-development@life-sciences with nf-core/atacseq or deeptools computeMatrix."
            ),
            "enable_with": "conda install -c bioconda deeptools  +  --bam aligned.bam",
        })

    # ── Build outputs ─────────────────────────────────────────────
    md_report = build_markdown(
        args, metrics, skipped, plots,
        blacklist_overlap=blacklist_overlap,
        frip_bedtools=frip_bedtools_result,
        chrom_warnings=chrom_warnings,
    )
    md_path = out_dir / "atac_qc_report.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md_report)
    log.info(f"Markdown report: {md_path}")

    summary = {
        "tool": TOOL_NAME,
        "version": VERSION,
        "python_version": sys.version,
        "bedtools_available": has_bedtools,
        "deeptools_available": has_deeptools,
        "generated": datetime.now().isoformat(),
        "inputs": vars(args),
        "metrics": metrics,
        "blacklist_overlap": blacklist_overlap,
        "frip_bedtools": frip_bedtools_result,
        "chrom_warnings": chrom_warnings,
        "skipped_metrics": skipped,
        "plots": plots,
    }
    json_path = out_dir / "atac_qc_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log.info(f"JSON summary: {json_path}")

    print(f"\nATAC QC complete. Outputs in: {out_dir}")
    if "fragments" in metrics:
        print(f"  Fragments loaded: {metrics['fragments'].get('n_fragments')}")
    if "frip" in metrics and metrics["frip"].get("frip") is not None:
        print(f"  FRiP (local): {metrics['frip']['frip']:.4f}")
    if frip_bedtools_result and frip_bedtools_result.get("frip_bedtools") is not None:
        print(f"  FRiP (bedtools): {frip_bedtools_result['frip_bedtools']:.4f}")
    if blacklist_overlap and blacklist_overlap.get("blacklist_fraction") is not None:
        print(f"  Blacklist fraction: {blacklist_overlap['blacklist_fraction']:.4f}")
    if chrom_warnings:
        print(f"  WARNING: {len(chrom_warnings)} chromosome naming mismatch(es) detected.")
    print(f"  bedtools available: {has_bedtools}")
    print(f"  Skipped metrics: {len(skipped)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
