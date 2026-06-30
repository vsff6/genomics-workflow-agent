"""
wgs_vcf_qc_local.py - Local WGS/WES BAM/CRAM/VCF QC summarization.

Computes alignment and variant QC metrics. Uses pysam when available.
Runs samtools flagstat/idxstats/stats when samtools is installed.
Runs bcftools stats when bcftools is installed and cross-validates with local parser.
Never makes clinical claims. Reports every skipped metric and why.
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
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

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

VERSION = "2.0.0"
TOOL_NAME = "wgs_vcf_qc_local.py"
MAX_READS_SCAN = 500_000
MAX_VARIANTS_SCAN = 1_000_000
SAMTOOLS_DEPTH_SIZE_LIMIT = 10 * 1024 * 1024  # 10 MB


def setup_logging(output_dir: Path, verbose: bool) -> logging.Logger:
    log_path = output_dir / "wgs_vcf_qc.log"
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
    """Return the first line of --version output for a CLI tool."""
    try:
        result = subprocess.run([cmd, "--version"], capture_output=True, text=True, timeout=5)
        output = (result.stdout + result.stderr).strip()
        return output.split("\n")[0] if output else "unknown"
    except Exception:
        return "unknown"


def parse_flagstat_output(stdout: str) -> dict:
    """Parse samtools flagstat text output into a structured dict.

    Each line like '500000 + 0 in total (QC-passed reads + QC-failed reads)'
    is mapped to a key derived from the description.
    """
    parsed = {}
    for line in stdout.splitlines():
        line = line.strip()
        m = re.match(r"^(\d+) \+ (\d+) (.+)$", line)
        if m:
            qc_pass = int(m.group(1))
            qc_fail = int(m.group(2))
            desc = m.group(3).strip()
            # Key: text before first '(' with spaces → underscores
            key = desc.split("(")[0].strip().replace(" ", "_").replace("/", "_")
            parsed[key] = {"qc_passed": qc_pass, "qc_failed": qc_fail, "description": desc}

    # Derive mapping rate
    total_val = parsed.get("in_total", {}).get("qc_passed", 0)
    mapped_val = parsed.get("mapped", {}).get("qc_passed", 0)
    if total_val > 0:
        parsed["_pct_mapped"] = round(100.0 * mapped_val / total_val, 2)

    # Derive duplicate rate
    dup_val = parsed.get("duplicates", {}).get("qc_passed", 0)
    if total_val > 0:
        parsed["_pct_duplicate"] = round(100.0 * dup_val / total_val, 2)

    return parsed


def parse_idxstats_output(stdout: str) -> list:
    """Parse samtools idxstats output into a list of per-chromosome records."""
    records = []
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            records.append({
                "chrom": parts[0],
                "length": int(parts[1]) if parts[1].isdigit() else parts[1],
                "mapped": int(parts[2]),
                "unmapped": int(parts[3].strip()),
            })
        except ValueError:
            continue
    return records


def parse_samtools_stats_summary(stdout: str) -> dict:
    """Extract SN (summary numbers) section from samtools stats output."""
    summary = {}
    for line in stdout.splitlines():
        if not line.startswith("SN"):
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            key = parts[1].rstrip(":").strip()
            val_str = parts[2].strip() if len(parts) > 2 else ""
            try:
                val = float(val_str)
                summary[key] = int(val) if val == int(val) else val
            except ValueError:
                summary[key] = val_str
    return summary


def parse_bcftools_stats(stdout: str) -> dict:
    """Parse bcftools stats output, extracting SN (summary numbers) section."""
    summary = {}
    for line in stdout.splitlines():
        if not line.startswith("SN"):
            continue
        parts = line.split("\t")
        if len(parts) >= 4:
            key = parts[2].rstrip(":").strip()
            val_str = parts[3].strip()
            try:
                val = float(val_str)
                summary[key] = int(val) if val == int(val) else val
            except ValueError:
                summary[key] = val_str
    return summary


def compare_vcf_stats(local: dict, bcftools: dict) -> list:
    """Compare local parser metrics with bcftools stats; return list of discrepancy warnings."""
    warnings = []
    comparisons = [
        ("n_snps", "number of SNPs"),
        ("n_indels", "number of indels"),
    ]
    for local_key, bc_key in comparisons:
        lv = local.get(local_key)
        bv = bcftools.get(bc_key)
        if lv is not None and bv is not None:
            if lv != bv:
                warnings.append({
                    "metric": local_key,
                    "local_value": lv,
                    "bcftools_value": bv,
                    "warning": (
                        f"Local parser reports {lv} {local_key}, bcftools reports {bv} {bc_key}. "
                        "Discrepancy may be due to multiallelic handling or scan cap. "
                        "Trust bcftools for authoritative counts."
                    ),
                })

    # Ti/Tv comparison (local uses "ti_tv_ratio", bcftools uses "ts/tv")
    ltv = local.get("ti_tv_ratio")
    btv = bcftools.get("ts/tv")
    if ltv is not None and btv is not None:
        if abs(float(ltv) - float(btv)) > 0.05:
            warnings.append({
                "metric": "ti_tv_ratio",
                "local_value": ltv,
                "bcftools_value": btv,
                "warning": (
                    f"Ti/Tv discrepancy: local={ltv}, bcftools={btv}. "
                    "Likely caused by different multiallelic treatment or scan cap. "
                    "Trust bcftools for authoritative Ti/Tv."
                ),
            })
    return warnings


def run_samtools_flagstat(bam_path: Path, log) -> dict:
    """Run samtools flagstat and return structured result with command metadata."""
    cmd = ["samtools", "flagstat", str(bam_path)]
    log.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        out = {
            "command": " ".join(cmd),
            "return_code": result.returncode,
            "tool_version": get_tool_version("samtools"),
            "stderr_snippet": result.stderr[:500] if result.stderr else "",
        }
        if result.returncode == 0:
            out["raw_stdout"] = result.stdout
            out["parsed"] = parse_flagstat_output(result.stdout)
        else:
            out["error"] = result.stderr[:500]
        return out
    except subprocess.TimeoutExpired:
        return {"command": " ".join(cmd), "error": "samtools flagstat timed out (120s)"}
    except Exception as e:
        return {"command": " ".join(cmd), "error": str(e)}


def run_samtools_idxstats(bam_path: Path, log) -> dict:
    """Run samtools idxstats if a BAM index (.bai or .csi) exists."""
    bai_candidates = [
        Path(str(bam_path) + ".bai"),
        bam_path.with_suffix(".bai"),
        Path(str(bam_path) + ".csi"),
        bam_path.with_suffix(".csi"),
    ]
    has_index = any(p.exists() for p in bai_candidates)
    if not has_index:
        return {
            "skipped": True,
            "reason": "No BAM index found (.bai or .csi). Run: samtools index <bam>",
            "missing_biological_conclusion": (
                "Cannot assess per-chromosome read distribution without an index. "
                "Cannot verify that all expected chromosomes were sequenced."
            ),
            "enable_with": "samtools index <bam>",
        }

    cmd = ["samtools", "idxstats", str(bam_path)]
    log.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        out = {
            "command": " ".join(cmd),
            "return_code": result.returncode,
            "stderr_snippet": result.stderr[:500] if result.stderr else "",
        }
        if result.returncode == 0:
            out["raw_stdout"] = result.stdout
            out["parsed"] = parse_idxstats_output(result.stdout)
        else:
            out["error"] = result.stderr[:500]
        return out
    except subprocess.TimeoutExpired:
        return {"command": " ".join(cmd), "error": "samtools idxstats timed out (60s)"}
    except Exception as e:
        return {"command": " ".join(cmd), "error": str(e)}


def run_samtools_stats(bam_path: Path, log) -> dict:
    """Run samtools stats and extract summary numbers."""
    cmd = ["samtools", "stats", str(bam_path)]
    log.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        out = {
            "command": " ".join(cmd),
            "return_code": result.returncode,
            "stderr_snippet": result.stderr[:500] if result.stderr else "",
        }
        if result.returncode == 0:
            out["parsed_summary"] = parse_samtools_stats_summary(result.stdout)
        else:
            out["error"] = result.stderr[:500]
        return out
    except subprocess.TimeoutExpired:
        return {"command": " ".join(cmd), "error": "samtools stats timed out (300s)"}
    except Exception as e:
        return {"command": " ".join(cmd), "error": str(e)}


def run_samtools_depth(bam_path: Path, intervals: str, log) -> dict:
    """Run samtools depth only for tiny files (<10 MB) or when --intervals is provided."""
    try:
        file_size = bam_path.stat().st_size
    except Exception:
        file_size = 0

    if not intervals and file_size > SAMTOOLS_DEPTH_SIZE_LIMIT:
        size_mb = file_size // (1024 * 1024)
        return {
            "skipped": True,
            "reason": (
                f"File size {size_mb} MB exceeds safety threshold for samtools depth "
                f"({SAMTOOLS_DEPTH_SIZE_LIMIT // (1024 * 1024)} MB). "
                "Provide --intervals BED to restrict depth calculation."
            ),
            "missing_biological_conclusion": (
                "Cannot assess per-base coverage depth without restricting to intervals. "
                "Coverage uniformity and low-coverage regions cannot be identified."
            ),
            "enable_with": "--intervals /path/to/targets.bed",
        }

    cmd = ["samtools", "depth"]
    if intervals:
        cmd += ["-b", intervals]
    cmd.append(str(bam_path))
    log.info(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        out = {
            "command": " ".join(cmd),
            "return_code": result.returncode,
            "stderr_snippet": result.stderr[:500] if result.stderr else "",
        }
        if result.returncode == 0:
            depths = []
            for line in result.stdout.splitlines():
                parts = line.split("\t")
                if len(parts) >= 3:
                    try:
                        depths.append(int(parts[2]))
                    except ValueError:
                        pass
            if depths:
                depths_sorted = sorted(depths)
                n = len(depths_sorted)
                out["depth_summary"] = {
                    "n_positions": n,
                    "mean": round(sum(depths_sorted) / n, 2),
                    "median": depths_sorted[n // 2],
                    "min": depths_sorted[0],
                    "max": depths_sorted[-1],
                }
            else:
                out["depth_summary"] = {"n_positions": 0}
        else:
            out["error"] = result.stderr[:500]
        return out
    except subprocess.TimeoutExpired:
        return {"command": " ".join(cmd), "error": "samtools depth timed out (300s)"}
    except Exception as e:
        return {"command": " ".join(cmd), "error": str(e)}


def run_bcftools_stats(vcf_path: Path, log) -> dict:
    """Run bcftools stats and return structured result with command metadata."""
    cmd = ["bcftools", "stats", str(vcf_path)]
    log.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        out = {
            "command": " ".join(cmd),
            "return_code": result.returncode,
            "tool_version": get_tool_version("bcftools"),
            "stderr_snippet": result.stderr[:500] if result.stderr else "",
        }
        if result.returncode == 0:
            out["parsed"] = parse_bcftools_stats(result.stdout)
        else:
            out["error"] = result.stderr[:500]
        return out
    except subprocess.TimeoutExpired:
        return {"command": " ".join(cmd), "error": "bcftools stats timed out (120s)"}
    except Exception as e:
        return {"command": " ".join(cmd), "error": str(e)}


def bam_qc_pysam(bam_path: Path, file_type: str, log) -> dict:
    """Compute BAM/CRAM QC metrics using pysam."""
    result = {"file_type": file_type, "path": str(bam_path)}

    mode = "rb" if file_type == "BAM" else "rc"
    try:
        n_total = 0
        n_mapped = 0
        n_unmapped = 0
        n_duplicate = 0
        n_supplementary = 0
        n_secondary = 0
        mapq_sum = 0

        with pysam.AlignmentFile(str(bam_path), mode) as f:
            header = f.header.to_dict()
            sq = header.get("SQ", [])
            result["n_references"] = len(sq)
            result["references_sample"] = [s["SN"] for s in sq[:5]]
            rg = header.get("RG", [])
            result["read_groups"] = [{"ID": r.get("ID"), "SM": r.get("SM"), "LB": r.get("LB")} for r in rg]

            for i, read in enumerate(f.fetch(until_eof=True)):
                if i >= MAX_READS_SCAN:
                    log.warning(f"Scan capped at {MAX_READS_SCAN} reads.")
                    break
                n_total += 1
                if not read.is_unmapped:
                    n_mapped += 1
                    if read.mapping_quality is not None:
                        mapq_sum += read.mapping_quality
                else:
                    n_unmapped += 1
                if read.is_duplicate:
                    n_duplicate += 1
                if read.is_supplementary:
                    n_supplementary += 1
                if read.is_secondary:
                    n_secondary += 1

        result["n_reads_scanned"] = n_total
        result["n_mapped"] = n_mapped
        result["n_unmapped"] = n_unmapped
        result["n_duplicate"] = n_duplicate
        result["n_supplementary"] = n_supplementary
        result["n_secondary"] = n_secondary
        if n_total > 0:
            result["pct_mapped"] = round(100 * n_mapped / n_total, 2)
            result["pct_duplicate"] = round(100 * n_duplicate / n_total, 2)
        if n_mapped > 0:
            result["mean_mapq"] = round(mapq_sum / n_mapped, 1)
        result["note"] = (
            f"Scanned first {min(n_total, MAX_READS_SCAN)} reads. "
            "For full stats run: samtools flagstat"
        )
    except Exception as e:
        result["error"] = str(e)
        log.error(f"pysam error on {bam_path}: {e}")

    return result


def parse_vcf(vcf_path: Path, log) -> dict:
    """Parse VCF for variant counts and quality metrics."""
    result = {
        "path": str(vcf_path),
        "n_snps": 0,
        "n_indels": 0,
        "n_other": 0,
        "n_filtered": 0,
        "n_pass": 0,
        "n_multiallelic": 0,
        "quals": [],
        "depths": [],
        "alt_afs": [],
        "transitions": 0,
        "transversions": 0,
        "n_het": 0,
        "n_hom_ref": 0,
        "n_hom_alt": 0,
    }

    TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    n_scanned = 0
    samples = []

    try:
        with opener(vcf_path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    cols = line.rstrip().split("\t")
                    samples = cols[9:]
                    result["n_samples"] = len(samples)
                    result["samples"] = samples[:5]
                    continue
                if not line.strip():
                    continue
                if n_scanned >= MAX_VARIANTS_SCAN:
                    log.warning(f"VCF scan capped at {MAX_VARIANTS_SCAN} variants.")
                    break

                parts = line.rstrip().split("\t")
                if len(parts) < 8:
                    continue

                ref = parts[3]
                alt_field = parts[4]
                qual = parts[5]
                filt = parts[6]
                info = parts[7]

                alts = [a for a in alt_field.split(",") if a != "."]
                if len(alts) > 1:
                    result["n_multiallelic"] += 1

                if filt in ("PASS", "."):
                    result["n_pass"] += 1
                else:
                    result["n_filtered"] += 1

                try:
                    result["quals"].append(float(qual))
                except ValueError:
                    pass

                for field in info.split(";"):
                    if field.startswith("DP="):
                        try:
                            result["depths"].append(int(field[3:]))
                        except ValueError:
                            pass
                        break

                for field in info.split(";"):
                    if field.startswith("AF="):
                        try:
                            result["alt_afs"].append(float(field[3:].split(",")[0]))
                        except ValueError:
                            pass
                        break

                main_alt = alts[0] if alts else ""
                if len(ref) == 1 and len(main_alt) == 1 and main_alt not in (".", "*"):
                    result["n_snps"] += 1
                    pair = (ref.upper(), main_alt.upper())
                    if pair in TRANSITIONS:
                        result["transitions"] += 1
                    else:
                        result["transversions"] += 1
                elif main_alt not in (".", "*", ""):
                    result["n_indels"] += 1
                else:
                    result["n_other"] += 1

                if len(parts) > 9:
                    fmt = parts[8].split(":")
                    gt_idx = fmt.index("GT") if "GT" in fmt else None
                    if gt_idx is not None:
                        gt_field = parts[9].split(":")[gt_idx] if len(parts[9].split(":")) > gt_idx else ""
                        gt = gt_field.replace("|", "/")
                        if gt == "0/0":
                            result["n_hom_ref"] += 1
                        elif gt in ("0/1", "1/0"):
                            result["n_het"] += 1
                        elif gt in ("1/1",):
                            result["n_hom_alt"] += 1

                n_scanned += 1

    except Exception as e:
        result["error"] = str(e)
        log.error(f"VCF parse error: {e}")
        return result

    result["n_variants_scanned"] = n_scanned
    result["n_total"] = result["n_snps"] + result["n_indels"] + result["n_other"]

    if result["transversions"] > 0:
        result["ti_tv_ratio"] = round(result["transitions"] / result["transversions"], 3)
    else:
        result["ti_tv_ratio"] = None

    if (result["n_het"] + result["n_hom_alt"]) > 0:
        result["het_hom_ratio"] = round(result["n_het"] / (result["n_het"] + result["n_hom_alt"]), 3)
    else:
        result["het_hom_ratio"] = None

    def dist_summary(vals):
        if not vals:
            return {}
        vals_s = sorted(vals)
        n = len(vals_s)
        return {
            "n": n,
            "min": round(vals_s[0], 2),
            "median": round(vals_s[n // 2], 2),
            "mean": round(sum(vals_s) / n, 2),
            "max": round(vals_s[-1], 2),
        }

    result["qual_summary"] = dist_summary(result.pop("quals", []))
    result["depth_summary"] = dist_summary(result.pop("depths", []))
    result["af_summary"] = dist_summary(result.pop("alt_afs", []))

    result["note"] = (
        "Ti/Tv interpretation: WGS ~2.0-2.1, WES ~2.8-3.0. "
        "Het/hom ratio interpretation varies by population and ploidy. "
        "No clinical interpretation provided."
    )

    return result


def make_variant_plots(vcf_metrics: dict, output_dir: Path, log) -> list:
    if not HAS_MATPLOTLIB:
        return []

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    generated = []

    fig, ax = plt.subplots(figsize=(6, 4))
    labels = ["SNPs", "Indels", "Other"]
    vals = [vcf_metrics.get("n_snps", 0), vcf_metrics.get("n_indels", 0), vcf_metrics.get("n_other", 0)]
    ax.bar(labels, vals, color=["#2196F3", "#4CAF50", "#FF9800"])
    ax.set_title("Variant Type Counts")
    ax.set_ylabel("Count")
    path = plots_dir / "variant_types.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    generated.append(str(path))

    return generated


def build_markdown(args, bam_metrics: dict, vcf_metrics: dict, skipped: list, plots: list,
                   samtools_qc: dict = None, bcftools_metrics: dict = None,
                   vcf_discrepancies: list = None) -> str:
    lines = [
        "# WGS/WES QC Report (Local)",
        f"\nGenerated: {datetime.now().isoformat()}",
        f"\n**Tool**: `{TOOL_NAME}` v{VERSION}",
        "\n**IMPORTANT: No clinical claims are made in this report.**",
        "\n## Inputs\n",
        f"- BAM: `{args.bam or 'not provided'}`",
        f"- CRAM: `{args.cram or 'not provided'}`",
        f"- VCF: `{args.vcf or 'not provided'}`",
        f"- Genome FASTA: `{args.genome_fasta or 'not provided'}`",
    ]

    if bam_metrics:
        lines.append("\n## Alignment QC (pysam)\n")
        for k, v in bam_metrics.items():
            if k not in ("error", "stdout") and not isinstance(v, list):
                lines.append(f"- **{k}**: {v}")
        if "error" in bam_metrics:
            lines.append(f"\n**Error**: {bam_metrics['error']}")

    if samtools_qc:
        lines.append("\n## Alignment QC (samtools)\n")

        flagstat = samtools_qc.get("flagstat", {})
        if flagstat and not flagstat.get("skipped"):
            lines.append(f"### samtools flagstat\n")
            lines.append(f"- **Command**: `{flagstat.get('command', '')}`")
            lines.append(f"- **Return code**: {flagstat.get('return_code', '?')}")
            lines.append(f"- **Tool version**: {flagstat.get('tool_version', 'unknown')}")
            parsed = flagstat.get("parsed", {})
            pct_mapped = parsed.get("_pct_mapped")
            pct_dup = parsed.get("_pct_duplicate")
            if pct_mapped is not None:
                lines.append(f"- **Mapping rate**: {pct_mapped}%")
            if pct_dup is not None:
                lines.append(f"- **Duplicate rate**: {pct_dup}%")
            if "error" in flagstat:
                lines.append(f"- **Error**: {flagstat['error']}")

        idxstats = samtools_qc.get("idxstats", {})
        if idxstats:
            if idxstats.get("skipped"):
                lines.append(f"\n### samtools idxstats\n")
                lines.append(f"- Skipped: {idxstats['reason']}")
                lines.append(f"  - _Missing biological conclusion_: {idxstats.get('missing_biological_conclusion', '')}")
            else:
                lines.append(f"\n### samtools idxstats\n")
                lines.append(f"- **Command**: `{idxstats.get('command', '')}`")
                records = idxstats.get("parsed", [])
                if records:
                    total_mapped = sum(r["mapped"] for r in records)
                    n_chroms = len([r for r in records if r["chrom"] != "*"])
                    lines.append(f"- **Chromosomes**: {n_chroms}")
                    lines.append(f"- **Total mapped (from idxstats)**: {total_mapped}")

        depth = samtools_qc.get("depth", {})
        if depth:
            if depth.get("skipped"):
                lines.append(f"\n### samtools depth\n")
                lines.append(f"- Skipped: {depth['reason']}")
                lines.append(f"  - _Missing biological conclusion_: {depth.get('missing_biological_conclusion', '')}")
            elif "depth_summary" in depth:
                lines.append(f"\n### samtools depth\n")
                lines.append(f"- **Command**: `{depth.get('command', '')}`")
                ds = depth["depth_summary"]
                lines.append(f"- **Mean depth**: {ds.get('mean')}")
                lines.append(f"- **Median depth**: {ds.get('median')}")
                lines.append(f"- **Positions assessed**: {ds.get('n_positions')}")

    if vcf_metrics:
        lines.append("\n## Variant QC (Local Parser)\n")
        ti_tv = vcf_metrics.get("ti_tv_ratio")
        het_hom = vcf_metrics.get("het_hom_ratio")
        lines.append(f"- **Total variants scanned**: {vcf_metrics.get('n_variants_scanned')}")
        lines.append(f"- **SNPs**: {vcf_metrics.get('n_snps')}")
        lines.append(f"- **Indels**: {vcf_metrics.get('n_indels')}")
        lines.append(f"- **PASS variants**: {vcf_metrics.get('n_pass')}")
        lines.append(f"- **Filtered variants**: {vcf_metrics.get('n_filtered')}")
        lines.append(f"- **Multiallelic sites**: {vcf_metrics.get('n_multiallelic')}")
        lines.append(f"- **Ti/Tv ratio**: {ti_tv} (WGS expected ~2.0-2.1; WES ~2.8-3.0)")
        lines.append(f"- **Het/hom ratio**: {het_hom} (varies by population and ploidy)")
        lines.append(f"- **Samples**: {vcf_metrics.get('n_samples')}")

        qs = vcf_metrics.get("qual_summary", {})
        if qs:
            lines.append(f"- **QUAL**: median={qs.get('median')}, mean={qs.get('mean')}, min={qs.get('min')}, max={qs.get('max')}")

        ds = vcf_metrics.get("depth_summary", {})
        if ds:
            lines.append(f"- **Depth (DP in INFO)**: median={ds.get('median')}, mean={ds.get('mean')}, min={ds.get('min')}, max={ds.get('max')}")

        if vcf_metrics.get("note"):
            lines.append(f"\n_{vcf_metrics['note']}_")

    if bcftools_metrics:
        lines.append("\n## Variant QC (bcftools stats)\n")
        if bcftools_metrics.get("skipped"):
            lines.append(f"- Skipped: {bcftools_metrics['reason']}")
            lines.append(f"  - _Missing biological conclusion_: {bcftools_metrics.get('missing_biological_conclusion', '')}")
        else:
            lines.append(f"- **Command**: `{bcftools_metrics.get('command', '')}`")
            lines.append(f"- **Return code**: {bcftools_metrics.get('return_code', '?')}")
            lines.append(f"- **Tool version**: {bcftools_metrics.get('tool_version', 'unknown')}")
            parsed = bcftools_metrics.get("parsed", {})
            for k, v in sorted(parsed.items()):
                lines.append(f"- **{k}**: {v}")
            if "error" in bcftools_metrics:
                lines.append(f"- **Error**: {bcftools_metrics['error']}")

    if vcf_discrepancies:
        lines.append("\n## Local vs bcftools Discrepancies\n")
        for d in vcf_discrepancies:
            lines.append(f"- **{d['metric']}**: local={d['local_value']}, bcftools={d['bcftools_value']}")
            lines.append(f"  - {d['warning']}")

    lines.append("\n## Skipped Metrics\n")
    for s in skipped:
        lines.append(f"- **{s['metric']}**: {s['reason']}")
        if "missing_biological_conclusion" in s:
            lines.append(f"  - _Missing biological conclusion_: {s['missing_biological_conclusion']}")
        if "enable_with" in s:
            lines.append(f"  - Enable with: `{s['enable_with']}`")

    lines.append("\n## Plots\n")
    for p in plots:
        lines.append(f"- `{p}`")
    if not plots:
        lines.append("- No plots generated.")

    lines.append("\n## Limitations and Caveats\n")
    lines.append("- Coverage calculation requires samtools depth or mosdepth.")
    lines.append("- Contamination estimation requires VerifyBamID or similar.")
    lines.append("- Ancestry and population inference not performed.")
    lines.append("- No clinical interpretation is provided. Medical decisions must never be based on this output.")
    lines.append("- Variant annotation must use a dedicated tool (VEP, ANNOVAR) with a validated annotation source.")
    lines.append("- For production WGS/WES analysis, use nf-core/sarek via nextflow-development@life-sciences.")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Local WGS/WES BAM/CRAM/VCF QC summarization tool.",
    )
    parser.add_argument("--bam", help="BAM file")
    parser.add_argument("--cram", help="CRAM file")
    parser.add_argument("--vcf", help="VCF file (gzipped or plain)")
    parser.add_argument("--genome-fasta", help="Reference genome FASTA (for CRAM decoding)")
    parser.add_argument("--annotation-vcf", help="Annotation VCF (optional)")
    parser.add_argument("--known-sites-vcf", help="Known sites VCF (optional)")
    parser.add_argument("--intervals", help="Target intervals BED (WES)")
    parser.add_argument("--output-dir", default="reports/wgs_qc")
    parser.add_argument("--json", action="store_true", help="Write JSON summary (default: always written)")
    parser.add_argument("--markdown", action="store_true", help="Write Markdown report (default: always written)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logging(out_dir, args.verbose)

    has_samtools = tool_available("samtools")
    has_bcftools = tool_available("bcftools")

    log.info(f"{TOOL_NAME} v{VERSION}")
    log.info(f"pysam available: {HAS_PYSAM}")
    log.info(f"samtools available: {has_samtools}")
    log.info(f"bcftools available: {has_bcftools}")

    bam_metrics = {}
    vcf_metrics = {}
    samtools_qc = {}
    bcftools_metrics = {}
    vcf_discrepancies = []
    skipped = []
    plots = []

    # ── BAM/CRAM QC ──────────────────────────────────────────────
    bam_path = Path(args.bam) if args.bam else (Path(args.cram) if args.cram else None)
    file_type = "BAM" if args.bam else ("CRAM" if args.cram else None)

    if bam_path and bam_path.exists():
        if HAS_PYSAM:
            log.info(f"Running pysam QC on {bam_path}")
            bam_metrics = bam_qc_pysam(bam_path, file_type, log)
        elif has_samtools:
            log.info("pysam not available. Using samtools flagstat as primary BAM QC.")
            flagstat_result = run_samtools_flagstat(bam_path, log)
            bam_metrics = {
                "source": "samtools flagstat",
                "command": flagstat_result.get("command"),
                "return_code": flagstat_result.get("return_code"),
                "tool_version": flagstat_result.get("tool_version"),
            }
            if "parsed" in flagstat_result:
                parsed = flagstat_result["parsed"]
                bam_metrics["pct_mapped"] = parsed.get("_pct_mapped")
                bam_metrics["pct_duplicate"] = parsed.get("_pct_duplicate")
            if "error" in flagstat_result:
                bam_metrics["error"] = flagstat_result["error"]
        else:
            skipped.append({
                "metric": "BAM/CRAM alignment QC",
                "reason": "pysam and samtools not available. Install pysam: pip install pysam",
                "missing_biological_conclusion": (
                    "Cannot assess mapping rate, duplicate rate, or read-level quality. "
                    "Cannot determine whether alignment quality is sufficient for variant calling."
                ),
                "enable_with": "pip install pysam  OR  conda install -c bioconda samtools",
            })


        if has_samtools:
            log.info("Running samtools extended QC suite...")
            samtools_qc["samtools_version"] = get_tool_version("samtools")
            samtools_qc["flagstat"] = run_samtools_flagstat(bam_path, log)
            samtools_qc["idxstats"] = run_samtools_idxstats(bam_path, log)
            samtools_qc["stats"] = run_samtools_stats(bam_path, log)
            samtools_qc["depth"] = run_samtools_depth(bam_path, args.intervals, log)
        else:
            skipped.append({
                "metric": "samtools extended QC (flagstat/idxstats/stats/depth)",
                "reason": "samtools not installed",
                "missing_biological_conclusion": (
                    "Cannot produce authoritative per-chromosome read counts (idxstats), "
                    "insert size distribution, or per-base coverage depth. "
                    "Install samtools for structured alignment QC."
                ),
                "enable_with": "conda install -c bioconda samtools",
            })

    elif bam_path:
        skipped.append({
            "metric": "BAM/CRAM QC",
            "reason": f"File not found: {bam_path}",
            "missing_biological_conclusion": (
                "Cannot assess alignment quality. Verify BAM/CRAM path and index (.bai/.crai)."
            ),
        })
    else:
        skipped.append({
            "metric": "BAM/CRAM QC",
            "reason": "No BAM or CRAM provided",
            "missing_biological_conclusion": (
                "Cannot assess mapping rate, duplicate rate, insert size, or read-level quality metrics. "
                "Provide a BAM or CRAM file with --bam or --cram."
            ),
            "enable_with": "--bam /path/to/aligned.bam or --cram /path/to/aligned.cram",
        })

    # ── VCF QC ──────────────────────────────────────────────────
    if args.vcf:
        vcf_path = Path(args.vcf)
        if vcf_path.exists():
            log.info(f"Parsing VCF (local parser): {vcf_path}")
            vcf_metrics = parse_vcf(vcf_path, log)
            plots = make_variant_plots(vcf_metrics, out_dir, log)

            if has_bcftools:
                log.info("Running bcftools stats for cross-validation...")
                bcftools_metrics = run_bcftools_stats(vcf_path, log)
                if "parsed" in bcftools_metrics and vcf_metrics:
                    vcf_discrepancies = compare_vcf_stats(vcf_metrics, bcftools_metrics["parsed"])
                    if vcf_discrepancies:
                        for d in vcf_discrepancies:
                            log.warning(d["warning"])
            else:
                bcftools_metrics = {
                    "skipped": True,
                    "reason": "bcftools not installed - local parser used only",
                    "missing_biological_conclusion": (
                        "Cannot cross-validate variant counts, Ti/Tv, or indel rates with bcftools. "
                        "Local parser results are provided but not independently validated. "
                        "Install bcftools for authoritative VCF statistics."
                    ),
                    "enable_with": "conda install -c bioconda bcftools",
                }
        else:
            skipped.append({
                "metric": "VCF QC",
                "reason": f"File not found: {vcf_path}",
                "missing_biological_conclusion": "Cannot assess variant counts or quality. Verify VCF path.",
            })
    else:
        skipped.append({
            "metric": "VCF QC",
            "reason": "No VCF provided",
            "missing_biological_conclusion": (
                "Cannot compute variant counts, Ti/Tv ratio, het/hom ratio, or quality distributions. "
                "Provide a VCF file with --vcf."
            ),
            "enable_with": "--vcf /path/to/variants.vcf.gz",
        })

    # ── Always-skipped metrics ────────────────────────────────────
    skipped.append({
        "metric": "Coverage statistics",
        "reason": "Requires mosdepth or samtools depth with full BAM. Not computed in local tool.",
        "missing_biological_conclusion": (
            "Cannot determine whether sequencing depth is sufficient for confident variant calling. "
            "Cannot identify under-covered regions where variants may be missed. "
            "Cannot assess GC-bias or target enrichment uniformity (WES). "
            "Recommended minimum: 30x for WGS germline, 100x for WES, 500x+ for somatic low-frequency variants."
        ),
        "enable_with": "mosdepth --quantize 0:5:10:30: output.prefix input.bam",
    })

    if not args.annotation_vcf:
        skipped.append({
            "metric": "Variant annotation (consequences)",
            "reason": "No annotation VCF or annotation tool configured.",
            "missing_biological_conclusion": (
                "Cannot determine predicted functional consequences (synonymous, missense, frameshift, splice). "
                "Cannot assess population allele frequency. "
                "NOTE: annotation alone does not constitute clinical interpretation."
            ),
            "enable_with": "VEP (Ensembl Variant Effect Predictor) or ANNOVAR with appropriate databases",
        })


    md = build_markdown(
        args, bam_metrics, vcf_metrics, skipped, plots,
        samtools_qc=samtools_qc,
        bcftools_metrics=bcftools_metrics,
        vcf_discrepancies=vcf_discrepancies,
    )
    md_path = out_dir / "wgs_vcf_qc_report.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md)
    log.info(f"Markdown report: {md_path}")

    summary = {
        "tool": TOOL_NAME,
        "version": VERSION,
        "python_version": sys.version,
        "pysam_available": HAS_PYSAM,
        "samtools_available": has_samtools,
        "bcftools_available": has_bcftools,
        "generated": datetime.now().isoformat(),
        "inputs": vars(args),
        "bam_metrics": bam_metrics,
        "vcf_metrics": vcf_metrics,
        "samtools_qc": samtools_qc,
        "bcftools_metrics": bcftools_metrics,
        "vcf_discrepancies": vcf_discrepancies,
        "skipped_metrics": skipped,
        "plots": plots,
        "clinical_disclaimer": "No clinical claims are made. Do not use for medical decisions.",
    }
    json_path = out_dir / "wgs_vcf_qc_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log.info(f"JSON summary: {json_path}")

    print(f"\nWGS/VCF QC complete. Outputs in: {out_dir}")
    if vcf_metrics:
        print(f"  Variants: {vcf_metrics.get('n_total')} (SNPs: {vcf_metrics.get('n_snps')}, Indels: {vcf_metrics.get('n_indels')})")
        print(f"  Ti/Tv: {vcf_metrics.get('ti_tv_ratio')}")
    if vcf_discrepancies:
        print(f"  WARNING: {len(vcf_discrepancies)} local vs bcftools discrepancy(ies) - see JSON.")
    print(f"  samtools available: {has_samtools}  bcftools available: {has_bcftools}")
    print(f"  Skipped metrics: {len(skipped)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
