"""
reference_validator.py - Validate local reference files before analysis.

Checks existence, compression, indexing, chromosome naming, and checksums.
Never downloads large files. Reports all issues clearly.
"""

import argparse
import gzip
import hashlib
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

VERSION = "1.0.0"
TOOL_NAME = "reference_validator.py"

KNOWN_HUMAN_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
KNOWN_HUMAN_CHROMS_NOCHR = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}
KNOWN_MOUSE_CHROMS = {f"chr{i}" for i in range(1, 20)} | {"chrX", "chrY", "chrM"}


def setup_logging(output_dir: Path, verbose: bool) -> logging.Logger:
    log_path = output_dir / "reference_validator.log"
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stderr)],
    )
    return logging.getLogger(TOOL_NAME)


def file_md5(path: Path, max_bytes: int = 50 * 1024 * 1024) -> str:
    """Compute MD5 of first max_bytes of file (fast approximation for large files)."""
    h = hashlib.md5()
    try:
        with open(path, "rb") as f:
            data = f.read(max_bytes)
            h.update(data)
        return h.hexdigest() + ("_partial" if path.stat().st_size > max_bytes else "")
    except Exception as e:
        return f"error:{e}"


def is_gzipped(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def safe_head(path: Path, n: int = 10) -> list:
    gz = is_gzipped(path)
    opener = gzip.open if gz else open
    lines = []
    try:
        with opener(path, "rt", errors="replace") as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                lines.append(line.rstrip())
    except Exception as e:
        lines = [f"ERROR: {e}"]
    return lines


def detect_chrom_style(chroms: list) -> str:
    """Detect chromosome naming style from a list of chrom names."""
    ucsc = sum(1 for c in chroms if c.startswith("chr"))
    ensembl = sum(1 for c in chroms if not c.startswith("chr") and c in KNOWN_HUMAN_CHROMS_NOCHR)
    if ucsc > ensembl:
        return "UCSC (chr prefix)"
    elif ensembl > 0:
        return "Ensembl (no chr prefix)"
    else:
        return "unknown"


def validate_fasta(path: Path) -> dict:
    result = {"path": str(path), "type": "genome FASTA", "status": "ok", "issues": [], "warnings": []}

    if not path.exists():
        result["status"] = "MISSING"
        result["issues"].append("File does not exist.")
        return result

    result["size_mb"] = round(path.stat().st_size / (1024 * 1024), 1)
    result["compressed"] = is_gzipped(path)
    if result["compressed"]:
        result["warnings"].append("Genome FASTA is gzipped. Some tools require uncompressed FASTA with .fai index.")

    # Check index
    fai = Path(str(path) + ".fai")
    dict_file = path.with_suffix(".dict")
    result["has_fai"] = fai.exists()
    result["has_dict"] = dict_file.exists()
    if not result["has_fai"]:
        result["warnings"].append("No .fai index found. Run: samtools faidx <fasta>")
    if not result["has_dict"]:
        result["warnings"].append("No .dict file found. Run: samtools dict <fasta> > <fasta>.dict")

    # Sample chromosomes
    chroms = []
    try:
        opener = gzip.open if result["compressed"] else open
        with opener(path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    chroms.append(line[1:].split()[0])
                if len(chroms) >= 30:
                    break
    except Exception as e:
        result["issues"].append(f"Error reading FASTA: {e}")

    result["chromosomes_sample"] = chroms[:10]
    result["chrom_style"] = detect_chrom_style(chroms)
    result["n_chroms_in_sample"] = len(chroms)

    result["partial_md5"] = file_md5(path)
    return result


def validate_gtf(path: Path) -> dict:
    result = {"path": str(path), "type": "GTF/GFF", "status": "ok", "issues": [], "warnings": []}

    if not path.exists():
        result["status"] = "MISSING"
        result["issues"].append("File does not exist.")
        return result

    result["size_mb"] = round(path.stat().st_size / (1024 * 1024), 1)
    result["compressed"] = is_gzipped(path)

    lines = safe_head(path, 20)
    meta = [l for l in lines if l.startswith("#")]
    data = [l for l in lines if not l.startswith("#") and l.strip()]

    result["meta_lines_in_sample"] = len(meta)
    result["meta_sample"] = meta[:5]

    # Try to detect chrom style and gene_id format
    chroms = []
    gene_ids = []
    for l in data:
        parts = l.split("\t")
        if len(parts) >= 9:
            chroms.append(parts[0])
            # Try to extract gene_id
            attrs = parts[8]
            if "gene_id" in attrs:
                try:
                    gid = attrs.split("gene_id")[1].split('"')[1]
                    gene_ids.append(gid)
                except Exception:
                    pass

    result["chrom_style"] = detect_chrom_style(chroms)
    result["chroms_sample"] = list(set(chroms))[:5]
    result["gene_id_sample"] = gene_ids[:5]

    if gene_ids:
        if gene_ids[0].startswith("ENSG") or gene_ids[0].startswith("ENSMUSG"):
            result["gene_id_style"] = "Ensembl stable IDs"
        elif gene_ids[0].startswith("NM_") or gene_ids[0].startswith("NR_"):
            result["gene_id_style"] = "RefSeq"
        else:
            result["gene_id_style"] = "possibly gene symbols or other"
        result["warnings"].append(
            f"Gene ID style '{result['gene_id_style']}' - ensure consistency with count matrix feature names."
        )

    # Check tabix index for compressed files
    if result["compressed"]:
        tbi = Path(str(path) + ".tbi")
        csi = Path(str(path) + ".csi")
        result["has_tabix"] = tbi.exists() or csi.exists()
        if not result["has_tabix"]:
            result["warnings"].append("Compressed GTF has no tabix index. Run: tabix -p gff <file>")

    result["partial_md5"] = file_md5(path)
    return result


def validate_bed(path: Path, label: str = "BED") -> dict:
    result = {"path": str(path), "type": label, "status": "ok", "issues": [], "warnings": []}

    if not path.exists():
        result["status"] = "MISSING"
        result["issues"].append("File does not exist.")
        return result

    result["size_mb"] = round(path.stat().st_size / (1024 * 1024), 1)
    result["compressed"] = is_gzipped(path)

    lines = safe_head(path, 10)
    data = [l for l in lines if not l.startswith("#") and not l.startswith("track") and l.strip()]

    chroms = [l.split("\t")[0] for l in data if "\t" in l]
    result["chrom_style"] = detect_chrom_style(chroms)
    result["chroms_sample"] = chroms[:5]

    if result["compressed"]:
        tbi = Path(str(path) + ".tbi")
        result["has_tabix"] = tbi.exists()

    result["partial_md5"] = file_md5(path)
    return result


def validate_vcf(path: Path, label: str = "VCF") -> dict:
    result = {"path": str(path), "type": label, "status": "ok", "issues": [], "warnings": []}

    if not path.exists():
        result["status"] = "MISSING"
        result["issues"].append("File does not exist.")
        return result

    result["size_mb"] = round(path.stat().st_size / (1024 * 1024), 1)
    result["compressed"] = is_gzipped(path)

    lines = safe_head(path, 20)
    meta = [l for l in lines if l.startswith("##")]
    header = next((l for l in lines if l.startswith("#CHROM")), None)
    result["n_meta_lines"] = len(meta)
    result["has_header"] = header is not None

    # Check for contig lines
    contig_lines = [l for l in meta if l.startswith("##contig")]
    result["n_contig_lines"] = len(contig_lines)
    if contig_lines:
        chrom = contig_lines[0].split("ID=")[-1].split(",")[0].rstrip(">")
        result["chrom_style"] = detect_chrom_style([chrom])

    if result["compressed"]:
        tbi = Path(str(path) + ".tbi")
        result["has_tabix"] = tbi.exists()
        if not result["has_tabix"]:
            result["warnings"].append("Compressed VCF has no tabix index. Run: tabix -p vcf <file>")

    result["partial_md5"] = file_md5(path)
    return result


def validate_chrom_sizes(path: Path) -> dict:
    result = {"path": str(path), "type": "chromosome sizes", "status": "ok", "issues": [], "warnings": []}

    if not path.exists():
        result["status"] = "MISSING"
        result["issues"].append("File does not exist.")
        return result

    lines = safe_head(path, 30)
    chroms = []
    for l in lines:
        if l and "\t" in l:
            chroms.append(l.split("\t")[0])

    result["chrom_style"] = detect_chrom_style(chroms)
    result["chroms_sample"] = chroms[:10]
    result["n_chroms_in_sample"] = len(chroms)
    return result


def validate_marker_list(path: Path) -> dict:
    result = {"path": str(path), "type": "marker gene list", "status": "ok", "issues": [], "warnings": []}

    if not path.exists():
        result["status"] = "MISSING"
        result["issues"].append("File does not exist.")
        return result

    lines = safe_head(path, 10)
    result["first_lines"] = lines[:5]
    return result


def check_chrom_compatibility(results: list) -> list:
    """Warn if different files have inconsistent chromosome naming."""
    styles = {}
    for r in results:
        style = r.get("chrom_style")
        if style and style != "unknown":
            styles[r["type"]] = style

    conflicts = []
    style_set = set(styles.values())
    if len(style_set) > 1:
        conflicts.append(
            f"WARNING: Inconsistent chromosome naming styles detected! "
            f"Files have: {styles}. This will cause tool failures. "
            f"Ensure all files use the same naming convention (UCSC 'chr' prefix or Ensembl no-prefix)."
        )
    return conflicts


def build_markdown(results: list, conflicts: list) -> str:
    lines = [
        "# Reference File Validation Report",
        f"\nGenerated: {datetime.now().isoformat()}",
        f"\n**Tool**: `{TOOL_NAME}` v{VERSION}",
        "\n## Summary\n",
    ]

    # Table
    lines.append("| File | Type | Status | Size (MB) | Chrom Style | Issues |")
    lines.append("|------|------|--------|-----------|-------------|--------|")
    for r in results:
        fname = Path(r["path"]).name if r["path"] else "N/A"
        status = r.get("status", "ok")
        size = r.get("size_mb", "N/A")
        chrom = r.get("chrom_style", "N/A")
        issues = "; ".join(r.get("issues", []))[:80]
        lines.append(f"| {fname} | {r['type']} | {status} | {size} | {chrom} | {issues} |")

    if conflicts:
        lines.append("\n## Chromosome Naming Conflicts\n")
        for c in conflicts:
            lines.append(f"- **{c}**")

    lines.append("\n## Details\n")
    for r in results:
        lines.append(f"### {r['type']}: `{Path(r['path']).name}`\n")
        for k, v in r.items():
            if k not in ("type", "path", "issues", "warnings", "partial_md5") and not isinstance(v, list):
                lines.append(f"- **{k}**: {v}")
        if r.get("partial_md5"):
            lines.append(f"- **MD5 (partial)**: `{r['partial_md5']}`")
        if r.get("issues"):
            for issue in r["issues"]:
                lines.append(f"- **ISSUE**: {issue}")
        if r.get("warnings"):
            for w in r["warnings"]:
                lines.append(f"- **WARNING**: {w}")
        lines.append("")

    lines.append("\n## Rules\n")
    lines.append("- Tools must not automatically download large reference files.")
    lines.append("- If a reference is missing, skip dependent metrics and report the limitation.")
    lines.append("- Never invent genome build or annotation version.")
    lines.append("- Chromosome naming must be consistent across all reference files and input data.")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Validate local reference files before genomics analysis.",
    )
    parser.add_argument("--genome-fasta", help="Genome FASTA")
    parser.add_argument("--gtf", help="GTF/GFF annotation")
    parser.add_argument("--chrom-sizes", help="Chromosome sizes file")
    parser.add_argument("--blacklist", help="Blacklist BED")
    parser.add_argument("--known-sites-vcf", help="Known sites VCF")
    parser.add_argument("--annotation-vcf", help="Annotation VCF")
    parser.add_argument("--marker-list", help="Marker gene list")
    parser.add_argument("--output-dir", default="reports/reference_check")
    parser.add_argument("--json", action="store_true", help="Write JSON summary (default: always written)")
    parser.add_argument("--markdown", action="store_true", help="Write Markdown report (default: always written)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logging(out_dir, args.verbose)

    log.info(f"{TOOL_NAME} v{VERSION}")

    results = []

    if args.genome_fasta:
        results.append(validate_fasta(Path(args.genome_fasta)))
    if args.gtf:
        results.append(validate_gtf(Path(args.gtf)))
    if args.chrom_sizes:
        results.append(validate_chrom_sizes(Path(args.chrom_sizes)))
    if args.blacklist:
        results.append(validate_bed(Path(args.blacklist), "blacklist BED"))
    if args.known_sites_vcf:
        results.append(validate_vcf(Path(args.known_sites_vcf), "known-sites VCF"))
    if args.annotation_vcf:
        results.append(validate_vcf(Path(args.annotation_vcf), "annotation VCF"))
    if args.marker_list:
        results.append(validate_marker_list(Path(args.marker_list)))

    if not results:
        log.warning("No reference files provided. Pass at least one --genome-fasta, --gtf, etc.")
        results = [{"path": "none", "type": "none", "status": "no input", "issues": [], "warnings": ["No reference files provided."]}]

    conflicts = check_chrom_compatibility(results)

    md = build_markdown(results, conflicts)
    md_path = out_dir / "reference_validation.md"
    with open(md_path, "w") as f:
        f.write(md)
    log.info(f"Markdown report: {md_path}")

    summary = {
        "tool": TOOL_NAME,
        "version": VERSION,
        "python_version": sys.version,
        "generated": datetime.now().isoformat(),
        "results": results,
        "conflicts": conflicts,
    }
    json_path = out_dir / "reference_validation.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log.info(f"JSON summary: {json_path}")

    n_issues = sum(len(r.get("issues", [])) for r in results)
    n_warnings = sum(len(r.get("warnings", [])) for r in results)
    print(f"\nReference validation complete. Outputs in: {out_dir}")
    print(f"  Files checked: {len(results)}")
    print(f"  Issues: {n_issues}")
    print(f"  Warnings: {n_warnings}")
    print(f"  Chromosome conflicts: {len(conflicts)}")

    return 1 if n_issues > 0 or conflicts else 0


if __name__ == "__main__":
    sys.exit(main())
