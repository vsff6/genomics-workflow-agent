"""Detects genomic file types and returns metadata. Never reads large files."""

from __future__ import annotations

import gzip
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from genomics_workflow_agent.safety.guardrails import LARGE_FILE_BYTES
from genomics_workflow_agent.tools.files import classify_file, discover_files

SAMPLE_LINES = 5


def _is_gzipped(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def _safe_head(path: Path, n: int = SAMPLE_LINES) -> list[str]:
    """Read first n non-empty lines without loading full file."""
    lines: list[str] = []
    try:
        compressed = _is_gzipped(path)
        opener = gzip.open if compressed else open
        mode = "rt"
        with opener(path, mode, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line:
                    lines.append(line[:300])
                if len(lines) >= n:
                    break
    except Exception as e:
        lines = [f"[read error: {e}]"]
    return lines


def _detect_fastq_encoding(head_lines: list[str]) -> str:
    """Guess Phred encoding from quality line (line 4 of FASTQ)."""
    if len(head_lines) >= 4:
        qual = head_lines[3]
        scores = [ord(c) - 33 for c in qual if 33 <= ord(c) <= 126]
        if scores:
            max_q = max(scores)
            if max_q > 40:
                return "Phred+33 (Illumina 1.8+, likely)"
            return "Phred+33 (likely)"
    return "unknown"


def _infer_assay(file_type: str, path: Path) -> str:
    name = path.name.lower()
    if file_type == "fastq":
        if "atac" in name:
            return "ATAC-seq (inferred from filename)"
        if "rna" in name or "mrna" in name or "cdna" in name:
            return "RNA-seq (inferred from filename)"
        if "wgs" in name or "wes" in name or "dna" in name:
            return "WGS/WES (inferred from filename)"
        if "16s" in name or "its" in name or "amplicon" in name or "microbiome" in name:
            return "Amplicon/Microbiome (inferred from filename)"
        return "unknown - check experimental records"
    if file_type == "alignment":
        return "Sequencing alignment (BAM/CRAM/SAM)"
    if file_type == "vcf":
        return "Variant calls (VCF/BCF)"
    if file_type == "h5":
        if "h5ad" in name:
            return "Single-cell (AnnData h5ad)"
        return "HDF5 (10x or similar)"
    if file_type == "annotation":
        return "Gene annotation (GTF/GFF)"
    if file_type == "fasta":
        return "Genome/Transcriptome sequence"
    return "unknown"


def inspect_file(path: str | Path) -> dict[str, Any]:
    """
    Inspect a single genomic file without loading large content.
    Returns structured metadata dict.
    """
    path = Path(path)
    result: dict[str, Any] = {
        "path": str(path),
        "name": path.name,
        "exists": path.exists(),
        "inspected_at": datetime.now(timezone.utc).isoformat(),
        "type": None,
        "compressed": None,
        "size_bytes": None,
        "size_mb": None,
        "safe_to_read": None,
        "head_lines": None,
        "inferred_assay": None,
        "warnings": [],
        "errors": [],
    }

    if not path.exists():
        result["errors"].append(f"File not found: {path}")
        return result

    if path.is_dir():
        # Check for 10x directory structure
        barcodes = path / "barcodes.tsv.gz"
        matrix = path / "matrix.mtx.gz"
        features = path / "features.tsv.gz"
        if barcodes.exists() or matrix.exists():
            result["type"] = "10x_directory"
            result["inferred_assay"] = "Single-cell (10x Genomics)"
            contents = [p.name for p in path.iterdir()]
            result["directory_contents"] = contents[:30]
            return result
        result["type"] = "directory"
        result["errors"].append("Path is a directory but not a recognized 10x structure")
        return result

    size = path.stat().st_size
    result["size_bytes"] = size
    result["size_mb"] = round(size / (1024 * 1024), 2)
    result["safe_to_read"] = size < LARGE_FILE_BYTES
    result["compressed"] = _is_gzipped(path)
    result["type"] = classify_file(path)

    if result["safe_to_read"] and result["type"] not in ("alignment", "h5"):
        result["head_lines"] = _safe_head(path)
    elif not result["safe_to_read"]:
        result["warnings"].append(
            f"File is {result['size_mb']} MB - content not loaded. Use CLI tools for analysis."
        )

    # FASTQ-specific
    if result["type"] == "fastq" and result.get("head_lines"):
        if len(result["head_lines"]) >= 1 and result["head_lines"][0].startswith("@"):
            result["fastq_valid_header"] = True
            result["fastq_encoding_guess"] = _detect_fastq_encoding(result["head_lines"])
        else:
            result["fastq_valid_header"] = False
            result["warnings"].append("First line does not start with '@' - may not be a valid FASTQ")

    result["inferred_assay"] = _infer_assay(result["type"], path)

    return result


def inspect_directory(
    input_path: str | Path, recursive: bool = False
) -> dict[str, Any]:
    """
    Inspect all files in a directory and return a summary.
    """
    input_path = Path(input_path)
    files = discover_files(input_path, recursive=recursive)

    type_counts: dict[str, int] = {}
    total_size = 0
    file_details = []

    for f in files:
        t = f["type"]
        type_counts[t] = type_counts.get(t, 0) + 1
        total_size += f["size_bytes"]
        file_details.append(f)

    # Infer likely workflow
    workflow_guess = _infer_workflow(type_counts, files)

    return {
        "input_path": str(input_path),
        "inspected_at": datetime.now(timezone.utc).isoformat(),
        "total_files": len(files),
        "total_size_mb": round(total_size / (1024 * 1024), 2),
        "file_type_counts": type_counts,
        "workflow_guess": workflow_guess,
        "files": file_details,
        "warnings": _directory_warnings(files),
    }


def _infer_workflow(type_counts: dict, files: list) -> str:
    names = " ".join(f["name"].lower() for f in files)
    if type_counts.get("fastq", 0) > 0:
        if "atac" in names:
            return "atacseq"
        if "16s" in names or "its" in names or "amplicon" in names or "microbiome" in names:
            return "amplicon"
        if "wgs" in names or "wes" in names:
            return "variant-qc"
        if "rna" in names or "mrna" in names:
            return "rnaseq"
        return "fastq-qc"
    if type_counts.get("alignment", 0) > 0 or type_counts.get("vcf", 0) > 0:
        return "variant-qc"
    if type_counts.get("h5", 0) > 0:
        return "scrna-qc"
    return "unknown"


def _directory_warnings(files: list) -> list[str]:
    warnings = []
    large = [f for f in files if f["size_mb"] > 50]
    if large:
        warnings.append(
            f"{len(large)} file(s) exceed 50 MB - content not read. Use CLI tools for analysis."
        )
    unknown = [f for f in files if f["type"] == "unknown"]
    if unknown:
        warnings.append(f"{len(unknown)} file(s) have unrecognized types: {[f['name'] for f in unknown[:5]]}")
    return warnings
