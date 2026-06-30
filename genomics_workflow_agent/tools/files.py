"""File discovery and validation utilities."""

from __future__ import annotations

import os
import re
from pathlib import Path

FASTQ_EXTS = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
FASTA_EXTS = {".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz"}
BAM_EXTS = {".bam", ".cram", ".sam"}
VCF_EXTS = {".vcf", ".vcf.gz", ".bcf"}
BED_EXTS = {".bed", ".bed.gz"}
GTF_EXTS = {".gtf", ".gff", ".gff3", ".gtf.gz", ".gff.gz"}
COUNT_EXTS = {".csv", ".tsv", ".txt", ".mtx"}
H5_EXTS = {".h5", ".h5ad", ".loom"}


def _match_ext(path: Path, ext_set: set) -> bool:
    name = path.name.lower()
    return any(name.endswith(e) for e in ext_set)


def classify_file(path: Path) -> str:
    """Return a broad file type label for a genomics file."""
    if _match_ext(path, FASTQ_EXTS):
        return "fastq"
    if _match_ext(path, FASTA_EXTS):
        return "fasta"
    if _match_ext(path, BAM_EXTS):
        return "alignment"
    if _match_ext(path, VCF_EXTS):
        return "vcf"
    if _match_ext(path, BED_EXTS):
        return "bed"
    if _match_ext(path, GTF_EXTS):
        return "annotation"
    if _match_ext(path, H5_EXTS):
        return "h5"
    if _match_ext(path, COUNT_EXTS):
        return "tabular"
    return "unknown"


def discover_files(input_path: str | Path, recursive: bool = False) -> list[dict]:
    """
    Discover genomics-relevant files under input_path.
    Returns list of {path, type, size_bytes, size_mb}.
    Never reads file content.
    """
    input_path = Path(input_path)
    results = []

    if input_path.is_file():
        candidates = [input_path]
    elif input_path.is_dir():
        if recursive:
            candidates = [p for p in input_path.rglob("*") if p.is_file()]
        else:
            candidates = [p for p in input_path.iterdir() if p.is_file()]
    else:
        return results

    for p in sorted(candidates):
        file_type = classify_file(p)
        size = p.stat().st_size
        results.append(
            {
                "path": str(p),
                "name": p.name,
                "type": file_type,
                "size_bytes": size,
                "size_mb": round(size / (1024 * 1024), 2),
            }
        )

    return results


def is_r1(path: Path) -> bool:
    return bool(re.search(r"[_.]R?1[_.]|[_.]R?1$|[_.]R?1\.", path.name))


def is_r2(path: Path) -> bool:
    return bool(re.search(r"[_.]R?2[_.]|[_.]R?2$|[_.]R?2\.", path.name))


def sample_name(path: Path) -> str:
    name = path.name
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return re.sub(r"[_.]R?[12]$", "", name)


def find_fastq_pairs(directory: Path) -> list[dict]:
    """
    Find FASTQ files and detect paired-end structure.
    Returns list of {sample, r1, r2, paired}.
    """
    fastqs = sorted(p for p in directory.iterdir() if p.is_file() and _match_ext(p, FASTQ_EXTS))
    r1_files = [f for f in fastqs if is_r1(f)]
    r2_by_name = {sample_name(f): f for f in fastqs if is_r2(f)}
    seen = set()
    pairs = []
    for r1 in r1_files:
        sname = sample_name(r1)
        r2 = r2_by_name.get(sname)
        pairs.append({"sample": sname, "r1": str(r1), "r2": str(r2) if r2 else None, "paired": r2 is not None})
        seen.add(sname)
    # Single-end (no R1/R2 suffix)
    for f in fastqs:
        if not is_r1(f) and not is_r2(f):
            sname = sample_name(f)
            if sname not in seen:
                pairs.append({"sample": sname, "r1": str(f), "r2": None, "paired": False})
                seen.add(sname)
    return pairs
