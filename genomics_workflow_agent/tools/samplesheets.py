"""Samplesheet generators — thin wrappers over nfcore_launcher.py logic."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from genomics_workflow_agent.tools.files import find_fastq_pairs

_log = logging.getLogger(__name__)


def build_rnaseq_samplesheet(input_dir: Path, output_path: Path) -> dict:
    pairs = find_fastq_pairs(input_dir)
    if not pairs:
        return {"created": False, "reason": f"No FASTQ files found in {input_dir}", "rows": 0}
    rows = [
        {
            "sample": p["sample"],
            "fastq_1": p["r1"],
            "fastq_2": p["r2"] or "",
            "strandedness": "auto",
        }
        for p in pairs
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        writer.writerows(rows)
    return {
        "created": True,
        "path": str(output_path),
        "rows": len(rows),
        "warnings": [
            "strandedness set to 'auto' — confirm from library preparation records",
            "verify sample names, pairing, and file paths before execution",
        ],
    }


def build_atacseq_samplesheet(input_dir: Path, output_path: Path) -> dict:
    pairs = find_fastq_pairs(input_dir)
    if not pairs:
        return {"created": False, "reason": f"No FASTQ files found in {input_dir}", "rows": 0}
    rows = [
        {
            "sample": p["sample"],
            "fastq_1": p["r1"],
            "fastq_2": p["r2"] or "",
            "replicate": "1",
        }
        for p in pairs
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "replicate"])
        writer.writeheader()
        writer.writerows(rows)
    return {
        "created": True,
        "path": str(output_path),
        "rows": len(rows),
        "warnings": [
            "replicate numbers are placeholders — must reflect actual experimental design",
            "blacklist BED must be provided separately",
            "genome build and GTF must be confirmed before execution",
        ],
    }


def build_sarek_samplesheet(input_dir: Path, output_path: Path) -> dict:
    pairs = find_fastq_pairs(input_dir)
    if not pairs:
        return {"created": False, "reason": f"No FASTQ files found in {input_dir}", "rows": 0}
    rows = [
        {
            "patient": "PATIENT_ID",
            "sex": "unknown",
            "status": "0",
            "sample": p["sample"],
            "lane": "L001",
            "fastq_1": p["r1"],
            "fastq_2": p["r2"] or "",
        }
        for p in pairs
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2"]
        )
        writer.writeheader()
        writer.writerows(rows)
    return {
        "created": True,
        "path": str(output_path),
        "rows": len(rows),
        "warnings": [
            "PATIENT_ID placeholders must be replaced with real patient identifiers",
            "tumor/normal status (0/1) cannot be inferred from filenames — set manually",
            "sex must be confirmed from sample metadata",
            "do not execute without manual review",
        ],
    }


def build_amplicon_samplesheet(input_dir: Path, output_path: Path) -> dict:
    """
    Build a minimal amplicon samplesheet (compatible with nf-core/ampliseq format).
    Columns: sampleID, forwardReads, reverseReads, run
    """
    pairs = find_fastq_pairs(input_dir)
    if not pairs:
        return {"created": False, "reason": f"No FASTQ files found in {input_dir}", "rows": 0}
    rows = [
        {
            "sampleID": p["sample"],
            "forwardReads": p["r1"],
            "reverseReads": p["r2"] or "",
            "run": "1",
        }
        for p in pairs
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sampleID", "forwardReads", "reverseReads", "run"])
        writer.writeheader()
        writer.writerows(rows)
    return {
        "created": True,
        "path": str(output_path),
        "rows": len(rows),
        "warnings": [
            "run numbers are placeholders — update to match actual sequencing run IDs",
            "primer sequences must be provided separately (--FW_primer / --RV_primer)",
            "taxonomy database path must be set before execution",
        ],
    }
