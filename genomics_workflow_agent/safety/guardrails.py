"""Safety guardrails: no clinical claims, no silent failures, no large-file loading."""

from __future__ import annotations

CLINICAL_CLAIM_PATTERNS = [
    "pathogenic",
    "likely pathogenic",
    "benign",
    "likely benign",
    "variant of uncertain significance",
    "vus",
    "disease-causing",
    "clinically significant",
    "diagnostic",
    "clinical interpretation",
    "medically actionable",
    "recommend treatment",
    "recommend surgery",
    "recommend medication",
    "diagnose",
    "prognosis",
]

# Files above this size are never read into memory directly.
LARGE_FILE_BYTES = 50 * 1024 * 1024  # 50 MB


def assert_no_clinical_claims(text: str) -> None:
    """Raise ValueError if text contains clinical claim language."""
    lower = text.lower()
    for pattern in CLINICAL_CLAIM_PATTERNS:
        if pattern in lower:
            raise ValueError(
                f"Clinical claim detected in output: '{pattern}'. "
                "This framework does not make clinical interpretations. "
                "Consult a certified clinical laboratory."
            )


def check_file_size(path: str) -> dict:
    """Return file size metadata without reading file content."""
    import os

    try:
        size = os.path.getsize(path)
        return {
            "path": path,
            "size_bytes": size,
            "size_mb": round(size / (1024 * 1024), 2),
            "safe_to_read": size < LARGE_FILE_BYTES,
        }
    except OSError as e:
        return {"path": path, "error": str(e), "safe_to_read": False}


def warn_large_file(path: str, size_bytes: int) -> str:
    size_mb = size_bytes / (1024 * 1024)
    return (
        f"WARNING: {path} is {size_mb:.1f} MB. "
        "Content will NOT be loaded into model context. "
        "Use dedicated CLI tools (samtools, bcftools, fastqc, etc.) for analysis."
    )


DISCLAIMER = (
    "This output is for research purposes only. "
    "It does not constitute clinical advice and must not be used for diagnostic, "
    "therapeutic, or medical decision-making without review by a qualified clinician."
)
