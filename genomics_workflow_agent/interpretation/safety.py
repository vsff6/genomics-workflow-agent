"""Safety checks: scan interpretation output for forbidden clinical terms."""
from __future__ import annotations

import json

FORBIDDEN_TERMS: list[str] = [
    "pathogenic",
    "benign",
    "likely pathogenic",
    "likely benign",
    "diagnostic",
    "diagnosis",
    "disease-causing",
    "treatment",
    "therapy recommendation",
    "clinical action",
    "medical action",
]


def scan_text(text: str) -> list[str]:
    """Return any forbidden clinical terms found in text (case-insensitive)."""
    lower = text.lower()
    return [t for t in FORBIDDEN_TERMS if t in lower]


def scan_result_dict(result: dict) -> list[str]:
    """Flatten a result dict to JSON text and scan for forbidden clinical terms."""
    return scan_text(json.dumps(result, default=str))


def raise_if_clinical_claims(result: dict) -> None:
    """Raise ValueError if any forbidden clinical terms appear in the result."""
    found = scan_result_dict(result)
    if found:
        raise ValueError(
            f"Forbidden clinical terms found in interpretation output: {found}. "
            "Remove these terms before returning the result."
        )
