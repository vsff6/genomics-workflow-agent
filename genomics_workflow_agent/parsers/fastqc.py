from __future__ import annotations

import io
import zipfile
from pathlib import Path
from typing import Any

TRACKED_MODULES = [
    "Basic Statistics",
    "Per base sequence quality",
    "Per sequence quality scores",
    "Per sequence GC content",
    "Sequence Length Distribution",
    "Adapter Content",
    "Overrepresented sequences",
]


def parse_fastqc_txt(text: str, source: str = "", sample: str = "") -> dict[str, Any]:
    """Parse the content of a fastqc_data.txt file into a structured dict."""
    modules: dict[str, dict] = {}
    summary: dict[str, list[str]] = {"pass": [], "warn": [], "fail": []}
    errors: list[str] = []

    current_module: str | None = None
    current_rows: list[list[str]] = []
    current_status: str | None = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip()

        if line.startswith(">>") and not line.startswith(">>END_MODULE"):
            parts = line.lstrip(">").split("\t", 1)
            current_module = parts[0].strip() if parts else None
            current_status = parts[1].strip().lower() if len(parts) > 1 else "unknown"
            current_rows = []
            continue

        if line == ">>END_MODULE":
            if current_module:
                modules[current_module] = {
                    "status": current_status,
                    "rows": current_rows[:],
                }
                if current_status in summary:
                    summary[current_status].append(current_module)
            current_module = None
            current_status = None
            current_rows = []
            continue

        if current_module and line and not line.startswith("#"):
            current_rows.append(line.split("\t"))

    return {
        "sample": sample,
        "source": source,
        "modules": modules,
        "summary": summary,
        "errors": errors,
        "parse_ok": len(errors) == 0,
    }


def parse_fastqc_zip(zip_path: str | Path) -> dict[str, Any]:
    """Parse a FastQC .zip output file and extract QC data."""
    zip_path = Path(zip_path)
    sample = zip_path.stem.replace("_fastqc", "")

    if not zip_path.exists():
        return _error_result(sample, str(zip_path), f"File not found: {zip_path}")

    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            names = zf.namelist()
            data_file = next(
                (n for n in names if n.endswith("fastqc_data.txt")), None
            )
            if data_file is None:
                return _error_result(
                    sample, str(zip_path),
                    "fastqc_data.txt not found inside zip"
                )
            text = zf.read(data_file).decode("utf-8", errors="replace")
    except zipfile.BadZipFile as e:
        return _error_result(sample, str(zip_path), f"BadZipFile: {e}")
    except Exception as e:
        return _error_result(sample, str(zip_path), f"Unexpected error: {e}")

    result = parse_fastqc_txt(text, source=str(zip_path), sample=sample)
    return result


def parse_fastqc_txt_file(txt_path: str | Path) -> dict[str, Any]:
    """Parse a FastQC fastqc_data.txt file directly."""
    txt_path = Path(txt_path)
    sample = txt_path.parent.name.replace("_fastqc", "")

    if not txt_path.exists():
        return _error_result(sample, str(txt_path), f"File not found: {txt_path}")

    try:
        text = txt_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return _error_result(sample, str(txt_path), f"Could not read file: {e}")

    return parse_fastqc_txt(text, source=str(txt_path), sample=sample)


def parse_fastqc_dir(fastqc_outdir: str | Path) -> list[dict[str, Any]]:
    """Parse all FastQC zip files found in a directory."""
    outdir = Path(fastqc_outdir)
    if not outdir.exists():
        return []
    results = []
    for zip_path in sorted(outdir.glob("*_fastqc.zip")):
        results.append(parse_fastqc_zip(zip_path))
    return results


def module_status(parsed: dict, module_name: str) -> str | None:
    """Return the status string for a named module, or None if absent."""
    return parsed.get("modules", {}).get(module_name, {}).get("status")


def _error_result(sample: str, source: str, error: str) -> dict[str, Any]:
    return {
        "sample": sample,
        "source": source,
        "modules": {},
        "summary": {"pass": [], "warn": [], "fail": []},
        "errors": [error],
        "parse_ok": False,
    }
