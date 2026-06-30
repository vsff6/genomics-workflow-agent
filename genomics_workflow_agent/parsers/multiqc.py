from __future__ import annotations

import csv
import io
import json
from pathlib import Path
from typing import Any


def parse_multiqc_fastqc_tsv(multiqc_data_dir: str | Path) -> dict[str, Any]:
    """
    Parse multiqc_data/multiqc_fastqc.txt (tab-separated, written by MultiQC
    after aggregating FastQC results) into per-sample status dicts.
    """
    data_dir = Path(multiqc_data_dir)
    tsv_path = data_dir / "multiqc_fastqc.txt"
    errors: list[str] = []
    samples: dict[str, dict] = {}

    if not tsv_path.exists():
        return {
            "source": str(tsv_path),
            "samples": {},
            "errors": [f"multiqc_fastqc.txt not found in {data_dir}"],
            "parse_ok": False,
        }

    try:
        text = tsv_path.read_text(encoding="utf-8", errors="replace")
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        for row in reader:
            sample = row.get("Sample", "").strip()
            if sample:
                samples[sample] = dict(row)
    except Exception as e:
        errors.append(f"Error parsing {tsv_path}: {e}")

    return {
        "source": str(tsv_path),
        "samples": samples,
        "errors": errors,
        "parse_ok": len(errors) == 0 and bool(samples),
    }


def parse_multiqc_general_stats(multiqc_data_dir: str | Path) -> dict[str, Any]:
    """Parse multiqc_data/multiqc_general_stats.txt for cross-tool summary stats."""
    data_dir = Path(multiqc_data_dir)
    tsv_path = data_dir / "multiqc_general_stats.txt"
    samples: dict[str, dict] = {}
    errors: list[str] = []

    if not tsv_path.exists():
        return {
            "source": str(tsv_path),
            "samples": {},
            "errors": [f"multiqc_general_stats.txt not found in {data_dir}"],
            "parse_ok": False,
        }

    try:
        text = tsv_path.read_text(encoding="utf-8", errors="replace")
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        for row in reader:
            sample = row.get("Sample", "").strip()
            if sample:
                samples[sample] = dict(row)
    except Exception as e:
        errors.append(f"Error parsing {tsv_path}: {e}")

    return {
        "source": str(tsv_path),
        "samples": samples,
        "errors": errors,
        "parse_ok": len(errors) == 0,
    }


def parse_multiqc_json(multiqc_dir: str | Path) -> dict[str, Any]:
    """Parse multiqc_data.json if present."""
    out_dir = Path(multiqc_dir)
    json_path = out_dir / "multiqc_data" / "multiqc_data.json"
    if not json_path.exists():
        json_path = out_dir / "multiqc_data.json"
    if not json_path.exists():
        return {
            "source": str(json_path),
            "data": {},
            "errors": ["multiqc_data.json not found"],
            "parse_ok": False,
        }
    try:
        data = json.loads(json_path.read_text(encoding="utf-8"))
        return {"source": str(json_path), "data": data, "errors": [], "parse_ok": True}
    except Exception as e:
        return {
            "source": str(json_path),
            "data": {},
            "errors": [f"JSON parse error: {e}"],
            "parse_ok": False,
        }


def parse_multiqc_output(multiqc_dir: str | Path) -> dict[str, Any]:
    """
    Parse all available MultiQC output files from a MultiQC output directory.
    Missing files are noted in errors rather than raised.
    """
    out_dir = Path(multiqc_dir)
    data_dir = out_dir / "multiqc_data"

    result: dict[str, Any] = {
        "multiqc_dir": str(out_dir),
        "html_report": None,
        "fastqc_per_sample": {},
        "general_stats": {},
        "raw_json": {},
        "errors": [],
        "parse_ok": False,
    }

    html = out_dir / "multiqc_report.html"
    if html.exists():
        result["html_report"] = str(html)

    if data_dir.exists():
        fastqc_parsed = parse_multiqc_fastqc_tsv(data_dir)
        result["fastqc_per_sample"] = fastqc_parsed.get("samples", {})
        result["errors"].extend(fastqc_parsed.get("errors", []))

        stats_parsed = parse_multiqc_general_stats(data_dir)
        result["general_stats"] = stats_parsed.get("samples", {})
        result["errors"].extend(stats_parsed.get("errors", []))
    else:
        result["errors"].append(f"multiqc_data/ directory not found in {out_dir}")

    json_parsed = parse_multiqc_json(out_dir)
    if json_parsed["parse_ok"]:
        result["raw_json"] = json_parsed["data"]
    else:
        result["errors"].extend(json_parsed.get("errors", []))

    result["parse_ok"] = html.exists()
    return result
