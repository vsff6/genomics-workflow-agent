"""Parse mosdepth summary outputs."""
from __future__ import annotations

from pathlib import Path
from typing import Any


def parse_mosdepth_summary(text: str, sample: str = "", source: str = "") -> dict[str, Any]:
    """Parse mosdepth *.mosdepth.summary.txt output."""
    result: dict[str, Any] = {
        "sample": sample,
        "source": source,
        "parse_ok": False,
        "errors": [],
        "mean_coverage": None,
        "regions": [],
    }

    if not text or not text.strip():
        result["errors"].append("Empty mosdepth summary output")
        return result

    parsed_any = False
    for line in text.strip().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        # Skip header row
        if parts[0].lower() in ("chrom", "chr"):
            continue
        try:
            chrom = parts[0]
            length = int(parts[1])
            bases = int(parts[2])
            mean = float(parts[3])
        except (ValueError, IndexError):
            continue

        parsed_any = True
        result["regions"].append({
            "chrom": chrom,
            "length": length,
            "bases": bases,
            "mean": mean,
        })

        if chrom.lower() in ("total", "all"):
            result["mean_coverage"] = mean

    if parsed_any:
        result["parse_ok"] = True
    else:
        result["errors"].append("No recognizable mosdepth summary lines found")

    return result


def parse_mosdepth_summary_file(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    sample = (
        path.stem
        .replace(".mosdepth.summary", "")
        .replace("_mosdepth", "")
    )
    try:
        return parse_mosdepth_summary(
            path.read_text(encoding="utf-8"), sample=sample, source=str(path)
        )
    except Exception as e:
        return {"sample": sample, "source": str(path), "parse_ok": False, "errors": [str(e)]}


def parse_mosdepth_dir(qc_dir: str | Path) -> list[dict[str, Any]]:
    """Parse all mosdepth summary files from a directory."""
    qc_dir = Path(qc_dir)
    if not qc_dir.exists():
        return []
    return [
        parse_mosdepth_summary_file(p)
        for p in sorted(qc_dir.glob("*.mosdepth.summary.txt"))
    ]
