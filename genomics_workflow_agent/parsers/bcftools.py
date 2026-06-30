"""Parse bcftools stats output."""
from __future__ import annotations

from pathlib import Path
from typing import Any


def parse_bcftools_stats(text: str, sample: str = "", source: str = "") -> dict[str, Any]:
    """Parse bcftools stats text output (SN summary lines)."""
    result: dict[str, Any] = {
        "sample": sample,
        "source": source,
        "parse_ok": False,
        "errors": [],
        "n_samples": None,
        "n_records": None,
        "n_snps": None,
        "n_indels": None,
        "n_multiallelic": None,
        "ts_tv": None,
        "raw_sn": {},
    }

    if not text or not text.strip():
        result["errors"].append("Empty bcftools stats output")
        return result

    # bcftools stats SN lines: SN\t0\tkey:\tvalue
    parsed_any = False
    for line in text.strip().splitlines():
        if line.startswith("#"):
            continue
        if not line.startswith("SN\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        key = parts[2].rstrip(":").strip()
        val_str = parts[3].strip()
        parsed_any = True
        try:
            result["raw_sn"][key] = float(val_str) if "." in val_str else int(val_str)
        except ValueError:
            result["raw_sn"][key] = val_str

    if parsed_any:
        result["parse_ok"] = True
        sn = result["raw_sn"]
        result["n_samples"] = sn.get("number of samples")
        result["n_records"] = sn.get("number of records")
        result["n_snps"] = sn.get("number of SNPs")
        result["n_indels"] = sn.get("number of indels")
        result["n_multiallelic"] = sn.get("number of multiallelic sites")
        result["ts_tv"] = sn.get("Ts/Tv")
    else:
        result["errors"].append("No SN: lines found in bcftools stats output")

    return result


def parse_bcftools_stats_file(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    sample = path.stem.replace("_bcftools_stats", "")
    try:
        return parse_bcftools_stats(
            path.read_text(encoding="utf-8"), sample=sample, source=str(path)
        )
    except Exception as e:
        return {"sample": sample, "source": str(path), "parse_ok": False, "errors": [str(e)]}


def parse_bcftools_dir(qc_dir: str | Path) -> list[dict[str, Any]]:
    """Parse all bcftools stats files from a directory."""
    qc_dir = Path(qc_dir)
    if not qc_dir.exists():
        return []
    return [
        parse_bcftools_stats_file(p)
        for p in sorted(qc_dir.glob("*_bcftools_stats.txt"))
    ]
