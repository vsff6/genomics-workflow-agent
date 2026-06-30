"""Parse samtools flagstat, idxstats, and stats outputs."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any


def parse_flagstat(text: str, sample: str = "", source: str = "") -> dict[str, Any]:
    """Parse samtools flagstat text output."""
    result: dict[str, Any] = {
        "sample": sample,
        "source": source,
        "parse_ok": False,
        "errors": [],
        "total_reads": None,
        "mapped_reads": None,
        "mapped_pct": None,
        "paired_reads": None,
        "properly_paired": None,
        "properly_paired_pct": None,
        "duplicates": None,
        "secondary": None,
        "supplementary": None,
        "raw": {},
    }

    if not text or not text.strip():
        result["errors"].append("Empty flagstat output")
        return result

    parsed_any = False
    for line in text.strip().splitlines():
        m = re.match(r"^(\d+)\s+\+\s+(\d+)\s+(.+)", line)
        if not m:
            continue
        qc_pass = int(m.group(1))
        description = m.group(3).strip()
        parsed_any = True

        if description.startswith("in total"):
            result["total_reads"] = qc_pass
            result["raw"]["total"] = line
        elif description.startswith("mapped ("):
            result["mapped_reads"] = qc_pass
            pct_m = re.search(r"([\d.]+)%", line)
            if pct_m:
                result["mapped_pct"] = float(pct_m.group(1))
            result["raw"]["mapped"] = line
        elif description.startswith("paired in sequencing"):
            result["paired_reads"] = qc_pass
            result["raw"]["paired"] = line
        elif description.startswith("properly paired"):
            result["properly_paired"] = qc_pass
            pct_m = re.search(r"([\d.]+)%", line)
            if pct_m:
                result["properly_paired_pct"] = float(pct_m.group(1))
            result["raw"]["properly_paired"] = line
        elif description.startswith("duplicates"):
            result["duplicates"] = qc_pass
            result["raw"]["duplicates"] = line
        elif description.startswith("secondary"):
            result["secondary"] = qc_pass
            result["raw"]["secondary"] = line
        elif description.startswith("supplementary"):
            result["supplementary"] = qc_pass
            result["raw"]["supplementary"] = line

    if parsed_any:
        result["parse_ok"] = True
    else:
        result["errors"].append("No recognizable flagstat lines found")

    return result


def parse_idxstats(text: str, sample: str = "", source: str = "") -> dict[str, Any]:
    """Parse samtools idxstats text output."""
    result: dict[str, Any] = {
        "sample": sample,
        "source": source,
        "parse_ok": False,
        "errors": [],
        "contigs": [],
        "total_mapped": 0,
        "total_unmapped": 0,
        "zero_read_contigs": [],
    }

    if not text or not text.strip():
        result["errors"].append("Empty idxstats output")
        return result

    parsed_any = False
    for line in text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            contig = parts[0]
            length = int(parts[1])
            mapped = int(parts[2])
            unmapped = int(parts[3])
        except (ValueError, IndexError):
            continue

        parsed_any = True
        result["contigs"].append({
            "name": contig,
            "length": length,
            "mapped": mapped,
            "unmapped": unmapped,
        })
        result["total_mapped"] += mapped
        result["total_unmapped"] += unmapped

        if mapped == 0 and contig != "*":
            result["zero_read_contigs"].append(contig)

    if parsed_any:
        result["parse_ok"] = True
    else:
        result["errors"].append("No recognizable idxstats lines found")

    return result


def parse_stats(text: str, sample: str = "", source: str = "") -> dict[str, Any]:
    """Parse samtools stats output (SN: summary lines)."""
    result: dict[str, Any] = {
        "sample": sample,
        "source": source,
        "parse_ok": False,
        "errors": [],
        "summary": {},
        "raw_sn": {},
    }

    if not text or not text.strip():
        result["errors"].append("Empty stats output")
        return result

    parsed_any = False
    for line in text.strip().splitlines():
        if not line.startswith("SN\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        key = parts[1].rstrip(":").strip()
        val_str = parts[2].strip()
        parsed_any = True
        try:
            result["raw_sn"][key] = float(val_str) if "." in val_str else int(val_str)
        except ValueError:
            result["raw_sn"][key] = val_str

    if parsed_any:
        result["parse_ok"] = True
        sn = result["raw_sn"]
        result["summary"] = {
            "total_reads": sn.get("raw total sequences"),
            "mapped_reads": sn.get("reads mapped"),
            "error_rate": sn.get("error rate"),
            "average_length": sn.get("average length"),
            "average_quality": sn.get("average quality"),
            "insert_size_average": sn.get("insert size average"),
        }
    else:
        result["errors"].append("No SN: lines found in samtools stats output")

    return result


def parse_flagstat_file(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    sample = path.stem.replace("_flagstat", "")
    try:
        return parse_flagstat(path.read_text(encoding="utf-8"), sample=sample, source=str(path))
    except Exception as e:
        return {"sample": sample, "source": str(path), "parse_ok": False, "errors": [str(e)]}


def parse_idxstats_file(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    sample = path.stem.replace("_idxstats", "")
    try:
        return parse_idxstats(path.read_text(encoding="utf-8"), sample=sample, source=str(path))
    except Exception as e:
        return {"sample": sample, "source": str(path), "parse_ok": False, "errors": [str(e)]}


def parse_stats_file(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    sample = path.stem.replace("_stats", "")
    try:
        return parse_stats(path.read_text(encoding="utf-8"), sample=sample, source=str(path))
    except Exception as e:
        return {"sample": sample, "source": str(path), "parse_ok": False, "errors": [str(e)]}


def parse_variant_qc_dir(qc_dir: str | Path) -> dict[str, Any]:
    """Parse all samtools outputs from a variant_qc/ directory."""
    qc_dir = Path(qc_dir)
    result: dict[str, Any] = {
        "flagstat": [],
        "idxstats": [],
        "stats": [],
        "parse_ok": qc_dir.exists(),
    }
    if not qc_dir.exists():
        return result
    for p in sorted(qc_dir.glob("*_flagstat.txt")):
        result["flagstat"].append(parse_flagstat_file(p))
    for p in sorted(qc_dir.glob("*_idxstats.txt")):
        result["idxstats"].append(parse_idxstats_file(p))
    for p in sorted(qc_dir.glob("*_stats.txt")):
        result["stats"].append(parse_stats_file(p))
    return result
