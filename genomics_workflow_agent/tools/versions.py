"""Checks which external bioinformatics tools are available on the system."""

from __future__ import annotations

import shutil
import subprocess

EXTERNAL_TOOLS = [
    "fastqc",
    "multiqc",
    "fastp",
    "cutadapt",
    "nextflow",
    "samtools",
    "bcftools",
    "bedtools",
    "qiime",
    "Rscript",
    "docker",
    "singularity",
    "conda",
]


def get_version(cmd: str, version_flag: str = "--version") -> str | None:
    """Return first line of version output, or None if tool not found."""
    path = shutil.which(cmd)
    if path is None:
        return None
    try:
        result = subprocess.run(
            [cmd, version_flag], capture_output=True, text=True, timeout=5
        )
        output = (result.stdout + result.stderr).strip()
        return output.split("\n")[0][:120] if output else "found (version unknown)"
    except Exception:
        return "found (version check failed)"


def check_tools(tool_list: list[str] | None = None) -> dict:
    """Return machine-readable tool availability dict."""
    tools = tool_list or EXTERNAL_TOOLS
    result: dict[str, dict] = {}
    for tool in tools:
        path = shutil.which(tool)
        version = get_version(tool) if path else None
        result[tool] = {
            "available": path is not None,
            "path": path,
            "version": version,
        }
    return result


def available_tools(tool_list: list[str] | None = None) -> set[str]:
    """Return set of available tool names."""
    return {t for t, info in check_tools(tool_list).items() if info["available"]}
