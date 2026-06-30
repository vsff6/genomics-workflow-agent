"""
check_environment.py - Pre-flight environment check for the genomics-agent workspace.

Checks:
  - Python version
  - Required Python packages
  - Optional Python packages
  - External CLI tools (samtools, bcftools, bedtools, nextflow, etc.)
  - Output directory writability
  - Official Claude Life Sciences skill availability hint

Writes JSON and Markdown outputs. Missing optional tools are warnings, not failures.
Required packages missing cause a non-zero exit code.
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

VERSION = "1.0.0"
TOOL_NAME = "check_environment.py"

# Packages required for local tools to function at all
REQUIRED_PACKAGES = [
    ("numpy", ">=1.24", "MAD-based QC thresholds"),
    ("pandas", ">=2.0", "Data loading and QC metrics CSV"),
    ("scipy", ">=1.11", "Statistical utilities"),
    ("matplotlib", ">=3.7", "QC distribution plots"),
    ("seaborn", ">=0.12", "QC violin and scatter plots"),
    ("h5py", ">=3.9", "h5ad and 10x HDF5 file reading"),
    ("yaml", None, "YAML configuration (pyyaml)"),
    ("jsonschema", None, "JSON validation"),
]

# Optional packages - missing = warning only
OPTIONAL_PACKAGES = [
    ("scanpy", ">=1.9", "scRNA-seq QC (local fallback); prefer single-cell-rna-qc@life-sciences"),
    ("anndata", ">=0.10", "AnnData object support for h5ad files"),
    ("scvi", None, "scVI/scANVI modeling (scvi-tools); prefer scvi-tools@life-sciences"),
    ("pysam", ">=0.22", "BAM/CRAM alignment QC; Linux/macOS only (no Windows wheel)"),
    ("biopython", None, "FASTA/FASTQ parsing utilities"),
    ("pyarrow", ">=14.0", "Parquet and Arrow format support"),
    ("sklearn", None, "scikit-learn for downstream ML utilities"),
]

# External CLI tools - all optional
EXTERNAL_TOOLS = [
    ("samtools", "BAM/CRAM manipulation and flagstat"),
    ("bcftools", "VCF stats and filtering"),
    ("bedtools", "BED interval operations (FRiP, blacklist fraction)"),
    ("tabix", "VCF/BED indexing"),
    ("fastqc", "FASTQ read quality control"),
    ("multiqc", "Aggregate QC report generation"),
    ("fastp", "FASTQ trimming and QC"),
    ("macs3", "ATAC/ChIP-seq peak calling"),
    ("deeptools", "TSS enrichment, coverage, heatmaps"),
    ("nextflow", "nf-core pipeline orchestration"),
    ("mosdepth", "WGS/WES coverage statistics"),
    ("docker", "Container runtime for nf-core profiles"),
    ("singularity", "HPC container runtime for nf-core"),
    ("conda", "Conda package manager for environment management"),
]


def setup_logging(output_dir: Path, verbose: bool) -> logging.Logger:
    log_path = output_dir / "check_environment.log"
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path, encoding="utf-8"), logging.StreamHandler(sys.stderr)],
    )
    return logging.getLogger(TOOL_NAME)


def check_python_version() -> dict:
    major, minor = sys.version_info.major, sys.version_info.minor
    ok = major == 3 and minor >= 11
    return {
        "name": "Python version",
        "version": f"{major}.{minor}.{sys.version_info.micro}",
        "executable": sys.executable,
        "status": "ok" if ok else "warning",
        "note": None if ok else f"Python 3.11+ recommended; found {major}.{minor}",
    }


def check_python_package(import_name: str, min_version_hint: str, purpose: str, required: bool) -> dict:
    result = {
        "name": import_name,
        "purpose": purpose,
        "min_version": min_version_hint,
        "required": required,
    }
    try:
        mod = __import__(import_name)
        ver = None
        for attr in ("__version__", "version", "VERSION"):
            if hasattr(mod, attr):
                ver = str(getattr(mod, attr))
                break
        if ver is None:
            try:
                from importlib.metadata import version as meta_ver
                ver = meta_ver(import_name)
            except Exception:
                ver = "installed (version unknown)"
        result["status"] = "ok"
        result["version"] = ver
    except ImportError as e:
        result["status"] = "MISSING" if required else "warning"
        result["version"] = None
        result["note"] = str(e)
    return result


def check_external_tool(name: str, purpose: str) -> dict:
    path = shutil.which(name)
    result = {"name": name, "purpose": purpose, "required": False}
    if path:
        result["status"] = "ok"
        result["path"] = path
        try:
            r = subprocess.run([name, "--version"], capture_output=True, text=True, timeout=5)
            first_line = (r.stdout or r.stderr or "").strip().splitlines()
            result["version"] = first_line[0][:120] if first_line else "version unknown"
        except Exception:
            result["version"] = "version unknown"
    else:
        result["status"] = "warning"
        result["path"] = None
        result["version"] = None
        result["note"] = f"{name} not found in PATH. Optional - see README for install instructions."
    return result


def check_output_dir(path: Path) -> dict:
    result = {"name": "Output directory writability", "path": str(path)}
    try:
        path.mkdir(parents=True, exist_ok=True)
        test_file = path / ".write_test"
        test_file.write_text("ok")
        test_file.unlink()
        result["status"] = "ok"
        result["note"] = "Directory is writable."
    except Exception as e:
        result["status"] = "FAIL"
        result["note"] = f"Cannot write to output directory: {e}"
    return result


def check_official_skills() -> dict:
    """
    Official Claude Life Sciences skills cannot be checked programmatically from a script.
    This check always returns a 'manual' status with install instructions.
    """
    return {
        "name": "Official Claude Life Sciences skills",
        "status": "manual",
        "note": (
            "Cannot verify from a script. In an interactive Claude Code terminal, run:\n"
            "  /plugin list\n"
            "If missing, install with:\n"
            "  /plugin marketplace add anthropics/life-sciences\n"
            "  /plugin install single-cell-rna-qc@life-sciences\n"
            "  /plugin install scvi-tools@life-sciences\n"
            "  /plugin install nextflow-development@life-sciences\n"
            "  /plugin install 10x-genomics@life-sciences\n"
            "  /plugin install pubmed@life-sciences"
        ),
        "skills": [
            "single-cell-rna-qc@life-sciences",
            "scvi-tools@life-sciences",
            "nextflow-development@life-sciences",
            "10x-genomics@life-sciences",
            "pubmed@life-sciences",
        ],
    }


def build_markdown(python_check, required_checks, optional_checks, tool_checks, dir_check, skill_check) -> str:
    lines = [
        "# Environment Check Report",
        f"\nGenerated: {datetime.now().isoformat()}",
        f"\n**Tool**: `{TOOL_NAME}` v{VERSION}",
        f"\n**Python**: {python_check['version']} - {python_check['executable']}",
        "",
    ]

    def status_icon(s):
        return {"ok": "[OK]", "warning": "[WARN]", "MISSING": "[MISSING]", "FAIL": "[FAIL]", "manual": "[MANUAL]"}.get(s, s)

    lines.append("\n## Python Version\n")
    icon = status_icon(python_check["status"])
    lines.append(f"- {icon} Python {python_check['version']} ({python_check['executable']})")
    if python_check.get("note"):
        lines.append(f"  - Note: {python_check['note']}")

    lines.append("\n## Required Python Packages\n")
    for c in required_checks:
        icon = status_icon(c["status"])
        ver = c.get("version", "not installed")
        lines.append(f"- {icon} **{c['name']}** {ver or 'NOT INSTALLED'} - {c['purpose']}")
        if c.get("note"):
            lines.append(f"  - Note: {c['note']}")

    lines.append("\n## Optional Python Packages\n")
    for c in optional_checks:
        icon = status_icon(c["status"])
        ver = c.get("version", "not installed")
        lines.append(f"- {icon} {c['name']} {ver or 'not installed'} - {c['purpose']}")
        if c.get("note"):
            lines.append(f"  - Note: {c['note']}")

    lines.append("\n## External CLI Tools\n")
    lines.append("All external tools are optional. Missing tools are warnings, not failures.\n")
    for c in tool_checks:
        icon = status_icon(c["status"])
        ver = c.get("version") or "not found"
        lines.append(f"- {icon} **{c['name']}** - {c['purpose']}")
        if c.get("path"):
            lines.append(f"  - Path: `{c['path']}`")
            lines.append(f"  - Version: {ver}")
        else:
            lines.append(f"  - {c.get('note', 'Not found in PATH')}")

    lines.append("\n## Output Directory\n")
    icon = status_icon(dir_check["status"])
    lines.append(f"- {icon} `{dir_check['path']}` - {dir_check['note']}")

    lines.append("\n## Official Claude Life Sciences Skills\n")
    lines.append(f"- [{status_icon(skill_check['status'])}] {skill_check['note']}")
    lines.append("\n**Skills to install:**")
    for s in skill_check.get("skills", []):
        lines.append(f"  - `{s}`")

    lines.append("\n## Summary\n")
    all_checks = required_checks + optional_checks + tool_checks + [dir_check]
    n_ok = sum(1 for c in all_checks if c["status"] == "ok")
    n_warn = sum(1 for c in all_checks if c["status"] == "warning")
    n_fail = sum(1 for c in all_checks if c["status"] in ("MISSING", "FAIL"))
    lines.append(f"- OK: {n_ok}")
    lines.append(f"- Warnings (optional, non-blocking): {n_warn}")
    lines.append(f"- Missing/Failures (blocking): {n_fail}")
    if n_fail == 0:
        lines.append("\n**Environment is ready for local tool execution.**")
    else:
        lines.append(f"\n**{n_fail} required package(s) missing. Install before running local tools.**")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Pre-flight environment check for the genomics-agent workspace.",
    )
    parser.add_argument("--output-dir", default="reports/env_check", help="Directory for JSON and Markdown outputs")
    parser.add_argument("--json", action="store_true", help="Write JSON output (default: always written)")
    parser.add_argument("--markdown", action="store_true", help="Write Markdown output (default: always written)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logging(out_dir, args.verbose)
    log.info(f"{TOOL_NAME} v{VERSION}")

    python_check = check_python_version()
    log.info(f"Python {python_check['version']} - {python_check['status']}")

    required_checks = []
    for pkg, ver, purpose in REQUIRED_PACKAGES:
        r = check_python_package(pkg, ver, purpose, required=True)
        log.info(f"  [required] {pkg}: {r['status']} {r.get('version', '')}")
        required_checks.append(r)

    optional_checks = []
    for pkg, ver, purpose in OPTIONAL_PACKAGES:
        r = check_python_package(pkg, ver, purpose, required=False)
        log.info(f"  [optional] {pkg}: {r['status']} {r.get('version', '')}")
        optional_checks.append(r)

    tool_checks = []
    for tool, purpose in EXTERNAL_TOOLS:
        r = check_external_tool(tool, purpose)
        log.info(f"  [cli] {tool}: {r['status']} {r.get('version', '')}")
        tool_checks.append(r)

    dir_check = check_output_dir(out_dir)
    log.info(f"Output dir: {dir_check['status']}")

    skill_check = check_official_skills()

    n_fail = sum(1 for c in required_checks if c["status"] in ("MISSING", "FAIL"))
    n_warn = sum(1 for c in optional_checks + tool_checks if c["status"] == "warning")
    n_ok = sum(1 for c in required_checks + optional_checks + tool_checks if c["status"] == "ok")

    md = build_markdown(python_check, required_checks, optional_checks, tool_checks, dir_check, skill_check)
    md_path = out_dir / "environment_check.md"
    md_path.write_text(md, encoding="utf-8")
    log.info(f"Markdown: {md_path}")

    summary = {
        "tool": TOOL_NAME,
        "version": VERSION,
        "generated": datetime.now().isoformat(),
        "python": python_check,
        "required_packages": required_checks,
        "optional_packages": optional_checks,
        "external_tools": tool_checks,
        "output_directory": dir_check,
        "official_skills": skill_check,
        "summary": {
            "n_ok": n_ok,
            "n_warnings": n_warn,
            "n_failures": n_fail,
            "ready": n_fail == 0 and dir_check["status"] == "ok",
        },
    }
    json_path = out_dir / "environment_check.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    log.info(f"JSON: {json_path}")

    print(f"\nEnvironment check complete. Outputs: {out_dir}")
    print(f"  Python: {python_check['version']}")
    print(f"  Required packages OK: {sum(1 for c in required_checks if c['status'] == 'ok')}/{len(required_checks)}")
    print(f"  Optional packages OK: {sum(1 for c in optional_checks if c['status'] == 'ok')}/{len(optional_checks)}")
    print(f"  CLI tools found: {sum(1 for c in tool_checks if c['status'] == 'ok')}/{len(tool_checks)}")
    print(f"  Warnings (non-blocking): {n_warn}")
    print(f"  Failures (blocking): {n_fail}")
    if n_fail == 0:
        print("  Status: READY")
    else:
        print(f"  Status: {n_fail} REQUIRED PACKAGES MISSING - install before running tools")

    return 1 if n_fail > 0 or dir_check["status"] == "FAIL" else 0


if __name__ == "__main__":
    sys.exit(main())
