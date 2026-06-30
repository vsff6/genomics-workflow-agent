"""
nfcore_launcher.py - nf-core workflow planner and safe launcher.

Builds samplesheets, validates preflight requirements, constructs Nextflow
commands, and optionally executes them.

Default is dry-run: writes a plan and samplesheet but does NOT launch Nextflow
unless --run is explicitly provided by the user.

Supported workflows: nf-core/rnaseq, nf-core/sarek, nf-core/atacseq
"""

import argparse
import csv
import json
import logging
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

VERSION = "1.0.0"
TOOL_NAME = "nfcore_launcher.py"

SUPPORTED_WORKFLOWS = ["rnaseq", "sarek", "atacseq"]
FASTQ_EXTENSIONS = {".fastq.gz", ".fastq", ".fq.gz", ".fq"}


def setup_logging(output_dir: Path, verbose: bool) -> logging.Logger:
    log_path = output_dir / "nfcore_launcher.log"
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stderr)],
    )
    return logging.getLogger(TOOL_NAME)


def get_tool_version(cmd: str) -> str:
    try:
        result = subprocess.run([cmd, "--version"], capture_output=True, text=True, timeout=5)
        output = (result.stdout + result.stderr).strip()
        return output.split("\n")[0][:120] if output else "unknown"
    except Exception:
        return "unknown"


def detect_executors() -> dict:
    """Check availability of Nextflow and container/environment executors."""
    executors = {}
    for tool in ["nextflow", "docker", "singularity", "apptainer", "conda"]:
        path = shutil.which(tool)
        executors[tool] = {
            "available": path is not None,
            "path": path,
            "version": get_tool_version(tool) if path is not None else None,
        }
    return executors


# ──────────────────────────────────────────────────────────
# FASTQ discovery helpers
# ──────────────────────────────────────────────────────────

def find_fastq_files(input_dir: Path) -> list:
    """Return sorted list of FASTQ files in input_dir (non-recursive)."""
    if not input_dir.exists():
        return []
    return sorted(
        p for p in input_dir.iterdir()
        if p.is_file() and any(str(p.name).endswith(ext) for ext in FASTQ_EXTENSIONS)
    )


def _sample_name(path: Path) -> str:
    """Derive sample name by stripping FASTQ extension and paired-end suffix."""
    name = path.name
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return re.sub(r"[_.]R?[12]$", "", name)


def _is_r1(path: Path) -> bool:
    return bool(re.search(r"[_.]R?1[_.]|[_.]R?1$", path.name))


def _is_r2(path: Path) -> bool:
    return bool(re.search(r"[_.]R?2[_.]|[_.]R?2$", path.name))


def _pair_fastqs(fastq_files: list) -> tuple[list, dict]:
    """Return (r1_files, r2_by_sample_name) from a list of FASTQ paths."""
    r1_files = [f for f in fastq_files if _is_r1(f)]
    r2_by_name = {_sample_name(f): f for f in fastq_files if _is_r2(f)}
    return r1_files, r2_by_name


# ──────────────────────────────────────────────────────────
# Samplesheet builders (pure - no subprocess)
# ──────────────────────────────────────────────────────────

def build_rnaseq_samplesheet(input_dir: Path, output_path: Path, log) -> dict:
    """
    Build an nf-core/rnaseq samplesheet from FASTQ files in input_dir.

    Detects paired-end (*_R1*/*_R2* or *_1*/*_2*) and single-end FASTQ files.
    Strandedness defaults to 'auto' - must be confirmed from library prep records.
    """
    fastq_files = find_fastq_files(input_dir)
    if not fastq_files:
        return {"created": False, "reason": f"No FASTQ files found in {input_dir}", "rows": 0}

    r1_files, r2_by_name = _pair_fastqs(fastq_files)
    rows = []
    paired_samples = set()

    for r1 in r1_files:
        sample = _sample_name(r1)
        r2 = r2_by_name.get(sample)
        rows.append({
            "sample": sample,
            "fastq_1": str(r1.resolve()),
            "fastq_2": str(r2.resolve()) if r2 else "",
            "strandedness": "auto",
        })
        paired_samples.add(sample)

    # Single-end: files with no R1/R2 suffix that were not already paired
    for f in fastq_files:
        if not _is_r1(f) and not _is_r2(f):
            sample = _sample_name(f)
            if sample not in paired_samples:
                rows.append({
                    "sample": sample,
                    "fastq_1": str(f.resolve()),
                    "fastq_2": "",
                    "strandedness": "auto",
                })

    if not rows:
        return {"created": False, "reason": "No R1 FASTQ files detected", "rows": 0}

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        writer.writerows(rows)

    log.info(f"rnaseq samplesheet: {output_path} ({len(rows)} rows)")
    return {
        "created": True,
        "path": str(output_path),
        "rows": len(rows),
        "warnings": [
            "strandedness set to 'auto' - confirm from library preparation records before running",
            "verify sample names, pairing, and file paths before execution",
        ],
    }


def build_sarek_samplesheet(input_dir: Path, output_path: Path, log) -> dict:
    """
    Build a conservative nf-core/sarek samplesheet draft from FASTQ files.

    patient, sex, status (tumor/normal), pedigree, ancestry, and sequencing type
    CANNOT be inferred from filenames. All fields use safe placeholders that MUST
    be manually reviewed and corrected before any pipeline is launched.
    """
    fastq_files = find_fastq_files(input_dir)
    if not fastq_files:
        return {"created": False, "reason": f"No FASTQ files found in {input_dir}", "rows": 0}

    r1_files, r2_by_name = _pair_fastqs(fastq_files)
    rows = []

    for r1 in r1_files:
        sample = _sample_name(r1)
        r2 = r2_by_name.get(sample)
        rows.append({
            "patient": "PATIENT_ID",
            "sex": "unknown",
            "status": "0",
            "sample": sample,
            "lane": "L001",
            "fastq_1": str(r1.resolve()),
            "fastq_2": str(r2.resolve()) if r2 else "",
        })

    if not rows:
        return {"created": False, "reason": "No R1 FASTQ files detected", "rows": 0}

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2"]
        )
        writer.writeheader()
        writer.writerows(rows)

    log.info(f"sarek samplesheet draft: {output_path} ({len(rows)} rows)")
    return {
        "created": True,
        "path": str(output_path),
        "rows": len(rows),
        "warnings": [
            "PATIENT_ID placeholders must be replaced with real patient identifiers",
            "tumor/normal status (status=0/1) cannot be inferred from filenames - set manually",
            "sex (XX/XY/unknown) must be confirmed from sample metadata",
            "pedigree, germline/somatic context, and disease status cannot be inferred automatically",
            "target intervals for WES must be added manually via --wes_intervals or pipeline config",
            "do not execute without manual review by the analyst responsible for this study",
        ],
    }


def build_atacseq_samplesheet(input_dir: Path, output_path: Path, log) -> dict:
    """
    Build a conservative nf-core/atacseq samplesheet draft.

    Replicate structure, control samples, and experimental design CANNOT be
    inferred from filenames. All replicate fields use placeholder values.
    """
    fastq_files = find_fastq_files(input_dir)
    if not fastq_files:
        return {"created": False, "reason": f"No FASTQ files found in {input_dir}", "rows": 0}

    r1_files, r2_by_name = _pair_fastqs(fastq_files)
    rows = []

    for r1 in r1_files:
        sample = _sample_name(r1)
        r2 = r2_by_name.get(sample)
        rows.append({
            "sample": sample,
            "fastq_1": str(r1.resolve()),
            "fastq_2": str(r2.resolve()) if r2 else "",
            "replicate": "1",
        })

    if not rows:
        return {"created": False, "reason": "No R1 FASTQ files detected", "rows": 0}

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "replicate"])
        writer.writeheader()
        writer.writerows(rows)

    log.info(f"atacseq samplesheet draft: {output_path} ({len(rows)} rows)")
    return {
        "created": True,
        "path": str(output_path),
        "rows": len(rows),
        "warnings": [
            "replicate numbers are placeholders - must reflect the actual experimental design",
            "control/input samples must be identified and annotated separately",
            "blacklist BED must be provided via --blacklist or pipeline config",
            "genome build and matching GTF/FASTA must be confirmed before execution",
            "peak-calling strategy (MACS3 settings) is tissue- and protocol-dependent",
        ],
    }


# ──────────────────────────────────────────────────────────
# Command builder (pure - no subprocess)
# ──────────────────────────────────────────────────────────

def build_nfcore_command(workflow: str, args, samplesheet_path: str) -> list:
    """Construct the nextflow run command as a list of strings."""
    cmd = ["nextflow", "run", f"nf-core/{workflow}"]

    input_path = samplesheet_path or args.samplesheet
    if input_path:
        cmd += ["--input", input_path]

    pipeline_outdir = str(Path(args.output_dir) / "pipeline_output")
    cmd += ["--outdir", pipeline_outdir]

    if args.genome:
        cmd += ["--genome", args.genome]
    if args.fasta:
        cmd += ["--fasta", args.fasta]
    if args.gtf:
        cmd += ["--gtf", args.gtf]
    if args.bed:
        cmd += ["--bed", args.bed]

    cmd += ["-profile", args.profile or "docker"]

    if args.max_cpus:
        cmd += ["--max_cpus", str(args.max_cpus)]
    if args.max_memory:
        cmd += ["--max_memory", args.max_memory]
    if args.max_time:
        cmd += ["--max_time", args.max_time]
    if args.resume:
        cmd += ["-resume"]
    if args.extra_args:
        cmd += shlex.split(args.extra_args)

    return cmd


# ──────────────────────────────────────────────────────────
# Biological caveats (required in every plan)
# ──────────────────────────────────────────────────────────

BIOLOGICAL_CAVEATS = {
    "rnaseq": [
        "Strandedness must match library preparation protocol - confirm before running.",
        "Genome build and GTF annotation version must match each other and all downstream tools.",
        "Batch, tissue, condition, and replicate metadata are required for meaningful differential "
        "expression analysis - quantification output alone is not interpretable.",
        "Successful quantification is not biological validation.",
        "Salmon/STAR count estimates require normalization and experimental context to interpret.",
    ],
    "sarek": [
        "No clinical interpretation is provided. Variant calls must not be used for medical decisions.",
        "Tumor/normal pairing, pedigree structure, germline/somatic context, and target intervals "
        "(WES) must be confirmed by the analyst before running.",
        "Reference build and known-sites VCF compatibility (BQSR) must be verified.",
        "Sequencing type (WGS vs WES), tumor purity, ploidy, and contamination are not automatically determined.",
        "Successful variant calling is not clinical validation - pathogenicity requires expert review "
        "with validated annotation sources.",
        "Ancestry and population context affect variant interpretation and must not be assumed.",
    ],
    "atacseq": [
        "Genome build, blacklist BED, and GTF annotation must be consistent across all pipeline inputs.",
        "Replicate structure and control/input sample assignment must reflect the experimental design.",
        "Peak-calling strategy (MACS3 settings, effective genome size, q-value) is tissue- and "
        "protocol-dependent - defaults may not be appropriate.",
        "FRiP and TSS enrichment thresholds are not universal - interpret in context of tissue, "
        "protocol, and cell type.",
        "Successful peak calling is not biological validation.",
    ],
}


# ──────────────────────────────────────────────────────────
# MultiQC output parser (no external tools required)
# ──────────────────────────────────────────────────────────

def parse_multiqc_output(output_dir: Path, log) -> dict:
    """
    Summarize MultiQC outputs in a pipeline output directory.
    Does not require MultiQC to be installed - reads existing files only.
    """
    result = {
        "searched_path": str(output_dir),
        "multiqc_report_html": None,
        "multiqc_data_dir": None,
        "general_stats_present": False,
        "multiqc_data_json_present": False,
        "summary": {},
        "skipped": [],
    }

    if not output_dir.exists():
        result["skipped"].append(f"Output directory not found: {output_dir}")
        return result

    html_candidates = list(output_dir.rglob("multiqc_report.html"))
    if html_candidates:
        result["multiqc_report_html"] = str(html_candidates[0])
        log.info(f"MultiQC HTML report: {html_candidates[0]}")
    else:
        result["skipped"].append(
            "multiqc_report.html not found - pipeline may not have completed or MultiQC was not run"
        )

    data_dirs = list(output_dir.rglob("multiqc_data"))
    if data_dirs:
        data_dir = data_dirs[0]
        result["multiqc_data_dir"] = str(data_dir)

        stats_file = data_dir / "multiqc_general_stats.txt"
        if stats_file.exists():
            result["general_stats_present"] = True
            try:
                lines = stats_file.read_text(errors="replace").splitlines()
                result["summary"]["general_stats_rows"] = max(0, len(lines) - 1)
                result["summary"]["general_stats_header"] = lines[0][:200] if lines else ""
            except Exception as e:
                result["skipped"].append(f"Could not parse multiqc_general_stats.txt: {e}")

        json_file = data_dir / "multiqc_data.json"
        if json_file.exists():
            result["multiqc_data_json_present"] = True
            try:
                data = json.loads(json_file.read_text(errors="replace"))
                result["summary"]["plot_types"] = list(data.get("report_plot_types", {}).keys())[:10]
                result["summary"]["saved_raw_data"] = list(data.get("report_saved_raw_data", {}).keys())[:10]
            except Exception as e:
                result["skipped"].append(f"Could not parse multiqc_data.json: {e}")

    return result


# ──────────────────────────────────────────────────────────
# Markdown report builder
# ──────────────────────────────────────────────────────────

def build_markdown(
    workflow: str, plan: dict, caveats: list, executors: dict,
    samplesheet_result: dict, command: list, dry_run: bool,
) -> str:
    mode_label = "DRY RUN (planning only - no workflow was launched)" if dry_run else "RUN"
    lines = [
        f"# nf-core/{workflow} Plan",
        f"\nGenerated: {datetime.now().isoformat()}",
        f"\n**Tool**: `{TOOL_NAME}` v{VERSION}",
        f"\n**Mode**: {mode_label}",
        "\n---\n",
        "## Command\n",
        "```bash",
        " ".join(command) if command else "(no command generated)",
        "```\n",
        "## Preflight Check\n",
        "### Nextflow and Executors\n",
    ]

    for tool, info in executors.items():
        icon = "[OK]" if info["available"] else "[MISSING]"
        ver = info.get("version") or "not found"
        lines.append(f"- {icon} **{tool}**: {ver}")

    blockers = plan.get("blockers", [])
    warnings = plan.get("warnings", [])

    if blockers:
        lines.append("\n### Blockers (must resolve before --run)\n")
        for b in blockers:
            lines.append(f"- **[BLOCKER]** {b}")

    if warnings:
        lines.append("\n### Warnings\n")
        for w in warnings:
            lines.append(f"- [WARN] {w}")

    missing = plan.get("missing_requirements", [])
    if missing:
        lines.append("\n### Missing Requirements\n")
        for m in missing:
            lines.append(f"- {m}")

    lines.append("\n## Samplesheet\n")
    if samplesheet_result.get("created"):
        lines.append(f"- **Path**: `{samplesheet_result['path']}`")
        lines.append(f"- **Rows**: {samplesheet_result['rows']}")
        ss_warnings = samplesheet_result.get("warnings", [])
        if ss_warnings:
            lines.append("\n**Samplesheet warnings - review before executing:**\n")
            for w in ss_warnings:
                lines.append(f"- [WARN] {w}")
    else:
        lines.append(f"- Not created: {samplesheet_result.get('reason', 'unknown')}")

    refs = plan.get("references", {})
    if refs:
        lines.append("\n## Reference Inputs\n")
        for k, v in refs.items():
            lines.append(f"- **{k}**: `{v}`")

    lines.append("\n## Biological and Experimental Caveats\n")
    lines.append(
        "> These caveats apply regardless of whether the pipeline completes successfully. "
        "Successful pipeline execution is not biological or clinical validation.\n"
    )
    for caveat in caveats:
        lines.append(f"- {caveat}")

    run_result = plan.get("run_result")
    if run_result:
        lines.append("\n## Execution Result\n")
        if not run_result.get("executed", False):
            lines.append(f"- **Not executed**: {run_result.get('reason', 'unknown')}")
            for b in run_result.get("blockers", []):
                lines.append(f"  - Blocker: {b}")
        else:
            lines.append(f"- **Return code**: {run_result.get('return_code')}")
            lines.append(f"- **Runtime**: {run_result.get('runtime_s')} s")
            if run_result.get("error"):
                lines.append(f"- **Error**: {run_result['error']}")
            if run_result.get("stdout_snippet"):
                lines.append(f"\n**stdout (tail)**:\n```\n{run_result['stdout_snippet']}\n```")
            if run_result.get("stderr_snippet"):
                lines.append(f"\n**stderr (tail)**:\n```\n{run_result['stderr_snippet']}\n```")

    lines.append("\n## Limitations\n")
    lines.append("- This tool generates plans and samplesheet drafts. It does not replace Nextflow or nf-core.")
    lines.append("- For production execution, use `nextflow-development@life-sciences` or Nextflow directly.")
    lines.append("- Samplesheets generated from filenames require human review before any workflow is launched.")
    lines.append("- Reference files must be provided locally. Large references are never downloaded automatically.")
    lines.append("- Always run the `biology-interpretation-reviewer` agent after parsing pipeline QC outputs.")

    return "\n".join(lines)


# ──────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "nf-core workflow planner and safe launcher. "
            "Builds samplesheets, validates requirements, and writes a Nextflow command. "
            "Default is dry-run - does not execute unless --run is explicitly provided."
        ),
    )
    parser.add_argument("--workflow", choices=SUPPORTED_WORKFLOWS, required=True,
                        help="nf-core workflow to plan: rnaseq, sarek, or atacseq")
    parser.add_argument("--input-dir", help="Directory containing input FASTQ files")
    parser.add_argument("--samplesheet", help="Pre-built samplesheet CSV (skips auto-generation)")
    parser.add_argument("--output-dir", default="reports/nfcore",
                        help="Directory for plan JSON, Markdown, and commands.sh")
    parser.add_argument("--genome", help="nf-core genome key (e.g. GRCh38, GRCm39)")
    parser.add_argument("--fasta", help="Local genome FASTA path")
    parser.add_argument("--gtf", help="GTF annotation path")
    parser.add_argument("--bed", help="Target intervals or blacklist BED path")
    parser.add_argument("--profile", default="docker",
                        help="Nextflow executor profile (docker, singularity, conda, etc.)")
    parser.add_argument("--dry-run", action="store_true", default=True,
                        help="Plan only - do not execute (default)")
    parser.add_argument("--run", action="store_true", default=False,
                        help="Execute the workflow after preflight checks pass")
    parser.add_argument("--max-cpus", type=int, help="Override max CPUs")
    parser.add_argument("--max-memory", help="Override max memory (e.g. '32.GB')")
    parser.add_argument("--max-time", help="Override max time (e.g. '48.h')")
    parser.add_argument("--resume", action="store_true", help="Pass -resume to Nextflow")
    parser.add_argument("--extra-args", help="Additional Nextflow arguments (quoted string)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    # --run overrides the --dry-run default
    if args.run:
        args.dry_run = False

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logging(out_dir, args.verbose)

    log.info(f"{TOOL_NAME} v{VERSION}")
    log.info(f"Workflow: nf-core/{args.workflow}")
    log.info(f"Mode: {'DRY RUN' if args.dry_run else 'RUN'}")

    executors = detect_executors()
    for tool, info in executors.items():
        status = "ok" if info["available"] else "missing"
        log.info(f"  [{status}] {tool}: {info.get('version') or 'not found'}")

    blockers = []
    warnings_list = []
    plan = {
        "blockers": blockers,
        "warnings": warnings_list,
        "references": {},
        "missing_requirements": [],
    }

    if not executors["nextflow"]["available"]:
        blockers.append(
            "nextflow not found in PATH. "
            "Install: https://www.nextflow.io/docs/latest/getstarted.html"
        )

    container_available = any(
        executors[t]["available"] for t in ["docker", "singularity", "apptainer"]
    )
    if not container_available and not executors["conda"]["available"]:
        warnings_list.append(
            "No container runtime (docker/singularity/apptainer) or conda found. "
            "nf-core pipelines require at least one executor."
        )
    elif not container_available:
        warnings_list.append(
            "No container runtime found. Will use -profile conda, but nf-core recommends "
            "docker or singularity for reproducibility."
        )

    if not args.genome and not args.fasta:
        plan["missing_requirements"].append(
            f"nf-core/{args.workflow} requires either --genome (nf-core key) or --fasta (local FASTA)"
        )
        warnings_list.append("No reference genome provided - command will be incomplete")

    if args.genome:
        plan["references"]["genome_key"] = args.genome
    if args.fasta:
        plan["references"]["fasta"] = args.fasta
        if not Path(args.fasta).exists():
            blockers.append(f"FASTA not found: {args.fasta}")
    if args.gtf:
        plan["references"]["gtf"] = args.gtf
        if not Path(args.gtf).exists():
            warnings_list.append(f"GTF not found: {args.gtf}")
    if args.bed:
        plan["references"]["bed"] = args.bed
        if not Path(args.bed).exists():
            warnings_list.append(f"BED not found: {args.bed}")

    samplesheet_path = None
    samplesheet_result = {}

    if args.samplesheet:
        ss_path = Path(args.samplesheet)
        if ss_path.exists():
            samplesheet_path = str(ss_path)
            samplesheet_result = {
                "created": True, "path": samplesheet_path,
                "rows": "unknown (pre-built)", "warnings": [],
            }
            log.info(f"Using provided samplesheet: {ss_path}")
        else:
            blockers.append(f"Samplesheet not found: {args.samplesheet}")
            samplesheet_result = {"created": False, "reason": f"File not found: {args.samplesheet}"}

    elif args.input_dir:
        input_dir = Path(args.input_dir)
        ss_output = out_dir / f"samplesheet_{args.workflow}.csv"
        builders = {
            "rnaseq": build_rnaseq_samplesheet,
            "sarek": build_sarek_samplesheet,
            "atacseq": build_atacseq_samplesheet,
        }
        samplesheet_result = builders[args.workflow](input_dir, ss_output, log)
        if samplesheet_result.get("created"):
            samplesheet_path = samplesheet_result["path"]
        else:
            warnings_list.append(
                f"Samplesheet auto-generation failed: {samplesheet_result.get('reason', 'unknown')}"
            )
    else:
        warnings_list.append(
            "No --input-dir or --samplesheet provided. Command will lack --input argument."
        )
        samplesheet_result = {
            "created": False,
            "reason": "No input directory or samplesheet provided",
        }

    command = build_nfcore_command(args.workflow, args, samplesheet_path)

    commands_path = out_dir / "commands.sh"
    with open(commands_path, "w") as fh:
        fh.write("#!/usr/bin/env bash\n")
        fh.write(f"# nf-core/{args.workflow} - generated by {TOOL_NAME} v{VERSION}\n")
        fh.write(f"# Generated: {datetime.now().isoformat()}\n")
        fh.write(f"# Mode: {'DRY RUN (not executed)' if args.dry_run else 'RUN'}\n")
        fh.write("#\n# Review all parameters, samplesheet, and references before running.\n\n")
        fh.write(" ".join(command) + "\n")
    log.info(f"commands.sh: {commands_path}")

    run_result = None
    if args.run:
        if blockers:
            log.error(f"Cannot run: {len(blockers)} blocker(s) present.")
            for b in blockers:
                log.error(f"  BLOCKER: {b}")
            run_result = {"executed": False, "reason": "Preflight blockers present", "blockers": blockers}
        else:
            log.info(f"Executing: {' '.join(command)}")
            t0 = time.time()
            try:
                result = subprocess.run(command, capture_output=True, text=True, timeout=86400)
                runtime = round(time.time() - t0, 1)
                run_result = {
                    "executed": True,
                    "return_code": result.returncode,
                    "runtime_s": runtime,
                    "stdout_snippet": result.stdout[-500:] if result.stdout else "",
                    "stderr_snippet": result.stderr[-500:] if result.stderr else "",
                }
                log.info(f"Workflow exited with code {result.returncode} in {runtime}s")
            except subprocess.TimeoutExpired:
                run_result = {"executed": True, "error": "Workflow timed out (24 h limit)"}
            except Exception as e:
                run_result = {"executed": True, "error": str(e)}

        plan["run_result"] = run_result

    caveats = BIOLOGICAL_CAVEATS.get(args.workflow, [])

    # pipeline_output may not exist on a dry run - parse_multiqc_output handles the missing dir
    pipeline_output = Path(args.output_dir) / "pipeline_output"
    multiqc_result = parse_multiqc_output(pipeline_output, log)

    md = build_markdown(args.workflow, plan, caveats, executors, samplesheet_result, command, args.dry_run)
    md_path = out_dir / "nfcore_plan.md"
    md_path.write_text(md, encoding="utf-8")
    log.info(f"Markdown: {md_path}")

    summary = {
        "tool": TOOL_NAME,
        "version": VERSION,
        "generated": datetime.now().isoformat(),
        "workflow": f"nf-core/{args.workflow}",
        "mode": "dry_run" if args.dry_run else "run",
        "executors": executors,
        "command": " ".join(command),
        "samplesheet": samplesheet_result,
        "references": plan.get("references", {}),
        "blockers": blockers,
        "warnings": warnings_list,
        "missing_requirements": plan.get("missing_requirements", []),
        "biological_caveats": caveats,
        "multiqc": multiqc_result,
        "run_result": run_result,
    }
    json_path = out_dir / "nfcore_plan.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    log.info(f"JSON: {json_path}")

    print(f"\nnf-core/{args.workflow} plan complete. Outputs in: {out_dir}")
    print(f"  Mode: {'DRY RUN' if args.dry_run else 'RUN'}")
    print(f"  Nextflow: {'found' if executors['nextflow']['available'] else 'NOT FOUND'}")
    if samplesheet_result.get("created"):
        print(f"  Samplesheet: {samplesheet_result['path']} ({samplesheet_result['rows']} rows)")
    else:
        print(f"  Samplesheet: not created - {samplesheet_result.get('reason', 'unknown')}")
    print(f"  Blockers: {len(blockers)}  Warnings: {len(warnings_list)}")
    print(f"  commands.sh: {commands_path}")
    if blockers:
        print("\n  BLOCKERS (resolve before --run):")
        for b in blockers:
            print(f"    - {b}")

    return 1 if (args.run and blockers) else 0


if __name__ == "__main__":
    sys.exit(main())
