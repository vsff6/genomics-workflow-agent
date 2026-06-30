"""RNA-seq workflow - nf-core/rnaseq planning and execution."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomics_workflow_agent.tools.nextflow import (
    build_nextflow_cmd,
    run_nextflow,
)
from genomics_workflow_agent.tools.samplesheets import build_rnaseq_samplesheet
from genomics_workflow_agent.tools.versions import check_tools

REQUIRED_TOOLS = ["nextflow"]
OPTIONAL_TOOLS = ["docker", "singularity", "conda", "fastqc", "multiqc"]

BIOLOGICAL_CAVEATS = [
    "Strandedness must match library preparation protocol - always confirm from records, not filenames.",
    "Genome build and GTF annotation must be consistent. Mismatched builds invalidate gene-level counts.",
    "Batch effects, tissue, condition, and replicate metadata are required for meaningful differential expression.",
    "Salmon/STAR quantification is not normalization. TPM and raw counts require appropriate statistical models.",
    "Successful alignment or quantification is not biological validation.",
    "Intronic reads in smart-seq or similar protocols may inflate counts - check protocol carefully.",
]


def plan(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    genome: str | None = None,
    fasta: str | None = None,
    gtf: str | None = None,
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
    dry_run: bool = True,
) -> dict[str, Any]:
    """Return an nf-core/rnaseq execution plan."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    pipeline_outdir = str(output_dir / "pipeline_output")

    tools = check_tools(REQUIRED_TOOLS + OPTIONAL_TOOLS)
    warnings: list[str] = []
    blockers: list[str] = []
    missing_requirements: list[str] = []

    ss_path = output_dir / "samplesheet_rnaseq.csv"
    samplesheet_result = build_rnaseq_samplesheet(input_dir, ss_path)

    if not tools["nextflow"]["available"]:
        blockers.append("nextflow not found in PATH. Install: https://www.nextflow.io")

    if not genome and not fasta:
        missing_requirements.append("--genome (nf-core key, e.g. GRCh38) or --fasta (local path) required")
        warnings.append("No reference genome specified - command will be incomplete")

    if fasta and not Path(fasta).exists():
        blockers.append(f"FASTA not found: {fasta}")
    if gtf and not Path(gtf).exists():
        warnings.append(f"GTF not found: {gtf}")

    container_ok = tools["docker"]["available"] or tools["singularity"]["available"]
    if not container_ok and not tools["conda"]["available"]:
        warnings.append("No container runtime (docker/singularity) or conda found")

    cmd = build_nextflow_cmd(
        "rnaseq",
        input_path=str(ss_path) if samplesheet_result.get("created") else None,
        outdir=pipeline_outdir,
        genome=genome,
        fasta=fasta,
        gtf=gtf,
        profile=profile,
        resume=resume,
        extra_args=extra_args,
    )

    steps = [{
        "name": "nextflow_rnaseq",
        "description": "Run nf-core/rnaseq via Nextflow",
        "command": cmd,
        "output_dir": pipeline_outdir,
        "expected_outputs": [pipeline_outdir],
        "dry_run": dry_run,
        "required_tools": ["nextflow"],
    }]

    return {
        "workflow": "rnaseq",
        "pipeline": "nf-core/rnaseq",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "pipeline_outdir": pipeline_outdir,
        "dry_run": dry_run,
        "samplesheet": samplesheet_result,
        "command": cmd,
        "command_str": " ".join(cmd),
        "steps": steps,
        "tool_status": tools,
        "blockers": blockers,
        "warnings": warnings,
        "missing_requirements": missing_requirements,
        "biological_caveats": BIOLOGICAL_CAVEATS,
        "next_actions": [
            "Review and correct samplesheet (strandedness, sample names, pairing)",
            "Confirm genome build and annotation version",
            "Run with -profile test first to validate environment",
            "After pipeline: run biology-interpretation-reviewer on MultiQC output",
        ],
        "limitations": [
            "This tool generates a plan and samplesheet. It does not replace nf-core/rnaseq.",
            "Large reference files are never downloaded automatically.",
            "Differential expression analysis is not included in this workflow plan.",
        ],
    }


def execute(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    genome: str | None = None,
    fasta: str | None = None,
    gtf: str | None = None,
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
    provenance_dir: Path | None = None,
    timeout: int = 86400,
) -> dict[str, Any]:
    """Execute nf-core/rnaseq. Requires nextflow in PATH."""
    output_dir = Path(output_dir)
    pipeline_outdir = output_dir / "pipeline_output"

    wf_plan = plan(
        input_dir, output_dir,
        genome=genome, fasta=fasta, gtf=gtf, profile=profile,
        resume=resume, extra_args=extra_args, dry_run=False,
    )

    if wf_plan["blockers"]:
        wf_plan["step_results"] = [{
            "label": "nextflow_rnaseq",
            "status": "skipped",
            "error": f"Blocked: {wf_plan['blockers']}",
            "executed": False,
        }]
        return wf_plan

    record = run_nextflow(
        "rnaseq",
        wf_plan["command"],
        outdir=pipeline_outdir,
        dry_run=False,
        provenance_dir=provenance_dir,
        timeout=timeout,
    )

    wf_plan["step_results"] = [record]
    wf_plan["dry_run"] = False
    if record.get("pipeline_output_validation"):
        wf_plan["pipeline_output_validation"] = record["pipeline_output_validation"]
    return wf_plan
