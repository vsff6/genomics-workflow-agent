"""ATAC-seq workflow — nf-core/atacseq planning and execution."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomics_workflow_agent.tools.nextflow import build_nextflow_cmd, run_nextflow
from genomics_workflow_agent.tools.samplesheets import build_atacseq_samplesheet
from genomics_workflow_agent.tools.versions import check_tools

REQUIRED_TOOLS = ["nextflow"]
OPTIONAL_TOOLS = ["docker", "singularity", "conda", "bedtools", "samtools"]

BIOLOGICAL_CAVEATS = [
    "FRiP (Fraction of Reads in Peaks) thresholds are not universal — interpret in context of tissue and protocol.",
    "TSS enrichment scores vary by cell type and genome build — do not apply generic cutoffs.",
    "Nucleosome periodicity in insert-size distribution depends on chromatin accessibility, not just library quality.",
    "Low FRiP may indicate technical failure OR global chromatin remodeling (e.g. differentiation, activation).",
    "Blacklist regions must match the genome build exactly — mismatched blacklists remove valid signal.",
    "Peak calling parameters (MACS3 q-value, effective genome size) are tissue- and protocol-dependent.",
    "Replicate concordance must be assessed before merging peak sets.",
]


def plan(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    genome: str | None = None,
    fasta: str | None = None,
    gtf: str | None = None,
    blacklist: str | None = None,
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
    dry_run: bool = True,
) -> dict[str, Any]:
    """Return an nf-core/atacseq execution plan."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    pipeline_outdir = str(output_dir / "pipeline_output")

    tools = check_tools(REQUIRED_TOOLS + OPTIONAL_TOOLS)
    warnings: list[str] = []
    blockers: list[str] = []
    missing_requirements: list[str] = []

    ss_path = output_dir / "samplesheet_atacseq.csv"
    samplesheet_result = build_atacseq_samplesheet(input_dir, ss_path)

    if not tools["nextflow"]["available"]:
        blockers.append("nextflow not found in PATH. Install: https://www.nextflow.io")

    if not genome and not fasta:
        missing_requirements.append("--genome (nf-core key) or --fasta (local path) required")
        warnings.append("No reference genome specified")

    if not blacklist:
        warnings.append(
            "No blacklist BED provided. Strongly recommended for ATAC-seq. "
            "Download: https://github.com/Bowtie-project/genome_blacklists"
        )

    if fasta and not Path(fasta).exists():
        blockers.append(f"FASTA not found: {fasta}")

    container_ok = tools["docker"]["available"] or tools["singularity"]["available"]
    if not container_ok and not tools["conda"]["available"]:
        warnings.append("No container runtime or conda found")

    cmd = build_nextflow_cmd(
        "atacseq",
        input_path=str(ss_path) if samplesheet_result.get("created") else None,
        outdir=pipeline_outdir,
        genome=genome,
        fasta=fasta,
        gtf=gtf,
        blacklist=blacklist,
        profile=profile,
        resume=resume,
        extra_args=extra_args,
    )

    steps = [{
        "name": "nextflow_atacseq",
        "description": "Run nf-core/atacseq via Nextflow",
        "command": cmd,
        "output_dir": pipeline_outdir,
        "expected_outputs": [pipeline_outdir],
        "dry_run": dry_run,
        "required_tools": ["nextflow"],
    }]

    return {
        "workflow": "atacseq",
        "pipeline": "nf-core/atacseq",
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
            "Review and correct samplesheet (replicate structure, control samples)",
            "Provide blacklist BED matching your genome build",
            "After pipeline: review FRiP, TSS enrichment, and insert-size plots",
            "Invoke biology-interpretation-reviewer on QC outputs",
        ],
        "limitations": [
            "This tool generates a plan and samplesheet only. It does not replace nf-core/atacseq.",
            "FRiP and TSS enrichment require completed alignment — not computed here.",
            "Peak-calling parameters depend on experimental context and are not set automatically.",
        ],
    }


def execute(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    genome: str | None = None,
    fasta: str | None = None,
    gtf: str | None = None,
    blacklist: str | None = None,
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
    provenance_dir: Path | None = None,
    timeout: int = 86400,
) -> dict[str, Any]:
    """Execute nf-core/atacseq. Requires nextflow in PATH."""
    output_dir = Path(output_dir)
    pipeline_outdir = output_dir / "pipeline_output"

    wf_plan = plan(
        input_dir, output_dir,
        genome=genome, fasta=fasta, gtf=gtf, blacklist=blacklist,
        profile=profile, resume=resume, extra_args=extra_args, dry_run=False,
    )

    if wf_plan["blockers"]:
        wf_plan["step_results"] = [{
            "label": "nextflow_atacseq", "status": "skipped",
            "error": f"Blocked: {wf_plan['blockers']}", "executed": False,
        }]
        return wf_plan

    record = run_nextflow(
        "atacseq", wf_plan["command"],
        outdir=pipeline_outdir, dry_run=False,
        provenance_dir=provenance_dir, timeout=timeout,
    )

    wf_plan["step_results"] = [record]
    wf_plan["dry_run"] = False
    if record.get("pipeline_output_validation"):
        wf_plan["pipeline_output_validation"] = record["pipeline_output_validation"]
    return wf_plan
