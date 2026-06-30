from __future__ import annotations

from pathlib import Path
from typing import Any

from genomics_workflow_agent.inspect.inspector import inspect_directory

WORKFLOW_NAMES = ["fastq-qc", "rnaseq", "atacseq", "amplicon", "variant-qc", "auto"]


def resolve_workflow(workflow: str, input_dir: str | Path) -> str:
    """If workflow is 'auto', infer from input files. Otherwise validate and return."""
    if workflow != "auto":
        if workflow not in WORKFLOW_NAMES:
            raise ValueError(f"Unknown workflow '{workflow}'. Choose from: {WORKFLOW_NAMES}")
        return workflow

    inspection = inspect_directory(input_dir)
    guess = inspection.get("workflow_guess", "fastq-qc")
    return guess if guess != "unknown" else "fastq-qc"


def build_plan(
    workflow: str,
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    dry_run: bool = True,
    extra_params: dict | None = None,
) -> dict[str, Any]:
    """
    Build a workflow execution plan (always dry-run in planning mode).

    Parameters
    ----------
    workflow : workflow name (or 'auto' to infer)
    input_dir : path to input data
    output_dir : where to write outputs
    dry_run : if True, commands are planned but not executed
    extra_params : workflow-specific parameters
    """
    params = extra_params or {}
    resolved = resolve_workflow(workflow, input_dir)
    return _dispatch_plan(resolved, input_dir, output_dir, dry_run=dry_run, params=params)


def execute_plan(
    workflow: str,
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    provenance_dir: Path | None = None,
    extra_params: dict | None = None,
) -> dict[str, Any]:
    """
    Execute a workflow and return annotated results.

    Calls the workflow-specific execute() function which:
    - Creates output directories
    - Runs external tools
    - Validates expected outputs
    - Captures provenance
    """
    params = extra_params or {}
    resolved = resolve_workflow(workflow, input_dir)
    return _dispatch_execute(resolved, input_dir, output_dir,
                             provenance_dir=provenance_dir, params=params)



def _dispatch_plan(
    workflow: str,
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    dry_run: bool,
    params: dict,
) -> dict[str, Any]:
    if workflow == "fastq-qc":
        from genomics_workflow_agent.workflows.fastq_qc import plan
        return plan(input_dir, output_dir, dry_run=dry_run,
                    **_pick(params, ["trim", "trimmer"]))

    if workflow == "rnaseq":
        from genomics_workflow_agent.workflows.rnaseq import plan
        return plan(input_dir, output_dir, dry_run=dry_run,
                    **_pick(params, ["genome", "fasta", "gtf", "profile", "resume"]))

    if workflow == "atacseq":
        from genomics_workflow_agent.workflows.atacseq import plan
        return plan(input_dir, output_dir, dry_run=dry_run,
                    **_pick(params, ["genome", "fasta", "gtf", "blacklist", "profile", "resume"]))

    if workflow == "amplicon":
        from genomics_workflow_agent.workflows.amplicon import plan
        return plan(input_dir, output_dir, dry_run=dry_run,
                    **_pick(params, ["primer_fw", "primer_rv", "taxonomy_db",
                                     "taxonomy_db_path", "denoiser", "profile", "resume"]))

    if workflow == "variant-qc":
        from genomics_workflow_agent.workflows.variant_qc import plan
        return plan(input_dir, output_dir, dry_run=dry_run,
                    **_pick(params, ["genome", "fasta", "known_sites", "profile", "resume"]))

    raise ValueError(f"Unhandled workflow: {workflow}")


def _dispatch_execute(
    workflow: str,
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    provenance_dir: Path | None,
    params: dict,
) -> dict[str, Any]:
    if workflow == "fastq-qc":
        from genomics_workflow_agent.workflows.fastq_qc import execute
        return execute(input_dir, output_dir, provenance_dir=provenance_dir,
                       **_pick(params, ["trim", "trimmer"]))

    if workflow == "rnaseq":
        from genomics_workflow_agent.workflows.rnaseq import execute
        return execute(input_dir, output_dir, provenance_dir=provenance_dir,
                       **_pick(params, ["genome", "fasta", "gtf", "profile", "resume"]))

    if workflow == "atacseq":
        from genomics_workflow_agent.workflows.atacseq import execute
        return execute(input_dir, output_dir, provenance_dir=provenance_dir,
                       **_pick(params, ["genome", "fasta", "gtf", "blacklist", "profile", "resume"]))

    if workflow == "amplicon":
        from genomics_workflow_agent.workflows.amplicon import execute
        return execute(input_dir, output_dir, provenance_dir=provenance_dir,
                       **_pick(params, ["primer_fw", "primer_rv", "taxonomy_db",
                                        "taxonomy_db_path", "denoiser", "profile", "resume"]))

    if workflow == "variant-qc":
        from genomics_workflow_agent.workflows.variant_qc import execute
        return execute(input_dir, output_dir, provenance_dir=provenance_dir,
                       **_pick(params, ["genome", "fasta", "known_sites", "profile", "resume"]))

    raise ValueError(f"Unhandled workflow: {workflow}")


def _pick(params: dict, allowed: list[str]) -> dict:
    return {k: v for k, v in params.items() if k in allowed and v is not None}
