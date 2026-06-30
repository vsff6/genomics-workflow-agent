from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from genomics_workflow_agent.agent.decision_engine import (
    CLINICAL_DISCLAIMER,
    LIMITATIONS,
    evaluate_fastqc_results,
)
from genomics_workflow_agent.agent.state import AgentState
from genomics_workflow_agent.parsers.fastqc import parse_fastqc_dir
from genomics_workflow_agent.parsers.multiqc import parse_multiqc_output


def run_fastq_agent(
    input_path: str | Path,
    outdir: str | Path,
    *,
    execute: bool = False,
    auto_trim: bool = False,
    trim_tool: str = "fastp",
    max_file_mb: float = 50.0,
) -> AgentState:
    """
    Observation → Decision → Action agent for FASTQ QC.

    Dry-run (execute=False):
      Inspects inputs and builds a plan. Does not call any external tool.
      Reports what would be observed, and what cannot yet be decided.

    Execute (execute=True):
      Runs FastQC and MultiQC. Parses outputs. Evaluates results.
      If auto_trim=True and trimming is recommended, also runs trimming.
    """
    if auto_trim and not execute:
        raise ValueError("--auto-trim requires --execute. Trimming cannot run in dry-run mode.")

    input_path = Path(input_path)
    outdir = Path(outdir)
    provenance_dir = outdir / "provenance"
    outdir.mkdir(parents=True, exist_ok=True)

    state = AgentState(
        input_path=str(input_path),
        workflow="fastq-qc",
    )

    state.limitations.extend(LIMITATIONS)

    # Step 1: inspect inputs
    state.input_summary = _inspect_inputs(input_path, max_file_mb, state)

    # Step 2: build plan
    plan = _build_qc_plan(input_path, outdir, trim_tool, state)
    state.planned_steps = plan.get("steps", [])

    if not execute:
        state.warnings.append(
            "Dry-run mode: no external tools were called. "
            "Observations, decisions, and trimming recommendations require --execute."
        )
        _add_dry_run_observations(state)
        return state

    # Step 3: execute FastQC/MultiQC
    qc_result = _run_qc(input_path, outdir, provenance_dir, state)
    state.executed_steps = qc_result.get("step_results", [])
    for p in (provenance_dir / f for f in [] if provenance_dir.exists()):
        state.provenance_paths.append(str(p))

    # Step 4: parse outputs
    fastqc_outdir = outdir / "fastqc"
    multiqc_outdir = outdir / "multiqc"
    parsed_fastqc = parse_fastqc_dir(fastqc_outdir)
    multiqc_data = parse_multiqc_output(multiqc_outdir)

    if not parsed_fastqc and not multiqc_data.get("parse_ok"):
        state.warnings.append(
            "FastQC/MultiQC outputs not found or could not be parsed. "
            "This may mean FastQC/MultiQC are not installed, or execution failed."
        )

    # Step 5: evaluate
    engine_result = evaluate_fastqc_results(parsed_fastqc, execute_allowed=execute)
    state.observations = engine_result["observations"]
    state.decisions = engine_result["decisions"]
    state.recommended_actions = engine_result["recommended_actions"]
    state.warnings.extend(engine_result["warnings"])

    # Step 6: auto-trim if requested and evidence supports it
    if auto_trim:
        _maybe_auto_trim(
            input_path, outdir, provenance_dir,
            trim_tool, state, engine_result["decisions"],
        )

    return state


def _inspect_inputs(input_path: Path, max_file_mb: float, state: AgentState) -> dict[str, Any]:
    from genomics_workflow_agent.inspect.inspector import inspect_file, inspect_directory

    try:
        if input_path.is_file():
            result = inspect_file(input_path)
        else:
            result = inspect_directory(input_path)
    except Exception as e:
        state.warnings.append(f"Input inspection failed: {e}")
        return {"error": str(e)}

    large_files = [
        f["name"] for f in result.get("files", [])
        if f.get("size_mb", 0) > max_file_mb
    ]
    if large_files:
        state.warnings.append(
            f"{len(large_files)} file(s) exceed {max_file_mb} MB. "
            "Large files are not loaded into memory — only shell commands are used."
        )
    return result


def _build_qc_plan(
    input_path: Path, outdir: Path, trim_tool: str, state: AgentState
) -> dict[str, Any]:
    from genomics_workflow_agent.workflows.fastq_qc import plan as fastq_plan

    try:
        return fastq_plan(input_path, outdir, dry_run=True)
    except Exception as e:
        state.warnings.append(f"Plan generation failed: {e}")
        return {"steps": [], "skipped_steps": [], "warnings": [str(e)]}


def _run_qc(
    input_path: Path, outdir: Path, provenance_dir: Path, state: AgentState
) -> dict[str, Any]:
    from genomics_workflow_agent.workflows.fastq_qc import execute as fastq_execute

    try:
        result = fastq_execute(
            input_path, outdir,
            provenance_dir=provenance_dir,
        )
    except Exception as e:
        state.warnings.append(f"FastQC/MultiQC execution failed: {e}")
        return {"step_results": []}

    failed = [r for r in result.get("step_results", []) if r.get("status") == "failed"]
    if failed:
        state.warnings.append(
            f"{len(failed)} step(s) failed during QC execution: "
            f"{[r.get('label') for r in failed]}"
        )
    return result


def _add_dry_run_observations(state: AgentState) -> None:
    from genomics_workflow_agent.agent.state import Observation, RecommendedAction

    state.observations.append(Observation(
        source="fastq_agent",
        sample="all",
        category="dry_run",
        status="missing",
        severity="info",
        message=(
            "Dry-run mode: FastQC/MultiQC were not executed. "
            "No QC outputs exist to parse. Decisions cannot be made without observed data."
        ),
        suggested_action="Re-run with --execute to generate QC outputs and enable decisions.",
    ))
    state.recommended_actions.append(RecommendedAction(
        action="Re-run with --execute to generate FastQC/MultiQC outputs",
        priority="high",
        reason="Dry-run mode cannot produce observations or decisions",
        requires_execute=True,
        requires_external_tool="fastqc, multiqc",
    ))


def _maybe_auto_trim(
    input_path: Path,
    outdir: Path,
    provenance_dir: Path,
    trim_tool: str,
    state: AgentState,
    decisions: list,
) -> None:
    from genomics_workflow_agent.agent.state import Decision, RecommendedAction
    from genomics_workflow_agent.workflows.fastq_qc import execute as fastq_execute

    trim_decisions = [d for d in decisions if d.decision_type == "trim"]

    if not trim_decisions:
        state.warnings.append(
            "auto-trim: trimming was requested but the decision engine found no trimming triggers. "
            "Trimming was NOT run."
        )
        state.recommended_actions.append(RecommendedAction(
            action="No trimming run — decision engine found no adapter or quality failures",
            priority="low",
            reason="auto-trim requires a positive trimming decision from QC evaluation",
            requires_execute=False,
        ))
        return

    trim_dir = outdir / "trimmed"
    state.warnings.append(
        f"auto-trim: trimming decision triggered. Running {trim_tool}. "
        "Original FASTQ files will NOT be modified."
    )

    try:
        trim_result = fastq_execute(
            input_path,
            outdir,
            trim=trim_tool,
            provenance_dir=provenance_dir,
        )
        trim_steps = trim_result.get("step_results", [])
        state.executed_steps.extend(trim_steps)

        for d in trim_decisions:
            d.executed = True

        trim_failures = [r for r in trim_steps if r.get("status") == "failed"]
        if trim_failures:
            state.warnings.append(
                f"auto-trim: {len(trim_failures)} trim step(s) failed: "
                f"{[r.get('label') for r in trim_failures]}"
            )
        else:
            state.recommended_actions.append(RecommendedAction(
                action=f"Review trimmed outputs in {trim_dir}",
                priority="medium",
                reason=f"Trimming completed with {trim_tool}",
                requires_execute=False,
            ))
    except Exception as e:
        state.warnings.append(f"auto-trim: execution failed: {e}")


def write_agent_report_json(state: AgentState, out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = state.to_dict()
    payload["generated_at"] = datetime.now(timezone.utc).isoformat()
    payload["clinical_disclaimer"] = CLINICAL_DISCLAIMER
    out_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    return out_path


def write_agent_report_md(state: AgentState, out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = [
        "# Agentic FASTQ QC Report",
        f"\n**Generated**: {datetime.now(timezone.utc).isoformat()}",
        f"\n**Input**: `{state.input_path}`",
        f"\n**Workflow**: `{state.workflow}`",
        "\n---\n",
        "## Input Summary\n",
    ]

    summary = state.input_summary
    lines.append(f"- Total files: **{summary.get('total_files', '?')}**")
    lines.append(f"- Total size: **{summary.get('total_size_mb', '?')} MB**")
    lines.append(f"- Workflow guess: `{summary.get('workflow_guess', 'unknown')}`")

    if state.planned_steps:
        lines.append(f"\n## Planned Steps ({len(state.planned_steps)})\n")
        for step in state.planned_steps:
            cmd = step.get("command", [])
            lines.append(f"- **{step.get('name', '?')}**: `{' '.join(str(c) for c in cmd)}`")

    if state.executed_steps:
        lines.append(f"\n## Executed Steps ({len(state.executed_steps)})\n")
        for step in state.executed_steps:
            status = step.get("status", "?")
            icon = {"succeeded": "✓", "failed": "✗", "skipped": "○", "error": "!"}.get(status, "?")
            lines.append(f"- {icon} **{step.get('label', '?')}** — {status}")
            if step.get("error"):
                lines.append(f"  - Error: {step['error']}")

    if state.observations:
        lines.append(f"\n## Observations ({len(state.observations)})\n")
        for obs in state.observations:
            sev = {"critical": "🔴", "warning": "🟡", "info": "🟢"}.get(obs.severity, "•")
            lines.append(f"### {sev} {obs.sample} — {obs.category}\n")
            lines.append(f"**Status**: `{obs.status}` | **Severity**: {obs.severity}\n")
            lines.append(f"{obs.message}\n")
            if obs.suggested_action:
                lines.append(f"> **Suggested action**: {obs.suggested_action}\n")

    if state.decisions:
        lines.append(f"\n## Decisions ({len(state.decisions)})\n")
        for dec in state.decisions:
            lines.append(f"### {dec.action}\n")
            lines.append(f"- **Type**: {dec.decision_type}")
            lines.append(f"- **Confidence**: {dec.confidence}")
            lines.append(f"- **Executed**: {'Yes' if dec.executed else 'No'}")
            lines.append(f"- **Reason**: {dec.reason}\n")
            if dec.evidence:
                lines.append("**Evidence**:")
                for e in dec.evidence:
                    lines.append(f"  - {e}")
            if dec.safety_notes:
                lines.append("\n**Safety notes**:")
                for n in dec.safety_notes:
                    lines.append(f"  - {n}")
            lines.append("")

    if state.recommended_actions:
        lines.append(f"\n## Recommended Actions ({len(state.recommended_actions)})\n")
        for i, act in enumerate(state.recommended_actions, 1):
            lines.append(f"{i}. **[{act.priority.upper()}]** {act.action}")
            lines.append(f"   - Reason: {act.reason}")
            if act.command_preview:
                lines.append(f"   - Preview: `{act.command_preview}`")
            if act.requires_execute:
                lines.append(f"   - Requires: `--execute` + `{act.requires_external_tool}`")

    if state.warnings:
        lines.append(f"\n## Warnings ({len(state.warnings)})\n")
        for w in state.warnings:
            lines.append(f"- {w}")

    lines.append("\n## Limitations\n")
    for lim in state.limitations:
        lines.append(f"- {lim}")

    lines.append(f"\n## Disclaimer\n\n> {CLINICAL_DISCLAIMER}\n")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path
