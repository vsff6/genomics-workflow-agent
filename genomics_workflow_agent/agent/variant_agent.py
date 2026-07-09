from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from genomics_workflow_agent.agent.variant_decision_engine import (
    CLINICAL_DISCLAIMER,
    LIMITATIONS,
    evaluate_variant_qc_results,
)
from genomics_workflow_agent.agent.state import AgentState
from genomics_workflow_agent.parsers.samtools import parse_variant_qc_dir
from genomics_workflow_agent.parsers.bcftools import parse_bcftools_dir
from genomics_workflow_agent.parsers.mosdepth import parse_mosdepth_dir


def run_variant_agent(
    input_path: str | Path,
    outdir: str | Path,
    *,
    execute: bool = False,
    max_file_mb: float = 50.0,
) -> AgentState:
    """Run the variant QC agent loop. Pass execute=True to call real tools."""
    input_path = Path(input_path)
    outdir = Path(outdir)
    provenance_dir = outdir / "provenance"
    outdir.mkdir(parents=True, exist_ok=True)

    state = AgentState(
        input_path=str(input_path),
        workflow="variant-qc",
    )
    state.limitations.extend(LIMITATIONS)

    state.input_summary = _inspect_inputs(input_path, max_file_mb, state)

    plan = _build_variant_plan(input_path, outdir, state)
    state.planned_steps = plan.get("steps", [])

    if not execute:
        state.warnings.append(
            "Dry-run mode: no external tools were called. "
            "Samtools/bcftools/mosdepth outputs do not exist yet. "
            "Observations and decisions require --execute."
        )
        _add_dry_run_observations(state)
        return state

    # Execute variant QC steps
    exec_result = _run_variant_qc(input_path, outdir, provenance_dir, state)
    state.executed_steps = exec_result.get("step_results", [])

    if provenance_dir.exists():
        state.provenance_paths = [
            str(p) for p in sorted(provenance_dir.glob("*.json"))
        ]

    # Parse outputs
    qc_outdir = outdir / "variant_qc"
    samtools_results = parse_variant_qc_dir(qc_outdir)
    bcftools_results = parse_bcftools_dir(qc_outdir)
    mosdepth_results = parse_mosdepth_dir(qc_outdir)

    if (
        not samtools_results.get("flagstat")
        and not bcftools_results
        and not mosdepth_results
    ):
        state.warnings.append(
            "No samtools/bcftools/mosdepth outputs could be parsed. "
            "This may mean the tools are not installed, or execution failed."
        )

    # Evaluate
    engine_result = evaluate_variant_qc_results(
        samtools_results,
        bcftools_results,
        mosdepth_results,
        execute_allowed=execute,
    )
    state.observations = engine_result["observations"]
    state.decisions = engine_result["decisions"]
    state.recommended_actions = engine_result["recommended_actions"]
    state.warnings.extend(engine_result["warnings"])

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
            "Large files are not loaded into memory - only shell commands are used."
        )
    return result


def _build_variant_plan(
    input_path: Path, outdir: Path, state: AgentState
) -> dict[str, Any]:
    from genomics_workflow_agent.workflows.variant_qc import plan as variant_plan

    try:
        return variant_plan(input_path, outdir, dry_run=True)
    except Exception as e:
        state.warnings.append(f"Plan generation failed: {e}")
        return {"steps": [], "skipped_steps": [], "warnings": [str(e)]}


def _run_variant_qc(
    input_path: Path, outdir: Path, provenance_dir: Path, state: AgentState
) -> dict[str, Any]:
    from genomics_workflow_agent.workflows.variant_qc import execute as variant_execute

    try:
        result = variant_execute(
            input_path, outdir, provenance_dir=provenance_dir
        )
    except Exception as e:
        state.warnings.append(f"Variant QC execution failed: {e}")
        return {"step_results": []}

    failed = [r for r in result.get("step_results", []) if r.get("status") == "failed"]
    if failed:
        state.warnings.append(
            f"{len(failed)} step(s) failed during variant QC execution: "
            f"{[r.get('label') for r in failed]}"
        )
    return result


def _add_dry_run_observations(state: AgentState) -> None:
    from genomics_workflow_agent.agent.state import Observation, RecommendedAction

    state.observations.append(Observation(
        source="variant_agent",
        sample="all",
        category="dry_run",
        status="missing",
        severity="info",
        message=(
            "Dry-run mode: samtools/bcftools/mosdepth were not executed. "
            "No QC outputs exist to parse. Decisions cannot be made without observed data."
        ),
        suggested_action="Re-run with --execute to generate QC outputs and enable decisions.",
    ))
    state.recommended_actions.append(RecommendedAction(
        action="Re-run with --execute to generate samtools/bcftools/mosdepth outputs",
        priority="high",
        reason="Dry-run mode cannot produce observations or decisions",
        requires_execute=True,
        requires_external_tool="samtools, bcftools",
    ))


def write_variant_agent_report_json(state: AgentState, out_path: str | Path) -> Path:
    from genomics_workflow_agent.interpretation import generate_interpretation

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = state.to_dict()
    payload["generated_at"] = datetime.now(timezone.utc).isoformat()
    payload["clinical_disclaimer"] = CLINICAL_DISCLAIMER
    payload["biological_interpretation"] = generate_interpretation(
        workflow="variant-qc",
        observations=[o.to_dict() for o in state.observations],
        decisions=[d.to_dict() for d in state.decisions],
    )
    out_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    return out_path


def write_variant_agent_report_md(state: AgentState, out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = [
        "# Agentic Variant QC Report",
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
            icon = {"succeeded": "[OK]", "failed": "[FAIL]", "skipped": "[SKIP]", "error": "[ERR]"}.get(status, "[?]")
            lines.append(f"- {icon} **{step.get('label', '?')}** - {status}")
            if step.get("error"):
                lines.append(f"  - Error: {step['error']}")

    if state.observations:
        lines.append(f"\n## Observations ({len(state.observations)})\n")
        for obs in state.observations:
            lines.append(f"### {obs.sample} - {obs.category}\n")
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

    # Biological interpretation section - rendered from structured JSON
    from genomics_workflow_agent.interpretation import generate_interpretation, render_interpretation_md

    interp = generate_interpretation(
        workflow="variant-qc",
        observations=[o.to_dict() for o in state.observations],
        decisions=[d.to_dict() for d in state.decisions],
    )
    lines.append("\n---\n")
    lines.append(render_interpretation_md(interp))

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path
