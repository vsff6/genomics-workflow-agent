"""Public Python API. All functions return a JSON-serializable dict."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def inspect_inputs(
    input_path: str | Path,
    max_file_mb: float = 50.0,
) -> dict[str, Any]:
    """Inspect a file or directory without running any external tool."""
    input_path = Path(input_path)

    if not input_path.exists():
        return _failure(
            workflow="inspect",
            input_path=str(input_path),
            outdir="",
            error=f"Input path does not exist: {input_path}",
        )

    try:
        from genomics_workflow_agent.inspect.inspector import inspect_file, inspect_directory

        if input_path.is_file():
            raw = inspect_file(input_path)
        else:
            raw = inspect_directory(input_path)

        warnings = list(raw.get("warnings", []))
        large = [
            f["name"] for f in raw.get("files", [])
            if f.get("size_mb", 0) > max_file_mb
        ]
        if large:
            warnings.append(
                f"{len(large)} file(s) exceed {max_file_mb} MB and will not be "
                "loaded into memory - operations are delegated to CLI tools."
            )

        return {
            "status": "success",
            "workflow": "inspect",
            "input_path": str(input_path),
            "outdir": "",
            "summary": raw,
            "warnings": warnings,
            "errors": [],
            "paths": {},
            "provenance_paths": [],
        }
    except Exception as e:
        return _failure(
            workflow="inspect",
            input_path=str(input_path),
            outdir="",
            error=str(e),
        )


def plan_workflow(
    input_path: str | Path,
    workflow: str = "auto",
    outdir: str | Path | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Build a dry-run workflow plan. No external tools are called."""
    input_path = Path(input_path)
    resolved_outdir = Path(outdir) if outdir else Path("results")

    if not input_path.exists():
        return _failure(
            workflow=workflow,
            input_path=str(input_path),
            outdir=str(resolved_outdir),
            error=f"Input path does not exist: {input_path}",
        )

    supported = {"auto", "fastq-qc", "rnaseq", "atacseq", "amplicon", "variant-qc"}
    if workflow not in supported:
        return _failure(
            workflow=workflow,
            input_path=str(input_path),
            outdir=str(resolved_outdir),
            error=f"Unsupported workflow: '{workflow}'. Supported: {sorted(supported)}",
        )

    try:
        from genomics_workflow_agent.workflows.planner import build_plan

        plan = build_plan(workflow, str(input_path), resolved_outdir, dry_run=True, extra_params=kwargs)

        resolved_outdir.mkdir(parents=True, exist_ok=True)
        json_path = resolved_outdir / "plan.json"
        json_path.write_text(json.dumps(plan, indent=2, default=str), encoding="utf-8")

        return {
            "status": "dry_run",
            "workflow": plan.get("workflow", workflow),
            "input_path": str(input_path),
            "outdir": str(resolved_outdir),
            "summary": {
                "steps": len(plan.get("steps", [])),
                "skipped_steps": len(plan.get("skipped_steps", [])),
                "blockers": len(plan.get("blockers", [])),
            },
            "warnings": list(plan.get("warnings", [])),
            "errors": list(plan.get("blockers", [])),
            "paths": {"plan_json": str(json_path)},
            "provenance_paths": [],
        }
    except Exception as e:
        return _failure(
            workflow=workflow,
            input_path=str(input_path),
            outdir=str(resolved_outdir),
            error=str(e),
        )


def run_workflow(
    input_path: str | Path,
    workflow: str = "auto",
    outdir: str | Path = "results",
    execute: bool = False,
    **kwargs: Any,
) -> dict[str, Any]:
    """Plan or execute a workflow. Dry-run by default (execute=False)."""
    input_path = Path(input_path)
    outdir = Path(outdir)

    if not input_path.exists():
        return _failure(
            workflow=workflow,
            input_path=str(input_path),
            outdir=str(outdir),
            error=f"Input path does not exist: {input_path}",
        )

    supported = {"auto", "fastq-qc", "rnaseq", "atacseq", "amplicon", "variant-qc"}
    if workflow not in supported:
        return _failure(
            workflow=workflow,
            input_path=str(input_path),
            outdir=str(outdir),
            error=f"Unsupported workflow: '{workflow}'. Supported: {sorted(supported)}",
        )

    try:
        from genomics_workflow_agent.workflows.planner import build_plan, execute_plan
        from genomics_workflow_agent.reports.json_report import write_json_report

        outdir.mkdir(parents=True, exist_ok=True)
        provenance_dir = outdir / "provenance"

        if not execute:
            result = build_plan(
                workflow, str(input_path), outdir, dry_run=True, extra_params=kwargs
            )
            status = "dry_run"
        else:
            result = execute_plan(
                workflow, str(input_path), outdir,
                provenance_dir=provenance_dir,
                extra_params=kwargs,
            )
            step_results = result.get("step_results", [])
            failed = [r for r in step_results if r.get("status") == "failed"]
            status = "failed" if failed else "success"

        json_path = write_json_report(result, outdir / "run_report.json")

        provenance_paths = (
            [str(p) for p in sorted(provenance_dir.glob("*.json"))]
            if provenance_dir.exists() else []
        )

        return {
            "status": status,
            "workflow": result.get("workflow", workflow),
            "input_path": str(input_path),
            "outdir": str(outdir),
            "summary": {
                "steps": len(result.get("steps", [])),
                "step_results": len(result.get("step_results", [])),
                "skipped_steps": len(result.get("skipped_steps", [])),
            },
            "warnings": list(result.get("warnings", [])),
            "errors": list(result.get("blockers", [])),
            "paths": {"run_report_json": str(json_path)},
            "provenance_paths": provenance_paths,
        }
    except Exception as e:
        return _failure(
            workflow=workflow,
            input_path=str(input_path),
            outdir=str(outdir),
            error=str(e),
        )


def run_fastq_qc_agent(
    input_path: str | Path,
    outdir: str | Path = "results_agent",
    execute: bool = False,
    auto_trim: bool = False,
    trim_tool: str = "fastp",
    **kwargs: Any,
) -> dict[str, Any]:
    """Run the FASTQ QC agent. Dry-run by default; auto_trim requires execute=True."""
    input_path = Path(input_path)
    outdir = Path(outdir)

    if not input_path.exists():
        return _failure(
            workflow="fastq-qc",
            input_path=str(input_path),
            outdir=str(outdir),
            error=f"Input path does not exist: {input_path}",
        )

    if auto_trim and not execute:
        return _failure(
            workflow="fastq-qc",
            input_path=str(input_path),
            outdir=str(outdir),
            error="auto_trim=True requires execute=True. Trimming cannot run in dry-run mode.",
        )

    try:
        from genomics_workflow_agent.agent.fastq_agent import (
            run_fastq_agent,
            write_agent_report_json,
            write_agent_report_md,
        )

        state = run_fastq_agent(
            input_path,
            outdir,
            execute=execute,
            auto_trim=auto_trim,
            trim_tool=trim_tool,
            max_file_mb=float(kwargs.get("max_file_mb", 50.0)),
        )

        json_path = write_agent_report_json(state, outdir / "agent_report.json")
        md_path = write_agent_report_md(state, outdir / "agent_report.md")

        failed_steps = [s for s in state.executed_steps if s.get("status") == "failed"]
        status = "dry_run" if not execute else ("partial" if failed_steps else "success")

        return {
            "status": status,
            "workflow": "fastq-qc",
            "input_path": str(input_path),
            "outdir": str(outdir),
            "summary": {
                "planned_steps": len(state.planned_steps),
                "executed_steps": len(state.executed_steps),
                "observations": len(state.observations),
                "decisions": len(state.decisions),
                "recommended_actions": len(state.recommended_actions),
            },
            "warnings": list(state.warnings),
            "errors": [s.get("label", "?") for s in failed_steps],
            "paths": {
                "agent_report_json": str(json_path),
                "agent_report_md": str(md_path),
            },
            "provenance_paths": list(state.provenance_paths),
            "observations": [o.to_dict() for o in state.observations],
            "decisions": [d.to_dict() for d in state.decisions],
            "recommended_actions": [r.to_dict() for r in state.recommended_actions],
        }
    except Exception as e:
        return _failure(
            workflow="fastq-qc",
            input_path=str(input_path),
            outdir=str(outdir),
            error=str(e),
        )


def run_variant_qc_agent(
    input_path: str | Path,
    outdir: str | Path = "results_variant_agent",
    execute: bool = False,
    **kwargs: Any,
) -> dict[str, Any]:
    """Run the variant QC agent. Dry-run by default."""
    input_path = Path(input_path)
    outdir = Path(outdir)

    if not input_path.exists():
        return _failure(
            workflow="variant-qc",
            input_path=str(input_path),
            outdir=str(outdir),
            error=f"Input path does not exist: {input_path}",
        )

    try:
        from genomics_workflow_agent.agent.variant_agent import (
            run_variant_agent,
            write_variant_agent_report_json,
            write_variant_agent_report_md,
        )

        state = run_variant_agent(
            input_path,
            outdir,
            execute=execute,
            max_file_mb=float(kwargs.get("max_file_mb", 50.0)),
        )

        json_path = write_variant_agent_report_json(state, outdir / "variant_agent_report.json")
        md_path = write_variant_agent_report_md(state, outdir / "variant_agent_report.md")

        failed_steps = [s for s in state.executed_steps if s.get("status") == "failed"]
        status = "dry_run" if not execute else ("partial" if failed_steps else "success")

        return {
            "status": status,
            "workflow": "variant-qc",
            "input_path": str(input_path),
            "outdir": str(outdir),
            "summary": {
                "planned_steps": len(state.planned_steps),
                "executed_steps": len(state.executed_steps),
                "observations": len(state.observations),
                "decisions": len(state.decisions),
                "recommended_actions": len(state.recommended_actions),
            },
            "warnings": list(state.warnings),
            "errors": [s.get("label", "?") for s in failed_steps],
            "paths": {
                "variant_agent_report_json": str(json_path),
                "variant_agent_report_md": str(md_path),
            },
            "provenance_paths": list(state.provenance_paths),
            "observations": [o.to_dict() for o in state.observations],
            "decisions": [d.to_dict() for d in state.decisions],
            "recommended_actions": [r.to_dict() for r in state.recommended_actions],
        }
    except Exception as e:
        return _failure(
            workflow="variant-qc",
            input_path=str(input_path),
            outdir=str(outdir),
            error=str(e),
        )


def write_report(
    results_path: str | Path,
    outdir: str | Path | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Aggregate all JSON reports in a results directory into one final report."""
    results_path = Path(results_path)

    if not results_path.exists():
        return _failure(
            workflow="report",
            input_path=str(results_path),
            outdir=str(outdir or results_path),
            error=f"Results directory does not exist: {results_path}",
        )

    try:
        from genomics_workflow_agent.reports.json_report import write_json_report
        from genomics_workflow_agent.reports.markdown import write_markdown_report

        out = Path(outdir) if outdir else results_path

        json_files = [
            f for f in sorted(results_path.rglob("*.json"))
            if "provenance" not in f.parts
        ]

        sections: dict[str, Any] = {}
        parse_warnings: list[str] = []
        for jf in json_files:
            try:
                sections[jf.stem] = json.loads(jf.read_text(encoding="utf-8"))
            except Exception as e:
                parse_warnings.append(f"Could not read {jf}: {e}")

        aggregate = {
            "workflow": "aggregate-report",
            "results_dir": str(results_path),
            "sections": sections,
            "clinical_disclaimer": (
                "This report is for research purposes only. "
                "It does not constitute clinical advice."
            ),
        }

        json_path = write_json_report(aggregate, out / "final_report.json")
        md_path = write_markdown_report(aggregate, out / "final_report.md", title="Final Genomics Report")

        return {
            "status": "success",
            "workflow": "report",
            "input_path": str(results_path),
            "outdir": str(out),
            "summary": {"sections": len(sections), "json_files_found": len(json_files)},
            "warnings": parse_warnings,
            "errors": [],
            "paths": {
                "final_report_json": str(json_path),
                "final_report_md": str(md_path),
            },
            "provenance_paths": [],
        }
    except Exception as e:
        return _failure(
            workflow="report",
            input_path=str(results_path),
            outdir=str(outdir or results_path),
            error=str(e),
        )


def _failure(
    workflow: str,
    input_path: str,
    outdir: str,
    error: str,
) -> dict[str, Any]:
    return {
        "status": "failed",
        "workflow": workflow,
        "input_path": input_path,
        "outdir": outdir,
        "summary": {},
        "warnings": [],
        "errors": [error],
        "paths": {},
        "provenance_paths": [],
    }
