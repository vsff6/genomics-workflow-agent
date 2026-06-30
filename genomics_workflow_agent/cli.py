from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out", "-o", default="results",
                        help="Output directory (default: results/)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")


def cmd_inspect(args: argparse.Namespace) -> int:
    from genomics_workflow_agent.inspect.inspector import inspect_directory, inspect_file
    from genomics_workflow_agent.reports.json_report import write_json_report
    from genomics_workflow_agent.reports.markdown import write_inspection_report

    input_path = Path(args.input)
    out_dir = Path(args.out)

    if input_path.is_file():
        result = inspect_file(input_path)
        print(json.dumps(result, indent=2, default=str))
    else:
        result = inspect_directory(input_path, recursive=getattr(args, "recursive", False))
        out_dir.mkdir(parents=True, exist_ok=True)
        json_path = write_json_report(result, out_dir / "inspection.json")
        md_path = write_inspection_report(result, out_dir / "inspection.md")
        _print_inspection_summary(result, json_path, md_path)

    return 0


def _print_inspection_summary(result: dict, json_path: Path, md_path: Path) -> None:
    print(f"\n=== File Inspection: {result['input_path']} ===")
    print(f"  Total files : {result['total_files']}")
    print(f"  Total size  : {result['total_size_mb']} MB")
    print(f"  File types  : {result['file_type_counts']}")
    print(f"  Workflow    : {result['workflow_guess']}")
    for w in result.get("warnings", []):
        print(f"  [WARN] {w}")
    print(f"\n  JSON : {json_path}")
    print(f"  MD   : {md_path}")


def cmd_plan(args: argparse.Namespace) -> int:
    from genomics_workflow_agent.workflows.planner import build_plan
    from genomics_workflow_agent.reports.json_report import write_json_report
    from genomics_workflow_agent.reports.markdown import write_markdown_report

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    extra = _collect_extra_params(args)
    plan = build_plan(args.workflow, args.input, out_dir, dry_run=True, extra_params=extra)

    json_path = write_json_report(plan, out_dir / "plan.json")
    md_path = write_markdown_report(plan, out_dir / "plan.md")
    _print_plan_summary(plan, json_path, md_path)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    from genomics_workflow_agent.workflows.planner import build_plan, execute_plan
    from genomics_workflow_agent.reports.json_report import write_json_report
    from genomics_workflow_agent.reports.markdown import write_markdown_report
    from genomics_workflow_agent.tools.runner import execution_summary

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    provenance_dir = out_dir / "provenance"
    dry_run = not args.execute

    if dry_run:
        print("\n[DRY RUN] Commands will be shown but NOT executed.")
        print("  Pass --execute to launch external tools.\n")

    extra = _collect_extra_params(args)

    if dry_run:
        result = build_plan(args.workflow, args.input, out_dir, dry_run=True, extra_params=extra)
        _print_dry_run_steps(result)
        json_path = write_json_report(result, out_dir / "run_report.json")
        md_path = write_markdown_report(result, out_dir / "run_report.md")
    else:
        print(f"[EXECUTE] Launching workflow: {args.workflow}")
        result = execute_plan(
            args.workflow, args.input, out_dir,
            provenance_dir=provenance_dir,
            extra_params=extra,
        )
        _print_execution_summary(result)
        json_path = write_json_report(result, out_dir / "run_report.json")
        md_path = write_markdown_report(result, out_dir / "run_report.md")

    print(f"\nOutputs:")
    print(f"  JSON      : {json_path}")
    print(f"  Markdown  : {md_path}")
    if not dry_run:
        print(f"  Provenance: {provenance_dir}")

    step_results = result.get("step_results", [])
    if step_results:
        summary = execution_summary(step_results)
        if summary["overall_status"] == "failed":
            print(f"\n[WARN] {len(summary['failed_steps'])} step(s) failed: {summary['failed_steps']}")
            return 1

    return 0


def _print_dry_run_steps(result: dict) -> None:
    steps = result.get("steps", [])
    skipped = result.get("skipped_steps", [])
    blockers = result.get("blockers", [])
    print(f"=== DRY RUN: {result.get('workflow', '?')} ===")
    print(f"  Steps : {len(steps)}  Skipped: {len(skipped)}  Blockers: {len(blockers)}")
    for b in blockers:
        print(f"  [BLOCKER] {b}")
    for step in steps:
        cmd = step.get("command")
        if cmd:
            print(f"\n  [PLANNED] {step['name']}")
            print(f"    $ {' '.join(str(c) for c in cmd)}")
            for out in step.get("expected_outputs", [])[:3]:
                print(f"    -> {out}")


def _print_execution_summary(result: dict) -> None:
    step_results = result.get("step_results", [])
    print(f"\n=== Execution: {result.get('workflow', '?')} ===")
    for r in step_results:
        status = r.get("status", "?")
        label = r.get("label", "?")
        icon = {"succeeded": "[OK]", "failed": "[FAIL]", "error": "[ERR]",
                "planned": "[DRY]", "skipped": "[SKIP]"}.get(status, "[?]")
        print(f"  {icon} {label}")
        if r.get("error"):
            print(f"       ERROR: {r['error']}")
        val = r.get("output_validation")
        if val and val.get("missing"):
            print(f"       MISSING outputs: {val['missing']}")


def cmd_report(args: argparse.Namespace) -> int:
    results_dir = Path(getattr(args, "results", args.out))
    if not results_dir.exists():
        print(f"Error: results directory not found: {results_dir}", file=sys.stderr)
        return 1

    json_files = [f for f in sorted(results_dir.rglob("*.json"))
                  if "provenance" not in f.parts]
    if not json_files:
        print(f"No JSON reports found in {results_dir}", file=sys.stderr)
        return 1

    from genomics_workflow_agent.reports.json_report import write_json_report
    from genomics_workflow_agent.reports.markdown import write_markdown_report

    sections = {}
    for jf in json_files:
        try:
            data = json.loads(jf.read_text(encoding="utf-8"))
            sections[jf.stem] = data
        except Exception as e:
            print(f"  [WARN] Could not read {jf}: {e}")

    aggregate = {
        "workflow": "aggregate-report",
        "results_dir": str(results_dir),
        "sections": sections,
        "biological_caveats": [
            "This report aggregates outputs from multiple workflow steps.",
            "Each section must be reviewed independently in biological context.",
            "Successful pipeline execution is not biological or clinical validation.",
        ],
        "clinical_disclaimer": (
            "This report is for research purposes only. "
            "It does not constitute clinical advice."
        ),
    }

    json_path = write_json_report(aggregate, results_dir / "final_report.json")
    md_path = write_markdown_report(aggregate, results_dir / "final_report.md",
                                    title="Final Genomics Report")

    print(f"\nFinal report:")
    print(f"  JSON     : {json_path}")
    print(f"  Markdown : {md_path}")
    return 0


def _collect_extra_params(args: argparse.Namespace) -> dict:
    params: dict = {}
    for key in [
        "genome", "fasta", "gtf", "blacklist", "profile",
        "primer_fw", "primer_rv", "taxonomy_db", "taxonomy_db_path",
        "denoiser", "known_sites", "trimmer", "resume",
    ]:
        val = getattr(args, key, None)
        if val is not None:
            params[key] = val

    trim_val = getattr(args, "trim", None)
    if trim_val is not None:
        params["trim"] = trim_val

    return params


def _print_plan_summary(plan: dict, json_path: Path, md_path: Path) -> None:
    workflow = plan.get("workflow", "?")
    steps = plan.get("steps", [])
    skipped = plan.get("skipped_steps", [])
    blockers = plan.get("blockers", [])
    warnings = plan.get("warnings", [])

    print(f"\n=== Plan: {workflow} ===")
    print(f"  Steps planned : {len(steps)}")
    print(f"  Steps skipped : {len(skipped)}")
    print(f"  Blockers      : {len(blockers)}")
    print(f"  Warnings      : {len(warnings)}")

    for b in blockers:
        print(f"  [BLOCKER] {b}")
    for w in warnings[:5]:
        print(f"  [WARN] {w}")

    cmd_str = plan.get("command_str")
    if cmd_str:
        print(f"\n  Command (dry-run):\n    {cmd_str}")

    print(f"\n  JSON : {json_path}")
    print(f"  MD   : {md_path}")


def cmd_agent(args: argparse.Namespace) -> int:
    from genomics_workflow_agent.agent.fastq_agent import (
        run_fastq_agent,
        write_agent_report_json,
        write_agent_report_md,
    )

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    execute = getattr(args, "execute", False)
    auto_trim = getattr(args, "auto_trim", False)

    if auto_trim and not execute:
        print("[ERROR] --auto-trim requires --execute.", file=sys.stderr)
        return 1

    mode = "EXECUTE" if execute else "DRY RUN"
    trim_note = f" + AUTO-TRIM ({args.trim_tool})" if auto_trim else ""
    print(f"\n[AGENT] FASTQ QC Agent — {mode}{trim_note}")
    print(f"  Input : {args.input}")
    print(f"  Output: {out_dir}\n")

    try:
        state = run_fastq_agent(
            args.input,
            out_dir,
            execute=execute,
            auto_trim=auto_trim,
            trim_tool=getattr(args, "trim_tool", "fastp"),
            max_file_mb=getattr(args, "max_file_mb", 50.0),
        )
    except ValueError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"[ERROR] Agent failed: {e}", file=sys.stderr)
        return 1

    _print_agent_summary(state)

    json_path = write_agent_report_json(state, out_dir / "agent_report.json")
    md_path = write_agent_report_md(state, out_dir / "agent_report.md")
    state.report_paths.extend([str(json_path), str(md_path)])

    print(f"\nReports:")
    print(f"  JSON     : {json_path}")
    print(f"  Markdown : {md_path}")

    failed_steps = [s for s in state.executed_steps if s.get("status") == "failed"]
    return 1 if failed_steps else 0


def _print_agent_summary(state) -> None:
    print(f"  Steps planned  : {len(state.planned_steps)}")
    print(f"  Steps executed : {len(state.executed_steps)}")
    print(f"  Observations   : {len(state.observations)}")
    print(f"  Decisions      : {len(state.decisions)}")
    print(f"  Actions        : {len(state.recommended_actions)}")

    if state.decisions:
        print("\nDecisions:")
        for d in state.decisions:
            icon = {"trim": "[TRIM]", "review": "[REVIEW]", "accept": "[OK]",
                    "flag": "[FLAG]", "skip": "[SKIP]"}.get(d.decision_type, "[?]")
            print(f"  {icon} {d.action}")

    if state.warnings:
        print(f"\nWarnings ({len(state.warnings)}):")
        for w in state.warnings[:5]:
            print(f"  [WARN] {w}")
        if len(state.warnings) > 5:
            print(f"  ... and {len(state.warnings) - 5} more (see report)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genomics_workflow_agent",
        description="Agentic Genomics Workflow Framework — inspect, plan, and execute reproducible NGS pipelines.",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 0.3.0")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    p = sub.add_parser("inspect", help="Inspect input files and detect data types")
    p.add_argument("--input", "-i", required=True, help="Input file or directory")
    p.add_argument("--recursive", "-r", action="store_true")
    _add_common_args(p)

    p = sub.add_parser("plan", help="Build a workflow execution plan (always dry-run)")
    p.add_argument("--input", "-i", required=True)
    p.add_argument("--workflow", "-w", default="auto",
                   choices=["auto", "fastq-qc", "rnaseq", "atacseq", "amplicon", "variant-qc"])
    _add_workflow_args(p)
    _add_common_args(p)

    p = sub.add_parser("run", help="Run a workflow (dry-run by default; --execute to run for real)")
    p.add_argument("--input", "-i", required=True)
    p.add_argument("--workflow", "-w", default="auto",
                   choices=["auto", "fastq-qc", "rnaseq", "atacseq", "amplicon", "variant-qc"])
    p.add_argument("--execute", action="store_true",
                   help="Actually execute commands. WARNING: launches external tools.")
    _add_workflow_args(p)
    _add_common_args(p)

    p = sub.add_parser("report", help="Aggregate results into a final report")
    p.add_argument("--results", required=True, help="Results directory")
    _add_common_args(p)

    p = sub.add_parser(
        "agent",
        help="Run the agentic FASTQ QC loop (observe, decide, act)",
    )
    p.add_argument("--input", "-i", required=True, help="Input directory of FASTQ files")
    p.add_argument("--workflow", "-w", default="fastq-qc",
                   choices=["fastq-qc"],
                   help="Agentic workflow (currently fastq-qc only)")
    p.add_argument("--execute", action="store_true",
                   help="Execute FastQC/MultiQC. Without this, only a dry-run plan is produced.")
    p.add_argument("--auto-trim", dest="auto_trim", action="store_true",
                   help="Automatically run trimming if QC evidence recommends it (requires --execute)")
    p.add_argument("--trim-tool", dest="trim_tool", default="fastp",
                   choices=["fastp", "cutadapt"],
                   help="Trimmer to use when --auto-trim is active (default: fastp)")
    p.add_argument("--max-file-mb", dest="max_file_mb", type=float, default=50.0,
                   help="Warn about input files larger than this size in MB (default: 50)")
    p.add_argument("--json", dest="write_json", action="store_true", default=True,
                   help="Write agent_report.json (default: on)")
    p.add_argument("--markdown", dest="write_markdown", action="store_true", default=True,
                   help="Write agent_report.md (default: on)")
    _add_common_args(p)

    return parser


def _add_workflow_args(parser: argparse.ArgumentParser) -> None:
    g = parser.add_argument_group("workflow parameters")
    g.add_argument("--genome", help="nf-core genome key (e.g. GRCh38)")
    g.add_argument("--fasta", help="Local genome FASTA path")
    g.add_argument("--gtf", help="GTF annotation path")
    g.add_argument("--blacklist", help="Blacklist BED (ATAC-seq)")
    g.add_argument("--profile", default="docker", help="Nextflow profile (default: docker)")
    g.add_argument("--resume", action="store_true", help="Pass -resume to Nextflow")
    g.add_argument("--primer-fw", dest="primer_fw", help="Forward primer sequence (amplicon)")
    g.add_argument("--primer-rv", dest="primer_rv", help="Reverse primer sequence (amplicon)")
    g.add_argument("--taxonomy-db", dest="taxonomy_db", default="SILVA",
                   choices=["SILVA", "GTDB", "UNITE", "Greengenes2", "custom"])
    g.add_argument("--taxonomy-db-path", dest="taxonomy_db_path",
                   help="Local taxonomy classifier path")
    g.add_argument("--denoiser", default="dada2", choices=["dada2", "deblur"])
    g.add_argument("--known-sites", dest="known_sites", help="Known sites VCF for BQSR")
    g.add_argument("--trim", nargs="?", const="fastp", default=None, metavar="TRIMMER",
                   help="Enable trimming. Optionally specify trimmer: fastp (default) or cutadapt")
    g.add_argument("--trimmer", default="fastp", choices=["fastp", "cutadapt"],
                   help="Trimmer when --trim is used as a flag (default: fastp)")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 0

    handlers = {
        "inspect": cmd_inspect,
        "plan": cmd_plan,
        "run": cmd_run,
        "report": cmd_report,
        "agent": cmd_agent,
    }

    handler = handlers.get(args.command)
    if handler is None:
        print(f"Unknown command: {args.command}", file=sys.stderr)
        return 1

    return handler(args)


if __name__ == "__main__":
    sys.exit(main())
