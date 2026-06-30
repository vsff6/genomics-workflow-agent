"""
Example: consume an agent_report.json produced by run_fastq_qc_agent() or
run_variant_qc_agent() from another Python project.

This script reads the structured result dict and extracts observations and
decisions without knowing the internal state classes.
"""
import json
from pathlib import Path


def load_agent_report(json_path: str | Path) -> dict:
    return json.loads(Path(json_path).read_text(encoding="utf-8"))


def summarize_report(report: dict) -> None:
    print(f"Workflow  : {report.get('workflow', '?')}")
    print(f"Input     : {report.get('input_path', '?')}")
    print(f"Generated : {report.get('generated_at', '?')}")
    print()

    observations = report.get("observations", [])
    print(f"Observations ({len(observations)}):")
    for obs in observations:
        print(f"  [{obs.get('severity', '?').upper()}] {obs.get('sample', '?')} "
              f"/ {obs.get('category', '?')}: {obs.get('status', '?')}")
        print(f"    {obs.get('message', '')[:100]}")

    decisions = report.get("decisions", [])
    print(f"\nDecisions ({len(decisions)}):")
    for dec in decisions:
        executed = "EXECUTED" if dec.get("executed") else "not executed"
        print(f"  [{dec.get('decision_type', '?').upper()}] {dec.get('action', '?')} "
              f"({executed}, confidence={dec.get('confidence', '?')})")

    actions = report.get("recommended_actions", [])
    print(f"\nRecommended actions ({len(actions)}):")
    for act in actions:
        print(f"  [{act.get('priority', '?').upper()}] {act.get('action', '?')}")

    warnings = report.get("warnings", [])
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for w in warnings[:5]:
            print(f"  [WARN] {w}")

    disclaimer = report.get("clinical_disclaimer", "")
    if disclaimer:
        print(f"\nDisclaimer: {disclaimer[:200]}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python consume_agent_report.py <agent_report.json>")
        print()
        print("Example: run the fastq agent first:")
        print("  python -m genomics_workflow_agent agent --input examples/minimal_fastq "
              "--workflow fastq-qc --out results_api/demo/")
        print("  python examples/api_usage/consume_agent_report.py "
              "results_api/demo/agent_report.json")
        sys.exit(0)

    report_path = sys.argv[1]
    if not Path(report_path).exists():
        print(f"File not found: {report_path}", file=sys.stderr)
        sys.exit(1)

    report = load_agent_report(report_path)
    summarize_report(report)
