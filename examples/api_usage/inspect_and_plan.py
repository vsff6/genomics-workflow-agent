"""
Example: inspect inputs and build a dry-run workflow plan using the public API.

Requires no external bioinformatics tools.
"""
import json
from pathlib import Path

from genomics_workflow_agent import inspect_inputs, plan_workflow

INPUT_DIR = Path(__file__).parent.parent / "minimal_fastq"

if not INPUT_DIR.exists():
    INPUT_DIR = Path(".")

result = inspect_inputs(INPUT_DIR)
print("=== Inspection ===")
print(f"Status         : {result['status']}")
print(f"Total files    : {result['summary'].get('total_files', 0)}")
print(f"Workflow guess : {result['summary'].get('workflow_guess', 'unknown')}")
if result["warnings"]:
    for w in result["warnings"]:
        print(f"[WARN] {w}")

plan = plan_workflow(INPUT_DIR, workflow="fastq-qc", outdir="results_api/plan")
print("\n=== Plan ===")
print(f"Status         : {plan['status']}")
print(f"Steps planned  : {plan['summary'].get('steps', 0)}")
print(f"Steps skipped  : {plan['summary'].get('skipped_steps', 0)}")
if plan["errors"]:
    for e in plan["errors"]:
        print(f"[BLOCKER] {e}")

print("\nFull plan JSON (first 500 chars):")
print(json.dumps(plan, indent=2)[:500])
