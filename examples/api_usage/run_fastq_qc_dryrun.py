"""
Example: run the FASTQ QC workflow in dry-run mode using the public API.

No external tools (FastQC, MultiQC, fastp) are required.
Shows what steps would be executed and why.
"""
import json
from pathlib import Path

from genomics_workflow_agent import run_workflow

INPUT_DIR = Path(__file__).parent.parent / "minimal_fastq"
if not INPUT_DIR.exists():
    INPUT_DIR = Path(".")

result = run_workflow(
    input_path=INPUT_DIR,
    workflow="fastq-qc",
    outdir="results_api/fastq_qc_dryrun",
    execute=False,
)

print("=== FASTQ QC Dry-Run ===")
print(f"Status      : {result['status']}")
print(f"Steps       : {result['summary'].get('steps', 0)}")
print(f"Skipped     : {result['summary'].get('skipped_steps', 0)}")

if result["warnings"]:
    print(f"\nWarnings ({len(result['warnings'])}):")
    for w in result["warnings"][:5]:
        print(f"  [WARN] {w}")

if result["errors"]:
    print(f"\nBlockers ({len(result['errors'])}):")
    for e in result["errors"]:
        print(f"  [BLOCK] {e}")

if result["paths"].get("run_report_json"):
    print(f"\nReport: {result['paths']['run_report_json']}")
