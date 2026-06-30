"""
Example: run the agentic FASTQ QC loop in dry-run mode using the public API.

No external tools (FastQC, MultiQC, fastp) are required.
The agent explains what observations and decisions would require --execute.
"""
import json
from pathlib import Path

from genomics_workflow_agent import run_fastq_qc_agent

INPUT_DIR = Path(__file__).parent.parent / "minimal_fastq"
if not INPUT_DIR.exists():
    INPUT_DIR = Path(".")

result = run_fastq_qc_agent(
    input_path=INPUT_DIR,
    outdir="results_api/fastq_agent_dryrun",
    execute=False,
)

print("=== FASTQ QC Agent (dry-run) ===")
print(f"Status              : {result['status']}")
print(f"Observations        : {result['summary'].get('observations', 0)}")
print(f"Decisions           : {result['summary'].get('decisions', 0)}")
print(f"Recommended actions : {result['summary'].get('recommended_actions', 0)}")

if result["observations"]:
    print(f"\nObservations:")
    for obs in result["observations"]:
        print(f"  [{obs['severity'].upper()}] {obs['sample']} / {obs['category']}: {obs['status']}")
        print(f"    {obs['message'][:120]}")

if result["recommended_actions"]:
    print(f"\nRecommended actions:")
    for act in result["recommended_actions"]:
        print(f"  [{act['priority'].upper()}] {act['action']}")

print(f"\nReport JSON : {result['paths'].get('agent_report_json', 'not written')}")
print(f"Report MD   : {result['paths'].get('agent_report_md', 'not written')}")
