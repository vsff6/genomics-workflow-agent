"""
Example: run the agentic variant QC loop in dry-run mode using the public API.

No external tools (samtools, bcftools, mosdepth) are required.
The agent explains what would need to run and what it cannot decide without data.
"""
import json
from pathlib import Path

from genomics_workflow_agent import run_variant_qc_agent

INPUT_DIR = Path(__file__).parent.parent / "wgs" if Path("examples/wgs").exists() else Path(".")

result = run_variant_qc_agent(
    input_path=INPUT_DIR,
    outdir="results_api/variant_agent_dryrun",
    execute=False,
)

print("=== Variant QC Agent (dry-run) ===")
print(f"Status              : {result['status']}")
print(f"Observations        : {result['summary'].get('observations', 0)}")
print(f"Decisions           : {result['summary'].get('decisions', 0)}")
print(f"Recommended actions : {result['summary'].get('recommended_actions', 0)}")

if result["observations"]:
    print(f"\nObservations:")
    for obs in result["observations"]:
        print(f"  [{obs['severity'].upper()}] {obs['sample']} / {obs['category']}: {obs['status']}")
        print(f"    {obs['message'][:120]}")

if result["warnings"]:
    print(f"\nWarnings ({len(result['warnings'])}):")
    for w in result["warnings"][:3]:
        print(f"  [WARN] {w}")

print(f"\nReport JSON : {result['paths'].get('variant_agent_report_json', 'not written')}")
print(f"Report MD   : {result['paths'].get('variant_agent_report_md', 'not written')}")
