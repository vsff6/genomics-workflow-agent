# Demo: run the FASTQ QC agent (dry-run) on the bundled examples.
# No real genomics tools are called in dry-run mode.
# No large data files or references are required.
# Usage: .\scripts\docker_demo_fastq_agent.ps1

Write-Host "=== FASTQ QC Agent - dry-run demo ===" -ForegroundColor Cyan

docker compose run --rm genomics-agent `
    python -m genomics_workflow_agent agent `
        --input examples/ `
        --workflow fastq-qc `
        --out results_agent_smoke/

Write-Host ""
Write-Host "Reports written to results_agent_smoke/" -ForegroundColor Green
