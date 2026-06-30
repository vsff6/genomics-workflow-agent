# Demo: run the variant QC agent (dry-run) on the bundled examples.
# No real genomics tools are called in dry-run mode.
# No large data files or references are required.
# Usage: .\scripts\docker_demo_variant_agent.ps1

Write-Host "=== Variant QC Agent - dry-run demo ===" -ForegroundColor Cyan

docker compose run --rm genomics-agent `
    python -m genomics_workflow_agent agent `
        --input examples/ `
        --workflow variant-qc `
        --out results_variant_agent_smoke/

Write-Host ""
Write-Host "Reports written to results_variant_agent_smoke/" -ForegroundColor Green
