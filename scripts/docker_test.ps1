# Run API import smoke test and full pytest suite inside the container.
# Usage: .\scripts\docker_test.ps1
# Requires Docker and docker compose on the host.

Write-Host "=== API import smoke test ===" -ForegroundColor Cyan
docker compose run --rm genomics-agent python -c @'
from genomics_workflow_agent import inspect_inputs, plan_workflow, run_workflow, run_fastq_qc_agent, run_variant_qc_agent
print("api ok")
'@

Write-Host ""
Write-Host "=== pytest ===" -ForegroundColor Cyan
docker compose run --rm genomics-agent python -m pytest tests/ -q
