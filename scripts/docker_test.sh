#!/usr/bin/env bash
# Run API import smoke test and full pytest suite inside the container.
# Run via: docker compose run --rm genomics-agent bash scripts/docker_test.sh
# Or from host:  bash scripts/docker_test.sh   (wraps docker compose)

set -euo pipefail

if [[ "${IN_CONTAINER:-}" == "1" ]]; then
    echo "=== API import smoke test ==="
    python -c "
from genomics_workflow_agent import inspect_inputs, plan_workflow, run_workflow, run_fastq_qc_agent, run_variant_qc_agent
print('api ok')
"

    echo ""
    echo "=== pytest ==="
    python -m pytest tests/ -q
else
    echo "Running tests inside container..."
    docker compose run --rm -e IN_CONTAINER=1 genomics-agent bash scripts/docker_test.sh
fi
