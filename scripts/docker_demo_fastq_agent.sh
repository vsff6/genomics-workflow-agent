#!/usr/bin/env bash
# Demo: run the FASTQ QC agent (dry-run) on the bundled examples.
# No real genomics tools are called in dry-run mode.
# No large data files or references are required.
#
# Run via: docker compose run --rm genomics-agent bash scripts/docker_demo_fastq_agent.sh
# Or from host:  bash scripts/docker_demo_fastq_agent.sh

set -euo pipefail

if [[ "${IN_CONTAINER:-}" == "1" ]]; then
    echo "=== FASTQ QC Agent - dry-run demo ==="
    python -m genomics_workflow_agent agent \
        --input examples/ \
        --workflow fastq-qc \
        --out results_agent_smoke/

    echo ""
    echo "Reports written to results_agent_smoke/"
    ls -lh results_agent_smoke/ 2>/dev/null || true
else
    docker compose run --rm -e IN_CONTAINER=1 genomics-agent bash scripts/docker_demo_fastq_agent.sh
fi
