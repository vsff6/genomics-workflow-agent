#!/usr/bin/env bash
# Check that all expected bioinformatics tools are available inside the container.
# Run via: docker compose run --rm genomics-agent bash scripts/docker_check_tools.sh

set -euo pipefail

echo "=== Tool version check ==="

check() {
    local tool="$1"
    shift
    if command -v "$tool" &>/dev/null; then
        echo -n "  $tool: "
        "$@" 2>&1 | head -1 || true
    else
        echo "  $tool: NOT FOUND"
    fi
}

check fastqc     fastqc --version
check multiqc    multiqc --version
check fastp      fastp --version
check cutadapt   cutadapt --version
check samtools   samtools --version
check bcftools   bcftools --version
check bedtools   bedtools --version
check mosdepth   mosdepth --version
check nextflow   nextflow -version
check python     python --version
check pip        pip --version

echo ""
echo "=== Python package check ==="
python -c "from genomics_workflow_agent import inspect_inputs, plan_workflow, run_workflow, run_fastq_qc_agent, run_variant_qc_agent; print('  api ok')"

echo ""
echo "=== Done ==="
