# Check that all expected bioinformatics tools are available inside the container.
# Run via PowerShell: .\scripts\docker_check_tools.ps1
# Requires Docker and docker compose to be installed on the host.

Write-Host "=== Tool version check (inside container) ===" -ForegroundColor Cyan

$script = @'
for tool in fastqc multiqc fastp cutadapt samtools bcftools bedtools mosdepth nextflow python pip; do
    if command -v "$tool" &>/dev/null; then
        printf "  %-12s " "$tool:"
        case "$tool" in
            nextflow) nextflow -version 2>&1 | head -1 ;;
            *)        "$tool" --version 2>&1 | head -1 ;;
        esac
    else
        echo "  $tool: NOT FOUND"
    fi
done
echo ""
echo "=== Python package check ==="
python -c "from genomics_workflow_agent import inspect_inputs, plan_workflow, run_workflow, run_fastq_qc_agent, run_variant_qc_agent; print('  api ok')"
'@

docker compose run --rm genomics-agent bash -c $script
