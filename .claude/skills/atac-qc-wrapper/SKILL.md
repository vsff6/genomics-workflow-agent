# Skill: atac-qc-wrapper

## Purpose
Wrapper for ATAC-seq and scATAC-seq QC. For full pipeline orchestration, prefers `nextflow-development@life-sciences` with `nf-core/atacseq`. For local QC summarization, uses `tools/atac_qc_local.py`. Clearly reports which metrics were skipped due to missing reference files.

## When to use
- Input is fragments file, peaks BED, or BAM from ATAC-seq experiment.
- User wants local ATAC QC metrics, insert-size plots, FRiP, or peak summaries.
- Before differential accessibility, motif analysis, or footprinting.

## Official tool to prefer
`nextflow-development@life-sciences` with `nf-core/atacseq` - for full ATAC pipeline (alignment, peak calling, QC). Use when Nextflow + Docker/Singularity are available and user needs the full pipeline.

## Required inputs
- Fragments file (`.tsv.gz` for scATAC) OR BAM (for bulk ATAC)

## Optional inputs
- Peaks BED (for FRiP)
- GTF/GFF (for TSS enrichment - never calculate without this)
- Blacklist BED (for blacklist fraction - never calculate without this)
- Genome FASTA (for reference validation)
- Chromosome sizes

## Workflow steps

1. **Assess which workflow is needed**: full pipeline or local QC summary?
2. **For full pipeline**: use `nextflow-pipeline-specialist` agent with `nf-core/atacseq`.
3. **For local QC**:
```bash
python tools/atac_qc_local.py \
  --fragments <path> \
  --peaks <path_or_omit> \
  --gtf <path_or_omit> \
  --blacklist <path_or_omit> \
  --output-dir <output_dir>
```
4. **Document skipped metrics**: for each optional input not provided, list what metric was skipped and why.
5. **Biological interpretation**: apply artifact-versus-biology review to all QC observations.

## Expected outputs
- ATAC QC metrics CSV
- Insert-size distribution plot (if data available)
- FRiP summary (only if peaks BED provided)
- TSS enrichment (only if GTF/GFF provided)
- Blacklist overlap (only if blacklist BED provided)
- Skipped-metrics report
- Markdown report section

## Never do
- Calculate TSS enrichment without GTF/GFF
- Calculate blacklist overlap without blacklist BED
- Download references automatically
- Treat FRiP thresholds as universal

## Biological reasoning checklist
- [ ] Insert-size periodicity consistent with nucleosomal banding?
- [ ] FRiP interpreted in context of tissue/cell type?
- [ ] TSS enrichment compared against expected value for protocol?
- [ ] Duplicate rate consistent with expected library complexity?
- [ ] Low-FRiP samples considered for global chromatin remodeling?
