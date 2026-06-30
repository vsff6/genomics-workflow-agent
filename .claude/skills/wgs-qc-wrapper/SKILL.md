# Skill: wgs-qc-wrapper

## Purpose
Wrapper for WGS/WES BAM/CRAM/VCF QC. For full pipeline orchestration, prefers `nextflow-development@life-sciences` with `nf-core/sarek`. For local QC summarization, uses `tools/wgs_vcf_qc_local.py`. Never makes clinical claims.

## When to use
- Input is BAM, CRAM, or VCF.
- User wants alignment QC, coverage summary, variant counts, or Ti/Tv ratio.
- Before variant interpretation or population analysis.

## Official tool to prefer
`nextflow-development@life-sciences` with `nf-core/sarek` - for full WGS/WES pipeline (alignment, BQSR, variant calling). Use when Nextflow + Docker/Singularity are available and user needs the full pipeline.

## Required inputs
- BAM or CRAM (for alignment QC) OR VCF (for variant QC)

## Optional inputs
- Genome FASTA (required for CRAM decoding)
- Known-sites VCF (for BQSR context)
- Annotation VCF (for consequence annotation)
- Target intervals (for WES coverage)
- Pedigree file (for family analysis)

## Workflow steps

1. **Assess which workflow is needed**: full pipeline or local QC summary?
2. **For full pipeline**: use `nextflow-pipeline-specialist` agent with `nf-core/sarek`.
3. **For local QC**:
```bash
python tools/wgs_vcf_qc_local.py \
  --bam <path_or_omit> \
  --vcf <path_or_omit> \
  --genome-fasta <path_or_omit> \
  --output-dir <output_dir>
```
4. **Document limitations**: list every metric skipped and why.
5. **Biological interpretation**: apply artifact-versus-biology review.
6. **Never make clinical claims** regardless of variant findings.

## Expected outputs
- Alignment QC (if BAM/CRAM)
- Coverage summary (if BAM/CRAM)
- Variant count summary (if VCF)
- Ti/Tv ratio, het/hom ratio
- Skipped-metrics report
- Markdown report section

## Never do
- Make clinical claims
- Interpret pathogenicity as fact
- Recommend medical action
- Annotate variants without annotation source
- Invent reference build, ancestry, or pedigree

## Biological reasoning checklist
- [ ] Ti/Tv ratio consistent with expected for variant type (WGS ~2.0-2.1, WES ~2.8-3.0)?
- [ ] Het/hom ratio consistent with expected ploidy and ancestry?
- [ ] Coverage uniformity considered (GC bias, repeat regions)?
- [ ] Allele balance reviewed for expected diploid genotypes?
- [ ] Duplicate rate consistent with library preparation method?
- [ ] Contamination estimation performed or noted as missing?
