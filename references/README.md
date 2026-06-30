# Reference Files

This directory is a placeholder. Reference files are **not** committed to the repository.  
Users must provide local reference files for reference-dependent metrics.

## Why references are not auto-downloaded

- Large reference files (genome FASTAs, annotation GTFs) are tens to hundreds of GB.
- Downloading them automatically during QC runs would be slow, expensive, and unpredictable.
- Reference file versions must be tracked explicitly for reproducibility.
- Using the wrong reference can silently produce incorrect results.

## Rules

- **Tools must not automatically download large reference files during ordinary QC runs.**
- **If a reference file is missing, skip dependent metrics and report the limitation.**
- **Never invent genome build or annotation version.**
- **Record full reference file paths and checksums where practical.**

## Expected reference files (examples)

Place files here or update config/paths accordingly:

| File | Description | Example Source |
|------|-------------|----------------|
| `GRCh38.primary_assembly.genome.fa` | Human genome FASTA | GENCODE, Ensembl, UCSC |
| `GRCh38.primary_assembly.genome.fa.fai` | FASTA index | samtools faidx |
| `gencode.v44.annotation.gtf` | Human gene annotation | GENCODE |
| `GRCh38-blacklist.v2.bed` | ENCODE blacklist regions | ENCODE portal |
| `GRCh38.chrom.sizes` | Chromosome sizes | UCSC genome browser |
| `homo_sapiens_snps.vcf.gz` | Known SNP sites | GATK resource bundle |
| `mm39.primary_assembly.genome.fa` | Mouse genome FASTA | GENCODE |
| `gencode.vM33.annotation.gtf` | Mouse gene annotation | GENCODE |

## Genome build guidance

| Build | Species | Note |
|-------|---------|------|
| GRCh38 / hg38 | Human | Current standard |
| GRCh37 / hg19 | Human | Legacy, still common for clinical |
| mm39 / GRCm39 | Mouse | Current standard |
| mm10 / GRCm38 | Mouse | Previous standard, still widely used |

## Chromosome naming conventions

| Style | Example | Used by |
|-------|---------|---------|
| UCSC | chr1, chrX, chrM | UCSC, 10x Genomics |
| Ensembl | 1, X, MT | Ensembl, GENCODE (non-UCSC) |

**Chromosome naming must be consistent across all reference files and input data.**  
Mismatch (chr1 vs 1) will silently drop all variants, peaks, or reads.

## Marker gene lists

For biological interpretation, create:
- `references/markers_human_pbmc.csv` - Known PBMC marker genes
- `references/markers_mouse_brain.csv` - Mouse brain cell type markers
- `references/markers_tumor_<type>.csv` - Tumor-specific markers

Format:
```csv
cell_type,gene,evidence
T cell,CD3D,canonical
B cell,CD19,canonical
NK cell,NCAM1,canonical
```

## Validation

Before analysis, run:
```bash
python tools/reference_validator.py \
  --genome-fasta references/GRCh38.fa \
  --gtf references/gencode.v44.annotation.gtf \
  --blacklist references/GRCh38-blacklist.v2.bed \
  --output-dir reports/reference_check
```
