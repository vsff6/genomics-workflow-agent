# Examples

## Toy data files

These minimal files are included for testing and demo purposes only. They are not biologically meaningful.

| File | Format | Purpose |
|------|--------|---------|
| `tiny_counts.csv` | Count matrix (genes × cells, CSV) | scRNA-seq QC fallback demo |
| `tiny_fragments.tsv` | ATAC fragments (TSV) | ATAC-seq QC demo |
| `tiny_peaks.bed` | Peak regions (BED) | ATAC FRiP and peak count demo |
| `tiny.vcf` | Variant calls (VCF) | WGS/VCF QC demo |
| `tiny.gtf` | Gene annotation (GTF) | Reference validation and TSS enrichment demo |
| `tiny.bed` | Blacklist regions (BED) | Reference validation demo |
| `tiny.fa` | FASTA sequence | Reference validation demo |
| `tiny.fastq` | FASTQ reads | File type detection demo |
| `rnaseq_sample_R1.fastq` | Synthetic FASTQ R1 | nf-core/rnaseq samplesheet builder demo |
| `rnaseq_sample_R2.fastq` | Synthetic FASTQ R2 | nf-core/rnaseq samplesheet builder demo |
| `atac_sample_R1.fastq` | Synthetic FASTQ R1 | nf-core/atacseq samplesheet builder demo |
| `atac_sample_R2.fastq` | Synthetic FASTQ R2 | nf-core/atacseq samplesheet builder demo |

The synthetic FASTQ files contain 3 reads of 52 bp (4 lines each). They have no biological content - they exist solely to exercise FASTQ detection and samplesheet generation logic. Do not use them as real sequencing data.

## End-to-end demo

```bash
bash examples/run_tiny_demo.sh
```

Runs from the repository root. Produces `reports/demo/final_report.md` and per-tool subdirectories under `reports/demo/`.

**Steps:**

1. Environment check (`check_environment.py`) - verifies Python packages and optional CLI tools
2. File inspection (`inspect_file.py`) - detects file type, size, and assay
3. Reference validation (`reference_validator.py`) - checks GTF and blacklist BED
4. scRNA-seq QC (`scrna_qc_local.py --recommend-only`) - local fallback on tiny_counts.csv
5. ATAC-seq QC (`atac_qc_local.py`) - fragments + peaks + GTF
6. WGS/VCF QC (`wgs_vcf_qc_local.py`) - variant counts, Ti/Tv, het/hom
7. nf-core/rnaseq plan (`nfcore_launcher.py --dry-run`) - samplesheet, preflight, command (dry run only - Nextflow not required)
8. Report assembly (`report_builder.py`) - combines all QC sections

**Requirements:** Python 3.11+ with numpy, pandas, scipy, matplotlib, seaborn, h5py. Run from the repository root with the conda environment active.

**Note:** This demo uses local fallback tools on toy data. For production analysis:
- Use `single-cell-rna-qc@life-sciences` for scRNA QC
- Use `nextflow-development@life-sciences` for full pipeline execution
- The nf-core planning step (step 7) generates a plan and samplesheet without launching Nextflow

**scRNA-seq raw FASTQ-to-count is not demonstrated here.** This demo operates on a pre-computed count matrix. There is no Cell Ranger, STARsolo, or nf-core/scrnaseq step.

See [`docs/DEMO_OUTPUTS.md`](../docs/DEMO_OUTPUTS.md) for example output excerpts.
