# Implementation Status

This document tracks what is implemented locally, what official Claude Life Sciences tools should be preferred, what external tools are required for missing functionality, and what biological conclusions cannot be drawn when a feature is unavailable.

> **scRNA-seq positioning**: The official Anthropic Life Sciences skill `single-cell-rna-qc@life-sciences` is the **primary and default** scRNA QC path. `tools/scrna_qc_local.py` is a **limited fallback** for when the official skill is unavailable or the input format is incompatible. It is not a replacement and must not be used in preference to the official skill when that skill is available.

---

## Feature Matrix

| Feature | Implemented Locally? | Official Claude LS Tool | External Tool | Biological Caveat if Skipped |
|---------|---------------------|------------------------|---------------|------------------------------|
| File type detection | Yes - `tools/inspect_file.py` | None | None | Cannot route to correct QC workflow without file type detection |
| scRNA-seq QC (primary path) | **Official-first wrapper implemented** - use `single-cell-rna-qc@life-sciences` when available. Full production scRNA QC should always prefer the official skill. | `single-cell-rna-qc@life-sciences` **(primary)** | None | Cannot identify low-quality cells or doublet-enriched populations |
| scRNA-seq QC (local fallback only) | **Local fallback implemented for limited compatible inputs** - `tools/scrna_qc_local.py` handles CSV/TSV count matrices when official skill is unavailable. Not for `.h5ad` or 10x `.h5` in production. | `single-cell-rna-qc@life-sciences` (use instead when available) | None | Local fallback cannot replicate full QC from official skill; do not use as substitute when official skill is accessible |
| scRNA-seq QC plots (local fallback) | Fallback only - `tools/scrna_qc_local.py` (requires matplotlib, scanpy) | `single-cell-rna-qc@life-sciences` | None | Cannot visually assess cell quality distributions without running QC |
| Doublet detection | **No** | `single-cell-rna-qc@life-sciences` (may include) | `scrublet`, `DoubletFinder`, `scDblFinder` | Cannot distinguish doublets from large/activated cells reliably |
| Ambient RNA estimation | **No** | `single-cell-rna-qc@life-sciences` (may include) | `SoupX`, `DecontX` | Cannot distinguish ambient contamination from true low-level expression |
| Cell type annotation | **No** | None | `celltypist`, `scType`, manual marker review | Cannot validate QC decisions without cell type context |
| Batch integration / correction | **No (local)** | `scvi-tools@life-sciences` (preferred) | `scVI`, `Harmony`, `Scanorama` | Cannot assess whether batch effects confound QC metrics |
| ATAC fragment counting | Yes - `tools/atac_qc_local.py` | None | None | Cannot assess per-cell or per-sample yield |
| ATAC insert-size distribution | Yes - `tools/atac_qc_local.py` (requires numpy, matplotlib) | None | `deeptools bamPEFragmentSize` | Cannot assess nucleosomal banding or over/under-digestion |
| ATAC FRiP - local Python overlap | Yes - `tools/atac_qc_local.py` (requires peaks BED) | None | None | Approximate; use bedtools FRiP for large datasets |
| ATAC FRiP - bedtools validation | **Yes (v2.0) - runs when bedtools installed** | None | `bedtools intersect` | Cannot cross-validate signal-to-noise without bedtools |
| ATAC blacklist fraction | **Yes (v2.0) - runs when bedtools installed + blacklist BED provided** | `nextflow-development@life-sciences` (nf-core/atacseq) | `bedtools intersect` | Cannot identify artifacts from problematic genomic regions |
| ATAC chromosome naming mismatch detection | **Yes (v2.0) - always runs when multiple inputs provided** | None | None | Silent zero-overlap from chr1 vs 1 naming will corrupt all intersect-based metrics |
| ATAC TSS enrichment - command documented | **Partial (v2.0) - deeptools command documented when BAM+deeptools present** | `nextflow-development@life-sciences` (nf-core/atacseq) | `deeptools computeMatrix` | Cannot assess whether open chromatin is enriched at gene regulatory elements |
| ATAC peak calling | **Not implemented** | `nextflow-development@life-sciences` (nf-core/atacseq) | `MACS3`, `Genrich` | Cannot define accessibility regions without peaks |
| VCF variant counting (SNPs, indels, Ti/Tv) | Yes - `tools/wgs_vcf_qc_local.py` | None | `bcftools stats` | Cannot assess basic variant call quality |
| VCF bcftools stats cross-validation | **Yes (v2.0) - runs when bcftools installed, discrepancies reported as warnings** | None | `bcftools stats` | Cannot independently validate local parser counts without bcftools |
| VCF allele balance distribution | Partial - `tools/wgs_vcf_qc_local.py` | None | `bcftools stats`, `GATK VariantEval` | Cannot identify systematic sequencing bias |
| BAM samtools flagstat (structured) | **Yes (v2.0) - structured JSON output with command+version recorded** | None | `samtools flagstat` | Cannot assess mapping rate or duplicate rate without samtools |
| BAM samtools idxstats | **Yes (v2.0) - runs when index (.bai/.csi) exists** | None | `samtools idxstats` | Cannot assess per-chromosome read distribution without index |
| BAM samtools stats | **Yes (v2.0) - runs when samtools installed** | None | `samtools stats` | Cannot assess insert size distribution or error rate without samtools stats |
| BAM samtools depth (safe) | **Yes (v2.0) - runs on tiny files (<10 MB) or with --intervals** | None | `samtools depth` | Cannot assess per-base coverage depth for large files without intervals |
| BAM alignment summary (pysam) | Yes - `tools/wgs_vcf_qc_local.py` (requires pysam) | None | `samtools flagstat` | Cannot assess mapping rate, duplicate rate, or insert size |
| BAM coverage statistics | **Not implemented** | `nextflow-development@life-sciences` (nf-core/sarek) | `mosdepth`, `samtools depth` | Cannot assess whether coverage is sufficient or uniform; cannot identify low-coverage regions that may mask variants |
| BQSR / variant recalibration | **Not implemented** | `nextflow-development@life-sciences` (nf-core/sarek) | `GATK BaseRecalibrator` | Cannot assess systematic base quality errors from sequencing chemistry |
| Variant annotation (consequences) | **Not implemented** | None | `VEP`, `ANNOVAR`, `SnpEff` | Cannot assess predicted functional impact; no consequence-level filtering possible |
| Contamination estimation | **Not implemented** | None | `VerifyBamID`, `Conpair` | Cannot detect sample mix-up or inter-sample contamination |
| Reference file validation | Yes - `tools/reference_validator.py` | None | None | Cannot detect chromosome naming conflicts or missing indexes before analysis |
| Report assembly | Yes - `tools/report_builder.py` | None | None | Manual assembly required if tool unavailable |
| Biological interpretation | Enforced in all agents/skills | `pubmed@life-sciences` (for context) | None | Without biological interpretation, QC thresholds are applied blind |
| Full RNA-seq alignment + quantification | **Not implemented locally** - use nf-core launcher to plan | `nextflow-development@life-sciences` (nf-core/rnaseq) | `STAR`, `HISAT2`, `Salmon`, `RSEM` | Cannot produce gene-level expression counts from raw FASTQ |
| Full WGS/WES pipeline | **Not implemented locally** - use nf-core launcher to plan | `nextflow-development@life-sciences` (nf-core/sarek) | Full GATK4 stack | Cannot call variants reproducibly from raw FASTQ |
| Full ATAC-seq pipeline | **Not implemented locally** - use nf-core launcher to plan | `nextflow-development@life-sciences` (nf-core/atacseq) | `Trim Galore`, `BWA/Bowtie2`, `MACS3`, `deeptools` | Cannot call peaks or compute full QC from raw FASTQ |
| nf-core samplesheet generation (rnaseq) | **Yes (v0.3) - `tools/nfcore_launcher.py`** detects R1/R2 pairs, writes CSV with strandedness=auto, warns that strandedness must be confirmed | `nextflow-development@life-sciences` | None | Without a valid samplesheet, no nf-core pipeline can launch |
| nf-core samplesheet generation (sarek) | **Yes (v0.3) - conservative draft** with PATIENT_ID placeholders; tumor/normal status and sex cannot be inferred from filenames | `nextflow-development@life-sciences` | None | Clinical metadata (status, patient, sex) must be added manually before execution |
| nf-core samplesheet generation (atacseq) | **Yes (v0.3) - draft** with replicate=1 placeholders; experimental design and blacklist BED must be confirmed manually | `nextflow-development@life-sciences` | None | Wrong replicates or missing blacklist will corrupt peak calls |
| nf-core preflight checks | **Yes (v0.3)** - checks nextflow, docker/singularity/apptainer/conda, reference file existence | `nextflow-development@life-sciences` | None | Cannot assess whether local environment is ready to run any nf-core pipeline |
| nf-core command builder | **Yes (v0.3)** - constructs nextflow run command; writes to commands.sh | `nextflow-development@life-sciences` | None | Command will be incomplete without reference genome or samplesheet |
| nf-core execution with provenance | **Yes (v0.3)** - `--run` executes after preflight pass; dry-run default; JSON provenance always written | `nextflow-development@life-sciences` | Nextflow | Execution cannot be attempted without Nextflow installed |
| MultiQC output parsing | **Yes (v0.3)** - finds multiqc_report.html, multiqc_general_stats.txt, multiqc_data.json; reports what is present and what is missing | `nextflow-development@life-sciences` | None | Cannot summarize QC across samples without MultiQC output |

---

## Skipped Metric → Missing Biological Conclusion

When a metric is skipped, the following biological conclusions **cannot** be made:

### TSS enrichment (ATAC-seq)
- **Cannot conclude**: whether transposase insertion is enriched at transcription start sites
- **Cannot conclude**: whether the library has acceptable signal-to-noise for accessibility analysis
- **Cannot conclude**: whether cells/nuclei were of sufficient quality for regulatory inference
- **Enable with**: `nextflow-development@life-sciences` with `nf-core/atacseq`, or `deeptools computeMatrix` + local GTF

### Blacklist fraction (ATAC-seq)
- **Cannot conclude**: what fraction of peaks or fragments overlap problematic genomic regions
- **Cannot conclude**: whether apparent accessibility signals are artifactual (satellite repeats, rDNA, centromeres)
- **Enable with**: `bedtools intersect -a fragments.bed -b blacklist.bed` with an ENCODE blacklist BED

### Coverage statistics (WGS/WES)
- **Cannot conclude**: whether sequencing depth is sufficient for confident variant calling
- **Cannot conclude**: whether coverage is uniform across the target region (exome) or genome
- **Cannot conclude**: which regions are under-covered and therefore have unreliable variant calls
- **Enable with**: `mosdepth --quantize 0:5:10:30: output.prefix input.bam`

### Contamination estimation (WGS/WES)
- **Cannot conclude**: whether the sample is a pure single-donor sample
- **Cannot conclude**: whether apparent heterozygosity is biological or contamination
- **Cannot conclude**: whether tumor purity estimates are reliable
- **Enable with**: `VerifyBamID2` with population VCF panel

### Doublet detection (scRNA-seq)
- **Cannot conclude**: whether high-count, high-gene cells are true doublets or large/activated cells
- **Cannot conclude**: the true number of singlet cells in the dataset
- **Enable with**: `scrublet`, `DoubletFinder`, or `scDblFinder`

### Variant annotation (WGS/WES)
- **Cannot conclude**: predicted functional consequences (synonymous, missense, frameshift, splice)
- **Cannot conclude**: population allele frequency from gnomAD/dbSNP
- **Cannot conclude**: ClinVar classification (NOTE: ClinVar classifications are not clinical diagnoses)
- **Enable with**: `VEP` or `ANNOVAR` with appropriate annotation databases

---

## Filter Decision Framework

For every recommended filter, tools in this workspace must state:

| Question | Required Answer | JSON field |
|----------|----------------|------------|
| What metric triggered it? | Specific metric name and bound (low/high) | `metric`, `direction`, `n_flagged_low/high` |
| What technical artifact could explain it? | At least one specific artifact | `technical_artifact` |
| What biological state could explain it? | At least one specific biological scenario | `biological_signal` |
| Evidence FOR filtering | Specific features that support removing these cells | `evidence_for_filtering` |
| Evidence AGAINST filtering | Specific features that argue for keeping these cells | `evidence_against_filtering` |
| Recommended action | Filter / Flag / Stratify / Preserve | `recommended_action_{low\|high}` |
| Validation / sensitivity analysis | How to verify the decision is correct | `validation_note` |

**Recommended action definitions:**
- **Filter**: Remove from downstream analysis. Use only when artifact evidence is strong and biological signal is unlikely.
- **Flag**: Keep in analysis but mark for sensitivity testing. Use when artifact and biology are both plausible.
- **Stratify**: Analyze separately (e.g., split high-mito cells as a separate subpopulation). Use when biological signal is plausible and the population is large enough.
- **Preserve**: Do not remove. Use when biological signal is likely and the population is rare or scientifically relevant.

---

## Known Limitations

1. No doublet detection
2. No ambient RNA estimation
3. No coverage calculation (requires mosdepth or samtools depth)
4. No TSS enrichment (requires GTF parsing + deeptools or nf-core)
5. No blacklist fraction (requires bedtools)
6. No variant annotation (requires VEP/ANNOVAR)
7. No contamination estimation (requires VerifyBamID)
8. No local alignment - use `tools/nfcore_launcher.py --workflow rnaseq` to plan and `nextflow-development@life-sciences` to run
9. No local variant calling - use `tools/nfcore_launcher.py --workflow sarek` to plan
10. No local peak calling - use `tools/nfcore_launcher.py --workflow atacseq` to plan
11. R-based workflows (Seurat, Signac, DESeq2, edgeR) not included

---

## Recommended Tool Stack by Assay

### scRNA-seq
```
Claude Code orchestration
  → single-cell-rna-qc@life-sciences  (preferred QC)
  → tools/scrna_qc_local.py           (fallback QC)
  → scrublet / scDblFinder            (doublet detection)
  → SoupX / DecontX                   (ambient RNA)
  → scvi-tools@life-sciences          (integration/modeling)
  → biology-interpretation-reviewer   (biological review)
```

### ATAC-seq
```
Claude Code orchestration
  → tools/nfcore_launcher.py --workflow atacseq  (samplesheet + plan)
  → nextflow-development@life-sciences           (nf-core/atacseq execution)
  → tools/atac_qc_local.py                      (local complementary QC)
  → deeptools                                    (TSS enrichment, coverage)
  → bedtools                                     (blacklist fraction)
  → biology-interpretation-reviewer              (biological review)
```

### WGS/WES
```
Claude Code orchestration
  → tools/nfcore_launcher.py --workflow sarek    (samplesheet + plan)
  → nextflow-development@life-sciences           (nf-core/sarek execution)
  → tools/wgs_vcf_qc_local.py                   (local complementary QC)
  → mosdepth                                     (coverage)
  → VEP / ANNOVAR                                (variant annotation)
  → VerifyBamID                                  (contamination)
  → biology-interpretation-reviewer              (biological review)
  [No clinical interpretation - ever]
```
