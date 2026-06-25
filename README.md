# Genomics Agent Workspace

> A reproducible Claude Code workspace for genomics QC orchestration with biological reasoning guardrails.

[![Tests](https://github.com/vfferreira/genomics-agent-workspace/actions/workflows/tests.yml/badge.svg)](https://github.com/vfferreira/genomics-agent-workspace/actions/workflows/tests.yml)

---

## Why this project exists

Production genomics pipelines are mature and well-served by nf-core, Scanpy, and Anthropic's Life Sciences skill suite. This workspace solves a different problem: **how to integrate Claude Code safely and reproducibly into those pipelines**.

The local tools here do not replace samtools, bcftools, bedtools, Scanpy, or nf-core. They handle the surrounding concerns that pipelines leave to the analyst:

- Detecting file types and routing to the correct workflow
- Running pre-flight environment checks with graceful degradation
- Parsing external tool output into structured, provenance-tracked JSON
- Validating reference files and detecting chromosome naming mismatches before analysis
- Enforcing biological reasoning (artifact vs. biology) on every QC decision
- Assembling timestamped, reproducible reports

When an official Anthropic Life Sciences skill is available for a task, the agents in this workspace use it. Local tools are fallbacks.

---

## Portfolio review

A reviewer should start here, in this order:

| What | Where |
|------|-------|
| Architecture and design rationale | [`docs/OVERVIEW.md`](docs/OVERVIEW.md) |
| End-to-end demo | [`examples/run_tiny_demo.sh`](examples/run_tiny_demo.sh) |
| Pre-flight environment check | [`tools/check_environment.py`](tools/check_environment.py) |
| WGS/VCF QC with samtools/bcftools | [`tools/wgs_vcf_qc_local.py`](tools/wgs_vcf_qc_local.py) |
| ATAC-seq QC with bedtools | [`tools/atac_qc_local.py`](tools/atac_qc_local.py) |
| Agent constraint definitions | [`.claude/agents/`](.claude/agents/) |
| Parser tests without external tools | [`tests/`](tests/) |
| CI configuration | [`.github/workflows/tests.yml`](.github/workflows/tests.yml) |

---

## Design constraints

**Claude Code should never load large genomic files into model context.** All large-file operations go through:

- Official Anthropic Life Sciences skills (preferred)
- Local CLI tools in `tools/` (fallback and wrapper layer)
- Nextflow/nf-core pipelines for full production workflows

Every QC decision must carry a biological interpretation — not just a threshold.

---

## Official Claude Life Sciences skills

Install the marketplace first:
```
/plugin marketplace add anthropics/life-sciences
```

Install individual skills:
```
/plugin install single-cell-rna-qc@life-sciences
/plugin install scvi-tools@life-sciences
/plugin install nextflow-development@life-sciences
/plugin install 10x-genomics@life-sciences
/plugin install pubmed@life-sciences
```

| Skill | Use when |
|-------|----------|
| `single-cell-rna-qc@life-sciences` | scRNA QC on `.h5ad` or 10x `.h5` — **primary path** |
| `scvi-tools@life-sciences` | Downstream modeling: integration, batch correction |
| `nextflow-development@life-sciences` | Full nf-core pipeline orchestration |
| `10x-genomics@life-sciences` | 10x Genomics Cloud data |
| `pubmed@life-sciences` | Biological context and marker interpretation |

> Local tools in `tools/` are fallbacks, not replacements for official skills.

---

## Repository structure

```
.claude/agents/     — Seven specialized agent configurations
.claude/skills/     — Six skill wrappers
.github/workflows/  — CI: pytest on Ubuntu + Windows, Python 3.11/3.12
CLAUDE.md           — Guardrails for Claude Code
docs/               — Implementation status, portfolio overview, roadmap
env/environment.yml — Conda environment spec
examples/           — Toy data + end-to-end demo script
references/README.md — Reference file guidance (files not committed)
templates/          — QC report template
tests/              — 279 tests (fixture-based; run without external tools)
tools/
  check_environment.py      — Pre-flight check
  inspect_file.py           — File type detection
  scrna_qc_local.py         — scRNA-seq QC fallback (CSV/TSV)
  atac_qc_local.py          — ATAC-seq QC (v2.0: bedtools, chrom-mismatch detection)
  wgs_vcf_qc_local.py       — WGS/VCF QC (v2.0: samtools, bcftools cross-validation)
  nfcore_launcher.py        — nf-core samplesheet builder, preflight, and safe launcher
  reference_validator.py    — Reference file validation
  report_builder.py         — Report assembly
```

---

## Installation

```bash
conda env create -f env/environment.yml
conda activate genomics-agent
```

Optional external tools (install separately; missing tools are warnings, not failures):

```bash
# Structured BAM QC and VCF cross-validation
conda install -c bioconda samtools bcftools

# ATAC-seq blacklist overlap and FRiP validation
conda install -c bioconda bedtools tabix

# QC and coverage tools
conda install -c bioconda fastqc multiqc fastp mosdepth

# ATAC-seq tools
conda install -c bioconda macs3 deeptools

# Pipeline orchestration
conda install -c bioconda nextflow
```

---

## Running the demo

```bash
bash examples/run_tiny_demo.sh
```

Runs all eight steps on toy data in `examples/` and writes `reports/demo/final_report.md`. See [`examples/README.md`](examples/README.md) for details.

---

## Running the tools individually

### Environment check
```bash
python tools/check_environment.py --output-dir reports/env_check
```

### File inspection
```bash
python tools/inspect_file.py --input examples/tiny_counts.csv --output-dir reports/inspect
```

### ATAC-seq QC
```bash
python tools/atac_qc_local.py \
  --fragments examples/tiny_fragments.tsv \
  --peaks examples/tiny_peaks.bed \
  --gtf examples/tiny.gtf \
  --output-dir reports/atac_tiny
```

### WGS/VCF QC
```bash
python tools/wgs_vcf_qc_local.py \
  --vcf examples/tiny.vcf \
  --output-dir reports/wgs_tiny
```

### scRNA-seq QC (fallback only)
```bash
# Prefer single-cell-rna-qc@life-sciences for .h5ad or 10x .h5 inputs.
# Use this only for CSV/TSV count matrices when the official skill is unavailable.
python tools/scrna_qc_local.py \
  --input examples/tiny_counts.csv \
  --species human --tissue unknown \
  --recommend-only \
  --output-dir reports/scrna_tiny
```

### Reference validation
```bash
python tools/reference_validator.py --gtf examples/tiny.gtf --output-dir reports/ref_check
```

---

## End-to-end sequencing support

`tools/nfcore_launcher.py` generates nf-core samplesheets, validates local requirements, and produces a Nextflow command without executing it by default:

```bash
# Plan only (dry-run default — no Nextflow required)
python tools/nfcore_launcher.py \
  --workflow rnaseq \
  --input-dir examples \
  --genome GRCh38 \
  --output-dir reports/nfcore_rnaseq \
  --dry-run

# Execute (only after reviewing the plan and resolving blockers)
python tools/nfcore_launcher.py \
  --workflow rnaseq \
  --input-dir /path/to/fastqs \
  --genome GRCh38 \
  --output-dir reports/nfcore_rnaseq \
  --run
```

Supported workflows: `rnaseq`, `sarek`, `atacseq`.

> **scRNA-seq raw FASTQ-to-count generation is not implemented in v0.3.** There is no nf-core/scrnaseq, STARsolo, or Cell Ranger integration. For scRNA-seq from raw FASTQ, use `nextflow-development@life-sciences` with nf-core/scrnaseq or run Cell Ranger/STARsolo directly. The local `scrna_qc_local.py` tool operates on count matrices only.

Every plan output includes a mandatory biological caveats section. Successful pipeline completion is never presented as biological or clinical validation. For sarek, patient ID, tumor/normal status, and sex are output as explicit placeholders — they cannot be inferred from filenames and must be set manually before running.

For production execution, prefer `nextflow-development@life-sciences`. The local launcher is a planning, provenance, and samplesheet layer only.

---

## Using Claude Code agents

```
Use the genomics-file-inspector agent to inspect the examples/ directory.
```

```
Use the atac-qc-wrapper skill on examples/tiny_fragments.tsv
with peaks from examples/tiny_peaks.bed.
```

```
Use the wgs-qc-wrapper skill on examples/tiny.vcf. Do not make clinical claims.
```

```
Use the biology-interpretation-reviewer agent on the outputs in reports/demo/.
Species: human. Tissue: PBMC.
```

```
Use the nextflow-pipeline-specialist agent to help me run nf-core/rnaseq.
Check if Nextflow and Docker are available first.
```

---

## Biological interpretation guardrails

Every QC report includes a structured biological interpretation section. For each warning, tools require:

1. What was observed?
2. What technical artifact could explain it?
3. What biological state could also explain it?
4. What metadata would help distinguish artifact from biology?
5. What validation should be done?
6. Should data be **filtered**, **flagged**, **stratified**, or **preserved**?

| Observation | Technical explanation | Biological explanation |
|-------------|----------------------|----------------------|
| High mito% | Dying cells, lysed membranes | Metabolic tissue, hypoxia, tumor stress, cardiac/muscle |
| High UMI counts | Doublets | Large cells, activated B/T cells, plasma cells |
| Low genes/cell | Empty droplets, dead cells | Mature RBCs, platelets, sparse mature cell types |
| Low FRiP (ATAC) | Poor transposition efficiency | Global chromatin remodeling |
| Low Ti/Tv | Calling artifact | Specific mutational signatures |
| Allele imbalance | Strand bias, PCR artifact | Mosaic variant, CNV, tumor purity effect |

High mitochondrial percentage, unusual variant distributions, low FRiP, and apparent outliers are not automatically artifacts. They require biological context.

---

## Reference files

Reference files are not included in this repository. See [`references/README.md`](references/README.md) for guidance on required files, genome builds, and chromosome naming conventions.

Tools never auto-download large references. Missing references produce a documented skip entry, not an error.

---

## Tests

```bash
pytest tests/ -v
```

Tests run without samtools, bcftools, or bedtools. Parser tests use fixture files in `tests/fixtures/`. Degradation tests verify that missing tools produce structured skip entries with `missing_biological_conclusion` fields.

| Test file | Coverage |
|-----------|---------|
| `test_inspect_file.py` | File types detected; JSON/Markdown always written |
| `test_check_environment.py` | Missing optional tools are warnings, not failures |
| `test_scrna_qc_local.py` | CSV input; MAD thresholds; biological filter notes |
| `test_atac_qc_local.py` | Fragments; FRiP; skipped metrics; degradation |
| `test_wgs_vcf_qc_local.py` | VCF parsed; Ti/Tv; no clinical claims in output |
| `test_reference_validator.py` | GTF; missing file; chromosome style detected |
| `test_report_builder.py` | All sections present; artifact table present |
| `test_claude_config.py` | Agent configs; official-skill-first routing verified |
| `test_nfcore_launcher.py` | Samplesheet builders; preflight; biological caveats; no-clinical-claims; no Nextflow required |
| `test_external_tools.py` | Parser functions; degradation; biological caveats |

---

## Privacy

Human genomic data is treated as privacy-sensitive by default. Local tools run entirely on-premises. Do not pass private genomic data to external services without explicit authorization. See [`SECURITY.md`](SECURITY.md).

---

## Known limitations

- Doublet detection requires scrublet, DoubletFinder, or scDblFinder
- Ambient RNA estimation requires SoupX or DecontX
- TSS enrichment requires deeptools + BAM (command documented in QC output when both present)
- Full coverage statistics require mosdepth or samtools depth with a full BAM
- Variant annotation requires VEP or ANNOVAR
- R-based workflows (Seurat, Signac, DESeq2) not included
- CRAM decoding requires a matching reference FASTA

See [`docs/ROADMAP.md`](docs/ROADMAP.md) and [`docs/IMPLEMENTATION_STATUS.md`](docs/IMPLEMENTATION_STATUS.md).
