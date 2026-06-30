# Agentic Genomics Workflow Framework

> Orchestrate reproducible NGS pipelines — inspect, plan, and execute genomics workflows with biological reasoning guardrails.

[![Tests](https://github.com/vfferreira/genomics-agent-workspace/actions/workflows/tests.yml/badge.svg)](https://github.com/vfferreira/genomics-agent-workspace/actions/workflows/tests.yml)

---

## What makes it agentic?

The FASTQ QC agent implements a real observation-decision-action loop. No LLM is called at runtime. All decisions are deterministic rules applied to parsed QC output.

**For FASTQ QC (`agent --workflow fastq-qc`):**

1. **Observe** — parses FastQC `.zip` files and MultiQC data files into structured per-sample results
2. **Decide** — applies deterministic rules:
   - Adapter Content WARN or FAIL: recommend trimming
   - Per base sequence quality FAIL: recommend trimming + sample review
   - Per sequence GC content FAIL: flag for biological review, never auto-filter
   - Overrepresented sequences WARN or FAIL: recommend contamination review
   - All tracked modules pass: record no trimming needed
3. **Act** — if `--execute --auto-trim` are both passed, runs trimming only when the decision engine has recommended it; original FASTQ files are never modified

All other subcommands (`run`, `plan`, `inspect`, `report`) are workflow runners: they build commands, execute them, validate outputs, and capture provenance. They do not parse outputs or make decisions.

---

## What this does

- Inspects input files and auto-detects data types (FASTQ, BAM, VCF, h5ad, BED, GTF, 10x directories)
- Infers the likely workflow from input files (`auto` mode)
- Generates samplesheets for nf-core pipelines (rnaseq, atacseq, sarek, ampliseq)
- Plans and executes FastQC/MultiQC, fastp/cutadapt trimming, samtools/bcftools QC, and Nextflow/nf-core pipelines when tools are installed and `--execute` is provided
- Validates that expected outputs exist after each step
- Captures provenance (command, exit code, stdout/stderr, timestamps) for every executed step
- Generates Markdown and JSON reports at every stage
- Makes deterministic agentic decisions for FASTQ QC: parses outputs, applies rules, recommends or runs trimming
- Delegates computation to FastQC, MultiQC, fastp, cutadapt, samtools, bcftools, mosdepth, Nextflow/nf-core

## What this does NOT do

- Does not reimplement FastQC, MultiQC, fastp, cutadapt, samtools, bcftools, mosdepth, Nextflow, nf-core, QIIME2, or DADA2
- Does not call an LLM at runtime — all decisions are deterministic rules
- Does not make agentic decisions for variant-qc, rnaseq, atacseq, or amplicon workflows yet
- Does not directly execute QIIME2 or DADA2 — amplicon workflow generates Nextflow/nf-core commands
- Does not make clinical claims — variant calls are never interpreted as diagnoses
- Does not load large genomic files into memory — file operations are delegated to CLI tools
- Does not download large reference files — references must be provided locally
- Does not silently fail — every skipped step and failed command is recorded

---

## Supported Workflows

| Workflow | CLI name | Plan | Execute | Agentic decisions | Tool backend | Status |
|----------|----------|------|---------|-------------------|--------------|--------|
| FASTQ QC | `fastq-qc` | Yes | Yes | Yes (FASTQ agent) | FastQC + MultiQC | Implemented |
| Trimming | `fastq-qc --trim` | Yes | Yes | Via FASTQ agent `--auto-trim` | fastp / cutadapt | Implemented |
| Variant QC | `variant-qc` | Yes | Yes | No | samtools + bcftools + mosdepth | Execution-capable |
| WGS full pipeline | `variant-qc` | Yes | Cautious | No | nf-core/sarek (Nextflow) | Planning / cautious execution |
| RNA-seq | `rnaseq` | Yes | Yes (needs Nextflow) | No | nf-core/rnaseq (Nextflow) | Execution-capable |
| ATAC-seq | `atacseq` | Yes | Yes (needs Nextflow) | No | nf-core/atacseq (Nextflow) | Execution-capable |
| Amplicon / Microbiome | `amplicon` | Yes | Yes (needs Nextflow) | No | nf-core/ampliseq (Nextflow) | Execution-capable |
| Auto-detect | `auto` | Yes | Yes | No | (inferred from input) | Implemented |
| scRNA-seq | — | No | No | No | External skill preferred | Not directly implemented |

**Notes:**
- "Execution-capable" means the framework builds correct commands and runs them via subprocess with provenance capture. It does not interpret the results.
- Nextflow workflows require Nextflow installed and a container runtime (Docker, Singularity) or conda.
- nf-core/sarek samplesheets contain placeholders that require human review before execution.
- scRNA-seq is out of scope for the Python package; use `single-cell-rna-qc@life-sciences` or `tools/scrna_qc_local.py` as a fallback.

---

## Installation

```bash
git clone <this-repo>
cd genomics-workflow-agent

# Minimal install (core package only)
pip install -e .

# Full install with scientific Python stack
pip install -e ".[full]"

# Or use the conda environment
conda env create -f env/environment.yml
conda activate genomics-agent
```

External tools must be installed separately:

```bash
conda install -c bioconda fastqc multiqc fastp cutadapt samtools bcftools bedtools nextflow
```

---

## Quick Start

```bash
# 1. Inspect input data and detect file types
python -m genomics_workflow_agent inspect --input data/

# 2. Build a plan (always dry-run)
python -m genomics_workflow_agent plan --input data/ --workflow fastq-qc --out results/

# 3. Execute (launches external tools when installed)
python -m genomics_workflow_agent run --input data/ --workflow fastq-qc --out results/ --execute

# 4. Agentic FASTQ QC: observe, decide, optionally act
python -m genomics_workflow_agent agent --input data/ --workflow fastq-qc --out results_agent/

# 5. Aggregate results into a final report
python -m genomics_workflow_agent report --results results/
```

---

## Agentic FASTQ QC

```bash
# Dry-run: inspect files, build plan, write reports (no external tools called)
python -m genomics_workflow_agent agent \
  --input data/fastqs/ \
  --workflow fastq-qc \
  --out results_agent/

# Execute: run FastQC/MultiQC, parse outputs, make decisions, write reports
python -m genomics_workflow_agent agent \
  --input data/fastqs/ \
  --workflow fastq-qc \
  --out results_agent/ \
  --execute

# Auto-trim: run QC, parse results, run trimming only if evidence recommends it
python -m genomics_workflow_agent agent \
  --input data/fastqs/ \
  --workflow fastq-qc \
  --out results_agent/ \
  --execute \
  --auto-trim \
  --trim-tool fastp
```

Output files:

- `results_agent/agent_report.json` — structured state: observations, decisions, recommended actions
- `results_agent/agent_report.md` — human-readable report with reasoning
- `results_agent/fastqc/` — FastQC HTML and zip outputs (if `--execute`)
- `results_agent/multiqc/` — MultiQC report (if `--execute`)
- `results_agent/trimmed/` — trimmed reads (only if `--auto-trim` and trimming was recommended)
- `results_agent/provenance/` — per-step provenance JSON files

---

## Other Example CLI Commands

### FASTQ QC (run subcommand)

```bash
# Plan only
python -m genomics_workflow_agent plan --input data/fastqs/ --workflow fastq-qc --out results/

# Execute FastQC + MultiQC
python -m genomics_workflow_agent run --input data/fastqs/ --workflow fastq-qc --out results/ --execute

# Execute with trimming
python -m genomics_workflow_agent run --input data/fastqs/ --workflow fastq-qc --out results/ --execute --trim fastp
```

### RNA-seq

```bash
# Plan nf-core/rnaseq (dry-run)
python -m genomics_workflow_agent plan \
  --input data/rnaseq/ --workflow rnaseq \
  --genome GRCh38 \
  --out results/rnaseq/

# Execute (requires Nextflow + Docker/Singularity/conda)
python -m genomics_workflow_agent run \
  --input data/rnaseq/ --workflow rnaseq \
  --genome GRCh38 \
  --out results/rnaseq/ \
  --execute
```

### ATAC-seq

```bash
python -m genomics_workflow_agent plan \
  --input data/atac/ --workflow atacseq \
  --genome GRCh38 \
  --blacklist /refs/hg38-blacklist.v2.bed \
  --out results/atacseq/
```

### Amplicon / Microbiome

```bash
# 16S V4 with SILVA taxonomy
python -m genomics_workflow_agent plan \
  --input data/16s/ --workflow amplicon \
  --primer-fw GTGYCAGCMGCCGCGGTAA \
  --primer-rv GGACTACNVGGGTWTCTAAT \
  --taxonomy-db SILVA \
  --out results/amplicon/
```

### WGS / Variant QC

```bash
python -m genomics_workflow_agent plan \
  --input data/wgs/ --workflow variant-qc \
  --out results/variantqc/

python -m genomics_workflow_agent run \
  --input data/wgs/ --workflow variant-qc \
  --out results/variantqc/ \
  --execute
```

---

## Package Structure

```
genomics_workflow_agent/
├── __init__.py
├── __main__.py
├── cli.py                    # CLI entry point (inspect/plan/run/report/agent)
├── agent/
│   ├── state.py              # AgentState, Observation, Decision, RecommendedAction
│   ├── decision_engine.py    # Deterministic FASTQ QC rules
│   └── fastq_agent.py        # Observation-decision-action orchestrator + report writers
├── parsers/
│   ├── fastqc.py             # FastQC zip/txt parser -> structured per-sample results
│   └── multiqc.py            # MultiQC TSV/JSON/HTML parser
├── inspect/
│   └── inspector.py          # File type detection, directory summary, workflow guess
├── workflows/
│   ├── planner.py            # Auto-detect + dispatch
│   ├── fastq_qc.py           # FastQC + MultiQC + fastp/cutadapt
│   ├── rnaseq.py             # nf-core/rnaseq planning and execution
│   ├── atacseq.py            # nf-core/atacseq planning and execution
│   ├── amplicon.py           # nf-core/ampliseq planning and execution
│   └── variant_qc.py         # samtools/bcftools/mosdepth + nf-core/sarek
├── tools/
│   ├── runner.py             # Safe subprocess runner with provenance capture
│   ├── nextflow.py           # Nextflow command builder and runner
│   ├── versions.py           # Tool availability checker
│   ├── files.py              # File discovery, FASTQ pair detection
│   └── samplesheets.py       # Samplesheet generators (rnaseq/atacseq/sarek/amplicon)
├── reports/
│   ├── markdown.py           # Markdown report generator
│   └── json_report.py        # JSON report with metadata
├── config/
│   └── schema.py             # Workflow defaults
└── safety/
    └── guardrails.py         # Clinical disclaimer, large-file protection
```

Legacy CLI tools (fallbacks, not part of the agentic layer):

```
tools/
├── scrna_qc_local.py         # scRNA-seq QC (fallback-only — prefer single-cell-rna-qc@life-sciences)
├── atac_qc_local.py          # ATAC-seq QC (fallback)
├── wgs_vcf_qc_local.py       # WGS/VCF QC (fallback)
├── nfcore_launcher.py        # nf-core workflow launcher
└── ...
```

---

## Amplicon Taxonomy Databases

The amplicon workflow supports planning for these taxonomy databases. Database files must be downloaded separately — the framework never downloads large references automatically.

| Database | Target |
|----------|--------|
| SILVA | 16S/18S/23S/28S rRNA |
| GTDB | Prokaryotic phylogeny |
| UNITE | Fungal ITS |
| Greengenes2 | 16S + ITS |
| custom | User-provided |

---

## Execution Model

```
inspect  ->  plan (always dry-run)  ->  run [--execute]  ->  report
                                        agent [--execute] [--auto-trim]
```

- **inspect**: Detects file types, reports sizes, guesses workflow
- **plan**: Generates samplesheets, builds commands, checks tool availability. Always dry-run.
- **run**: Executes steps with provenance capture. Dry-run by default. `--execute` required to launch tools.
- **agent**: FASTQ QC only. Dry-run by default. With `--execute`: runs FastQC/MultiQC, parses outputs, makes decisions. With `--execute --auto-trim`: also runs trimming if evidence recommends it.
- **report**: Aggregates JSON outputs from results directory

Default behavior is safe: no external tools are launched without `--execute`.

---

## Provenance

Every executed command is recorded in `results/provenance/`:

```json
{
  "label": "fastqc",
  "command": ["fastqc", "--outdir", "results/fastqc", "--threads", "4", "sample_R1.fastq.gz"],
  "command_str": "fastqc --outdir results/fastqc --threads 4 sample_R1.fastq.gz",
  "cwd": "/data/project",
  "dry_run": false,
  "started_at": "2026-06-30T10:00:00+00:00",
  "ended_at": "2026-06-30T10:01:32+00:00",
  "return_code": 0,
  "stdout_snippet": "...",
  "stderr_snippet": "...",
  "runtime_s": 92.3,
  "executed": true,
  "output_validation": {"all_present": true, "missing": []}
}
```

---

## Claude Agents and Skills

Seven specialized Claude Code agents are available in `.claude/agents/` for delegation:

| Agent | Task |
|-------|------|
| `genomics-file-inspector` | File inspection and routing |
| `scrna-qc-specialist` | Single-cell RNA-seq QC |
| `atac-qc-specialist` | ATAC-seq QC |
| `wgs-qc-variant-specialist` | WGS/WES BAM/VCF QC |
| `nextflow-pipeline-specialist` | nf-core pipeline orchestration |
| `biology-interpretation-reviewer` | Artifact vs. biology review |
| `single-cell-modeling-specialist` | scVI/scANVI/totalVI modeling |

These are Claude Code agent definitions, not Python code. They are used when running Claude Code interactively.

When official Anthropic Life Sciences skills are available, agents prefer them over local fallbacks:

| Official Skill | Use Case |
|---------------|----------|
| `single-cell-rna-qc@life-sciences` | scRNA-seq QC |
| `scvi-tools@life-sciences` | Downstream single-cell modeling |
| `nextflow-development@life-sciences` | nf-core orchestration |
| `10x-genomics@life-sciences` | 10x Cloud data access |
| `pubmed@life-sciences` | Biological background and marker interpretation |

---

## Safety and Limitations

- **No clinical claims** — enforced by tests. Variant QC and FASTQ agent outputs explicitly disclaim clinical interpretation.
- **No LLM at runtime** — all decisions are deterministic rules applied to parsed QC flags.
- **No large-file loading** — files above 50 MB are not read into memory. All large-file operations are delegated to CLI tools.
- **No silent failures** — every failed command is recorded in provenance with exit code and stderr.
- **No automatic reference downloads** — reference files must be provided locally.
- **Samplesheets require human review** — auto-generated samplesheets for nf-core pipelines contain placeholders that must be corrected before execution.
- **Dry-run by default** — expensive workflows are never launched without an explicit `--execute` flag.

---

## Tests

```bash
# Run the full test suite
python -m pytest tests/ -q

# Agent and execution tests only
python -m pytest tests/test_agent.py tests/test_execution.py -q

# With coverage
python -m pytest tests/ --cov=genomics_workflow_agent --cov=tools
```

The test suite covers: FastQC zip/txt parsing, MultiQC output parsing, decision engine rules (no-trim on all-pass, trim on adapter warn, review on quality fail, no auto-filter on GC fail), agent dry-run, agent execution with mocked tools, auto-trim gating, report writing, CLI agent subcommand, workflow planning, execution with mocked subprocess, provenance capture, output validation, failed command handling, clinical-claims guardrail.

Integration tests that require external tools (FastQC, MultiQC, fastp, samtools, bcftools, nextflow) are automatically skipped when the tool is not in PATH.

---

## Running the Demo

```bash
bash examples/run_tiny_demo.sh
```

This runs the end-to-end demo on toy data in `examples/`.

---

## Privacy

Human genomic data is treated as privacy-sensitive by default. Do not upload or transmit private genomic data without appropriate authorization. The framework never transmits data to external services.

---

## Future Work

- Agentic decision layer for variant-qc, rnaseq, atacseq, and amplicon workflows
- Per-sample MultiQC TSV parsing connected to the FASTQ agent decision engine
- Direct QIIME2/DADA2 execution (currently generates Nextflow/nf-core commands only)
- HTML report export
- Adaptive Nextflow pipeline decisions based on parsed pipeline outputs
- Benchmark datasets and integration test fixtures for CI
- R-based downstream workflow documentation (DESeq2, phyloseq)
