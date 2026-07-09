# Agentic Genomics Workflow Framework

> Orchestrate reproducible NGS pipelines - inspect, plan, and execute genomics workflows with biological reasoning guardrails.

## What makes it agentic?

Two workflows have a real observation-decision-action loop. No LLM is called at runtime. All decisions are deterministic rules applied to parsed QC output.

**FASTQ QC (`agent --workflow fastq-qc`):**

1. **Observe** - parses FastQC `.zip` files and MultiQC data files into structured per-sample results
2. **Decide** - applies deterministic rules:
   - Adapter Content WARN or FAIL: recommend trimming
   - Per base sequence quality FAIL: recommend trimming + sample review
   - Per sequence GC content FAIL: flag for biological review, never auto-filter
   - Overrepresented sequences WARN or FAIL: recommend contamination review
   - All tracked modules pass: record no trimming needed
3. **Act** - if `--execute --auto-trim` are both passed, runs trimming only when the decision engine has recommended it; original FASTQ files are never modified

**Variant QC (`agent --workflow variant-qc`):**

1. **Observe** - parses samtools flagstat/idxstats/stats, bcftools stats, and mosdepth summary outputs into structured per-sample results
2. **Decide** - applies deterministic rules:
   - Mapped reads below 80%: recommend alignment investigation
   - Zero VCF records: recommend variant calling review
   - Many zero-read contigs: flag for reference compatibility review
   - Low coverage (below 10x): warn; do not make diagnostic claims
   - No issues found: record accept decision
3. **Recommend** - never makes clinical claims; never uses terms like pathogenic, benign, or diagnostic value

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
- Makes deterministic agentic decisions for FASTQ QC and Variant QC: parses outputs, applies rules, recommends or runs actions
- Exposes a stable Python API (`from genomics_workflow_agent import inspect_inputs, plan_workflow, run_workflow, run_fastq_qc_agent, run_variant_qc_agent`) for use by other Python projects
- Delegates computation to FastQC, MultiQC, fastp, cutadapt, samtools, bcftools, mosdepth, Nextflow/nf-core

## What this does NOT do

- Does not reimplement FastQC, MultiQC, fastp, cutadapt, samtools, bcftools, mosdepth, Nextflow, nf-core, QIIME2, or DADA2
- Does not call an LLM at runtime - all decisions are deterministic rules
- Does not make agentic decisions for rnaseq, atacseq, or amplicon workflows yet
- Does not directly execute QIIME2 or DADA2 - amplicon workflow generates Nextflow/nf-core commands
- Does not make clinical claims - variant calls are never interpreted as diagnoses
- Does not load large genomic files into memory - file operations are delegated to CLI tools
- Does not download large reference files - references must be provided locally
- Does not silently fail - every skipped step and failed command is recorded

---

## Supported Workflows

| Workflow | CLI name | Plan | Execute | Agentic decisions | Tool backend | Status |
|----------|----------|------|---------|-------------------|--------------|--------|
| FASTQ QC | `fastq-qc` | Yes | Yes | Yes (FASTQ agent) | FastQC + MultiQC | Implemented + Agentic |
| Trimming | `fastq-qc --trim` | Yes | Yes | Via FASTQ agent `--auto-trim` | fastp / cutadapt | Implemented + Agentic |
| Variant QC | `variant-qc` | Yes | Yes | Yes (Variant agent) | samtools + bcftools + mosdepth | Implemented + Agentic |
| WGS full pipeline | `variant-qc` | Yes | Cautious | No | nf-core/sarek (Nextflow) | Planning / cautious execution |
| RNA-seq | `rnaseq` | Yes | Yes (needs Nextflow) | No | nf-core/rnaseq (Nextflow) | Execution-capable |
| ATAC-seq | `atacseq` | Yes | Yes (needs Nextflow) | No | nf-core/atacseq (Nextflow) | Execution-capable |
| Amplicon / Microbiome | `amplicon` | Yes | Yes (needs Nextflow) | No | nf-core/ampliseq (Nextflow) | Execution-capable |
| Auto-detect | `auto` | Yes | Yes | No | (inferred from input) | Implemented |
| scRNA-seq | - | No | No | No | External skill preferred | Not directly implemented |

**Notes:**
- "Execution-capable" means the framework builds correct commands and runs them via subprocess with provenance capture. It does not interpret the results.
- Nextflow workflows require Nextflow installed and a container runtime (Docker, Singularity) or conda.
- nf-core/sarek samplesheets contain placeholders that require human review before execution.
- scRNA-seq is out of scope for the Python package; use `single-cell-rna-qc@life-sciences` or `tools/scrna_qc_local.py` as a fallback.

---

## Using as a Python Library

```python
from genomics_workflow_agent import (
    inspect_inputs,
    plan_workflow,
    run_workflow,
    run_fastq_qc_agent,
    run_variant_qc_agent,
    write_report,
)

# Inspect what is in a directory
result = inspect_inputs("data/fastqs/")
print(result["summary"]["workflow_guess"])  # e.g. "fastq-qc"

# Build a dry-run plan (no tools executed)
plan = plan_workflow("data/fastqs/", workflow="fastq-qc", outdir="results/")
print(plan["summary"]["steps"])

# Run the FASTQ QC agent (dry-run by default)
agent_result = run_fastq_qc_agent("data/fastqs/", outdir="results_agent/")
for obs in agent_result["observations"]:
    print(obs["sample"], obs["category"], obs["status"])

# Run the variant QC agent (dry-run by default)
variant_result = run_variant_qc_agent("data/wgs/", outdir="results_variant_agent/")
for dec in variant_result["decisions"]:
    print(dec["action"], dec["decision_type"])
```

All public functions return a plain JSON-serializable dict. See `genomics_workflow_agent/schemas/` for field contracts. See `examples/api_usage/` for runnable examples.

**Agentic workflows (v0.3.0):**
- `fastq-qc` - full observe-decide-act loop (FastQC/MultiQC parsing, trimming decisions)
- `variant-qc` - observe-decide-recommend loop (samtools/bcftools/mosdepth parsing, alignment and VCF checks)
- `rnaseq`, `atacseq`, `amplicon` - execution-capable runners, not agentic interpreters yet

---

## Running with Docker

Docker provides a reproducible bioinformatics runtime with all tools pre-installed.
**Claude Code runs on the host**; Docker provides the tool environment.

```bash
# Build the image (first time or after Dockerfile changes)
docker compose build genomics-agent

# CLI help
docker compose run --rm genomics-agent python -m genomics_workflow_agent --help

# API import smoke test
docker compose run --rm genomics-agent python -c \
  "from genomics_workflow_agent import inspect_inputs, plan_workflow, run_workflow, run_fastq_qc_agent, run_variant_qc_agent; print('api ok')"

# Run tests
docker compose run --rm genomics-agent python -m pytest tests/ -q

# FASTQ QC agent - dry-run (no external tools called)
docker compose run --rm genomics-agent python -m genomics_workflow_agent agent \
  --input examples/ --workflow fastq-qc --out results_agent_smoke/

# Variant QC agent - dry-run (no external tools called)
docker compose run --rm genomics-agent python -m genomics_workflow_agent agent \
  --input examples/ --workflow variant-qc --out results_variant_agent_smoke/

# FASTQ QC agent - execute with real tools (FastQC/MultiQC run inside container)
# Put data under ./data/ first, then:
docker compose run --rm genomics-agent python -m genomics_workflow_agent agent \
  --input data/ --workflow fastq-qc --out results/fastq_agent --execute

# Variant QC agent - execute with real BAM/VCF files
docker compose run --rm genomics-agent python -m genomics_workflow_agent agent \
  --input data/ --workflow variant-qc --out results/variant_agent --execute
```

Helper scripts (run from the host):

```bash
# Check all tool versions inside the container
bash scripts/docker_check_tools.sh
# Or on Windows PowerShell:
.\scripts\docker_check_tools.ps1

# Run tests
bash scripts/docker_test.sh          # Linux/macOS
.\scripts\docker_test.ps1            # Windows

# Run demos
bash scripts/docker_demo_fastq_agent.sh
.\scripts\docker_demo_fastq_agent.ps1
```

**Nextflow / nf-core note:**
Nextflow is installed in the image. Running nf-core pipelines from inside Docker
requires extra setup depending on the execution profile:

- `docker` profile needs the host Docker socket mounted (`-v /var/run/docker.sock:/var/run/docker.sock`)
- `singularity` profile is preferred on HPC systems
- `conda` profile works but can be slower inside the container
- Do not run full nf-core/rnaseq, nf-core/atacseq, nf-core/ampliseq, or nf-core/sarek in demos; they require references and configuration
- Full pipeline execution requires explicit user approval and a configured reference genome

---

## Claude Code usage

`.claude/agents/` contains Claude Code agent definitions for specialized tasks.
These are orchestration definitions - they are not Python runtime code.

**Recommended setup:**
- Claude Code runs on the **host**
- Docker provides the reproducible bioinformatics tool runtime
- Claude Code issues `docker compose run --rm genomics-agent ...` for real tool execution
- Dry-run commands (`python -m genomics_workflow_agent ... ` without `--execute`) can run on the host directly

**Security:**
- Do not put Claude API keys, tokens, or credentials into Dockerfiles or docker-compose.yml
- Do not commit `.env` files
- Docker is for tool reproducibility, not for storing secrets

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

- `results_agent/agent_report.json` - structured state: observations, decisions, recommended actions
- `results_agent/agent_report.md` - human-readable report with reasoning
- `results_agent/fastqc/` - FastQC HTML and zip outputs (if `--execute`)
- `results_agent/multiqc/` - MultiQC report (if `--execute`)
- `results_agent/trimmed/` - trimmed reads (only if `--auto-trim` and trimming was recommended)
- `results_agent/provenance/` - per-step provenance JSON files

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
├── api.py                    # Public Python API (inspect_inputs, plan_workflow, run_workflow, ...)
├── agent/
│   ├── state.py              # AgentState, Observation, Decision, RecommendedAction
│   ├── decision_engine.py    # Deterministic FASTQ QC rules
│   ├── fastq_agent.py        # FASTQ observe-decide-act orchestrator + report writers
│   ├── variant_decision_engine.py  # Deterministic variant/BAM/VCF QC rules
│   └── variant_agent.py      # Variant QC observe-decide-recommend orchestrator + report writers
├── parsers/
│   ├── fastqc.py             # FastQC zip/txt parser -> structured per-sample results
│   ├── multiqc.py            # MultiQC TSV/JSON/HTML parser
│   ├── samtools.py           # samtools flagstat/idxstats/stats parser
│   ├── bcftools.py           # bcftools stats parser
│   └── mosdepth.py           # mosdepth summary parser
├── schemas/
│   ├── inspection_result.schema.json
│   ├── workflow_plan.schema.json
│   ├── run_result.schema.json
│   ├── agent_result.schema.json
│   └── provenance_record.schema.json
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
├── scrna_qc_local.py         # scRNA-seq QC (fallback-only - prefer single-cell-rna-qc@life-sciences)
├── atac_qc_local.py          # ATAC-seq QC (fallback)
├── wgs_vcf_qc_local.py       # WGS/VCF QC (fallback)
├── nfcore_launcher.py        # nf-core workflow launcher
└── ...
```

---

## Amplicon Taxonomy Databases

The amplicon workflow supports planning for these taxonomy databases. Database files must be downloaded separately - the framework never downloads large references automatically.

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

- **No clinical claims** - enforced by tests. Variant QC and FASTQ agent outputs explicitly disclaim clinical interpretation.
- **No LLM at runtime** - all decisions are deterministic rules applied to parsed QC flags.
- **No large-file loading** - files above 50 MB are not read into memory. All large-file operations are delegated to CLI tools.
- **No silent failures** - every failed command is recorded in provenance with exit code and stderr.
- **No automatic reference downloads** - reference files must be provided locally.
- **Samplesheets require human review** - auto-generated samplesheets for nf-core pipelines contain placeholders that must be corrected before execution.
- **Dry-run by default** - expensive workflows are never launched without an explicit `--execute` flag.

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

## Biological interpretation and hypothesis generation

The framework includes a **deterministic interpretation scaffold** that converts
QC observations into structured findings and testable hypotheses.

### What it is

- Rule-based, not LLM-based. No language model is called at Python runtime.
- Outputs are JSON-first: the Markdown report is a rendering of the structured JSON.
- Every finding separates **technical artifact** from **plausible biological signal**.
- All outputs are explicitly labelled as hypotheses for human review.
- **No clinical claims are made.** `clinical_claims_allowed: false` is always present.

### What it covers

| Workflow | Observations interpreted |
|----------|--------------------------|
| `fastq-qc` | Adapter content, per-base quality, GC content, overrepresented sequences |
| `variant-qc` | Mapping rate, VCF record count, contig coverage, mean depth |

### CLI usage

```bash
# Generate interpretation from an existing agent report
python -m genomics_workflow_agent interpret \
  --input results_agent/agent_report.json \
  --out results_agent/interpretation/

# Output:
#   interpretation_report.json
#   interpretation_report.md
```

The `agent` command now automatically includes `biological_interpretation` in
its JSON and Markdown reports. No additional flags needed.

### Python API

```python
from genomics_workflow_agent import generate_interpretation

result = generate_interpretation(
    workflow="fastq-qc",
    observations=[...],  # from run_fastq_qc_agent()
    decisions=[...],
)
# result is a plain JSON-serializable dict
print(result["clinical_claims_allowed"])  # False
print(result["findings"][0]["should_filter"])  # False (GC content never auto-filtered)
```

### Example JSON output (fragment)

```json
{
  "interpretation_version": "1.0",
  "workflow": "fastq-qc",
  "clinical_claims_allowed": false,
  "findings": [
    {
      "finding_id": "FASTQ_GC_F001",
      "sample": "sample_A",
      "observation": "Per sequence GC content FAIL detected by FastQC.",
      "technical_explanations": [
        "PCR amplification bias can distort GC distribution.",
        "Contamination from another species may shift the GC profile."
      ],
      "plausible_biological_explanations": [
        "The sequenced organism may have a naturally GC-rich or GC-poor genome.",
        "For metagenomics or environmental samples, community composition shifts the GC distribution."
      ],
      "recommended_action": "Do not filter samples based on GC content alone. Biological review required.",
      "should_filter": false,
      "should_preserve_until_review": true
    }
  ],
  "hypotheses": [
    {
      "hypothesis_id": "FASTQ_GC_F001_H1",
      "statement": "The GC-content deviation may reflect technical bias or a biologically distinct composition.",
      "clinical_claim": false,
      "interpretation_type": "ambiguous"
    }
  ]
}
```

### Limitations

- The Python layer does not consult literature or external databases.
- Interpretations are scaffolds: they enumerate possibilities, not conclusions.
- A Claude Code `hypothesis-generation-specialist` agent (`.claude/agents/`) can
  be invoked separately to extend findings with scientific context.

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

### Interpretation layer extensions
- **RNA-seq**: map nf-core/rnaseq summary metrics (mapping rate, strandedness, dedup rate, gene-body coverage) to interpretation findings; flag strand-specific bias and rRNA contamination.
- **ATAC-seq**: interpret peak counts, FRiP score, fragment size distribution (mono/di-nucleosomal), and TSS enrichment as technical or biological signals.
- **Amplicon / metagenomics**: interpret rarefaction curves, taxa assignment quality, alpha-diversity flags, and chimera rates with community-context hypotheses.
- **scRNA-seq**: extend the scrna-qc-specialist with interpretation scaffolds for doublet rate, ambient RNA, mitochondrial fraction, and per-cluster QC.
- **Hypothesis ranking**: score hypotheses by prior probability given assay type and metadata, to surface the most likely explanation first.
- **Literature links**: integrate PubMed search (via `pubmed@life-sciences`) to attach relevant references to biological explanations.

### Pipeline and tooling
- Direct QIIME2/DADA2 execution (currently generates Nextflow/nf-core commands only)
- Richer MultiQC TSV parsing integrated into the FASTQ agent decision engine
- HTML report export with interpretation section
- Adaptive Nextflow decisions based on parsed pipeline completion outputs
- Benchmark datasets and integration test fixtures for CI
- R-based downstream workflow documentation (DESeq2, phyloseq)
- Stable service/API wrapper if deployment as a microservice is needed
