# Genomics Agent Guardrails

## Bioinformatics runtime

Claude Code runs on the **host**. Do not assume FastQC, MultiQC, samtools, bcftools,
fastp, cutadapt, mosdepth, bedtools, or Nextflow are installed on the host.

Use Docker Compose when real bioinformatics tools must be executed:

```bash
# General pattern
docker compose run --rm genomics-agent <command>

# CLI help
docker compose run --rm genomics-agent python -m genomics_workflow_agent --help

# Run tests
docker compose run --rm genomics-agent python -m pytest tests/ -q

# FASTQ QC agent - dry-run (no tools needed)
docker compose run --rm genomics-agent python -m genomics_workflow_agent agent \
  --input examples/ --workflow fastq-qc --out results_agent_smoke/

# FASTQ QC agent - execute (FastQC + MultiQC run inside container)
docker compose run --rm genomics-agent python -m genomics_workflow_agent agent \
  --input examples/ --workflow fastq-qc --out results_agent_smoke/ --execute

# Variant QC agent - dry-run
docker compose run --rm genomics-agent python -m genomics_workflow_agent agent \
  --input examples/ --workflow variant-qc --out results_variant_agent_smoke/
```

Dry-run commands (no `--execute`) do not call any external bioinformatics tool and
can be run directly on the host using `python -m genomics_workflow_agent ...`.

## Execution rules

- Do not run expensive Nextflow/nf-core workflows unless the user explicitly approves.
- Do not download large references, genomes, indices, or taxonomy databases unless
  the user explicitly approves and the data directory is already configured.
- Do not commit or push unless the user explicitly asks.
- Do not add generated result folders (results/, results_*/) to Git.
- Use Docker for bioinformatics tools when possible.
- Use PowerShell-compatible commands on Windows. Avoid `rm -rf`, `&&` chaining, and
  other bash-only constructs unless running inside Linux/container.
- Check `.gitignore` before creating output directories to confirm they will be ignored.



Claude Code is an orchestrator, not a raw large-file parser.

## Data handling

- Never load large genomics files directly into model context.
- Always inspect file type, size, compression, schema, and metadata before analysis.
- Preserve raw data. Never overwrite original input files.
- Treat human genomic data as privacy-sensitive by default.
- Do not upload or transmit private genomic data unless explicitly instructed and appropriate.

## Reproducibility

- Always record commands, package versions, tool versions, reference files, genome build, annotation version, and parameters.
- Prefer reproducible scripts, plots, logs, JSON summaries, CSV outputs, and Markdown reports.
- Never automatically download large reference files during routine QC.
- If a reference file is missing, skip dependent metrics and report the limitation clearly.

## Biological reasoning

- Never apply fixed QC thresholds without plotting distributions and considering species, tissue, protocol, chemistry, batch, disease state, expected cell types, sample preparation, and genome build.
- Use robust outlier detection such as MAD-based thresholds where appropriate.
- Always distinguish technical artifact from plausible biological signal.
- High mitochondrial percentage, low RNA complexity, high counts, unusual accessibility, or unusual variant patterns are not automatically artifacts. They require biological context.
- Rare cell populations must not be removed simply because they are outliers.
- Do not interpret clusters, peaks, variants, or QC filters as final biological truth without validation.

## Clinical safety

- Never make clinical claims from genomic variants.
- Clearly separate exploratory research interpretation from validated clinical interpretation.
- Do not recommend medical action.

## Official Claude Life Sciences tools

When the following official skills/plugins are available, prefer them over local fallbacks:

| Skill | Use Case |
|-------|----------|
| `single-cell-rna-qc@life-sciences` | Standard scRNA-seq QC on `.h5ad` or 10x `.h5` |
| `scvi-tools@life-sciences` | Downstream single-cell modeling, integration, batch correction |
| `nextflow-development@life-sciences` | nf-core pipeline orchestration (rnaseq, sarek, atacseq) |
| `10x-genomics@life-sciences` | 10x Genomics Cloud data access |
| `pubmed@life-sciences` | Biological background and marker interpretation |

Local tools in `tools/` are fallbacks and wrappers - not replacements for official skills.

## Agent delegation

Use project subagents in `.claude/agents/` to delegate specialized tasks:
- File inspection → `genomics-file-inspector`
- scRNA QC → `scrna-qc-specialist`
- ATAC QC → `atac-qc-specialist`
- WGS/VCF QC → `wgs-qc-variant-specialist`
- Nextflow pipelines → `nextflow-pipeline-specialist`
- Biological review → `biology-interpretation-reviewer`
- Downstream modeling → `single-cell-modeling-specialist`

## Biological interpretation is mandatory

Every analysis must include a biological interpretation section. Every QC warning must address:
1. What was observed?
2. What technical artifact could explain it?
3. What biological state could explain it?
4. What metadata would help distinguish artifact from biology?
5. What validation should be done?
6. Should data be filtered, flagged, stratified, or preserved?
