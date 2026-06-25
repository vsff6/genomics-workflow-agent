---
name: nfcore-orchestration
description: Plans nf-core pipeline runs using tools/nfcore_launcher.py. Generates samplesheets, validates preflight requirements, and writes provenance records. Default is dry-run — never executes without explicit --run flag and passing preflight.
---

# nf-core Orchestration Skill

## What this skill does

Uses `tools/nfcore_launcher.py` to:
1. Generate nf-core-compatible samplesheets from a FASTQ directory
2. Run preflight checks (Nextflow, Docker/Singularity, references)
3. Construct the Nextflow command
4. Write a plan: `nfcore_plan.json`, `nfcore_plan.md`, `commands.sh`
5. Optionally execute the pipeline (`--run` only, after preflight passes)

This is the **local planning and provenance layer**. For production execution, prefer `nextflow-development@life-sciences`.

## When to invoke

- User provides a directory of FASTQ files and asks to run rnaseq, sarek, or atacseq
- User wants a samplesheet generated from FASTQ filenames
- User wants to check whether local requirements are met before launching a pipeline
- User wants a reproducible record of the pipeline command and parameters

## When NOT to invoke

- User has `.h5ad` / 10x `.h5` input → use `single-cell-rna-qc@life-sciences`
- User wants batch correction or downstream modeling → use `scvi-tools@life-sciences`
- Production pipeline with validated data → use `nextflow-development@life-sciences`

## Invocation pattern

```bash
# Dry-run plan (default — always start here)
python tools/nfcore_launcher.py \
  --workflow rnaseq \
  --input-dir /path/to/fastqs \
  --genome GRCh38 \
  --output-dir reports/nfcore_rnaseq \
  --dry-run

# Execute (only after reviewing the plan and resolving all blockers)
python tools/nfcore_launcher.py \
  --workflow rnaseq \
  --input-dir /path/to/fastqs \
  --genome GRCh38 \
  --output-dir reports/nfcore_rnaseq \
  --run
```

## Output files

| File | Contents |
|------|----------|
| `nfcore_plan.json` | Structured plan: executors, blockers, warnings, samplesheet result, caveats |
| `nfcore_plan.md` | Human-readable plan with biological caveats section |
| `commands.sh` | Exact Nextflow command — review before running |
| `samplesheet_{workflow}.csv` | Generated samplesheet — always review before running |
| `nfcore_launcher.log` | Execution log |

## Safety requirements

- **Default is dry-run** — never execute without `--run`
- **--run requires preflight pass** — tool exits non-zero if blockers exist
- **Never download large references** — provide `--genome` key or `--fasta` path
- **Samplesheets require human review** — especially sarek (tumor/normal, patient ID)
- **Always invoke `biology-interpretation-reviewer`** after parsing pipeline outputs

## Mandatory biological caveat

Every plan output must include the biological caveats section. Successful pipeline execution is not biological or clinical validation.

## Supported workflows

| Workflow | Key samplesheet field to review |
|----------|---------------------------------|
| `rnaseq` | `strandedness` (must confirm from library prep) |
| `sarek` | `status` (tumor/normal), `patient`, `sex` (cannot be inferred from filenames) |
| `atacseq` | `replicate` (must reflect experimental design), blacklist BED |
