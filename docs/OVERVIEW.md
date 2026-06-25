# Portfolio Overview

## What this project is

A reproducible Claude Code workspace for genomics QC orchestration. It demonstrates how to integrate an LLM safely into bioinformatics workflows — handling the surrounding engineering concerns that production pipelines leave to the analyst.

## What it is not

This project does not reimplement samtools, bcftools, bedtools, Scanpy, nf-core, or the Anthropic Life Sciences skill suite. Those tools exist and are good. The local tools here are fallbacks and wrappers — used when official tools are unavailable or when structured provenance is needed.

## The problems this project addresses

**1. Safe large-file handling**
Claude Code cannot safely read a 50 GB BAM or 1 M-cell h5ad file into model context. This workspace enforces a strict rule: inspect metadata, route to the right tool, never paste raw data into the model.

**2. Graceful degradation**
In a research environment, not every compute node has samtools installed. Every tool here detects what is available and produces structured skip entries for unavailable metrics — with `missing_biological_conclusion` and `enable_with` fields — rather than silently failing or producing incomplete results.

**3. Provenance tracking**
External tool outputs (samtools flagstat, bcftools stats, bedtools intersect) are captured as structured JSON with command string, return code, stderr snippet, and tool version. This ensures reproducibility even when the same analysis is repeated later.

**4. Biological reasoning enforcement**
Genomics QC thresholds are not universal. High mitochondrial percentage in heart tissue is not the same as in PBMC. The agent configurations enforce a structured artifact-vs-biology review for every QC warning, using Filter / Flag / Stratify / Preserve language rather than binary pass/fail.

**5. Clinical safety guardrails**
WGS tools include a mandatory clinical disclaimer and are tested to ensure no clinical claims appear in any output. This is checked as a test (`test_no_clinical_claims_in_report`).

## Architecture

```
Claude Code (orchestrator)
  │
  ├── .claude/agents/           — Specialist agents with enforced constraints
  │     ├── genomics-file-inspector
  │     ├── scrna-qc-specialist (official-skill-first)
  │     ├── atac-qc-specialist
  │     ├── wgs-qc-variant-specialist (no clinical claims)
  │     ├── nextflow-pipeline-specialist
  │     ├── biology-interpretation-reviewer
  │     └── single-cell-modeling-specialist
  │
  ├── .claude/skills/           — Skill wrappers that route to agents or tools
  │
  ├── tools/                    — Local CLI tools (fallback and orchestration layer)
  │     ├── check_environment.py    — Pre-flight check
  │     ├── inspect_file.py         — File type detection and routing
  │     ├── scrna_qc_local.py       — scRNA QC (CSV/TSV fallback)
  │     ├── atac_qc_local.py        — ATAC QC + bedtools integration
  │     ├── wgs_vcf_qc_local.py     — WGS/VCF QC + samtools/bcftools
  │     ├── nfcore_launcher.py      — nf-core samplesheet builder, preflight, safe launcher
  │     ├── reference_validator.py  — Reference file validation
  │     └── report_builder.py       — Report assembly
  │
  └── Official Anthropic Life Sciences skills (preferred)
        single-cell-rna-qc@life-sciences
        scvi-tools@life-sciences
        nextflow-development@life-sciences
```

## Key engineering decisions

| Decision | Rationale |
|----------|-----------|
| Pure parser functions (importable without running external tools) | Enables fixture-based CI tests that do not require samtools/bcftools/bedtools |
| `skipped_metrics` list with `missing_biological_conclusion` | Every unavailable metric documents what biological question cannot be answered |
| bcftools cross-validation against local parser | Named discrepancy warnings; does not silently override local results |
| samtools depth safety guard (10 MB limit without `--intervals`) | Prevents accidental slow execution on large BAM files |
| Chromosome naming mismatch detection before any intersect | A silent zero-overlap from chr1 vs 1 style mismatch corrupts all bedtools results |
| Official-skill-first in all agent configs | Local tools are fallbacks, not defaults; agents stop and delegate when official skills are available |
| nf-core launcher dry-run default | `tools/nfcore_launcher.py` never executes a pipeline unless `--run` is explicit and all preflight checks pass; prevents accidental large-scale execution |
| Sarek samplesheet placeholders | Patient ID, tumor/normal status, and sex cannot be inferred from filenames; output uses explicit `PATIENT_ID` placeholders to force human review before execution |

## What the tests demonstrate

- **Parser tests**: Pure functions tested with fixture files. `tests/test_external_tools.py` runs 40 tests in CI without samtools, bcftools, or bedtools.
- **Degradation tests**: Tools run with no inputs; JSON output verified to contain structured skip entries.
- **Biological caveat tests**: Every skipped metric verified to have `missing_biological_conclusion`.
- **Clinical disclaimer test**: WGS report verified to contain no clinical claims and to include the clinical disclaimer.
- **Config tests**: Agent and skill configs verified to contain official-skill-first routing.
- **nf-core launcher tests**: Samplesheet builders, preflight, biological caveats, and no-clinical-claims tested without Nextflow installed.

## Test count

279 tests, 0 failures. Runs on Ubuntu and Windows, Python 3.11 and 3.12.
