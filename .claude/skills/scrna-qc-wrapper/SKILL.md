# Skill: scrna-qc-wrapper

## Purpose
Wrapper for single-cell RNA-seq QC. The **primary path is always `single-cell-rna-qc@life-sciences`** when available and compatible. `tools/scrna_qc_local.py` is a **fallback only** - used exclusively when the official skill is unavailable or the input format is incompatible. This wrapper adds metadata validation, biological interpretation, provenance, and reporting layers around whichever path is taken.

## When to use
- Input is `.h5ad`, 10x directory, or count matrix (CSV/TSV/MTX).
- User wants QC metrics, filtering recommendations, or QC plots.
- Before clustering, annotation, integration, or modeling.

## Positioning - read before proceeding
`single-cell-rna-qc@life-sciences` is the Anthropic Life Sciences official scRNA QC skill. It is the **default and preferred** path. This local wrapper exists to:
1. Check official skill availability before doing anything else.
2. Add metadata, privacy, biological interpretation, and provenance layers around the official skill's output.
3. Fall back to `tools/scrna_qc_local.py` **only** when the official skill is unavailable or incompatible.

**Do not use `tools/scrna_qc_local.py` when `single-cell-rna-qc@life-sciences` is available.**
**Do not duplicate scRNA QC unnecessarily.** The local fallback does not replicate the full capability of the official skill.

## Required inputs
- Count matrix: `.h5ad`, 10x directory, CSV, or TSV

## Optional inputs
- `--species` (human / mouse / other)
- `--tissue` (UBERON term or free text)
- `--disease` (healthy / tumor / inflamed / other)
- `--protocol` (10x_v3 / 10x_v2 / dropseq / smartseq / other)
- `--chemistry` (v2 / v3 / v3.1 / other)
- `--genome-build` (GRCh38 / GRCh37 / mm10 / mm39 / other)
- `--sample-metadata` CSV with sample-level annotations
- `--mito-prefix` (default: MT- for human, mt- for mouse)
- `--mad-threshold` (default: 5.0)
- `--recommend-only` (do not apply filters, only recommend)

## Workflow steps

### Pre-flight
1. Confirm input file exists and is readable.
2. **Check if `single-cell-rna-qc@life-sciences` is available. This is the first and most important check.**
3. Collect or prompt for: species, tissue, disease state, protocol, chemistry, genome build.
4. If metadata is missing, proceed with labeled conservative assumptions.
5. For human data: confirm data handling is appropriate (no upload to external services unless authorized).

### Official skill path - use this when available
If `single-cell-rna-qc@life-sciences` is available and input is compatible (`.h5ad` or 10x `.h5`):
- Invoke `single-cell-rna-qc@life-sciences`.
- Do NOT also run `tools/scrna_qc_local.py`. Do not duplicate QC unnecessarily.
- Capture outputs (metrics CSV, plots, report section).
- Add metadata, provenance, and biological interpretation wrapper (Steps below).
- Record in report: "QC performed by official `single-cell-rna-qc@life-sciences` skill."

### Local fallback path - use only when official skill is unavailable
If `single-cell-rna-qc@life-sciences` is **not available** or the input format is incompatible:
- Record reason: "Official skill unavailable" or "Input format incompatible with official skill."
- Flag in report: "Local fallback used - `tools/scrna_qc_local.py`. Results may be less comprehensive than the official skill. Consider re-running with official skill when available."
```bash
python tools/scrna_qc_local.py \
  --input <path> \
  --species <species> \
  --tissue <tissue> \
  --protocol <protocol> \
  --chemistry <chemistry> \
  --genome-build <build> \
  --recommend-only \
  --output-dir <output_dir>
```

### Biological interpretation
After QC output is collected:
1. Build artifact-versus-biology table for all QC warnings.
2. Document biological justification for each proposed filter.
3. Note every assumption made.
4. List limitations.

### Report assembly
```bash
python tools/report_builder.py \
  --scrna-qc-dir <qc_output_dir> \
  --output-dir reports/ \
  --title "scRNA-seq QC Report"
```

## Expected outputs
- QC metrics CSV
- QC plots directory
- Summary JSON
- Markdown report section
- Artifact-versus-biology table
- Assumptions and limitations

## Failure modes
- Corrupt or truncated h5ad: inspect_file.py will report the issue
- Missing mitochondrial genes: report which mito prefix was tried; do not fail silently
- All cells flagged by a criterion: do not apply - report instead
- Missing sample metadata: proceed with assumptions, document them

## Reproducibility requirements
- Record input path, file hash where feasible
- Record species, tissue, protocol, genome build assumptions
- Record tool version, scanpy version, Python version
- Record all filtering parameters used

## Biological reasoning checklist
- [ ] Mitochondrial prefix confirmed for species?
- [ ] High mito% cells reviewed in tissue context?
- [ ] High UMI cells reviewed for doublet vs. large-cell biology?
- [ ] Low-gene cells reviewed for mature cell type biology?
- [ ] Rare populations identified before filtering?
- [ ] Before/after cell counts documented?
- [ ] Known marker sanity check performed?
