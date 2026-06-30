# Skill: single-cell-modeling-wrapper

## Purpose
Wrapper for downstream single-cell modeling. Prefers `scvi-tools@life-sciences`. Decides whether modeling is biologically justified. Documents assumptions and caveats. Does not apply models when simpler approaches suffice.

## When to use
- After QC is complete.
- User needs batch correction, integration, reference mapping, or deep generative modeling.
- Multiple batches or modalities are present.

## Official tool to prefer
`scvi-tools@life-sciences` - for scVI, scANVI, totalVI, PeakVI, MultiVI workflows.

## Required inputs
- Filtered `.h5ad` (QC complete)
- Batch covariate column name (if integration)

## Optional inputs
- Cell type labels (for scANVI)
- Reference dataset (for reference mapping)
- Protein panel (for totalVI / CITE-seq)
- Peak matrix (for PeakVI / MultiVI)

## Workflow steps

1. **Confirm QC is complete**: modeling requires filtered, QC-passed data.
2. **Decide whether modeling is justified**:
   - Multiple technical batches? → integration may be appropriate
   - Multiple biological conditions? → check whether batch correction is appropriate
   - Single modality, single batch? → standard normalization + PCA may suffice
3. **Check official skill availability**: `scvi-tools@life-sciences`
4. **Select model** (see model selection table in agent).
5. **Document assumptions**: batch covariate definitions, covariates included/excluded, model architecture, training parameters.
6. **Run model via official skill** if available.
7. **Validate with known biology**: check that known cell types cluster appropriately, marker genes are enriched correctly.
8. **Produce output** with interpretation caveats.

## Expected outputs
- Model choice rationale
- Preprocessing requirement checklist
- Commands/workflow (via official skill or script)
- Interpretation caveats
- Validation steps recommended

## Biological reasoning checklist
- [ ] Confirmed QC is complete before modeling?
- [ ] Justified that batch correction is appropriate (not biological differences)?
- [ ] Documented all batch covariates and their origin?
- [ ] Checked reference compatibility before reference mapping?
- [ ] Validated model output with known biology?
- [ ] Reported uncertainty estimates where available?
