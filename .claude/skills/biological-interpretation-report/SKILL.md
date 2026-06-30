# Skill: biological-interpretation-report

## Purpose
Produces a final biological interpretation section for any genomics analysis. Reviews outputs from official skills, nf-core workflows, and local QC scripts. Produces the standard artifact-versus-biology table. Lists limitations and next analyses. This skill is mandatory - never skip biological interpretation.

## When to use
- After any QC run, pipeline, or analysis step.
- Before reporting results to users.
- When any QC observation, filter, or threshold is proposed.

## Official tool to prefer
`pubmed@life-sciences` - for biological background, marker interpretation, and tissue/disease biology. Use only when evidence is directly relevant to the species, tissue, assay, and condition. Do not use PubMed to justify arbitrary thresholds.

## Required inputs
- QC output directory (JSON summaries, metrics CSVs, Markdown sections)

## Optional inputs
- Species, tissue, disease context
- Protocol and chemistry
- Genome build
- Known marker genes
- Expected cell types

## Workflow steps

1. **Read available QC outputs** from reports directory.
2. **Gather biological context**: species, tissue, disease, protocol, expected biology.
3. **For every observation**: apply artifact-versus-biology framework.
4. **Build standard table** (see below).
5. **Write plausibility notes** per observation.
6. **List validation suggestions** for ambiguous observations.
7. **List limitations** - what cannot be determined.
8. **List next analyses** - what would resolve ambiguities.
9. **Assemble report section** and pass to `report_builder.py`.

## Standard artifact-versus-biology table

| Observation | Possible Technical Explanation | Possible Biological Explanation | Evidence Supporting Artifact | Evidence Supporting Biology | Recommended Follow-up | Confidence |
|-------------|-------------------------------|--------------------------------|------------------------------|-----------------------------|-----------------------|------------|

Confidence scale:
- `low` - weak or ambiguous evidence
- `moderate` - some support, alternatives remain plausible
- `high` - strong specific evidence; use sparingly

## Expected outputs
- Artifact-versus-biology table (Markdown)
- Biological plausibility notes (per observation)
- Validation suggestions (per ambiguous observation)
- Limitations section
- Suggested next analyses section

## Failure modes
- No QC outputs available: request QC to be run first
- Context metadata missing: document assumptions and reduce confidence to `low`

## Reproducibility requirements
- Record biological context used in interpretation
- Record literature references if `pubmed@life-sciences` was used
- Record confidence levels with justification

## Biological reasoning checklist
- [ ] Every QC observation has at least one technical AND one biological explanation?
- [ ] No observation marked `high` confidence artifact without strong evidence?
- [ ] Rare populations explicitly addressed?
- [ ] Tissue-specific biology considered?
- [ ] Assay-specific artifacts considered?
- [ ] Clinical claims avoided?
- [ ] Validation steps suggested for ambiguous cases?
