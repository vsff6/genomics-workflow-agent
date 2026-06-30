# Changelog

All notable changes are documented here. Versioning follows [Semantic Versioning](https://semver.org/).

---

## [0.2.0] - 2025-06-24

### Added

- `tools/wgs_vcf_qc_local.py` v2.0: samtools flagstat/idxstats/stats/depth integration with structured JSON output; commands, return codes, and tool versions recorded for provenance
- `tools/wgs_vcf_qc_local.py` v2.0: bcftools stats cross-validation; discrepancies reported as named warnings, not silent overrides
- `tools/atac_qc_local.py` v2.0: bedtools intersect for blacklist overlap fraction and authoritative FRiP
- `tools/atac_qc_local.py` v2.0: chromosome naming mismatch detection (UCSC chr1 vs Ensembl 1) across fragments, peaks, blacklist, and GTF inputs
- `tools/atac_qc_local.py` v2.0: TSS enrichment deeptools command documented in skip entry when BAM and deeptools are both present
- `tests/test_external_tools.py`: 40 parser-only tests using fixture files; run without samtools, bcftools, or bedtools installed
- `tests/fixtures/`: samtools_flagstat.txt, samtools_idxstats.txt, bcftools_stats.txt
- GitHub presentation files: LICENSE, CITATION.cff, CONTRIBUTING.md, SECURITY.md, CHANGELOG.md
- Docs: docs/PORTFOLIO_OVERVIEW.md, docs/ROADMAP.md, docs/DEMO_OUTPUTS.md, examples/README.md

### Fixed

- `.gitignore`: removed `env/` directory pattern that prevented `env/environment.yml` from being tracked
- `run_samtools_idxstats`: removed dead code - stray duplicate dict key with broken string literal
- `atac_qc_local.py`: removed dead `cmd_count_total` variable (assigned but never used)

---

## [0.1.0] - 2024

### Added

- Seven local CLI tools: `check_environment.py`, `inspect_file.py`, `scrna_qc_local.py`, `atac_qc_local.py`, `wgs_vcf_qc_local.py`, `reference_validator.py`, `report_builder.py`
- Seven Claude Code agent configurations in `.claude/agents/`
- Six skill wrappers in `.claude/skills/`
- `CLAUDE.md`: data handling, reproducibility, biological reasoning, and clinical safety guardrails
- `env/environment.yml`: conda environment (Python 3.11, scverse stack)
- GitHub Actions CI: Ubuntu + Windows, Python 3.11/3.12
- End-to-end toy demo: `examples/run_tiny_demo.sh`
- All skipped metrics carry `missing_biological_conclusion` and `enable_with` fields
