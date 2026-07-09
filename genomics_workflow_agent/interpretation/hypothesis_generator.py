"""Deterministic rules mapping QC observations to biological interpretation scaffolds."""
from __future__ import annotations

from typing import Any

from genomics_workflow_agent.interpretation.evidence import (
    collect_obs_by_category,
    obs_message_list,
    obs_to_evidence_source,
)
from genomics_workflow_agent.interpretation.models import (
    INTERPRETATION_VERSION,
    Finding,
    Hypothesis,
    InterpretationResult,
)

_FASTQ_LIMITATIONS = [
    "FASTQ QC metrics are technical quality measures. They cannot establish biological causation.",
    "GC content deviations may reflect the organism's natural GC composition rather than contamination.",
    "Overrepresented sequences may be biologically abundant transcripts rather than contaminants.",
    "Adapter content reflects library preparation choices, not biological signal.",
    "These interpretations are scaffolds for human review, not conclusions.",
]

_VARIANT_LIMITATIONS = [
    "Variant QC metrics are technical quality measures. They cannot establish variant significance.",
    "Low mapping rate may reflect biological divergence from the reference, not necessarily failure.",
    "Zero VCF records cannot be interpreted as absence of biological variation.",
    "Coverage metrics are sensitive to capture design and reference choice.",
    "These interpretations are scaffolds for human review, not conclusions.",
]


def _fastq_gc_fail(obs_list: list[dict], finding_id: str, sample: str) -> tuple[Finding, Hypothesis]:
    finding = Finding(
        finding_id=finding_id,
        workflow="fastq-qc",
        sample=sample,
        observation="Per sequence GC content FAIL detected by FastQC.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "fastqc:gc_content",
        technical_explanations=[
            "PCR amplification bias can distort GC distribution.",
            "Contamination from another species or sample may shift the GC profile.",
            "Library preparation artifacts (e.g., AT/GC clamps, primer composition) can cause deviation.",
            "Reference genome mismatch if GC model was computed from a different species.",
        ],
        plausible_biological_explanations=[
            "The sequenced organism or strain may have a naturally GC-rich or GC-poor genome.",
            "For metagenomics or environmental samples, community composition shifts the GC distribution.",
            "AT-rich or GC-rich genomic regions (e.g., centromeres, promoters) can dominate targeted assays.",
            "Single-cell samples from GC-biased cell types may show skewed distributions.",
        ],
        metadata_needed=[
            "Expected GC content of the target organism or community.",
            "Library preparation protocol (PCR-free vs. PCR-amplified).",
            "Assay type (whole-genome, targeted panel, amplicon, single-cell).",
            "Whether multiple samples from the same experiment show the same pattern.",
        ],
        recommended_validation=[
            "Compare GC profile across all samples in the batch.",
            "Check the organism's published GC content against the observed distribution.",
            "Run BLAST or Kraken2 on a subset of reads to identify contaminating sequences.",
            "Inspect the FastQC HTML GC content plot for the shape of the deviation.",
        ],
        recommended_action=(
            "Do not filter samples or discard reads based on GC content alone. "
            "Biological review with organism/protocol metadata is required before any action."
        ),
        confidence="medium",
        should_filter=False,
        should_preserve_until_review=True,
    )

    hypothesis = Hypothesis(
        hypothesis_id=f"{finding_id}_H1",
        statement=(
            "The GC-content deviation may reflect either a technical bias "
            "(PCR amplification or contamination) or a biologically distinct composition "
            "(organism GC content, community structure, or targeted region composition). "
            "These cannot be distinguished from QC data alone."
        ),
        supporting_observations=obs_message_list(obs_list),
        alternative_explanations=[
            "Technical: PCR amplification bias, contamination, adapter sequence GC composition.",
            "Biological: organism-specific GC content, microbiome community shift, targeted assay design.",
        ],
        validation_steps=[
            "Obtain the known GC content of the reference organism.",
            "Compare GC profile across samples and batches.",
            "Run taxonomic classification on a subset of reads.",
        ],
        confidence="low",
        interpretation_type="ambiguous",
        clinical_claim=False,
    )

    return finding, hypothesis


def _fastq_adapter(obs_list: list[dict], finding_id: str, sample: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="fastq-qc",
        sample=sample,
        observation="Adapter Content WARN or FAIL detected by FastQC.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "fastqc:adapter_content",
        technical_explanations=[
            "Insert size shorter than read length causes adapter read-through.",
            "Incomplete demultiplexing or index contamination from adjacent libraries.",
            "Trimming was not applied upstream or the trimming was insufficient.",
        ],
        plausible_biological_explanations=[
            "For short-fragment libraries (ATAC-seq, ChIP-seq, small RNA), adapter content is expected "
            "and does not indicate a biological problem.",
        ],
        metadata_needed=[
            "Expected insert size distribution for this library type.",
            "Whether trimming was applied prior to this QC run.",
            "Library preparation protocol (ATAC-seq, short-read WGS, small RNA, etc.).",
        ],
        recommended_validation=[
            "Inspect the FastQC Adapter Content plot to identify which adapter sequence is present.",
            "Check insert size distribution if available.",
            "Confirm whether trimming is appropriate for the downstream alignment tool.",
        ],
        recommended_action=(
            "Run adapter trimming (fastp or cutadapt) with parameters matched to the library protocol. "
            "Preserve original FASTQ files before trimming."
        ),
        confidence="high",
        should_filter=False,
        should_preserve_until_review=False,
    )


def _fastq_per_base_quality(obs_list: list[dict], finding_id: str, sample: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="fastq-qc",
        sample=sample,
        observation="Per base sequence quality FAIL detected by FastQC.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "fastqc:per_base_quality",
        technical_explanations=[
            "3' end quality decay is a normal feature of sequencing-by-synthesis chemistry.",
            "Instrument or flowcell issues can cause systematic quality drops.",
            "Library complexity issues or degraded input material can reduce quality scores.",
            "Run-specific reagent degradation may affect later cycles.",
        ],
        plausible_biological_explanations=[
            "Quality decay at 3' ends is rarely biological; it is almost always a technical artifact.",
            "In rare cases, samples with unusual base composition may show systematic quality effects.",
        ],
        metadata_needed=[
            "Whether the quality drop is at the 3' end (expected) or systematic throughout.",
            "Run date and instrument/flowcell ID to identify batch effects.",
            "Comparison to other samples on the same run.",
        ],
        recommended_validation=[
            "Compare quality profiles across all samples on the same run.",
            "Check whether the quality drop is read-end-specific or systematic.",
            "Evaluate downstream alignment quality after trimming.",
        ],
        recommended_action=(
            "Apply quality-based trimming (fastp or cutadapt) with parameters appropriate for the "
            "downstream aligner. Review whether the quality drop is systematic or sample-specific."
        ),
        confidence="high",
        should_filter=False,
        should_preserve_until_review=False,
    )


def _fastq_overrepresented(obs_list: list[dict], finding_id: str, sample: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="fastq-qc",
        sample=sample,
        observation="Overrepresented sequences WARN or FAIL detected by FastQC.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "fastqc:overrepresented_sequences",
        technical_explanations=[
            "Adapter or primer sequences present in reads.",
            "rRNA contamination (common in RNA-seq without ribo-depletion).",
            "PCR duplicates inflating a small number of sequences.",
            "Index-hopping contamination from adjacent libraries.",
            "Low-complexity library producing repeated sequence content.",
        ],
        plausible_biological_explanations=[
            "Highly abundant transcripts may legitimately dominate RNA-seq libraries.",
            "Amplicon assays produce overrepresented sequences by design.",
            "Ribosomal or mitochondrial sequences are biologically overrepresented in some sample types.",
        ],
        metadata_needed=[
            "Assay type (RNA-seq, amplicon, ChIP-seq, etc.).",
            "Whether ribo-depletion or poly-A selection was applied (for RNA-seq).",
            "BLAST annotation of the overrepresented sequences from the FastQC report.",
        ],
        recommended_validation=[
            "Inspect the overrepresented sequence table in the FastQC HTML report.",
            "BLAST the flagged sequences to identify their source.",
            "Check whether the same sequences are overrepresented across samples or sample-specific.",
        ],
        recommended_action=(
            "Do not filter reads without identifying the source of overrepresentation. "
            "If adapter-derived: trim. If rRNA: consider ribo-depletion. "
            "If biologically abundant: preserve and annotate."
        ),
        confidence="medium",
        should_filter=False,
        should_preserve_until_review=True,
    )


def _fastq_all_pass(finding_id: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="fastq-qc",
        sample="all",
        observation="All tracked FastQC modules passed. No QC-derived biological warnings from FASTQ-level checks.",
        evidence_source="fastqc:overall_qc",
        technical_explanations=[
            "No adapter, quality, GC, or overrepresentation flags were triggered.",
        ],
        plausible_biological_explanations=[
            "FASTQ-level QC passing does not confirm biological validity of the data.",
            "Downstream analysis (alignment, variant calling, quantification) may reveal additional issues.",
        ],
        metadata_needed=[
            "Downstream QC metrics (alignment rate, duplication rate, coverage) for complete assessment.",
        ],
        recommended_validation=[
            "Proceed to alignment or downstream analysis.",
            "Review alignment-level QC after mapping.",
        ],
        recommended_action="Proceed to downstream analysis. No QC-based intervention required at FASTQ level.",
        confidence="high",
        should_filter=False,
        should_preserve_until_review=False,
    )


def _variant_low_mapping(obs_list: list[dict], finding_id: str, sample: str) -> tuple[Finding, Hypothesis]:
    finding = Finding(
        finding_id=finding_id,
        workflow="variant-qc",
        sample=sample,
        observation="Mapped read percentage below the expected threshold.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "samtools:alignment",
        technical_explanations=[
            "Reference genome mismatch (wrong species, build, or strain).",
            "Library preparation failure producing low-quality or adapter-only reads.",
            "Contamination from another species reducing effective mapping rate.",
            "Incorrect alignment parameters or aligner misconfiguration.",
        ],
        plausible_biological_explanations=[
            "The sample may be from an organism or strain significantly diverged from the reference.",
            "Metagenomic or mixed samples may have reads from multiple organisms, reducing single-reference mapping.",
            "Structural variation (large deletions, inversions, or translocations) can reduce local mapping.",
        ],
        metadata_needed=[
            "Expected organism and reference genome build.",
            "Library preparation protocol and quality metrics.",
            "Whether contamination screening was applied prior to alignment.",
            "Alignment log and parameters used.",
        ],
        recommended_validation=[
            "Check alignment log for warnings about unmapped reads.",
            "Run a contamination screen on a subset of unmapped reads (e.g., Kraken2 or FastQ Screen).",
            "Compare alignment rate to other samples on the same run.",
            "Verify the reference genome build matches the sample preparation.",
        ],
        recommended_action=(
            "Investigate alignment failure before proceeding to variant calling. "
            "Do not proceed to variant calling without understanding the cause of low mapping."
        ),
        confidence="high",
        should_filter=False,
        should_preserve_until_review=True,
    )

    hypothesis = Hypothesis(
        hypothesis_id=f"{finding_id}_H1",
        statement=(
            "Low mapping rate may reflect either a technical failure (reference mismatch, "
            "library quality) or a biological signal (divergent strain, metagenomic mixture, "
            "structural variation). These cannot be distinguished from mapping statistics alone."
        ),
        supporting_observations=obs_message_list(obs_list),
        alternative_explanations=[
            "Technical: wrong reference, aligner misconfiguration, library degradation.",
            "Biological: strain divergence from reference, mixed-species sample, large structural variants.",
        ],
        validation_steps=[
            "Classify unmapped reads using a taxonomic tool.",
            "Attempt re-alignment with a more permissive reference or alternative aligner.",
            "Compare to positive-control samples with known mapping rates.",
        ],
        confidence="low",
        interpretation_type="ambiguous",
        clinical_claim=False,
    )

    return finding, hypothesis


def _variant_zero_vcf(obs_list: list[dict], finding_id: str, sample: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="variant-qc",
        sample=sample,
        observation="VCF contains zero variant records.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "bcftools:vcf_content",
        technical_explanations=[
            "Incorrect variant caller invocation (wrong inputs, wrong mode, wrong reference).",
            "Empty or corrupt input BAM file passed to the variant caller.",
            "Wrong genomic intervals or target regions specified.",
            "Aggressive pre-filtering removed all variants before this QC stage.",
            "Variant caller output may have been written to a different location.",
        ],
        plausible_biological_explanations=[
            "Zero records is not evidence of no biological variation. "
            "A true absence of variation in a target region is possible but cannot be inferred from a failed call.",
        ],
        metadata_needed=[
            "Variant caller command and version used.",
            "Input BAM file integrity and read count.",
            "Target genomic intervals.",
            "Filtering thresholds applied post-calling.",
        ],
        recommended_validation=[
            "Re-run variant calling with explicit logging and check for error messages.",
            "Verify the input BAM is non-empty and properly indexed.",
            "Confirm the target intervals overlap with sequenced regions.",
            "Compare to expected variant count for the assay and sample type.",
        ],
        recommended_action=(
            "Investigate the variant calling pipeline for errors. "
            "Zero VCF records should be treated as a pipeline failure until proven otherwise. "
            "Do not interpret as biological absence of variation."
        ),
        confidence="high",
        should_filter=False,
        should_preserve_until_review=True,
    )


def _variant_zero_read_contigs(obs_list: list[dict], finding_id: str, sample: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="variant-qc",
        sample=sample,
        observation="Many reference contigs have zero mapped reads.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "samtools:contig_coverage",
        technical_explanations=[
            "Reference assembly includes contigs not targeted by the assay (decoy contigs, unplaced scaffolds).",
            "Contig naming mismatch between BAM header and reference (chr vs. numeric naming).",
            "Targeted capture panel does not cover all reference contigs by design.",
            "Low overall sequencing depth leaving many contigs uncovered.",
        ],
        plausible_biological_explanations=[
            "Some organisms or strains have genomic regions with no representation in standard references.",
            "Sample-specific structural deletions could result in missing coverage on some contigs.",
        ],
        metadata_needed=[
            "Reference genome assembly version and contig annotation.",
            "Whether decoy contigs are included in the reference.",
            "Capture kit design and targeted intervals.",
            "Expected coverage distribution for the assay.",
        ],
        recommended_validation=[
            "List the zero-read contigs and check whether they are decoy or unplaced scaffolds.",
            "Verify contig naming conventions match between BAM and reference.",
            "Confirm that targeted intervals map to covered contigs.",
        ],
        recommended_action=(
            "Review reference compatibility and capture design. "
            "Many zero-read contigs in decoy/unplaced sequences is expected in some assays. "
            "Contig naming mismatches should be corrected before variant calling."
        ),
        confidence="medium",
        should_filter=False,
        should_preserve_until_review=True,
    )


def _variant_low_coverage(obs_list: list[dict], finding_id: str, sample: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="variant-qc",
        sample=sample,
        observation="Mean sequencing coverage is below the recommended threshold.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "mosdepth:coverage",
        technical_explanations=[
            "Insufficient total sequencing depth for the target genome size.",
            "Uneven coverage distribution due to capture inefficiency or GC bias.",
            "Library complexity too low, producing many duplicate reads but low unique coverage.",
            "Sample preparation issues reducing usable read count.",
        ],
        plausible_biological_explanations=[
            "Some assays (targeted panels, amplicon) may intentionally have uneven coverage by design. "
            "Low overall mean may not reflect poor coverage of targeted regions.",
        ],
        metadata_needed=[
            "Expected coverage target for the assay.",
            "Coverage uniformity metrics (e.g., fraction of target bases at adequate depth).",
            "Duplicate rate and library complexity metrics.",
            "Capture panel design if applicable.",
        ],
        recommended_validation=[
            "Generate per-region coverage plots for target intervals.",
            "Compare coverage to the minimum required for the variant caller settings.",
            "Review library complexity and duplication rate.",
            "Check whether coverage uniformity, not only mean coverage, meets assay requirements.",
        ],
        recommended_action=(
            "Review coverage distribution before interpreting variant calls. "
            "Low coverage reduces sensitivity for heterozygous or low-frequency variants. "
            "Do not infer genotypes from low-coverage data without appropriate statistical models."
        ),
        confidence="high",
        should_filter=False,
        should_preserve_until_review=True,
    )


def _variant_no_issues(finding_id: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="variant-qc",
        sample="all",
        observation=(
            "No critical technical issues detected in tracked variant QC metrics "
            "(mapping rate, VCF record count, coverage)."
        ),
        evidence_source="variant_decision_engine:overall_qc",
        technical_explanations=[
            "Mapping rates, VCF record counts, and coverage metrics are within expected ranges.",
        ],
        plausible_biological_explanations=[
            "Technical QC passing does not confirm biological validity of variant calls.",
            "Ti/Tv ratio, het/hom ratio, and allele frequency distributions require expert review.",
            "Variant accuracy depends on caller settings, filtering, and annotation not assessed here.",
        ],
        metadata_needed=[
            "Variant annotation and filtering thresholds for the specific assay.",
            "Expected Ti/Tv and het/hom ratios for the organism and assay type.",
        ],
        recommended_validation=[
            "Proceed to variant annotation and functional filtering.",
            "Review Ti/Tv and het/hom ratios for the variant set.",
            "Compare to positive-control samples or expected population statistics.",
        ],
        recommended_action=(
            "Technical QC did not detect major red flags. "
            "Proceed to variant annotation and biological review with appropriate domain expertise."
        ),
        confidence="medium",
        should_filter=False,
        should_preserve_until_review=False,
    )


def generate_fastq_interpretation(
    observations: list[dict],
    decisions: list[dict],
) -> InterpretationResult:
    findings: list[Finding] = []
    hypotheses: list[Hypothesis] = []
    finding_counter = 0

    def next_id(prefix: str) -> str:
        nonlocal finding_counter
        finding_counter += 1
        return f"{prefix}_F{finding_counter:03d}"

    # Group observations by category and sample
    seen_categories: set[str] = set()

    gc_obs = collect_obs_by_category(observations, "gc_content", status="fail")
    if gc_obs:
        samples = {o.get("sample", "unknown") for o in gc_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in gc_obs if o.get("sample") == sample]
            fid = next_id("FASTQ_GC")
            f, h = _fastq_gc_fail(sample_obs, fid, sample)
            findings.append(f)
            hypotheses.append(h)
        seen_categories.add("gc_content")

    adapter_obs = collect_obs_by_category(observations, "adapter_content")
    adapter_obs = [o for o in adapter_obs if o.get("status") in ("warn", "fail")]
    if adapter_obs:
        samples = {o.get("sample", "unknown") for o in adapter_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in adapter_obs if o.get("sample") == sample]
            fid = next_id("FASTQ_ADAPTER")
            findings.append(_fastq_adapter(sample_obs, fid, sample))
        seen_categories.add("adapter_content")

    pbq_obs = collect_obs_by_category(observations, "per_base_quality")
    pbq_obs = [o for o in pbq_obs if o.get("status") in ("warn", "fail")]
    if pbq_obs:
        samples = {o.get("sample", "unknown") for o in pbq_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in pbq_obs if o.get("sample") == sample]
            fid = next_id("FASTQ_QUAL")
            findings.append(_fastq_per_base_quality(sample_obs, fid, sample))
        seen_categories.add("per_base_quality")

    overrep_obs = collect_obs_by_category(observations, "overrepresented_sequences")
    overrep_obs = [o for o in overrep_obs if o.get("status") in ("warn", "fail")]
    if overrep_obs:
        samples = {o.get("sample", "unknown") for o in overrep_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in overrep_obs if o.get("sample") == sample]
            fid = next_id("FASTQ_OVERREP")
            findings.append(_fastq_overrepresented(sample_obs, fid, sample))
        seen_categories.add("overrepresented_sequences")

    if not seen_categories:
        findings.append(_fastq_all_pass(next_id("FASTQ_PASS")))

    return InterpretationResult(
        interpretation_version=INTERPRETATION_VERSION,
        workflow="fastq-qc",
        scope=(
            "FASTQ-level QC interpretation based on FastQC module flags. "
            "Covers adapter content, per-base quality, GC content, and overrepresented sequences. "
            "Does not cover alignment, variant calling, or downstream biological analysis."
        ),
        limitations=_FASTQ_LIMITATIONS,
        findings=findings,
        hypotheses=hypotheses,
        validation_recommendations=[
            "Review all findings with knowledge of the organism, protocol, and experimental design.",
            "Do not filter samples or reads based on QC flags alone.",
            "Proceed to alignment-level QC after FASTQ-level remediation.",
            "Interpret findings in the context of the full batch, not individual samples.",
        ],
        safety_flags=[],
        clinical_claims_allowed=False,
    )


def generate_variant_interpretation(
    observations: list[dict],
    decisions: list[dict],
) -> InterpretationResult:
    findings: list[Finding] = []
    hypotheses: list[Hypothesis] = []
    finding_counter = 0

    def next_id(prefix: str) -> str:
        nonlocal finding_counter
        finding_counter += 1
        return f"{prefix}_F{finding_counter:03d}"

    low_map_obs = collect_obs_by_category(observations, "alignment")
    low_map_obs = [o for o in low_map_obs if o.get("status") in ("warn", "fail")]
    if low_map_obs:
        samples = {o.get("sample", "unknown") for o in low_map_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in low_map_obs if o.get("sample") == sample]
            fid = next_id("VAR_MAP")
            f, h = _variant_low_mapping(sample_obs, fid, sample)
            findings.append(f)
            hypotheses.append(h)

    zero_vcf_obs = collect_obs_by_category(observations, "vcf_content", status="fail")
    if zero_vcf_obs:
        samples = {o.get("sample", "unknown") for o in zero_vcf_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in zero_vcf_obs if o.get("sample") == sample]
            fid = next_id("VAR_VCF0")
            findings.append(_variant_zero_vcf(sample_obs, fid, sample))

    contig_obs = collect_obs_by_category(observations, "contig_coverage", status="warn")
    if contig_obs:
        samples = {o.get("sample", "unknown") for o in contig_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in contig_obs if o.get("sample") == sample]
            fid = next_id("VAR_CONTIG")
            findings.append(_variant_zero_read_contigs(sample_obs, fid, sample))

    cov_obs = collect_obs_by_category(observations, "coverage", status="warn")
    if cov_obs:
        samples = {o.get("sample", "unknown") for o in cov_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in cov_obs if o.get("sample") == sample]
            fid = next_id("VAR_COV")
            findings.append(_variant_low_coverage(sample_obs, fid, sample))

    if not findings:
        findings.append(_variant_no_issues(next_id("VAR_PASS")))

    return InterpretationResult(
        interpretation_version=INTERPRETATION_VERSION,
        workflow="variant-qc",
        scope=(
            "Variant QC interpretation based on samtools flagstat, bcftools stats, and mosdepth outputs. "
            "Covers alignment rate, VCF record count, contig coverage, and mean sequencing depth. "
            "Does not cover variant annotation, functional impact, population frequency, or clinical significance."
        ),
        limitations=_VARIANT_LIMITATIONS,
        findings=findings,
        hypotheses=hypotheses,
        validation_recommendations=[
            "Review all findings with knowledge of the reference genome, assay design, and sample context.",
            "Do not infer variant significance from QC metrics alone.",
            "Proceed to variant annotation after QC review.",
            "Compare metrics to positive-control samples from the same experiment where available.",
        ],
        safety_flags=[],
        clinical_claims_allowed=False,
    )


_RNASEQ_LIMITATIONS = [
    "RNA-seq QC metrics are technical quality measures. They cannot establish gene-level biological causation.",
    "rRNA fraction thresholds are protocol-dependent; poly-A selected libraries behave differently from ribo-depleted ones.",
    "Gene-body coverage skew reflects both RNA integrity and native transcript biology and cannot be disambiguated by QC alone.",
    "Intronic/intergenic mapping rates depend heavily on the reference annotation completeness.",
    "These interpretations are scaffolds for human review, not conclusions.",
]

_RRNA_WARN_THRESHOLD = 0.10
_RRNA_FAIL_THRESHOLD = 0.20
_INTRONIC_WARN_THRESHOLD = 0.15
_INTERGENIC_WARN_THRESHOLD = 0.10
_UNMAPPED_WARN_THRESHOLD = 0.30
_MULTIMAPPED_WARN_THRESHOLD = 0.25


def _rnaseq_high_rrna(obs_list: list[dict], finding_id: str, sample: str) -> tuple[Finding, Hypothesis]:
    finding = Finding(
        finding_id=finding_id,
        workflow="rna-seq-qc",
        sample=sample,
        observation="Elevated rRNA fraction detected in RNA-seq library.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "rnaseq:rrna_fraction",
        technical_explanations=[
            "Incomplete enzymatic rRNA depletion (e.g., RiboZero, RiboMinus) leaving residual ribosomal sequences.",
            "Poly-A selection failure or sub-optimal oligo-dT capture efficiency producing rRNA carry-over.",
            "RNA degradation prior to library preparation reducing the poly-A-tail-bearing fraction relative to rRNA.",
            "Cross-contamination during library preparation from high-rRNA samples processed in parallel.",
        ],
        plausible_biological_explanations=[
            "Certain highly metabolic tissue types (liver, heart, activated immune cells) have elevated "
            "baseline ribosomal activity, producing naturally higher rRNA abundance prior to depletion.",
            "Proliferating cell populations have high ribosomal biogenesis rates that can exceed "
            "standard depletion capacity.",
        ],
        metadata_needed=[
            "Library preparation protocol (ribo-depletion kit vs. poly-A selection) and lot numbers.",
            "RNA Integrity Number (RIN) or RNA Quality Number (RQN) from capillary electrophoresis.",
            "Tissue type and expected metabolic state of the source material.",
            "Whether multiple samples from the same batch show concordant rRNA fractions.",
        ],
        recommended_validation=[
            "Review ribo-depletion or poly-A selection efficiency metrics from the library QC step.",
            "Compare rRNA fraction across all samples in the batch to identify batch-level vs. sample-level issues.",
            "Calculate effective mRNA sequencing depth (total reads minus rRNA reads) to assess usable coverage.",
            "If rRNA fraction exceeds 0.20, evaluate whether remaining mRNA coverage meets minimum depth requirements.",
        ],
        recommended_action=(
            "Do not filter or discard the sample based on rRNA fraction alone. "
            "Assess whether effective mRNA coverage is sufficient for the downstream analysis. "
            "If coverage is adequate, proceed with the sample while flagging the rRNA metric in the report."
        ),
        confidence="medium",
        should_filter=False,
        should_preserve_until_review=True,
    )

    hypothesis = Hypothesis(
        hypothesis_id=f"{finding_id}_H1",
        statement=(
            "Elevated rRNA levels suggest either lower technical library preparation efficiency "
            "(incomplete depletion or poly-A capture failure) or a biologically high baseline "
            "ribosomal activity in the source tissue or cell type. "
            "Primary recommendation is technical validation via effective mRNA coverage analysis "
            "before drawing biological conclusions."
        ),
        supporting_observations=obs_message_list(obs_list),
        alternative_explanations=[
            "Technical: ribo-depletion kit failure, poly-A selection inefficiency, RNA degradation, batch contamination.",
            "Biological: elevated ribosomal biogenesis in metabolically active or proliferating cell populations.",
        ],
        validation_steps=[
            "Compare rRNA fraction to manufacturer-specified depletion efficiency benchmarks.",
            "Obtain RIN/RQN values from pre-library RNA QC to distinguish degradation from depletion failure.",
            "Check whether effective mRNA depth (total reads minus rRNA reads) meets analysis requirements.",
        ],
        confidence="medium",
        interpretation_type="ambiguous",
        clinical_claim=False,
    )

    return finding, hypothesis


def _rnaseq_gene_body_bias(obs_list: list[dict], finding_id: str, sample: str) -> tuple[Finding, Hypothesis]:
    finding = Finding(
        finding_id=finding_id,
        workflow="rna-seq-qc",
        sample=sample,
        observation="Severe gene-body coverage skew or 3'/5' bias detected in RNA-seq library.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "rnaseq:gene_body_coverage",
        technical_explanations=[
            "Input RNA degradation (low RIN/RQN) causing fragmentation biased towards the 3' end of transcripts.",
            "Incomplete reverse transcription elongation truncating cDNA synthesis before the 5' end.",
            "Magnetic bead-based poly-A selection introducing 3' enrichment bias in the captured fragment pool.",
            "Priming artifacts from random hexamer or oligo-dT primers producing non-uniform coverage.",
        ],
        plausible_biological_explanations=[
            "Differential transcript stability: some transcripts are inherently 3'-biased due to regulatory "
            "structures in 5' UTRs that reduce ribosome engagement and increase 5' end turnover.",
            "Alternative polyadenylation producing transcript isoforms with different effective 3' ends, "
            "creating apparent 3' bias when aligned to the canonical transcript model.",
            "Extensive 3' UTR length variation across expressed gene sets in specific cell states or "
            "developmental stages can shift the apparent coverage distribution.",
        ],
        metadata_needed=[
            "RNA Integrity Number (RIN) or RNA Quality Number (RQN) from capillary electrophoresis "
            "of the input RNA prior to library construction.",
            "Library preparation protocol (strand-specific, poly-A, ribo-depleted, total RNA).",
            "Read length and whether paired-end sequencing was used.",
            "Whether the bias is consistent across samples or sample-specific.",
        ],
        recommended_validation=[
            "Cross-reference with pre-library RNA integrity profiles (RIN/RQN) — a low RIN strongly "
            "implicates technical degradation as the primary cause.",
            "Inspect gene-body coverage plots (e.g., RSeQC or RNA-SeQC output) across multiple "
            "housekeeping genes to distinguish systematic vs. gene-specific bias.",
            "Check whether 5' or 3' bias is concordant across all samples or enriched in specific batches.",
            "For poly-A selected libraries, compare to a ribo-depleted control if available.",
        ],
        recommended_action=(
            "Do not filter the sample based on gene-body bias alone. "
            "If RIN/RQN indicates degradation, flag the sample and assess whether "
            "downstream quantification is robust to 3' bias (e.g., using 3'-end quantification methods). "
            "If RIN is adequate, investigate alternative polyadenylation as a biological signal."
        ),
        confidence="high",
        should_filter=False,
        should_preserve_until_review=True,
    )

    hypothesis = Hypothesis(
        hypothesis_id=f"{finding_id}_H1",
        statement=(
            "Observed gene-body coverage skew indicates either technical RNA degradation "
            "(detectable via pre-library RIN/RQN values) or native alternative polyadenylation "
            "patterns in the source material. Verification requires inspection of pre-library "
            "RNA integrity numbers (RIN) and comparison across the batch."
        ),
        supporting_observations=obs_message_list(obs_list),
        alternative_explanations=[
            "Technical: RNA degradation (low RIN), incomplete reverse transcription, poly-A selection bias.",
            "Biological: alternative polyadenylation, differential transcript stability, 3' UTR length variation.",
        ],
        validation_steps=[
            "Retrieve RIN/RQN values from pre-library RNA QC records.",
            "Plot gene-body coverage for a panel of housekeeping genes with known uniform coverage.",
            "Compare coverage profile to a positive-control sample with known good RNA integrity.",
        ],
        confidence="high",
        interpretation_type="ambiguous",
        clinical_claim=False,
    )

    return finding, hypothesis


def _rnaseq_high_intronic(obs_list: list[dict], finding_id: str, sample: str) -> tuple[Finding, Hypothesis]:
    finding = Finding(
        finding_id=finding_id,
        workflow="rna-seq-qc",
        sample=sample,
        observation="High intronic or intergenic mapping rate detected in RNA-seq library.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "rnaseq:intronic_mapping",
        technical_explanations=[
            "Genomic DNA (gDNA) contamination from missing or exhausted DNase digestion during RNA extraction.",
            "Incomplete DNase I digestion leaving residual gDNA co-purified with the RNA fraction.",
            "Reference annotation incompleteness causing exonic reads to be misclassified as intronic.",
            "Aligner soft-clipping artefacts at exon-intron junctions inflating intronic counts.",
        ],
        plausible_biological_explanations=[
            "High levels of unannotated transcription (enhancer RNAs, novel lncRNAs, unannotated exons) "
            "being captured and aligning to intronic or intergenic regions.",
            "Nascent pre-mRNA capture in nuclear or total RNA fractions where intron retention is a "
            "regulated biological event (e.g., neuronal differentiation, immune activation).",
            "Systematic intron retention events in specific cell states, which are functionally relevant "
            "post-transcriptional regulatory mechanisms.",
        ],
        metadata_needed=[
            "Whether DNase digestion was applied and the DNase I concentration/incubation conditions.",
            "Library type (total RNA, poly-A, ribo-depleted, nuclear RNA fraction).",
            "Reference annotation version and completeness for the target organism.",
            "Proportion of intergenic vs. intronic reads to distinguish gDNA from nascent transcription.",
        ],
        recommended_validation=[
            "Introduce or verify inline DNase digestion protocols and confirm complete DNA removal "
            "via gel or capillary electrophoresis of the RNA fraction.",
            "Inspect the ratio of intronic to intergenic reads: predominantly intronic may suggest "
            "intron retention, while predominantly intergenic suggests gDNA or annotation gaps.",
            "Evaluate using a specialized intron-retention analysis tool (e.g., IRFinder, rMATS) "
            "to distinguish regulated retention from technical noise.",
            "Compare to samples prepared without DNase digestion as a positive control if available.",
        ],
        recommended_action=(
            "Investigate the DNase digestion protocol before discarding the sample. "
            "If gDNA contamination is confirmed, the sample may require re-preparation. "
            "If intron retention is the primary signal, preserve the sample and apply appropriate "
            "intron-retention analysis tools downstream."
        ),
        confidence="medium",
        should_filter=False,
        should_preserve_until_review=True,
    )

    hypothesis = Hypothesis(
        hypothesis_id=f"{finding_id}_H1",
        statement=(
            "Elevated intronic or intergenic mapping suggests either unresolved genomic DNA presence "
            "(technical: absent or insufficient DNase digestion) or extensive intron retention / "
            "unannotated transcriptional activity (biological: regulated intron retention, nascent "
            "pre-mRNA capture, or novel transcripts). These cannot be distinguished without "
            "additional protocol and annotation metadata."
        ),
        supporting_observations=obs_message_list(obs_list),
        alternative_explanations=[
            "Technical: gDNA contamination from missing or exhausted DNase I digestion.",
            "Biological: regulated intron retention, nascent pre-mRNA, unannotated lncRNAs or enhancer RNAs.",
        ],
        validation_steps=[
            "Confirm DNase digestion conditions from the extraction protocol record.",
            "Compute the intronic-to-intergenic read ratio to separate retention from gDNA signal.",
            "Run IRFinder or equivalent tool to quantify intron retention events at the gene level.",
        ],
        confidence="medium",
        interpretation_type="ambiguous",
        clinical_claim=False,
    )

    return finding, hypothesis


def _rnaseq_low_mapping(obs_list: list[dict], finding_id: str, sample: str) -> tuple[Finding, Hypothesis]:
    finding = Finding(
        finding_id=finding_id,
        workflow="rna-seq-qc",
        sample=sample,
        observation="Low unique mapping rate or high multimapping rate detected in RNA-seq library.",
        evidence_source=obs_to_evidence_source(obs_list[0]) if obs_list else "rnaseq:mapping_rate",
        technical_explanations=[
            "Outdated or incomplete reference transcriptome annotation missing novel splice junctions or isoforms.",
            "Read lengths too short to uniquely resolve highly repetitive genomic regions or "
            "gene families with high sequence similarity.",
            "Sample cross-contamination introducing reads from a different species or sample type.",
            "Mismatched genome build between the library preparation reference and the alignment reference.",
            "Aligner parameter settings too stringent for the sequencing quality or species divergence.",
        ],
        plausible_biological_explanations=[
            "Non-model organism features with limited reference annotation leading to unmapped reads "
            "that represent genuine transcripts not yet in the reference.",
            "Hyper-activation of repetitive element or transposon landscapes (e.g., in stress, "
            "senescence, or specific developmental stages) producing reads that multi-map to "
            "transposon families.",
            "Highly homologous gene family expansions (e.g., olfactory receptors, immunoglobulin genes) "
            "where biological multi-mapping is expected.",
        ],
        metadata_needed=[
            "Reference genome and annotation version used for alignment.",
            "Alignment software and key parameter settings (mismatch tolerance, multi-map threshold).",
            "Whether the sample is from a model organism with a well-annotated reference.",
            "Taxonomic composition of unmapped reads (contamination screen).",
        ],
        recommended_validation=[
            "Extract unmapped reads and submit to taxonomic classification (e.g., Kraken2, FastQ Screen) "
            "to screen for exogenous contamination.",
            "Verify alignment reference genome build and annotation version match the library preparation.",
            "Re-align with updated annotation (e.g., latest Ensembl or GENCODE release) if annotation "
            "completeness may be limiting.",
            "For multimapper analysis, inspect the distribution of multi-mapping loci to identify "
            "whether transposable elements or gene families are the primary contributors.",
        ],
        recommended_action=(
            "Investigate the cause of low unique mapping before proceeding to quantification. "
            "Do not discard the sample without first screening unmapped reads for contamination "
            "and verifying the alignment reference."
        ),
        confidence="high",
        should_filter=False,
        should_preserve_until_review=True,
    )

    hypothesis = Hypothesis(
        hypothesis_id=f"{finding_id}_H1",
        statement=(
            "Low unique mapping rate points to either technical reference mismatch or misalignment "
            "(outdated annotation, parameter mismatch, contamination) or structural biological features "
            "(repetitive element activation, highly homologous gene families, non-model organism sequences). "
            "Taxonomic classification of unmapped reads is the most informative first validation step."
        ),
        supporting_observations=obs_message_list(obs_list),
        alternative_explanations=[
            "Technical: reference mismatch, outdated annotation, aligner misconfiguration, "
            "sample cross-contamination.",
            "Biological: transposon activation, gene family expansions, non-model organism features, "
            "novel splice junctions absent from the reference.",
        ],
        validation_steps=[
            "Run Kraken2 or FastQ Screen on a representative subset of unmapped reads.",
            "Re-align a test batch using the most recent genome/annotation release.",
            "Inspect multimapper loci to identify whether transposable elements are the primary source.",
        ],
        confidence="high",
        interpretation_type="ambiguous",
        clinical_claim=False,
    )

    return finding, hypothesis


def _rnaseq_all_pass(finding_id: str) -> Finding:
    return Finding(
        finding_id=finding_id,
        workflow="rna-seq-qc",
        sample="all",
        observation=(
            "No critical RNA-seq QC flags detected in tracked metrics "
            "(rRNA fraction, gene-body coverage, mapping rate, intronic/intergenic rate)."
        ),
        evidence_source="rnaseq_qc:overall",
        technical_explanations=[
            "All evaluated QC metrics are within acceptable ranges for the tracked modules.",
        ],
        plausible_biological_explanations=[
            "RNA-seq QC passing does not confirm biological accuracy of quantification or "
            "the validity of downstream differential expression analysis.",
            "Batch effects, normalization choices, and annotation completeness require separate evaluation.",
        ],
        metadata_needed=[
            "Downstream QC metrics (duplication rate, gene detection rate, normalization metrics) "
            "for a complete assessment.",
        ],
        recommended_validation=[
            "Proceed to alignment-level and quantification QC.",
            "Review gene detection rate and sample clustering before downstream analysis.",
        ],
        recommended_action=(
            "No RNA-seq QC intervention required based on tracked metrics. "
            "Proceed to quantification and downstream analysis QC."
        ),
        confidence="high",
        should_filter=False,
        should_preserve_until_review=False,
    )


def _obs_numeric(obs: dict, key: str) -> float | None:
    """Extract a numeric metric value from an observation's evidence dict."""
    val = obs.get("evidence", {}).get(key)
    try:
        return float(val) if val is not None else None
    except (TypeError, ValueError):
        return None


def generate_rnaseq_interpretation(
    observations: list[dict],
    decisions: list[dict],
) -> InterpretationResult:
    findings: list[Finding] = []
    hypotheses: list[Hypothesis] = []
    finding_counter = 0

    def next_id(prefix: str) -> str:
        nonlocal finding_counter
        finding_counter += 1
        return f"{prefix}_F{finding_counter:03d}"

    seen_categories: set[str] = set()

    # --- rRNA fraction ---
    rrna_obs = [
        o for o in observations
        if o.get("category") == "rrna_fraction" and o.get("status") in ("warn", "fail")
    ]
    # Also trigger on numeric evidence above threshold
    rrna_numeric = [
        o for o in observations
        if o.get("category") == "rrna_fraction" and o.get("status") not in ("warn", "fail")
        and (
            (_obs_numeric(o, "rrna_fraction") or 0.0) >= _RRNA_WARN_THRESHOLD
            or (_obs_numeric(o, "rrna_fraction_value") or 0.0) >= _RRNA_WARN_THRESHOLD
        )
    ]
    rrna_obs = rrna_obs or rrna_numeric
    if rrna_obs:
        samples = {o.get("sample", "unknown") for o in rrna_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in rrna_obs if o.get("sample") == sample]
            fid = next_id("RNASEQ_RRNA")
            f, h = _rnaseq_high_rrna(sample_obs, fid, sample)
            findings.append(f)
            hypotheses.append(h)
        seen_categories.add("rrna_fraction")

    # --- Gene-body coverage / 3' bias ---
    bias_obs = [
        o for o in observations
        if o.get("category") in ("gene_body_coverage", "three_prime_bias", "five_prime_bias")
        and o.get("status") in ("warn", "fail")
    ]
    if bias_obs:
        samples = {o.get("sample", "unknown") for o in bias_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in bias_obs if o.get("sample") == sample]
            fid = next_id("RNASEQ_BIAS")
            f, h = _rnaseq_gene_body_bias(sample_obs, fid, sample)
            findings.append(f)
            hypotheses.append(h)
        seen_categories.add("gene_body_coverage")

    # --- Intronic / intergenic mapping ---
    intronic_obs = [
        o for o in observations
        if o.get("category") in ("intronic_mapping", "intergenic_mapping")
        and o.get("status") in ("warn", "fail")
    ]
    intronic_numeric = [
        o for o in observations
        if o.get("category") in ("intronic_mapping", "intergenic_mapping")
        and o.get("status") not in ("warn", "fail")
        and (
            (_obs_numeric(o, "intronic_mapping_rate") or 0.0) > _INTRONIC_WARN_THRESHOLD
            or (_obs_numeric(o, "intergenic_mapping_rate") or 0.0) > _INTERGENIC_WARN_THRESHOLD
        )
    ]
    intronic_obs = intronic_obs or intronic_numeric
    if intronic_obs:
        samples = {o.get("sample", "unknown") for o in intronic_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in intronic_obs if o.get("sample") == sample]
            fid = next_id("RNASEQ_INTRONIC")
            f, h = _rnaseq_high_intronic(sample_obs, fid, sample)
            findings.append(f)
            hypotheses.append(h)
        seen_categories.add("intronic_mapping")

    # --- Low mapping / high multimappers ---
    mapping_obs = [
        o for o in observations
        if o.get("category") in ("mapping_rate", "multimapping_rate", "unmapped_rate")
        and o.get("status") in ("warn", "fail")
    ]
    mapping_numeric = [
        o for o in observations
        if o.get("category") in ("mapping_rate", "multimapping_rate", "unmapped_rate")
        and o.get("status") not in ("warn", "fail")
        and (
            (_obs_numeric(o, "unmapped_reads_rate") or 0.0) > _UNMAPPED_WARN_THRESHOLD
            or (_obs_numeric(o, "multimapped_reads_rate") or 0.0) > _MULTIMAPPED_WARN_THRESHOLD
        )
    ]
    mapping_obs = mapping_obs or mapping_numeric
    if mapping_obs:
        samples = {o.get("sample", "unknown") for o in mapping_obs}
        for sample in sorted(samples):
            sample_obs = [o for o in mapping_obs if o.get("sample") == sample]
            fid = next_id("RNASEQ_MAP")
            f, h = _rnaseq_low_mapping(sample_obs, fid, sample)
            findings.append(f)
            hypotheses.append(h)
        seen_categories.add("mapping_rate")

    if not seen_categories:
        findings.append(_rnaseq_all_pass(next_id("RNASEQ_PASS")))

    return InterpretationResult(
        interpretation_version=INTERPRETATION_VERSION,
        workflow="rna-seq-qc",
        scope=(
            "RNA-seq QC interpretation based on rRNA fraction, gene-body coverage, "
            "intronic/intergenic mapping rate, and unique mapping rate. "
            "Does not cover differential expression, normalization, batch correction, "
            "or downstream biological analysis."
        ),
        limitations=_RNASEQ_LIMITATIONS,
        findings=findings,
        hypotheses=hypotheses,
        validation_recommendations=[
            "Review all findings with knowledge of the library preparation protocol, "
            "RNA source material quality (RIN/RQN), and organism annotation completeness.",
            "Do not filter samples based on individual QC metrics without assessing "
            "effective mRNA coverage and downstream usability.",
            "Compare all metrics across the full batch before making per-sample decisions.",
            "Proceed to quantification-level QC (gene detection rate, sample clustering) "
            "after FASTQ and alignment QC are resolved.",
        ],
        safety_flags=[],
        clinical_claims_allowed=False,
    )


def generate_interpretation(
    workflow: str,
    observations: list[dict] | None = None,
    decisions: list[dict] | None = None,
    run_result: dict | None = None,
    metadata: dict | None = None,
) -> dict:
    """
    Generate a deterministic biological interpretation scaffold.

    No LLM is called. This function applies rule-based logic to map
    QC observations to structured findings and testable hypotheses.

    Returns a plain JSON-serializable dict.
    """
    obs = observations or []
    dec = decisions or []

    if workflow == "fastq-qc":
        result = generate_fastq_interpretation(obs, dec)
    elif workflow == "variant-qc":
        result = generate_variant_interpretation(obs, dec)
    elif workflow in ("rna-seq-qc", "rnaseq"):
        result = generate_rnaseq_interpretation(obs, dec)
    else:
        result = InterpretationResult(
            interpretation_version=INTERPRETATION_VERSION,
            workflow=workflow,
            scope=f"Interpretation not yet implemented for workflow: {workflow}",
            limitations=[
                f"The interpretation engine does not have rules for workflow '{workflow}'.",
                "Supported workflows: fastq-qc, variant-qc, rna-seq-qc.",
            ],
            findings=[],
            hypotheses=[],
            validation_recommendations=[
                "Use the biology-interpretation-reviewer Claude Code agent "
                "to review outputs manually for unsupported workflows."
            ],
            safety_flags=[f"No deterministic rules available for workflow: {workflow}"],
            clinical_claims_allowed=False,
        )

    return result.to_dict()
