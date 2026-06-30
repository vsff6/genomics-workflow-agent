"""Deterministic decision rules for variant/BAM/VCF QC observations."""
from __future__ import annotations

from typing import Any

from genomics_workflow_agent.agent.state import Decision, Observation, RecommendedAction

CLINICAL_DISCLAIMER = (
    "Variant QC metrics are technical quality measures, not clinical findings. "
    "No variant is assessed for clinical significance. "
    "No interpretation of variant impact is made here. "
    "Variant interpretation requires validated clinical annotation and expert review."
)

LIMITATIONS = [
    "Mapped read percentage thresholds depend on organism, library type, and reference completeness.",
    "Low variant counts may reflect caller settings, aggressive filtering, or intentionally low-VAF targets.",
    "Coverage metrics from mosdepth require BAM indexing and may reflect capture target design.",
    "Ti/Tv ratio varies by variant caller, target region, and filtering strategy.",
    "BAM completeness cannot be verified without a reference manifest.",
    "All decisions are deterministic rules applied to summary statistics. "
    "Biological and clinical interpretation requires domain expertise.",
]

LOW_MAPPING_THRESHOLD = 80.0
LOW_COVERAGE_THRESHOLD = 10.0
LOW_VARIANT_COUNT_THRESHOLD = 100


def _obs(
    sample: str,
    category: str,
    status: str,
    severity: str,
    message: str,
    evidence: dict | None = None,
    suggested_action: str = "",
) -> Observation:
    return Observation(
        source="variant_decision_engine",
        sample=sample,
        category=category,
        status=status,
        severity=severity,
        message=message,
        evidence=evidence or {},
        suggested_action=suggested_action,
    )


def evaluate_variant_qc_results(
    samtools_results: dict[str, Any],
    bcftools_results: list[dict[str, Any]],
    mosdepth_results: list[dict[str, Any]],
    execute_allowed: bool = False,
) -> dict[str, Any]:
    """
    Apply deterministic rules to parsed samtools/bcftools/mosdepth outputs.

    Returns observations, decisions, recommended_actions, warnings, limitations,
    and the non-clinical disclaimer.
    """
    observations: list[Observation] = []
    decisions: list[Decision] = []
    recommended_actions: list[RecommendedAction] = []
    warnings: list[str] = []

    no_data = (
        not samtools_results.get("flagstat")
        and not samtools_results.get("idxstats")
        and not samtools_results.get("stats")
        and not bcftools_results
        and not mosdepth_results
    )

    if no_data:
        observations.append(_obs(
            sample="all",
            category="qc_output",
            status="missing",
            severity="warning",
            message=(
                "No samtools/bcftools/mosdepth outputs found. "
                "QC interpretation is not possible in dry-run mode."
            ),
            suggested_action="Re-run with --execute to generate QC outputs.",
        ))
        recommended_actions.append(RecommendedAction(
            action="Run variant-qc with --execute to generate QC outputs",
            priority="high",
            reason="No QC data available for evaluation",
            requires_execute=True,
            requires_external_tool="samtools, bcftools",
        ))
        return _package(observations, decisions, recommended_actions, warnings)

    low_mapping_samples: list[str] = []
    zero_vcf_samples: list[str] = []

    # Flagstat evaluation
    for fs in samtools_results.get("flagstat", []):
        sample = fs.get("sample", "unknown")
        if not fs.get("parse_ok"):
            for err in fs.get("errors", []):
                observations.append(_obs(
                    sample=sample,
                    category="flagstat_parse_error",
                    status="error",
                    severity="warning",
                    message=f"Could not parse flagstat output: {err}",
                ))
            continue

        total = fs.get("total_reads")
        mapped_pct = fs.get("mapped_pct")

        if total is not None and total == 0:
            observations.append(_obs(
                sample=sample,
                category="alignment",
                status="fail",
                severity="critical",
                message=(
                    f"{sample}: BAM contains zero reads. "
                    "File may be empty, corrupt, or from the wrong sample."
                ),
                evidence={"total_reads": 0},
                suggested_action="Verify BAM file integrity. Check alignment pipeline logs.",
            ))
            low_mapping_samples.append(sample)
        elif mapped_pct is not None and mapped_pct < LOW_MAPPING_THRESHOLD:
            observations.append(_obs(
                sample=sample,
                category="alignment",
                status="warn",
                severity="critical",
                message=(
                    f"{sample}: Mapped reads = {mapped_pct:.1f}% "
                    f"(threshold: {LOW_MAPPING_THRESHOLD}%). "
                    "Possible causes: wrong reference genome, contamination, "
                    "low-quality library, or species mismatch. "
                    "Do not proceed to variant calling without investigation."
                ),
                evidence={"mapped_pct": mapped_pct, "threshold_pct": LOW_MAPPING_THRESHOLD},
                suggested_action=(
                    "Check reference/sample compatibility. "
                    "Review alignment log. Consider re-alignment."
                ),
            ))
            low_mapping_samples.append(sample)
        elif mapped_pct is not None:
            observations.append(_obs(
                sample=sample,
                category="alignment",
                status="pass",
                severity="info",
                message=f"{sample}: Alignment looks normal - {mapped_pct:.1f}% reads mapped.",
                evidence={"mapped_pct": mapped_pct},
            ))

    # Idxstats: check for many zero-read contigs
    for ix in samtools_results.get("idxstats", []):
        sample = ix.get("sample", "unknown")
        if not ix.get("parse_ok"):
            continue
        zero_contigs = ix.get("zero_read_contigs", [])
        total_contigs = len(ix.get("contigs", []))
        if total_contigs > 0 and len(zero_contigs) > total_contigs * 0.5:
            observations.append(_obs(
                sample=sample,
                category="contig_coverage",
                status="warn",
                severity="warning",
                message=(
                    f"{sample}: {len(zero_contigs)}/{total_contigs} contigs have zero mapped reads. "
                    "This may indicate reference/sample mismatch, decoy contigs, "
                    "or an unusual reference assembly."
                ),
                evidence={
                    "zero_read_contigs_count": len(zero_contigs),
                    "total_contigs": total_contigs,
                },
                suggested_action=(
                    "Review contig list. "
                    "Confirm reference build matches library preparation."
                ),
            ))

    # bcftools stats evaluation
    for bc in bcftools_results:
        sample = bc.get("sample", "unknown")
        if not bc.get("parse_ok"):
            observations.append(_obs(
                sample=sample,
                category="bcftools_unavailable",
                status="missing",
                severity="warning",
                message=(
                    f"{sample}: bcftools stats output could not be parsed. "
                    "VCF interpretation is not available."
                ),
                suggested_action="Re-run bcftools stats and ensure output is captured.",
            ))
            continue

        n_records = bc.get("n_records")
        n_snps = bc.get("n_snps")

        if n_records is not None and n_records == 0:
            observations.append(_obs(
                sample=sample,
                category="vcf_content",
                status="fail",
                severity="critical",
                message=(
                    f"{sample}: VCF contains zero variant records. "
                    "Possible causes: incorrect variant caller invocation, empty input BAM, "
                    "wrong reference, aggressive pre-filtering, "
                    "or sample with no variants in target region."
                ),
                evidence={"n_records": 0},
                suggested_action=(
                    "Verify variant caller command and input BAM. "
                    "Check whether target region is correct."
                ),
            ))
            zero_vcf_samples.append(sample)
        elif n_records is not None and n_records < LOW_VARIANT_COUNT_THRESHOLD:
            observations.append(_obs(
                sample=sample,
                category="vcf_content",
                status="warn",
                severity="warning",
                message=(
                    f"{sample}: VCF contains very few variants "
                    f"({n_records} records, {n_snps} SNPs). "
                    "This may reflect targeted panel design, aggressive filtering, "
                    "or data quality issues."
                ),
                evidence={"n_records": n_records, "n_snps": n_snps},
                suggested_action=(
                    "Review variant caller settings and filtering. "
                    "Compare to expected variant count for this assay."
                ),
            ))
        elif n_records is not None:
            observations.append(_obs(
                sample=sample,
                category="vcf_content",
                status="pass",
                severity="info",
                message=(
                    f"{sample}: VCF contains {n_records} records "
                    f"({n_snps} SNPs, {bc.get('n_indels')} indels, "
                    f"Ts/Tv={bc.get('ts_tv')})."
                ),
                evidence={
                    "n_records": n_records,
                    "n_snps": n_snps,
                    "n_indels": bc.get("n_indels"),
                    "ts_tv": bc.get("ts_tv"),
                },
            ))

    # mosdepth coverage evaluation
    for md in mosdepth_results:
        sample = md.get("sample", "unknown")
        if not md.get("parse_ok"):
            continue
        cov = md.get("mean_coverage")
        if cov is not None and cov < LOW_COVERAGE_THRESHOLD:
            observations.append(_obs(
                sample=sample,
                category="coverage",
                status="warn",
                severity="warning",
                message=(
                    f"{sample}: Mean coverage = {cov:.1f}x "
                    f"(below {LOW_COVERAGE_THRESHOLD}x). "
                    "Low coverage affects variant calling sensitivity. "
                    "Do not make diagnostic claims based on low-coverage data."
                ),
                evidence={"mean_coverage": cov, "threshold": LOW_COVERAGE_THRESHOLD},
                suggested_action=(
                    "Review coverage distribution. "
                    "Check library quality and target capture efficiency."
                ),
            ))
        elif cov is not None:
            observations.append(_obs(
                sample=sample,
                category="coverage",
                status="pass",
                severity="info",
                message=f"{sample}: Mean coverage = {cov:.1f}x.",
                evidence={"mean_coverage": cov},
            ))

    # Build decisions
    if low_mapping_samples:
        decisions.append(Decision(
            action="review_alignment",
            decision_type="review",
            reason="One or more samples have low or zero mapping rates.",
            evidence=[f"Low/zero mapping: {s}" for s in low_mapping_samples],
            confidence="high",
            execute_allowed=False,
            safety_notes=[
                "Do not proceed to variant calling without investigating mapping failures.",
            ],
        ))
        recommended_actions.append(RecommendedAction(
            action="Investigate alignment for flagged samples",
            priority="high",
            reason=f"Low mapping rate: {', '.join(low_mapping_samples)}",
            requires_execute=False,
        ))

    if zero_vcf_samples:
        decisions.append(Decision(
            action="review_vcf_output",
            decision_type="review",
            reason="One or more VCFs contain zero variant records.",
            evidence=[f"Zero records: {s}" for s in zero_vcf_samples],
            confidence="high",
            execute_allowed=False,
            safety_notes=[
                "A zero-variant VCF may indicate a pipeline problem. "
                "Do not interpret as a true negative without investigation.",
            ],
        ))
        recommended_actions.append(RecommendedAction(
            action="Investigate variant calling for zero-record VCF samples",
            priority="high",
            reason=f"Zero VCF records: {', '.join(zero_vcf_samples)}",
            requires_execute=False,
        ))

    has_any_parsed = (
        any(fs.get("parse_ok") for fs in samtools_results.get("flagstat", []))
        or any(bc.get("parse_ok") for bc in bcftools_results)
    )

    if not low_mapping_samples and not zero_vcf_samples and has_any_parsed:
        decisions.append(Decision(
            action="no_critical_issues",
            decision_type="accept",
            reason="No critical alignment or VCF issues found in tracked metrics.",
            evidence=["Mapping rates and VCF record counts within expected ranges"],
            confidence="medium",
            execute_allowed=False,
            safety_notes=[
                "QC passing does not guarantee variant accuracy or completeness.",
                "Ti/Tv ratio, het/hom ratio, and coverage uniformity require expert interpretation.",
                "Do not interpret variants clinically without validated annotation and expert review.",
            ],
        ))

    return _package(observations, decisions, recommended_actions, warnings)


def _package(
    observations: list[Observation],
    decisions: list[Decision],
    recommended_actions: list[RecommendedAction],
    warnings: list[str],
) -> dict[str, Any]:
    return {
        "observations": observations,
        "decisions": decisions,
        "recommended_actions": recommended_actions,
        "warnings": warnings,
        "limitations": LIMITATIONS,
        "clinical_disclaimer": CLINICAL_DISCLAIMER,
    }
