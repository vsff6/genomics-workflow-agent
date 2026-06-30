from __future__ import annotations

from typing import Any

from genomics_workflow_agent.agent.state import Decision, Observation, RecommendedAction

CLINICAL_DISCLAIMER = (
    "FastQC/MultiQC results are technical QC metrics, not biological or clinical findings. "
    "No variant, diagnosis, or medical interpretation is possible from QC data alone."
)

LIMITATIONS = [
    "FastQC flags are defined for generic short-read sequencing - thresholds may not apply "
    "to all library types (amplicon, long-read, targeted panels, single-cell).",
    "Per-sequence GC content failures may reflect the organism's GC content, not contamination.",
    "High duplication rates may be expected for low-complexity libraries or targeted sequencing.",
    "Adapter content warnings depend on insert size - very short inserts are expected in some "
    "protocols (e.g., ATAC-seq, ChIP-seq).",
    "All decisions here are deterministic rules applied to QC summary flags. "
    "Biological interpretation requires domain expertise and sample metadata.",
]


def _obs(sample: str, category: str, status: str, severity: str, message: str,
         evidence: dict | None = None, suggested_action: str = "") -> Observation:
    return Observation(
        source="decision_engine",
        sample=sample,
        category=category,
        status=status,
        severity=severity,
        message=message,
        evidence=evidence or {},
        suggested_action=suggested_action,
    )


def evaluate_fastqc_results(
    parsed_results: list[dict[str, Any]],
    execute_allowed: bool = False,
) -> dict[str, Any]:
    """
    Apply deterministic rules to a list of parse_fastqc_* outputs.

    Returns observations, decisions, recommended_actions, warnings, limitations.
    """
    observations: list[Observation] = []
    decisions: list[Decision] = []
    recommended_actions: list[RecommendedAction] = []
    warnings: list[str] = []

    if not parsed_results:
        observations.append(_obs(
            sample="all",
            category="fastqc_output",
            status="missing",
            severity="warning",
            message="No FastQC results found. QC interpretation is not possible.",
            suggested_action="Run FastQC on input files with --execute.",
        ))
        recommended_actions.append(RecommendedAction(
            action="Run FastQC to generate QC data",
            priority="high",
            reason="No FastQC outputs available for evaluation",
            requires_execute=True,
            requires_external_tool="fastqc",
        ))
        return _package(observations, decisions, recommended_actions, warnings)

    trim_evidence: list[str] = []
    review_evidence: list[str] = []
    gc_warnings: list[str] = []

    for result in parsed_results:
        sample = result.get("sample", "unknown")

        if not result.get("parse_ok") or result.get("errors"):
            for err in result.get("errors", []):
                observations.append(_obs(
                    sample=sample,
                    category="parse_error",
                    status="error",
                    severity="warning",
                    message=f"Could not parse FastQC output: {err}",
                ))
            continue

        modules = result.get("modules", {})

        # Per base sequence quality
        pbq_status = _module_status(modules, "Per base sequence quality")
        if pbq_status == "fail":
            msg = (
                f"{sample}: Per base sequence quality FAIL. "
                "Quality drops may indicate sequencing issues, flowcell problems, "
                "or normal 3' degradation. Trimming low-quality bases is often appropriate."
            )
            observations.append(_obs(
                sample=sample,
                category="per_base_quality",
                status="fail",
                severity="critical",
                message=msg,
                evidence={"module": "Per base sequence quality", "fastqc_status": "fail"},
                suggested_action="Consider quality trimming. Review whether this is systematic or sample-specific.",
            ))
            trim_evidence.append(f"{sample}: Per base sequence quality FAIL")
            review_evidence.append(f"{sample}: Per base sequence quality FAIL")
        elif pbq_status == "warn":
            observations.append(_obs(
                sample=sample,
                category="per_base_quality",
                status="warn",
                severity="warning",
                message=f"{sample}: Per base sequence quality WARN. Consider trimming if downstream tools are quality-sensitive.",
                evidence={"module": "Per base sequence quality", "fastqc_status": "warn"},
                suggested_action="Consider quality trimming, especially for alignment-based workflows.",
            ))
            trim_evidence.append(f"{sample}: Per base sequence quality WARN")

        # Adapter content
        adapter_status = _module_status(modules, "Adapter Content")
        if adapter_status in ("warn", "fail"):
            observations.append(_obs(
                sample=sample,
                category="adapter_content",
                status=adapter_status,
                severity="critical" if adapter_status == "fail" else "warning",
                message=(
                    f"{sample}: Adapter Content {adapter_status.upper()}. "
                    "Adapter sequences are present in reads. "
                    "This is expected when insert size is shorter than read length. "
                    "Adapter trimming is recommended before alignment or assembly."
                ),
                evidence={"module": "Adapter Content", "fastqc_status": adapter_status},
                suggested_action="Run adapter trimming (fastp or cutadapt) before alignment.",
            ))
            trim_evidence.append(f"{sample}: Adapter Content {adapter_status.upper()}")

        # Overrepresented sequences
        overrep_status = _module_status(modules, "Overrepresented sequences")
        if overrep_status in ("warn", "fail"):
            observations.append(_obs(
                sample=sample,
                category="overrepresented_sequences",
                status=overrep_status,
                severity="warning",
                message=(
                    f"{sample}: Overrepresented sequences {overrep_status.upper()}. "
                    "Possible causes: adapter contamination, rRNA, PCR amplification bias, "
                    "or genuinely abundant transcripts. Review overrepresented sequence annotations."
                ),
                evidence={"module": "Overrepresented sequences", "fastqc_status": overrep_status},
                suggested_action="Inspect overrepresented sequences in the FastQC HTML report. Do not filter without identifying the source.",
            ))
            review_evidence.append(f"{sample}: Overrepresented sequences {overrep_status.upper()}")

        # Per sequence GC content
        gc_status = _module_status(modules, "Per sequence GC content")
        if gc_status == "fail":
            observations.append(_obs(
                sample=sample,
                category="gc_content",
                status="fail",
                severity="warning",
                message=(
                    f"{sample}: Per sequence GC content FAIL. "
                    "GC distribution deviates from the expected model. "
                    "Possible causes: species with unusual GC content, contamination, "
                    "PCR bias, or library-specific patterns. "
                    "Do not filter samples based on GC content without biological justification."
                ),
                evidence={"module": "Per sequence GC content", "fastqc_status": "fail"},
                suggested_action="Review GC distribution. Check species GC content. Do not auto-filter.",
            ))
            gc_warnings.append(f"{sample}: GC content FAIL - requires biological context review")
        elif gc_status == "warn":
            observations.append(_obs(
                sample=sample,
                category="gc_content",
                status="warn",
                severity="info",
                message=f"{sample}: Per sequence GC content WARN. Monitor across samples for batch effects.",
                evidence={"module": "Per sequence GC content", "fastqc_status": "warn"},
            ))

        # Summary observation for this sample
        tracked_fail = any(s == "fail" for s in [pbq_status, adapter_status, overrep_status, gc_status] if s)
        tracked_warn = any(s in ("warn", "fail") for s in [adapter_status, overrep_status] if s)
        if not tracked_fail and not tracked_warn:
            observations.append(_obs(
                sample=sample,
                category="overall_qc",
                status="pass",
                severity="info",
                message=(
                    f"{sample}: No trimming triggers found in tracked FastQC modules "
                    "(Per base quality, Adapter Content, Overrepresented sequences, GC content)."
                ),
                evidence={"tracked_modules_checked": [
                    "Per base sequence quality", "Adapter Content",
                    "Overrepresented sequences", "Per sequence GC content",
                ]},
            ))

    # Decisions
    if trim_evidence:
        decisions.append(Decision(
            action="trim_reads",
            decision_type="trim",
            reason="FastQC identified quality or adapter issues that trimming can address.",
            evidence=trim_evidence,
            confidence="high",
            execute_allowed=execute_allowed,
            safety_notes=[
                "Trimming modifies reads - always preserve original FASTQ files.",
                "Trimming parameters should match library type and downstream tool requirements.",
                "Aggressive trimming can remove genuine signal in some library types.",
            ],
        ))
        recommended_actions.append(RecommendedAction(
            action="Run adapter/quality trimming",
            priority="high",
            reason="; ".join(trim_evidence),
            command_preview="fastp -i sample_R1.fastq.gz -o trimmed_R1.fastq.gz --thread 4",
            requires_execute=True,
            requires_external_tool="fastp or cutadapt",
        ))

    if review_evidence:
        decisions.append(Decision(
            action="review_samples",
            decision_type="review",
            reason="FastQC identified issues requiring manual or expert review before proceeding.",
            evidence=review_evidence,
            confidence="high",
            execute_allowed=False,
            safety_notes=["Do not remove samples or filter reads without biological justification."],
        ))
        recommended_actions.append(RecommendedAction(
            action="Review FastQC HTML reports for flagged samples",
            priority="medium",
            reason="; ".join(review_evidence),
            requires_execute=False,
        ))

    if gc_warnings:
        warnings.extend(gc_warnings)
        recommended_actions.append(RecommendedAction(
            action="Review GC content distribution with biological context",
            priority="medium",
            reason="GC content failures require biological interpretation, not automatic filtering.",
            requires_execute=False,
        ))

    if not trim_evidence and not review_evidence:
        decisions.append(Decision(
            action="no_trimming_required",
            decision_type="accept",
            reason="No FastQC modules triggered trimming or review rules.",
            evidence=["All evaluated samples: no adapter or quality failures"],
            confidence="high",
            execute_allowed=False,
            safety_notes=[
                "FastQC passing does not guarantee data quality for all downstream uses.",
                "Always validate final results in biological context.",
            ],
        ))
        recommended_actions.append(RecommendedAction(
            action="Proceed to alignment or downstream analysis",
            priority="low",
            reason="FastQC results do not indicate QC-based trimming is needed.",
            requires_execute=False,
        ))

    return _package(observations, decisions, recommended_actions, warnings)


def _module_status(modules: dict, name: str) -> str | None:
    return modules.get(name, {}).get("status")


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
