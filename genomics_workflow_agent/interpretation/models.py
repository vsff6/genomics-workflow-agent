"""Structured data models for biological interpretation scaffolds."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

INTERPRETATION_VERSION = "1.0"


@dataclass
class Finding:
    finding_id: str
    workflow: str
    sample: str
    observation: str
    evidence_source: str
    technical_explanations: list[str] = field(default_factory=list)
    plausible_biological_explanations: list[str] = field(default_factory=list)
    metadata_needed: list[str] = field(default_factory=list)
    recommended_validation: list[str] = field(default_factory=list)
    recommended_action: str = ""
    confidence: str = "medium"
    should_filter: bool = False
    should_preserve_until_review: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class Hypothesis:
    hypothesis_id: str
    statement: str
    supporting_observations: list[str] = field(default_factory=list)
    alternative_explanations: list[str] = field(default_factory=list)
    validation_steps: list[str] = field(default_factory=list)
    confidence: str = "low"
    interpretation_type: str = "ambiguous"  # "technical" | "biological" | "ambiguous"
    clinical_claim: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class InterpretationResult:
    interpretation_version: str
    workflow: str
    scope: str
    limitations: list[str] = field(default_factory=list)
    findings: list[Finding] = field(default_factory=list)
    hypotheses: list[Hypothesis] = field(default_factory=list)
    validation_recommendations: list[str] = field(default_factory=list)
    safety_flags: list[str] = field(default_factory=list)
    clinical_claims_allowed: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "interpretation_version": self.interpretation_version,
            "workflow": self.workflow,
            "scope": self.scope,
            "limitations": self.limitations,
            "findings": [f.to_dict() for f in self.findings],
            "hypotheses": [h.to_dict() for h in self.hypotheses],
            "validation_recommendations": self.validation_recommendations,
            "safety_flags": self.safety_flags,
            "clinical_claims_allowed": self.clinical_claims_allowed,
        }
