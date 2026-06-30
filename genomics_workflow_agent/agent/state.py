from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class Observation:
    source: str
    sample: str
    category: str
    status: str          # "pass" | "warn" | "fail" | "missing" | "error"
    severity: str        # "info" | "warning" | "critical"
    message: str
    evidence: dict[str, Any] = field(default_factory=dict)
    suggested_action: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "sample": self.sample,
            "category": self.category,
            "status": self.status,
            "severity": self.severity,
            "message": self.message,
            "evidence": self.evidence,
            "suggested_action": self.suggested_action,
        }


@dataclass
class Decision:
    action: str
    decision_type: str   # "trim" | "review" | "accept" | "flag" | "skip"
    reason: str
    evidence: list[str] = field(default_factory=list)
    confidence: str = "high"   # "high" | "medium" | "low"
    execute_allowed: bool = False
    executed: bool = False
    safety_notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "action": self.action,
            "decision_type": self.decision_type,
            "reason": self.reason,
            "evidence": self.evidence,
            "confidence": self.confidence,
            "execute_allowed": self.execute_allowed,
            "executed": self.executed,
            "safety_notes": self.safety_notes,
        }


@dataclass
class RecommendedAction:
    action: str
    priority: str        # "high" | "medium" | "low"
    reason: str
    command_preview: str = ""
    requires_execute: bool = False
    requires_external_tool: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "action": self.action,
            "priority": self.priority,
            "reason": self.reason,
            "command_preview": self.command_preview,
            "requires_execute": self.requires_execute,
            "requires_external_tool": self.requires_external_tool,
        }


@dataclass
class AgentState:
    input_path: str = ""
    workflow: str = ""
    input_summary: dict[str, Any] = field(default_factory=dict)
    planned_steps: list[dict] = field(default_factory=list)
    executed_steps: list[dict] = field(default_factory=list)
    observations: list[Observation] = field(default_factory=list)
    decisions: list[Decision] = field(default_factory=list)
    recommended_actions: list[RecommendedAction] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    provenance_paths: list[str] = field(default_factory=list)
    report_paths: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "input_path": self.input_path,
            "workflow": self.workflow,
            "input_summary": self.input_summary,
            "planned_steps": self.planned_steps,
            "executed_steps": self.executed_steps,
            "observations": [o.to_dict() for o in self.observations],
            "decisions": [d.to_dict() for d in self.decisions],
            "recommended_actions": [r.to_dict() for r in self.recommended_actions],
            "warnings": self.warnings,
            "limitations": self.limitations,
            "provenance_paths": self.provenance_paths,
            "report_paths": self.report_paths,
        }
