"""Helpers for linking observations to evidence sources."""
from __future__ import annotations

from typing import Any


def obs_to_evidence_source(obs: dict[str, Any]) -> str:
    source = obs.get("source", "unknown")
    category = obs.get("category", "unknown")
    return f"{source}:{category}"


def collect_obs_by_category(
    observations: list[dict[str, Any]],
    category: str,
    status: str | None = None,
) -> list[dict[str, Any]]:
    results = [o for o in observations if o.get("category") == category]
    if status is not None:
        results = [o for o in results if o.get("status") == status]
    return results


def obs_message_list(observations: list[dict[str, Any]]) -> list[str]:
    return [o.get("message", "") for o in observations if o.get("message")]
