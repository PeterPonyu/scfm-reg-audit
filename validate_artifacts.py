#!/usr/bin/env python3
"""Validate the public code-only capsule boundary."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FORBIDDEN_ROOTS = {"paper", "submission", "site", ".github", ".codex", ".omx", ".omc"}
PRIVATE_NEEDLES = ("/home/" + "zeyufu", "/Users/" + "zeyufu", "BEGIN OPENSSH PRIVATE KEY", "github_pat_", "ghp_")


class ValidationError(Exception):
    """The capsule violates its public release contract."""


def require(condition, message):
    """Raise a stable validation error instead of relying on ``assert``."""
    if not condition:
        raise ValidationError(message)


def bh(pvalues):
    """Benjamini-Hochberg step-up adjusted p-values in input order."""
    n = len(pvalues)
    order = sorted(range(n), key=pvalues.__getitem__)
    out = [0.0] * n
    running = 1.0
    for reverse_index in range(n - 1, -1, -1):
        index = order[reverse_index]
        running = min(running, pvalues[index] * n / (reverse_index + 1), 1.0)
        out[index] = running
    return out


def _walk_numbers(value):
    if isinstance(value, dict):
        for child in value.values():
            yield from _walk_numbers(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_numbers(child)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield value


def validate() -> None:
    present = {path.name for path in ROOT.iterdir()}
    require(not (present & FORBIDDEN_ROOTS), f"forbidden public roots: {sorted(present & FORBIDDEN_ROOTS)}")
    require((ROOT / "src").is_dir(), "src/ is missing")
    require((ROOT / "data/manifest/shared_genes.v2.json").is_file(), "frozen panel manifest is missing")

    for path in ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        if path.suffix == ".json":
            value = json.loads(path.read_text(), parse_constant=lambda token: (_ for _ in ()).throw(
                ValidationError(f"{path.relative_to(ROOT)} contains non-JSON constant {token}")))
            require(all(math.isfinite(number) for number in _walk_numbers(value)),
                    f"{path.relative_to(ROOT)} contains a non-finite number")
        if path.suffix in {".py", ".md", ".cff", ".example", ".json", ".yaml", ".yml"}:
            text = path.read_text(errors="replace")
            for needle in PRIVATE_NEEDLES:
                require(needle not in text, f"private marker {needle!r} in {path.relative_to(ROOT)}")


if __name__ == "__main__":
    validate()
    print("PASS code-capsule boundary")
