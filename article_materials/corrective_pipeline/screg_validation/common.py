"""Small contracts, content identities and atomic, verifiable receipts."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any


class ContractError(ValueError):
    """An input or execution contract failed before a scientific claim was made."""


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def sha256(path: Path | str) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(2 ** 20), b""):
            h.update(block)
    return h.hexdigest()


def identities(values, label="identity") -> list[str]:
    out = []
    for item in values:
        value = str(item).strip() if item is not None else ""
        if not value or value.casefold() in {"na", "nan", "unknown", "null", "none", "n/a"}:
            raise ContractError(f"missing {label}")
        out.append(value)
    return out


def unique(values, label="identity") -> list[str]:
    out = identities(values, label)
    if len(out) != len(set(out)):
        raise ContractError(f"duplicate {label}")
    return out


def write_json(path: Path | str, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    fd, name = tempfile.mkstemp(prefix=".receipt-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        os.replace(name, path)
    finally:
        if Path(name).exists():
            Path(name).unlink()


def load_json(path: Path | str):
    return json.loads(Path(path).read_text())


def authorize(plan: dict, path: Path | None, action: str, item: str) -> dict:
    """Check the local run record before a download, inference or analysis step.

    The record is written by the person running the pipeline; no command in this package
    creates it. It must name the current plan and list the permitted action and item.
    """
    if path is None:
        raise ContractError("execution requires a run record (see INPUT_CONTRACTS.md)")
    auth = load_json(path)
    if (auth.get("approved") is not True or auth.get("approved_by") != "user"
            or not str(auth.get("instruction", "")).strip()
            or auth.get("plan_sha256") != digest(plan)):
        raise ContractError("authorization absent, unapproved or bound to a different plan")
    if item not in auth.get("actions", {}).get(action, []):
        raise ContractError(f"action not authorized: {action}/{item}")
    return auth


def confined_output(root: Path, relative: str) -> Path:
    raw = Path(relative)
    if raw.is_absolute() or ".." in raw.parts:
        raise ContractError("output must be a relative path within the new validation run")
    destination = (root / raw).resolve()
    if not destination.is_relative_to(root.resolve()):
        raise ContractError("output escapes run directory")
    return destination


def cache_receipt(path: Path, context: dict) -> dict | None:
    if not path.exists():
        return None
    record = load_json(path)
    if record.get("context_sha256") != digest(context):
        raise ContractError("cache context drift")
    for relative, expected in record.get("outputs", {}).items():
        output = confined_output(path.parent, relative)
        if not output.is_file() or sha256(output) != expected:
            raise ContractError("cache output missing or hash drift")
    if not record.get("outputs"):
        raise ContractError("empty cache receipt")
    return record
