"""Query-only GPU readiness. Idle watts never imply insufficient readiness."""
from __future__ import annotations

import subprocess
import os
import xml.etree.ElementTree as ET

from .common import ContractError


def number(value):
    try:
        return float(str(value).split()[0])
    except (ValueError, TypeError, IndexError):
        return None


def parse_gpu(xml, policy, *, own_pid=None):
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as exc:
        raise ContractError("unreadable NVIDIA query") from exc
    devices = root.findall("gpu")
    index = int(policy.get("index", 0))
    if not 0 <= index < len(devices):
        raise ContractError("configured GPU absent")
    gpu = devices[index]
    temperature = number(gpu.findtext("temperature/gpu_temp"))
    slowdown = number(gpu.findtext("temperature/gpu_temp_slow_threshold"))
    limit_margin = number(gpu.findtext("temperature/gpu_temp_tlimit"))
    slow_margin = number(gpu.findtext("temperature/gpu_temp_slow_tlimit_threshold"))
    temperature_source = "absolute_driver_slowdown_threshold"
    if slowdown is None and None not in (temperature, limit_margin, slow_margin):
        # New NVIDIA XML expresses thresholds as offsets from the target limit.
        slowdown = temperature + limit_margin - slow_margin
        temperature_source = "driver_temperature_plus_TLimit_minus_slowdown_offset"
    explicit = policy.get("max_temperature_c")
    limits = [v for v in (slowdown - float(policy.get("slowdown_margin_c", 5)) if slowdown is not None else None,
                         float(explicit) if explicit is not None else None) if v is not None]
    blockers = []
    if temperature is None or not limits:
        blockers.append("temperature_or_driver_slowdown_limit_unavailable")
    elif temperature >= min(limits):
        blockers.append("temperature_not_recovered")
    for family in ("clocks_event_reasons", "clocks_throttle_reasons"):
        reasons = gpu.find(family)
        if reasons is not None:
            for reason in reasons:
                if "thermal_slowdown" in reason.tag and reason.text == "Active":
                    blockers.append(reason.tag)
    processes = gpu.findall("processes/process_info")
    if any(p.findtext("type") in {"C", "C+G"} and p.findtext("pid") != str(own_pid) for p in processes):
        blockers.append("another_compute_process_occupies_gpu")
    recovery = gpu.findtext("gpu_recovery_action")
    if recovery and recovery.strip().casefold() not in {"none", "n/a"}:
        blockers.append("driver_requests_recovery_action")
    return {"ready": not blockers, "blockers": blockers, "uuid": gpu.findtext("uuid"),
            "temperature_c": temperature, "driver_slowdown_c": slowdown,
            "threshold_source": temperature_source,
            "effective_admission_temperature_c": min(limits) if limits else None,
            "power_draw_w": number(gpu.findtext("gpu_power_readings/instant_power_draw", gpu.findtext("power_readings/power_draw"))),
            "power_limit_w": number(gpu.findtext("gpu_power_readings/current_power_limit", gpu.findtext("power_readings/power_limit"))),
            "recovery_action": recovery,
            "compute_process_count": sum(p.findtext("type") in {"C", "C+G"} for p in processes),
            "scope": "read-only admission; no power, temperature, clock, reset or reboot changes"}


def gpu_readiness(policy, *, command=subprocess.run):
    try:
        result = command(["nvidia-smi", "-q", "-x"], capture_output=True, text=True, timeout=10)
        if result.returncode:
            return {"ready": False, "blockers": ["driver_query_failed"],
                    "diagnostic": result.stderr.strip() or result.stdout.strip()}
        return parse_gpu(result.stdout, policy, own_pid=os.getpid())
    except (OSError, subprocess.TimeoutExpired, ContractError) as exc:
        return {"ready": False, "blockers": ["gpu_readiness_unknown"], "diagnostic": str(exc)}


def require_gpu(policy):
    state = gpu_readiness(policy)
    if not state["ready"]:
        raise ContractError("GPU work refused: " + ", ".join(state["blockers"]))
    return state
