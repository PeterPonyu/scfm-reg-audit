"""Explicitly authorized, length-pinned data transfer.

Uses the Python standard library unless a local download gate module is registered as a resource.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import re
import shutil
import urllib.request
from urllib.parse import urlsplit

from .common import ContractError, authorize, confined_output, digest, load_json, sha256, write_json


class StandardGate:
    """Standard-library transport; the proxy comes from SC_FETCH_PROXY or the usual environment variables."""

    @staticmethod
    def _req(url, method="GET", headers=None):
        return urllib.request.Request(url, method=method, headers={"User-Agent": "screg-validation", **(headers or {})})

    @staticmethod
    def _urlopen(request, timeout, proxy):
        handlers = [urllib.request.ProxyHandler({"http": proxy, "https": proxy})] if proxy else []
        return urllib.request.build_opener(*handlers).open(request, timeout=timeout)


def load_fetch_gate(path):
    if path is None:
        return StandardGate
    path = Path(path)
    if not path.is_file():
        raise ContractError("registered download gate module missing")
    spec = importlib.util.spec_from_file_location("download_gate", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fetch_registered(plan, resource, destination_root, authorization, *, gate):
    auth = authorize(plan, authorization, "download", resource["id"])
    maximum = auth.get("download_limits", {}).get(resource["id"])
    if type(maximum) is not int or not 0 < maximum <= 5 * 1024 ** 3:
        raise ContractError("explicit per-file approved byte limit missing or exceeds 5 GiB")
    if maximum > resource["max_bytes"]:
        raise ContractError("approved download budget exceeds registered maximum")
    proxy = os.environ.get("SC_FETCH_PROXY", "")
    if proxy and urlsplit(proxy).scheme not in {"http", "https"}:
        raise ContractError("SC_FETCH_PROXY must be an http or https proxy URL")
    target = confined_output(Path(destination_root), resource["filename"])
    expected_hash = resource.get("sha256")
    if target.exists():
        sidecar = target.with_name(target.name + ".download.json")
        if expected_hash and sha256(target) == expected_hash:
            return {"status": "reused_verified", "sha256": expected_hash, "bytes": target.stat().st_size}
        if sidecar.exists():
            record = load_json(sidecar)
            if record.get("resource_sha256") == digest(resource) and record.get("sha256") == sha256(target):
                return {**record, "status": "reused_verified"}
        raise ContractError("existing destination unverified; refuse overwrite")
    # Strict HEAD; no range-GET fallback when Content-Length is unknown.
    with gate._urlopen(gate._req(resource["url"], "HEAD"), 20, proxy) as head:
        if int(head.status) != 200:
            raise ContractError("required HEAD did not return 200")
        length_text = head.headers.get("Content-Length", "")
        if not length_text.isdigit() or int(length_text) <= 0:
            raise ContractError("Content-Length missing; no body downloaded")
        length = int(length_text)
        if length > maximum:
            raise ContractError("HEAD exceeds explicitly approved bytes")
        if resource.get("expected_bytes") is not None and length != resource["expected_bytes"]:
            raise ContractError("HEAD length differs from frozen resource")
        head_url = head.geturl()
        etag = head.headers.get("ETag")
        modified = head.headers.get("Last-Modified")
    target.parent.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(target.parent).free < length + resource.get("disk_reserve_bytes", 1024 ** 3):
        raise ContractError("insufficient disk space before download")
    partial = target.with_name(target.name + ".part")
    partial_receipt = target.with_name(target.name + ".part.json")
    context = {"resource_sha256": digest(resource), "length": length, "url": head_url,
               "etag": etag, "last_modified": modified}
    start = partial.stat().st_size if partial.exists() else 0
    if start:
        if (not partial_receipt.exists() or load_json(partial_receipt) != context
                or not (etag or modified) or start > length):
            raise ContractError("partial download identity changed or cannot support safe resume")
    write_json(partial_receipt, context)
    if start < length:
        headers = {"Accept-Encoding": "identity"}
        if start:
            headers.update({"Range": f"bytes={start}-", "If-Range": etag or modified})
        with gate._urlopen(gate._req(resource["url"], "GET", headers), 60, proxy) as body:
            if body.geturl() != head_url:
                raise ContractError("GET/HEAD resource redirect changed")
            if start:
                match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", body.headers.get("Content-Range", ""))
                if (int(body.status) != 206 or match is None or tuple(map(int, match.groups())) != (start, length-1, length)):
                    raise ContractError("server failed exact resume range; partial preserved")
            elif int(body.status) != 200:
                raise ContractError("initial body did not return 200")
            if body.headers.get("Content-Encoding", "identity") != "identity":
                raise ContractError("unexpected transfer encoding prevents byte-pinned validation")
            with partial.open("ab" if start else "wb") as handle:
                total = start
                while True:
                    block = body.read(min(2 ** 20, length - total + 1))
                    if not block:
                        break
                    total += len(block)
                    if total > length:
                        raise ContractError("body exceeds pinned HEAD length")
                    handle.write(block)
        if total != length:
            raise ContractError("truncated body; verified partial preserved for resume")
    actual = sha256(partial)
    if expected_hash and actual != expected_hash:
        raise ContractError("download hash differs from frozen resource")
    partial.replace(target)
    record = {"status": "downloaded_verified", "resource_sha256": digest(resource),
              "sha256": actual, "bytes": length, "head": context,
              "authorization_sha256": digest(auth)}
    write_json(target.with_name(target.name + ".download.json"), record)
    return record
