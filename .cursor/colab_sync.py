#!/usr/bin/env python3
"""将 Colab notebook 输出的同步标记写回本地 OmniScreen-AI 项目目录。"""
from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

SYNC_START = "__OMNISCREEN_SYNC_START__"
SYNC_END = "__OMNISCREEN_SYNC_END__"

DEFAULT_LOCAL_ROOT = Path(__file__).resolve().parents[1]


def extract_sync_payload(text: str) -> dict | None:
    payloads = extract_all_sync_payloads(text)
    if not payloads:
        return None
    merged: dict = {"files": {}}
    for p in payloads:
        merged["files"].update(p.get("files", {}))
    return merged


def extract_all_sync_payloads(text: str) -> list[dict]:
    payloads: list[dict] = []
    pos = 0
    while True:
        start = text.find(SYNC_START, pos)
        if start < 0:
            break
        end = text.find(SYNC_END, start)
        if end < 0:
            break
        try:
            payloads.append(json.loads(text[start + len(SYNC_START) : end].strip()))
        except json.JSONDecodeError:
            pass
        pos = end + len(SYNC_END)
    return payloads


def write_files_local(payload: dict, local_root: Path) -> list[str]:
    written: list[str] = []
    for rel, b64 in payload.get("files", {}).items():
        out = local_root / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(base64.b64decode(b64))
        written.append(str(out))
    return written


def sync_from_output(text: str, local_root: str | Path = DEFAULT_LOCAL_ROOT) -> list[str]:
    data = extract_sync_payload(text)
    if not data:
        return []
    return write_files_local(data, Path(local_root))


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] != "-":
        text = Path(sys.argv[1]).read_text(encoding="utf-8")
    else:
        text = sys.stdin.read()
    root = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_LOCAL_ROOT
    written: list[str] = []
    for payload in extract_all_sync_payloads(text):
        written.extend(write_files_local(payload, root))
    if not written:
        print("No sync payload found.", file=sys.stderr)
        sys.exit(1)
    for p in written:
        print(f"Wrote: {p}")
