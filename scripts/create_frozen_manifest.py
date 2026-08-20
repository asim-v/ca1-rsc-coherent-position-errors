#!/usr/bin/env python3
"""Create the write-once pre-confirmation provenance manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


FIXED_FILES = (
    "sub-M01_ses-20240318T100000_behavior+ecephys.nwb",
    "sub-M02_ses-20240318T100000_behavior+ecephys.nwb",
    "sub-M03_ses-20240624T100000_behavior+ecephys.nwb",
)

REPO_FILES = (
    "README.md",
    "requirements.txt",
    "docs/discovery_protocol.md",
    "docs/frozen_confirmation_spec.md",
    "docs/novelty_boundary.md",
    "scripts/coherent_error_core.py",
    "scripts/run_discovery.py",
    "scripts/run_confirmation.py",
    "scripts/calibrate_rotation_null.py",
    "scripts/create_frozen_manifest.py",
    "tests/test_core.py",
    "outputs/calibration/summary.json",
    "outputs/calibration/simulation_results.csv",
    "outputs/discovery/summary.json",
    "outputs/discovery/primary_record_scores.csv",
    "outputs/discovery/primary_mouse_scores.csv",
    "outputs/discovery/population_totals_record_scores.csv",
    "outputs/discovery/ca3_partial_record_scores.csv",
    "outputs/discovery_tracking_odd/summary.json",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record(path: Path, logical_path: str) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {"path": logical_path, "bytes": path.stat().st_size, "sha256": sha256(path)}


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, default=Path("outputs/provenance/frozen_manifest.json")
    )
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    output = (repo / args.output).resolve() if not args.output.is_absolute() else args.output.resolve()
    canonical = (repo / "outputs/provenance/frozen_manifest.json").resolve()
    if output != canonical:
        raise ValueError(f"Manifest path must be {canonical}")
    if output.exists():
        raise FileExistsError(output)
    if git(repo, "status", "--porcelain"):
        raise RuntimeError("Manifest creation requires a clean committed tree")
    calibration = json.loads((repo / "outputs/calibration/summary.json").read_text())
    if calibration.get("passed") is not True:
        raise RuntimeError("Rotation calibration did not pass")
    repo_records = [record(repo / relative, relative) for relative in REPO_FILES]
    raw_dir = args.raw_dir.resolve()
    raw_records = [record(raw_dir / name, name) for name in FIXED_FILES]
    payload = {
        "schema_version": "coherent-errors-confirmation-manifest-v1",
        "analysis_commit": git(repo, "rev-parse", "HEAD"),
        "python": sys.version,
        "confirmation_seed": 20260830,
        "null_draws": 9999,
        "fixed_inputs": raw_records,
        "repo_files": repo_records,
    }
    canonical_bytes = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["manifest_self_hash"] = hashlib.sha256(canonical_bytes).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "self_hash": payload["manifest_self_hash"]}, indent=2))


if __name__ == "__main__":
    main()

