#!/usr/bin/env python3
"""Independently recompute and optionally hash the saved confirmation outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


VARIANTS = ("primary", "population_totals", "ca3_partial", "tracking_odd")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def recompute_variant(output: Path, name: str) -> dict:
    records = pd.read_csv(output / f"{name}_record_scores.csv")
    mouse_table = pd.read_csv(output / f"{name}_mouse_scores.csv")
    arrays = np.load(output / f"{name}_nulls.npz")
    record_nulls = arrays["record_level"].astype(float)
    if record_nulls.shape != (24, 9999):
        raise RuntimeError(f"Unexpected null shape for {name}: {record_nulls.shape}")
    if len(records) != 24 or sorted(records.subject.unique()) != ["M01", "M02", "M03"]:
        raise RuntimeError(f"Unexpected record hierarchy for {name}")
    mouse_observed, mouse_nulls, mouse_rows = [], [], []
    for subject in ("M01", "M02", "M03"):
        index = np.flatnonzero(records.subject.to_numpy() == subject)
        observed = float(records.observed_fisher_z.to_numpy()[index].mean())
        null = record_nulls[index].mean(axis=0)
        mouse_observed.append(observed)
        mouse_nulls.append(null)
        mouse_rows.append(
            {
                "subject": subject,
                "observed_fisher_z": observed,
                "null_mean_fisher_z": float(null.mean()),
                "excess_fisher_z": observed - float(null.mean()),
            }
        )
    observed = float(np.mean(mouse_observed))
    null = np.mean(mouse_nulls, axis=0)
    p_value = float((1 + np.sum(null >= observed)) / 10000)
    saved_null = arrays["overall"].astype(float)
    if not np.allclose(saved_null, null, rtol=0, atol=1e-15):
        raise RuntimeError(f"Saved and recomputed overall null differ for {name}")
    expected_mice = pd.DataFrame(mouse_rows)
    if not np.allclose(
        mouse_table[["observed_fisher_z", "null_mean_fisher_z", "excess_fisher_z"]],
        expected_mice[["observed_fisher_z", "null_mean_fisher_z", "excess_fisher_z"]],
        rtol=1e-12,
        atol=1e-12,
    ):
        raise RuntimeError(f"Mouse table mismatch for {name}")
    return {
        "observed_fisher_z": observed,
        "null_mean_fisher_z": float(null.mean()),
        "excess_fisher_z": observed - float(null.mean()),
        "p_value": p_value,
        "positive_mouse_margins": int(np.sum(expected_mice.excess_fisher_z > 0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    output = repo / "outputs/confirmation"
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    recomputed = {name: recompute_variant(output, name) for name in VARIANTS}
    for name, values in recomputed.items():
        saved = summary["results"][name]
        for key, value in values.items():
            if isinstance(value, float):
                if not np.isclose(value, saved[key], rtol=1e-12, atol=1e-12):
                    raise RuntimeError(f"Summary mismatch: {name}.{key}")
            elif value != saved[key]:
                raise RuntimeError(f"Summary mismatch: {name}.{key}")
    if summary["primary_confirmed"] is not True:
        raise RuntimeError("Saved primary decision is not confirmed")
    report = {
        "verified": True,
        "recomputed": recomputed,
        "manifest_self_hash": summary["execution"]["manifest_self_hash"],
        "execution_commit": summary["execution"]["execution_commit"],
    }
    if args.write_manifest:
        manifest_path = output / "result_manifest.json"
        if manifest_path.exists():
            raise FileExistsError(manifest_path)
        files = []
        for path in sorted(output.iterdir()):
            if path.is_file() and path.name != manifest_path.name:
                files.append(
                    {
                        "path": path.name,
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )
        sentinel = repo / "outputs/.confirmation_unblinded"
        files.append(
            {
                "path": "../.confirmation_unblinded",
                "bytes": sentinel.stat().st_size,
                "sha256": sha256(sentinel),
            }
        )
        payload = {**report, "files": files}
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        payload["result_manifest_self_hash"] = hashlib.sha256(canonical).hexdigest()
        manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        report["result_manifest_self_hash"] = payload["result_manifest_self_hash"]
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
