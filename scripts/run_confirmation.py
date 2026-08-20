#!/usr/bin/env python3
"""One-shot held-out fixed-session confirmation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from coherent_error_core import ErrorRecord, hierarchy_weights
from run_discovery import process_session, sha256, summarize_variant


FIXED_FILES = (
    "sub-M01_ses-20240318T100000_behavior+ecephys.nwb",
    "sub-M02_ses-20240318T100000_behavior+ecephys.nwb",
    "sub-M03_ses-20240624T100000_behavior+ecephys.nwb",
)
SEEDS = {
    "primary": 20260830,
    "population_totals": 20260831,
    "ca3_partial": 20260832,
    "tracking_odd": 20260833,
    "literal_rsc": 20260834,
}
NULL_DRAWS = 9999


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def verify_manifest(repo: Path, raw_dir: Path, manifest_path: Path) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    claimed = manifest.pop("manifest_self_hash")
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
    observed = hashlib.sha256(canonical).hexdigest()
    manifest["manifest_self_hash"] = claimed
    if observed != claimed:
        raise RuntimeError("Manifest self-hash mismatch")
    if manifest["confirmation_seed"] != SEEDS["primary"] or manifest["null_draws"] != NULL_DRAWS:
        raise RuntimeError("Frozen seed/draw mismatch")
    for entry in manifest["repo_files"]:
        path = repo / entry["path"]
        if path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
            raise RuntimeError(f"Repository provenance mismatch: {path}")
    expected_names = set(FIXED_FILES)
    if {entry["path"] for entry in manifest["fixed_inputs"]} != expected_names:
        raise RuntimeError("Frozen input set mismatch")
    for entry in manifest["fixed_inputs"]:
        path = raw_dir / entry["path"]
        if path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
            raise RuntimeError(f"Raw input mismatch: {path}")
    subprocess.check_call(
        ["git", "merge-base", "--is-ancestor", manifest["analysis_commit"], "HEAD"], cwd=repo
    )
    calibration = json.loads((repo / "outputs/calibration/summary.json").read_text())
    if calibration != {
        "simulations_per_condition": 1000,
        "null_draws": 199,
        "seed": 20260826,
        "injected_correlation": 0.1,
        "null_rejection_rate": 0.05,
        "injected_power": 0.999,
        "passed": True,
    }:
        raise RuntimeError("Calibration summary does not match the frozen passing gate")
    return manifest


def subset_margin(
    records: list[ErrorRecord], observed: np.ndarray, null_bank: np.ndarray, keep: np.ndarray
) -> dict:
    chosen_records = [record for record, selected in zip(records, keep) if selected]
    weights = hierarchy_weights(chosen_records)
    selected_observed = observed[keep]
    selected_null = null_bank[keep]
    value = float(weights @ selected_observed)
    null = weights @ selected_null
    return {
        "observed_fisher_z": value,
        "null_mean_fisher_z": float(np.mean(null)),
        "excess_fisher_z": value - float(np.mean(null)),
    }


def influence_table(
    records: list[ErrorRecord], record_table: pd.DataFrame, record_nulls: np.ndarray
) -> pd.DataFrame:
    observed = record_table.observed_fisher_z.to_numpy(dtype=float)
    rows = []
    for direction in ("LtoR", "RtoL"):
        keep = np.asarray([record.direction == direction for record in records])
        rows.append({"kind": "direction", "omitted": "", "level": direction, **subset_margin(records, observed, record_nulls, keep)})
    for fold in (0, 1):
        keep = np.asarray([record.fold == fold for record in records])
        rows.append({"kind": "fold", "omitted": "", "level": str(fold), **subset_margin(records, observed, record_nulls, keep)})
    for subject in sorted({record.subject for record in records}):
        keep = np.asarray([record.subject != subject for record in records])
        rows.append({"kind": "leave_one_mouse_out", "omitted": subject, "level": "", **subset_margin(records, observed, record_nulls, keep)})
    return pd.DataFrame(rows)


def save_variant(
    output: Path,
    name: str,
    records: list[ErrorRecord],
    partial_ca3: bool,
    seed: int,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    summary, record_table, null, mouse_table, record_nulls = summarize_variant(
        records,
        NULL_DRAWS,
        seed,
        partial_ca3=partial_ca3,
        expected_records=24,
    )
    influence = influence_table(records, record_table, record_nulls)
    record_table.to_csv(output / f"{name}_record_scores.csv", index=False)
    mouse_table.to_csv(output / f"{name}_mouse_scores.csv", index=False)
    influence.to_csv(output / f"{name}_influence.csv", index=False)
    np.savez_compressed(
        output / f"{name}_nulls.npz", overall=null, record_level=record_nulls
    )
    return summary, mouse_table, influence, null, record_nulls


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    raw_dir = args.raw_dir.resolve()
    output = repo / "outputs/confirmation"
    sentinel = repo / "outputs/.confirmation_unblinded"
    manifest_path = repo / "outputs/provenance/frozen_manifest.json"
    if Path.cwd().resolve() != repo.resolve():
        raise RuntimeError("Confirmation must run from the repository root")
    if git(repo, "status", "--porcelain"):
        raise RuntimeError("Confirmation requires a clean committed tree")
    if output.exists() or sentinel.exists():
        raise FileExistsError("The one-shot confirmation has already been opened")
    manifest = verify_manifest(repo, raw_dir, manifest_path)
    missing = [name for name in FIXED_FILES if not (raw_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(missing)

    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel_payload = {
        "opened_utc": datetime.now(timezone.utc).isoformat(),
        "execution_commit": git(repo, "rev-parse", "HEAD"),
        "analysis_commit": manifest["analysis_commit"],
        "manifest_self_hash": manifest["manifest_self_hash"],
        "seed": SEEDS["primary"],
        "null_draws": NULL_DRAWS,
    }
    with sentinel.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(sentinel_payload, indent=2) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    output.mkdir(parents=True, exist_ok=False)

    even_records = {"primary": [], "population_totals": []}
    odd_records = {"primary": [], "population_totals": []}
    literal_records = {"primary": [], "population_totals": []}
    ledgers, equivalence = [], []
    literal_ledgers = []
    for name in FIXED_FILES:
        path = raw_dir / name
        records, rows, diagnostics = process_session(path, tracking_parity=0)
        for variant in even_records:
            even_records[variant].extend(records[variant])
        ledgers.extend(rows)
        equivalence.extend(diagnostics)

        records, _, _ = process_session(path, tracking_parity=1)
        for variant in odd_records:
            odd_records[variant].extend(records[variant])

        records, rows, _ = process_session(path, tracking_parity=0, literal_rsc=True)
        for variant in literal_records:
            literal_records[variant].extend(records[variant])
        literal_ledgers.extend(rows)

    eligibility = pd.DataFrame(ledgers)
    eligibility.to_csv(output / "eligibility.csv", index=False)
    pd.DataFrame(equivalence).to_csv(output / "continuous_vs_fixed_center.csv", index=False)
    pd.DataFrame(literal_ledgers).to_csv(output / "literal_rsc_eligibility.csv", index=False)
    if len(eligibility) != 6 or not eligibility.eligible.all():
        raise RuntimeError("All six held-out blocks must be eligible")
    max_equivalence_error = float(
        pd.DataFrame(equivalence).max_abs_residual_difference.max()
    )
    if max_equivalence_error > 1e-10:
        raise RuntimeError("Continuous/fixed-center residual equivalence failed")

    results = {}
    result_tables = {}
    specifications = (
        ("primary", even_records["primary"], False),
        ("population_totals", even_records["population_totals"], False),
        ("ca3_partial", even_records["primary"], True),
        ("tracking_odd", odd_records["primary"], False),
    )
    for name, records, partial_ca3 in specifications:
        summary, mice, influence, null, _ = save_variant(
            output, name, records, partial_ca3, SEEDS[name]
        )
        results[name] = summary
        result_tables[name] = (mice, influence, null)

    literal_evaluable = len(literal_records["primary"]) == 24
    if literal_evaluable:
        summary, mice, influence, null, _ = save_variant(
            output,
            "literal_rsc",
            literal_records["primary"],
            False,
            SEEDS["literal_rsc"],
        )
        results["literal_rsc"] = summary
        result_tables["literal_rsc"] = (mice, influence, null)
    else:
        results["literal_rsc"] = {
            "evaluable": False,
            "records": len(literal_records["primary"]),
        }

    primary_mice, primary_influence, _ = result_tables["primary"]
    direction_rows = primary_influence[primary_influence.kind == "direction"]
    fold_rows = primary_influence[primary_influence.kind == "fold"]
    lomo_rows = primary_influence[primary_influence.kind == "leave_one_mouse_out"]
    primary_pass = bool(
        results["primary"]["excess_fisher_z"] > 0
        and results["primary"]["p_value"] <= 0.05
        and np.all(primary_mice.excess_fisher_z > 0)
        and np.all(direction_rows.excess_fisher_z > 0)
        and np.all(fold_rows.excess_fisher_z > 0)
        and np.all(lomo_rows.excess_fisher_z > 0)
    )
    population_pass = bool(
        results["population_totals"]["excess_fisher_z"] > 0
        and results["population_totals"]["p_value"] <= 0.05
        and np.all(result_tables["population_totals"][0].excess_fisher_z > 0)
    )
    ca3_pass = bool(
        results["ca3_partial"]["excess_fisher_z"] > 0
        and results["ca3_partial"]["p_value"] <= 0.05
        and np.all(result_tables["ca3_partial"][0].excess_fisher_z > 0)
    )
    tracking_pass = bool(
        results["tracking_odd"]["excess_fisher_z"] > 0
        and np.sum(result_tables["tracking_odd"][0].excess_fisher_z > 0) >= 2
    )
    final = {
        "analysis": "one-shot held-out-session confirmation in three previously recorded mice",
        "execution": sentinel_payload,
        "results": results,
        "maximum_continuous_vs_fixed_residual_difference": max_equivalence_error,
        "primary_confirmed": primary_pass,
        "population_total_claim_available": bool(primary_pass and population_pass),
        "ca3_conditional_claim_available": bool(primary_pass and ca3_pass),
        "tracking_sensitivity_passed": tracking_pass,
        "narrow_candidate_contribution_supported": bool(
            primary_pass and population_pass and ca3_pass and tracking_pass
        ),
        "new_animal_replication": False,
    }
    (output / "summary.json").write_text(
        json.dumps(final, indent=2) + "\n", encoding="utf-8"
    )

    figure, axes = plt.subplots(1, 4, figsize=(15, 3.5), constrained_layout=True)
    for axis, name in zip(axes, ("primary", "population_totals", "ca3_partial", "tracking_odd")):
        null = result_tables[name][2]
        observed = results[name]["observed_fisher_z"]
        axis.hist(null, bins=45, color="#b8c7d9", edgecolor="none")
        axis.axvline(observed, color="#a43b3b", linewidth=2)
        axis.set(title=name.replace("_", " ").title(), xlabel="Hierarchical Fisher z")
    axes[0].set_ylabel("Rotation draws")
    figure.savefig(output / "confirmation_rotation_nulls.png", dpi=220)
    plt.close(figure)
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    main()

