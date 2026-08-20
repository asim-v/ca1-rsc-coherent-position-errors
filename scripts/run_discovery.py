#!/usr/bin/env python3
"""Run the exploratory 12-session CA1--RSC coherent-deviation analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO

from coherent_error_core import (
    DIRECTIONS,
    ErrorRecord,
    active_unit_masks,
    behavior_blocks,
    decoded_bin_table,
    hierarchy_weights,
    linearize_position,
    mouse_weights,
    fisher_z,
    record_correlation,
    residual_design,
    residualize,
    rotation_draws,
    table_to_record,
    traversal_observation,
    traversals,
)


DISCOVERY_SESSIONS = {
    "M01": ("20240312", "20240313", "20240314"),
    "M02": ("20240312", "20240313", "20240314"),
    "M03": ("20240621", "20240622", "20240623"),
    "M05": ("20240729", "20240730", "20240731"),
}


def file_identity(path: Path) -> tuple[str, str]:
    match = re.search(r"sub-(M\d+)_ses-(\d{8})", path.name)
    if match is None:
        raise ValueError(f"Cannot parse session identity: {path.name}")
    return match.group(1), match.group(2)


def expected_paths(raw_dir: Path) -> list[Path]:
    paths = sorted(raw_dir.glob("*.nwb"))
    selected = []
    for path in paths:
        subject, session = file_identity(path)
        if session in DISCOVERY_SESSIONS.get(subject, ()):
            selected.append(path)
    expected = {(mouse, session) for mouse, sessions in DISCOVERY_SESSIONS.items() for session in sessions}
    observed = {file_identity(path) for path in selected}
    if observed != expected or len(selected) != 12:
        raise RuntimeError(f"Discovery set mismatch: missing={expected-observed}, extra={observed-expected}")
    return selected


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def process_session(
    path: Path, tracking_parity: int = 0, literal_rsc: bool = False
) -> tuple[dict[str, list[ErrorRecord]], list[dict], list[dict]]:
    subject, session = file_identity(path)
    with NWBHDF5IO(str(path), "r", load_namespaces=True) as io:
        nwb = io.read()
        behavior = nwb.processing["behavior"]
        series = behavior["AnimalPosition"]["Position"]
        timestamps = np.asarray(series.timestamps[:]).squeeze().astype(float)
        position = linearize_position(np.asarray(series.data[:]).squeeze())
        speed = np.asarray(behavior["Speed"].data[:]).squeeze().astype(float)
        units = nwb.units.to_dataframe()
        populations = {
            "CA1": units[(units.cell_area == "CA1") & (units.cell_type == "Pyramidal Cell")],
            "RSC": units[
                (units.cell_area == "RSC")
                & (
                    (units.cell_type == "Pyramidal Cell")
                    if literal_rsc
                    else (units.cell_type != "Narrow Interneuron")
                )
            ],
            "CA3": units[(units.cell_area == "CA3") & (units.cell_type == "Pyramidal Cell")],
        }
        traversal_rows = traversals(position, timestamps)
        observations = []
        for row in traversal_rows:
            observation = traversal_observation(
                row, position, timestamps, speed, populations, frame_parity=tracking_parity
            )
            if observation is not None:
                observations.append(observation)

    masks = active_unit_masks(observations, populations, minimum_spikes=20)
    records = {"primary": [], "population_totals": []}
    ledgers: list[dict] = []
    equivalence: list[dict] = []
    for block in sorted({observation.block for observation in observations}):
        block_masks = {area: masks[(block, area)] for area in populations}
        counts = {area: int(np.sum(mask)) for area, mask in block_masks.items()}
        by_direction = {
            direction: sorted(
                [
                    observation
                    for observation in observations
                    if observation.block == block and observation.direction == direction
                ],
                key=lambda observation: observation.direction_trial,
            )
            for direction in DIRECTIONS
        }
        eligible = (
            counts["CA1"] >= 10
            and counts["RSC"] >= 5
            and all(len(by_direction[direction]) >= 16 for direction in DIRECTIONS)
        )
        ledgers.append(
            {
                "subject": subject,
                "session": session,
                "block": block,
                "behavior_blocks": len(behavior_blocks(timestamps)),
                "ca1_units": counts["CA1"],
                "rsc_units": counts["RSC"],
                "ca3_units": counts["CA3"],
                "ltr_traversals": len(by_direction["LtoR"]),
                "rtl_traversals": len(by_direction["RtoL"]),
                "eligible": eligible,
            }
        )
        if not eligible:
            continue
        decoder_masks = {area: mask for area, mask in block_masks.items() if np.sum(mask) > 0}
        for direction in DIRECTIONS:
            direction_observations = by_direction[direction]
            for fold in (0, 1):
                train = [
                    observation
                    for observation in direction_observations
                    if observation.direction_trial % 2 == fold
                ]
                test = [
                    observation
                    for observation in direction_observations
                    if observation.direction_trial % 2 != fold
                ]
                table = decoded_bin_table(train, test, decoder_masks)
                if len(table) < 50 or table.trial.nunique() < 7:
                    raise RuntimeError(
                        f"Sparse held-out table in {subject} {session} block {block} {direction} fold {fold}"
                    )
                records["primary"].append(
                    table_to_record(table, subject, session, block, direction, fold, False)
                )
                records["population_totals"].append(
                    table_to_record(table, subject, session, block, direction, fold, True)
                )
                design = residual_design(table, include_population_totals=False)
                bin_center = (table.position_bin.to_numpy(dtype=float) + 0.5) / 24
                for area in ("CA1", "RSC"):
                    decoded = table[f"decoded_{area}"].to_numpy(dtype=float)
                    continuous = residualize(decoded - table.coordinate.to_numpy(dtype=float), design)
                    fixed_center = residualize(decoded - bin_center, design)
                    equivalence.append(
                        {
                            "subject": subject,
                            "session": session,
                            "block": block,
                            "direction": direction,
                            "fold": fold,
                            "area": area,
                            "max_abs_residual_difference": float(
                                np.max(np.abs(continuous - fixed_center))
                            ),
                        }
                    )
    return records, ledgers, equivalence


def summarize_variant(
    records: list[ErrorRecord],
    null_draws: int,
    seed: int,
    partial_ca3: bool = False,
    expected_records: int = 48,
) -> tuple[dict, pd.DataFrame, np.ndarray, pd.DataFrame, np.ndarray]:
    if len(records) != expected_records:
        raise RuntimeError(f"Expected {expected_records} records, found {len(records)}")
    rng = np.random.default_rng(seed)
    record_nulls = np.empty((len(records), null_draws), dtype=float)
    observed_record_values = np.asarray(
        [fisher_z(record_correlation(record, partial_ca3=partial_ca3)) for record in records]
    )
    rows = []
    for index, record in enumerate(records):
        shifts = rng.integers(0, len(record.trials), size=null_draws)
        record_nulls[index] = rotation_draws(record, shifts, partial_ca3=partial_ca3)
        if not np.all(np.isfinite(record_nulls[index])):
            raise RuntimeError("Nonfinite rotation score")
        rows.append(
            {
                "subject": record.subject,
                "session": record.session,
                "block": record.block,
                "direction": record.direction,
                "fold": record.fold,
                "trials": len(record.trials),
                "rows": record.row_count,
                "observed_r": float(record_correlation(record, partial_ca3=partial_ca3)),
                "observed_fisher_z": observed_record_values[index],
                "null_mean_fisher_z": float(np.mean(record_nulls[index])),
            }
        )
    weights = hierarchy_weights(records)
    observed = float(weights @ observed_record_values)
    null = weights @ record_nulls
    p_value = float((1 + np.sum(null >= observed)) / (null_draws + 1))
    mouse_rows = []
    for subject in sorted({record.subject for record in records}):
        weights_mouse = mouse_weights(records, subject)
        mouse_observed = float(weights_mouse @ observed_record_values)
        mouse_null = weights_mouse @ record_nulls
        mouse_rows.append(
            {
                "subject": subject,
                "observed_fisher_z": mouse_observed,
                "null_mean_fisher_z": float(np.mean(mouse_null)),
                "excess_fisher_z": mouse_observed - float(np.mean(mouse_null)),
            }
        )
    summary = {
        "records": len(records),
        "subjects": len({record.subject for record in records}),
        "sessions": len({(record.subject, record.session) for record in records}),
        "null_draws": null_draws,
        "seed": seed,
        "observed_fisher_z": observed,
        "null_mean_fisher_z": float(np.mean(null)),
        "excess_fisher_z": observed - float(np.mean(null)),
        "p_value": p_value,
        "positive_mouse_margins": int(
            np.sum(pd.DataFrame(mouse_rows).excess_fisher_z.to_numpy() > 0)
        ),
    }
    return summary, pd.DataFrame(rows), null, pd.DataFrame(mouse_rows), record_nulls


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/discovery"))
    parser.add_argument("--null-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260820)
    parser.add_argument("--tracking-parity", type=int, choices=(0, 1), default=0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.null_draws < 999:
        raise ValueError("Discovery requires at least 999 rotation draws")
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.force:
        raise FileExistsError(f"Refusing to overwrite {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_records = {"primary": [], "population_totals": []}
    ledger_rows: list[dict] = []
    equivalence_rows: list[dict] = []
    inputs = []
    for path in expected_paths(args.raw_dir):
        print(f"Processing {path.name}", flush=True)
        records, ledgers, equivalence = process_session(path, tracking_parity=args.tracking_parity)
        for variant in all_records:
            all_records[variant].extend(records[variant])
        ledger_rows.extend(ledgers)
        equivalence_rows.extend(equivalence)
        inputs.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})

    summaries = {}
    null_arrays = {}
    variant_specs = (
        ("primary", "primary", False),
        ("population_totals", "population_totals", False),
        ("ca3_partial", "primary", True),
    )
    for offset, (variant, record_source, partial_ca3) in enumerate(variant_specs):
        summary, record_table, null, mouse_table, _ = summarize_variant(
            all_records[record_source], args.null_draws, args.seed + offset,
            partial_ca3=partial_ca3,
        )
        summaries[variant] = summary
        null_arrays[variant] = null
        record_table.to_csv(args.output_dir / f"{variant}_record_scores.csv", index=False)
        mouse_table.to_csv(args.output_dir / f"{variant}_mouse_scores.csv", index=False)

    pd.DataFrame(ledger_rows).to_csv(args.output_dir / "eligibility.csv", index=False)
    pd.DataFrame(equivalence_rows).to_csv(
        args.output_dir / "continuous_vs_fixed_center.csv", index=False
    )
    np.savez_compressed(args.output_dir / "rotation_nulls.npz", **null_arrays)
    result = {
        "analysis": "exploratory discovery; not the sealed confirmation",
        "estimator_version": "discovery-v0.2-temporal-residuals-complete-product-group",
        "tracking_frame_parity": args.tracking_parity,
        "inputs": inputs,
        "variants": summaries,
        "maximum_continuous_vs_fixed_residual_difference": float(
            pd.DataFrame(equivalence_rows).max_abs_residual_difference.max()
        ),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )

    figure, axes = plt.subplots(1, 3, figsize=(12.5, 3.6), constrained_layout=True)
    for axis, variant in zip(axes, ("primary", "population_totals", "ca3_partial")):
        null = null_arrays[variant]
        observed = summaries[variant]["observed_fisher_z"]
        axis.hist(null, bins=45, color="#b8c7d9", edgecolor="none")
        axis.axvline(observed, color="#a43b3b", linewidth=2)
        axis.set(
            title=variant.replace("_", " ").title(),
            xlabel="Hierarchical Fisher z",
            ylabel="Rotation draws",
        )
    figure.savefig(args.output_dir / "discovery_rotation_nulls.png", dpi=220)
    plt.close(figure)
    print(json.dumps(result["variants"], indent=2))


if __name__ == "__main__":
    main()
