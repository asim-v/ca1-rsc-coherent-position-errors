#!/usr/bin/env python3
"""Calibrate the complete-product cyclic-rotation test on synthetic records."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from coherent_error_core import ErrorRecord, fisher_z, hierarchy_weights


def synthetic_records(rng: np.random.Generator, correlation: float) -> list[ErrorRecord]:
    records = []
    for subject in ("M01", "M02", "M03"):
        for block in (1, 2):
            for direction in ("LtoR", "RtoL"):
                for fold in (0, 1):
                    left = rng.normal(size=(10, 16))
                    noise = rng.normal(size=(10, 16))
                    right = correlation * left + np.sqrt(1 - correlation**2) * noise
                    observed_r = float(np.corrcoef(left.ravel(), right.ravel())[0, 1])
                    records.append(
                        ErrorRecord(
                            subject=subject,
                            session=f"{subject}-heldout",
                            block=block,
                            direction=direction,
                            fold=fold,
                            trials=np.arange(10),
                            ca1=left,
                            rsc=right,
                            ca3=None,
                            observed_r=observed_r,
                            observed_z=fisher_z(observed_r),
                            row_count=left.size,
                        )
                    )
    return records


def one_test(
    rng: np.random.Generator, correlation: float, null_draws: int
) -> tuple[float, float]:
    records = synthetic_records(rng, correlation)
    weights = hierarchy_weights(records)
    observed = float(weights @ np.asarray([record.observed_z for record in records]))
    null_bank = np.empty((len(records), null_draws), dtype=float)
    for index, record in enumerate(records):
        shifts = rng.integers(0, len(record.trials), size=null_draws)
        left = record.ca1.ravel().astype(float)
        left -= np.mean(left)
        shifted = np.stack(
            [np.roll(record.rsc, int(shift), axis=0).ravel() for shift in shifts]
        ).astype(float)
        shifted -= np.mean(shifted, axis=1, keepdims=True)
        correlations = (shifted @ left) / (
            np.linalg.norm(shifted, axis=1) * np.linalg.norm(left)
        )
        null_bank[index] = np.arctanh(np.clip(correlations, -1 + 1e-12, 1 - 1e-12))
    null = weights @ null_bank
    p_value = float((1 + np.sum(null >= observed)) / (null_draws + 1))
    return p_value, observed - float(np.mean(null))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/calibration"))
    parser.add_argument("--simulations", type=int, default=1000)
    parser.add_argument("--null-draws", type=int, default=199)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--injected-correlation", type=float, default=0.10)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {args.output_dir}")
    if args.simulations < 1000 or args.null_draws < 199:
        raise ValueError("Calibration requires >=1000 simulations and >=199 null draws")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    rows = []
    for condition, correlation in (("null", 0.0), ("injected", args.injected_correlation)):
        for simulation in range(args.simulations):
            rng = np.random.default_rng(
                np.random.SeedSequence([args.seed, 0 if condition == "null" else 1, simulation])
            )
            p_value, margin = one_test(rng, correlation, args.null_draws)
            rows.append(
                {
                    "condition": condition,
                    "simulation": simulation,
                    "p_value": p_value,
                    "margin": margin,
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(args.output_dir / "simulation_results.csv", index=False)
    null_rate = float(np.mean(table.loc[table.condition == "null", "p_value"] <= 0.05))
    power = float(np.mean(table.loc[table.condition == "injected", "p_value"] <= 0.05))
    summary = {
        "simulations_per_condition": args.simulations,
        "null_draws": args.null_draws,
        "seed": args.seed,
        "injected_correlation": args.injected_correlation,
        "null_rejection_rate": null_rate,
        "injected_power": power,
        "passed": bool(0.025 <= null_rate <= 0.075 and power >= 0.80),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
