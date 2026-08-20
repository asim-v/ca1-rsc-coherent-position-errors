#!/usr/bin/env python3
"""Core utilities for cross-fitted CA1--RSC position-deviation analyses."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.special import logsumexp


DIRECTIONS = ("LtoR", "RtoL")
N_POSITION_BINS = 24
CENTRAL_BINS = tuple(range(4, 20))


def linearize_position(position: np.ndarray) -> np.ndarray:
    position = np.asarray(position).squeeze()
    if position.ndim == 1:
        return position.astype(float)
    valid = np.all(np.isfinite(position), axis=1)
    center = np.nanmean(position[valid], axis=0)
    _, _, vh = np.linalg.svd(position[valid] - center, full_matrices=False)
    return ((position - center) @ vh[0]).astype(float)


def behavior_blocks(timestamps: np.ndarray) -> list[tuple[int, int]]:
    delta = np.diff(timestamps)
    threshold = max(1.0, float(np.median(delta) * 10))
    splits = np.flatnonzero(delta > threshold)
    return [
        (int(a), int(b))
        for a, b in zip(np.r_[0, splits + 1], np.r_[splits, len(timestamps) - 1])
    ]


def endpoint_events(normalized: np.ndarray, low: float = 0.15, high: float = 0.85):
    events: list[tuple[int, str]] = []
    state = None
    for index, value in enumerate(normalized):
        new = "L" if value <= low else "R" if value >= high else None
        if new is not None and new != state:
            events.append((index, new))
            state = new
    return events


def traversals(position: np.ndarray, timestamps: np.ndarray) -> list[dict]:
    rows: list[dict] = []
    for block, (block_start, block_stop) in enumerate(behavior_blocks(timestamps), start=1):
        local = position[block_start : block_stop + 1]
        lo, hi = np.nanpercentile(local, [1, 99])
        if not np.isfinite(lo + hi) or hi <= lo:
            continue
        normalized = np.clip((local - lo) / (hi - lo), 0, 1)
        events = endpoint_events(normalized)
        global_trial = 0
        direction_trials = {direction: 0 for direction in DIRECTIONS}
        for (a, side_a), (b, side_b) in zip(events[:-1], events[1:]):
            if side_a == side_b:
                continue
            start, stop = block_start + a, block_start + b
            duration = float(timestamps[stop] - timestamps[start])
            coverage = float(np.nanmax(normalized[a : b + 1]) - np.nanmin(normalized[a : b + 1]))
            if not (0 < duration <= 120 and coverage >= 0.70):
                continue
            direction = side_a + "to" + side_b
            global_trial += 1
            direction_trials[direction] += 1
            rows.append(
                {
                    "block": block,
                    "trial": global_trial,
                    "direction_trial": direction_trials[direction],
                    "direction": direction,
                    "start": start,
                    "stop": stop,
                    "lo": float(lo),
                    "hi": float(hi),
                }
            )
    return rows


def pava(values: np.ndarray) -> np.ndarray:
    """Unweighted increasing pool-adjacent-violators fit."""
    y = np.asarray(values, dtype=float)
    if y.ndim != 1 or not len(y) or not np.all(np.isfinite(y)):
        raise ValueError("PAVA requires a finite one-dimensional vector")
    levels: list[float] = []
    weights: list[int] = []
    for value in y:
        levels.append(float(value))
        weights.append(1)
        while len(levels) >= 2 and levels[-2] > levels[-1]:
            weight = weights[-2] + weights[-1]
            level = (levels[-2] * weights[-2] + levels[-1] * weights[-1]) / weight
            levels[-2:] = [level]
            weights[-2:] = [weight]
    return np.concatenate([np.full(weight, level) for level, weight in zip(levels, weights)])


def monotone_tracking_coordinate(
    row: dict,
    position: np.ndarray,
    timestamps: np.ndarray,
    frame_parity: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit monotone position using one parity of tracking frames only."""
    indices = np.arange(row["start"], row["stop"] + 1)
    chosen = indices[frame_parity::2]
    if len(chosen) < 2:
        raise ValueError("Traversal has fewer than two selected tracking frames")
    raw = np.clip((position[chosen] - row["lo"]) / (row["hi"] - row["lo"]), 0, 1)
    sign = 1.0 if row["direction"] == "LtoR" else -1.0
    fitted = sign * pava(sign * raw)
    return timestamps[chosen], np.clip(fitted, 0, 1)


@dataclass
class TraversalObservation:
    block: int
    direction: str
    direction_trial: int
    centers: np.ndarray
    coordinate: np.ndarray
    speed: np.ndarray
    position_bin: np.ndarray
    counts: dict[str, np.ndarray]


def count_spikes(spike_times: Iterable[np.ndarray], starts: np.ndarray, stops: np.ndarray) -> np.ndarray:
    spike_times = list(spike_times)
    counts = np.zeros((len(starts), len(spike_times)), dtype=np.int16)
    for unit, values in enumerate(spike_times):
        values = np.asarray(values, dtype=float)
        left = np.searchsorted(values, starts, side="left")
        right = np.searchsorted(values, stops, side="left")
        counts[:, unit] = right - left
    return counts


def traversal_observation(
    row: dict,
    position: np.ndarray,
    timestamps: np.ndarray,
    speed: np.ndarray,
    populations: dict[str, pd.DataFrame],
    window_seconds: float = 0.250,
    speed_threshold: float = 2.5,
    frame_parity: int = 0,
) -> TraversalObservation | None:
    start_time = float(timestamps[row["start"]])
    stop_time = float(timestamps[row["stop"]])
    starts = np.arange(start_time, stop_time - window_seconds + 1e-12, window_seconds)
    if not len(starts):
        return None
    stops = starts + window_seconds
    centers = starts + window_seconds / 2
    track_t, track_x = monotone_tracking_coordinate(row, position, timestamps, frame_parity)
    coordinate = np.interp(centers, track_t, track_x)
    local_speed = np.interp(centers, timestamps, speed)
    keep = np.isfinite(coordinate) & np.isfinite(local_speed) & (local_speed > speed_threshold)
    if not np.any(keep):
        return None
    starts, stops, centers = starts[keep], stops[keep], centers[keep]
    coordinate, local_speed = coordinate[keep], local_speed[keep]
    position_bin = np.clip((coordinate * N_POSITION_BINS).astype(int), 0, N_POSITION_BINS - 1)
    counts = {
        area: count_spikes(frame.spike_times, starts, stops)
        for area, frame in populations.items()
    }
    return TraversalObservation(
        block=int(row["block"]),
        direction=str(row["direction"]),
        direction_trial=int(row["direction_trial"]),
        centers=centers,
        coordinate=coordinate,
        speed=local_speed,
        position_bin=position_bin,
        counts=counts,
    )


@dataclass
class ConditionalMultinomialDecoder:
    log_cell_probability: np.ndarray
    position_centers: np.ndarray

    def predict(self, counts: np.ndarray) -> np.ndarray:
        counts = np.asarray(counts, dtype=float)
        log_likelihood = counts @ self.log_cell_probability.T
        log_posterior = log_likelihood - logsumexp(log_likelihood, axis=1, keepdims=True)
        posterior = np.exp(log_posterior)
        return posterior @ self.position_centers


def fit_decoder(
    counts: np.ndarray,
    labels: np.ndarray,
    pseudocount: float = 0.5,
    smoothing_sigma: float = 1.0,
) -> ConditionalMultinomialDecoder:
    counts = np.asarray(counts, dtype=float)
    labels = np.asarray(labels, dtype=int)
    if counts.ndim != 2 or len(counts) != len(labels) or counts.shape[1] == 0:
        raise ValueError("Invalid decoder training arrays")
    spatial_counts = np.zeros((N_POSITION_BINS, counts.shape[1]), dtype=float)
    for position_bin in range(N_POSITION_BINS):
        spatial_counts[position_bin] = counts[labels == position_bin].sum(axis=0)
    spatial_counts = gaussian_filter1d(
        spatial_counts, smoothing_sigma, axis=0, mode="nearest"
    )
    spatial_counts += pseudocount
    probability = spatial_counts / spatial_counts.sum(axis=1, keepdims=True)
    return ConditionalMultinomialDecoder(
        log_cell_probability=np.log(probability),
        position_centers=(np.arange(N_POSITION_BINS) + 0.5) / N_POSITION_BINS,
    )


def active_unit_masks(
    observations: list[TraversalObservation],
    populations: dict[str, pd.DataFrame],
    minimum_spikes: int = 20,
) -> dict[tuple[int, str], np.ndarray]:
    masks: dict[tuple[int, str], np.ndarray] = {}
    blocks = sorted({observation.block for observation in observations})
    for block in blocks:
        selected = [observation for observation in observations if observation.block == block]
        for area, frame in populations.items():
            total = np.zeros(len(frame), dtype=np.int64)
            for observation in selected:
                total += observation.counts[area].sum(axis=0)
            masks[(block, area)] = total >= minimum_spikes
    return masks


def _stack_training(
    observations: list[TraversalObservation], area: str, mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.concatenate([observation.counts[area][:, mask] for observation in observations]),
        np.concatenate([observation.position_bin for observation in observations]),
    )


def decoded_bin_table(
    train: list[TraversalObservation],
    test: list[TraversalObservation],
    masks: dict[str, np.ndarray],
) -> pd.DataFrame:
    """Return one held-out row per 250-ms behavioral window.

    Covariates are removed at this temporal resolution before any traversal x
    position-bin averaging.  Zero-spike observations remain present.
    """
    decoders = {
        area: fit_decoder(*_stack_training(train, area, mask))
        for area, mask in masks.items()
    }
    rows: list[dict] = []
    for observation in test:
        decoded = {
            area: decoders[area].predict(observation.counts[area][:, mask])
            for area, mask in masks.items()
        }
        for temporal_index, position_bin in enumerate(observation.position_bin):
            if int(position_bin) not in CENTRAL_BINS:
                continue
            row = {
                "trial": observation.direction_trial,
                "position_bin": int(position_bin),
                "coordinate": float(observation.coordinate[temporal_index]),
                "speed": float(observation.speed[temporal_index]),
            }
            for area in masks:
                row[f"decoded_{area}"] = float(decoded[area][temporal_index])
                row[f"count_{area}"] = int(
                    np.sum(observation.counts[area][temporal_index, masks[area]])
                )
            rows.append(row)
    return pd.DataFrame(rows)


def residual_design(table: pd.DataFrame, include_population_totals: bool = False) -> np.ndarray:
    bins = table.position_bin.to_numpy(dtype=int)
    bin_effects = np.column_stack([bins == value for value in CENTRAL_BINS]).astype(float)
    centers = (bins + 0.5) / N_POSITION_BINS
    within_bin = table.coordinate.to_numpy(dtype=float) - centers
    speed = table.speed.to_numpy(dtype=float)
    speed = (speed - np.mean(speed)) / (np.std(speed, ddof=1) or 1.0)
    trials = table.trial.to_numpy(dtype=float)
    trials = (trials - np.mean(trials)) / (np.std(trials, ddof=1) or 1.0)
    columns = [
        bin_effects,
        within_bin[:, None],
        (within_bin**2)[:, None],
        (within_bin**3)[:, None],
        speed[:, None],
        (speed**2)[:, None],
        (speed**3)[:, None],
        trials[:, None],
        (trials**2)[:, None],
    ]
    if include_population_totals:
        count_columns = sorted(column for column in table if column.startswith("count_"))
        totals = np.column_stack([np.log1p(table[column].to_numpy(dtype=float)) for column in count_columns])
        totals -= np.mean(totals, axis=0, keepdims=True)
        totals /= np.where(np.std(totals, axis=0, ddof=1) > 0, np.std(totals, axis=0, ddof=1), 1)
        columns.append(totals)
    return np.column_stack(columns)


def residualize(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    beta, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ beta


@dataclass
class ErrorRecord:
    subject: str
    session: str
    block: int
    direction: str
    fold: int
    trials: np.ndarray
    ca1: np.ndarray
    rsc: np.ndarray
    ca3: np.ndarray | None
    observed_r: float
    observed_z: float
    row_count: int


def table_to_record(
    table: pd.DataFrame,
    subject: str,
    session: str,
    block: int,
    direction: str,
    fold: int,
    include_population_totals: bool = False,
) -> ErrorRecord:
    design = residual_design(table, include_population_totals=include_population_totals)
    residuals: dict[str, np.ndarray] = {}
    for area in ("CA1", "RSC", "CA3"):
        column = f"decoded_{area}"
        if column in table:
            signed_error = table[column].to_numpy(dtype=float) - table.coordinate.to_numpy(dtype=float)
            residuals[area] = residualize(signed_error, design)
    aggregate = table[["trial", "position_bin"]].copy()
    for area, values in residuals.items():
        aggregate[area] = values
    aggregate = aggregate.groupby(["trial", "position_bin"], as_index=False).mean()
    trial_values = np.sort(aggregate.trial.unique())
    shape = (len(trial_values), len(CENTRAL_BINS))
    matrices = {area: np.full(shape, np.nan) for area in residuals}
    trial_lookup = {value: index for index, value in enumerate(trial_values)}
    bin_lookup = {value: index for index, value in enumerate(CENTRAL_BINS)}
    for row_index, row in aggregate.reset_index(drop=True).iterrows():
        i = trial_lookup[int(row.trial)]
        j = bin_lookup[int(row.position_bin)]
        for area, values in residuals.items():
            matrices[area][i, j] = values[row_index]
    observed_r = paired_correlation(matrices["CA1"], matrices["RSC"])
    return ErrorRecord(
        subject=subject,
        session=session,
        block=block,
        direction=direction,
        fold=fold,
        trials=trial_values,
        ca1=matrices["CA1"],
        rsc=matrices["RSC"],
        ca3=matrices.get("CA3"),
        observed_r=observed_r,
        observed_z=fisher_z(observed_r),
        row_count=len(aggregate),
    )


def paired_correlation(left: np.ndarray, right: np.ndarray) -> float:
    keep = np.isfinite(left) & np.isfinite(right)
    if np.sum(keep) < 10:
        return np.nan
    x, y = left[keep], right[keep]
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def partial_correlation(
    left: np.ndarray, right: np.ndarray, covariate: np.ndarray
) -> float:
    """Correlation after removing a cubic function of a third decoded deviation."""
    keep = np.isfinite(left) & np.isfinite(right) & np.isfinite(covariate)
    if np.sum(keep) < 10:
        return np.nan
    x, y, z = left[keep], right[keep], covariate[keep]
    z = (z - np.mean(z)) / (np.std(z, ddof=1) or 1.0)
    design = np.column_stack([np.ones(len(z)), z, z**2, z**3])
    x = residualize(x, design)
    y = residualize(y, design)
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def fisher_z(value: float) -> float:
    return float(np.arctanh(np.clip(value, -1 + 1e-12, 1 - 1e-12)))


def record_correlation(record: ErrorRecord, partial_ca3: bool = False) -> float:
    if partial_ca3:
        if record.ca3 is None:
            return np.nan
        return partial_correlation(record.ca1, record.rsc, record.ca3)
    return paired_correlation(record.ca1, record.rsc)


def rotation_draws(
    record: ErrorRecord, shifts: np.ndarray, partial_ca3: bool = False
) -> np.ndarray:
    if len(shifts) == 0:
        return np.empty(0, dtype=float)
    values = np.empty(len(shifts), dtype=float)
    for draw, shift in enumerate(shifts):
        shifted = np.roll(record.rsc, int(shift), axis=0)
        correlation = (
            partial_correlation(record.ca1, shifted, record.ca3)
            if partial_ca3 and record.ca3 is not None
            else paired_correlation(record.ca1, shifted)
        )
        values[draw] = fisher_z(correlation)
    return values


def hierarchy_weights(records: list[ErrorRecord]) -> np.ndarray:
    """Equal fold -> direction -> block -> session-within-mouse -> mouse weights."""
    keys = pd.DataFrame(
        [
            {
                "subject": record.subject,
                "session": record.session,
                "block": record.block,
                "direction": record.direction,
                "fold": record.fold,
            }
            for record in records
        ]
    )
    mice = sorted(keys.subject.unique())
    weights = np.zeros(len(keys), dtype=float)
    for mouse in mice:
        mouse_rows = keys.subject == mouse
        sessions = sorted(keys.loc[mouse_rows, "session"].unique())
        for session in sessions:
            session_rows = mouse_rows & (keys.session == session)
            blocks = sorted(keys.loc[session_rows, "block"].unique())
            for block in blocks:
                block_rows = session_rows & (keys.block == block)
                directions = sorted(keys.loc[block_rows, "direction"].unique())
                for direction in directions:
                    rows = block_rows & (keys.direction == direction)
                    folds = int(np.sum(rows))
                    weights[rows] = (
                        1 / len(mice) / len(sessions) / len(blocks) / len(directions) / folds
                    )
    if not np.isclose(np.sum(weights), 1):
        raise RuntimeError(f"Hierarchy weights sum to {np.sum(weights)}")
    return weights


def mouse_weights(records: list[ErrorRecord], subject: str) -> np.ndarray:
    chosen = [record for record in records if record.subject == subject]
    local = hierarchy_weights(chosen)
    result = np.zeros(len(records), dtype=float)
    result[[index for index, record in enumerate(records) if record.subject == subject]] = local
    return result
