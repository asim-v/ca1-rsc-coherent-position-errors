from __future__ import annotations

import numpy as np
import pandas as pd

from coherent_error_core import (
    ErrorRecord,
    fit_decoder,
    hierarchy_weights,
    paired_correlation,
    pava,
    residual_design,
    residualize,
    rotation_draws,
)


def test_pava_is_monotone_and_pools_violations():
    result = pava(np.array([0.0, 0.4, 0.2, 0.8]))
    assert np.all(np.diff(result) >= 0)
    assert np.allclose(result, [0.0, 0.3, 0.3, 0.8])


def test_zero_spike_decoder_returns_uniform_prior_mean():
    decoder = fit_decoder(
        np.array([[2, 0], [0, 2], [1, 1]]), np.array([0, 12, 23])
    )
    prediction = decoder.predict(np.zeros((3, 2)))
    assert np.allclose(prediction, 0.5)


def test_continuous_and_fixed_center_errors_have_same_residual():
    rng = np.random.default_rng(4)
    bins = np.repeat(np.arange(4, 20), 8)
    centers = (bins + 0.5) / 24
    within = rng.uniform(-0.018, 0.018, len(bins))
    table = pd.DataFrame(
        {
            "trial": np.tile(np.arange(8), 16),
            "position_bin": bins,
            "coordinate": centers + within,
            "speed": rng.uniform(3, 35, len(bins)),
            "decoded_CA1": rng.uniform(0, 1, len(bins)),
        }
    )
    design = residual_design(table)
    decoded = table.decoded_CA1.to_numpy()
    continuous = residualize(decoded - table.coordinate.to_numpy(), design)
    fixed = residualize(decoded - centers, design)
    assert np.allclose(continuous, fixed, rtol=1e-10, atol=1e-10)


def test_zero_rotation_equals_observed_correlation():
    rng = np.random.default_rng(8)
    left = rng.normal(size=(9, 16))
    right = rng.normal(size=(9, 16))
    correlation = paired_correlation(left, right)
    record = ErrorRecord(
        subject="M01",
        session="S1",
        block=1,
        direction="LtoR",
        fold=0,
        trials=np.arange(9),
        ca1=left,
        rsc=right,
        ca3=None,
        observed_r=correlation,
        observed_z=float(np.arctanh(correlation)),
        row_count=left.size,
    )
    assert np.allclose(rotation_draws(record, np.array([0]))[0], record.observed_z)


def test_hierarchy_weights_equalize_mice():
    records = []
    for subject, sessions in (("M01", ("A", "B")), ("M02", ("C",))):
        for session in sessions:
            for direction in ("LtoR", "RtoL"):
                for fold in (0, 1):
                    matrix = np.arange(96, dtype=float).reshape(6, 16)
                    records.append(
                        ErrorRecord(
                            subject=subject,
                            session=session,
                            block=1,
                            direction=direction,
                            fold=fold,
                            trials=np.arange(6),
                            ca1=matrix,
                            rsc=matrix,
                            ca3=None,
                            observed_r=1.0,
                            observed_z=1.0,
                            row_count=96,
                        )
                    )
    weights = hierarchy_weights(records)
    assert np.isclose(weights.sum(), 1)
    assert np.isclose(weights[[r.subject == "M01" for r in records]].sum(), 0.5)
    assert np.isclose(weights[[r.subject == "M02" for r in records]].sum(), 0.5)

