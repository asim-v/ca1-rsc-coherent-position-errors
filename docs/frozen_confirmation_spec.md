# Frozen held-out confirmation specification v1.0

Status: **to be committed before any aligned score is calculated in the fixed
sessions**.

## Sealed inputs

Exactly three fixed-condition sessions are reserved:

- `sub-M01_ses-20240318T100000_behavior+ecephys.nwb`
- `sub-M02_ses-20240318T100000_behavior+ecephys.nwb`
- `sub-M03_ses-20240624T100000_behavior+ecephys.nwb`

Each contains two maze blocks. File hashing and interval/unit-count feasibility
checks do not constitute opening the outcome. No position-deviation correlation,
aligned scatter, or trial-rotation score may be calculated before the frozen code,
calibration, and manifest are committed.

These are held-out recordings, not held-out animals.

## Exact primary estimator

The estimator is `discovery-v0.2-temporal-residuals-complete-product-group` in
`scripts/coherent_error_core.py`:

- 250-ms nonoverlapping windows; speed >2.5 cm/s;
- block-level physical coordinates and endpoint-defined traversals;
- PAVA coordinate from even-indexed tracking frames;
- CA1 Pyramidal Cell and source-convention non-narrow RSC populations;
- block-active units with at least 20 locomotor spikes; minimum ten CA1/five RSC;
- at least 16 accepted traversals in each direction in every block;
- conditional-multinomial decoder, 24 bins, pseudocount 0.5, Gaussian sigma one
  bin, uniform prior;
- direction-specific odd/even traversal cross-fitting;
- score only true bins 4--19;
- retain every eligible behavioral window, including zero-spike windows;
- separately residualize each area's 250-ms signed errors against true-bin fixed
  effects, cubic within-bin continuous position, cubic speed, and quadratic trial
  order;
- average residuals within traversal x true bin;
- Pearson correlation per block x direction x fold, Fisher transformed;
- equal fold -> direction -> block -> session -> mouse aggregation.

All six blocks must remain eligible. There are no score-based exclusions and no
reduced-cohort replacement test.

## Primary trial-rotation reference

Use seed `20260830` and exactly 9,999 draws. In every block x direction x held-out
fold, draw one RSC traversal shift uniformly from the complete cyclic group
`0, ..., n_trials - 1`. Shifts are independent across strata and draws. Zero shifts
in any number of strata, including all strata, are allowed. Preserve the entire RSC
trial row, its missing-bin mask, and all position-bin columns. Aggregate each draw
through the exact primary hierarchy. Use

`p = (1 + count(T_null >= T_observed)) / 10000`.

## Frozen confirmation decision

The primary is confirmed only if every condition holds:

1. observed hierarchical Fisher z minus mean null is positive;
2. one-sided Monte Carlo `p <= 0.05`;
3. all three mouse-level observed-minus-null margins are positive;
4. all three leave-one-mouse-out margins are positive;
5. both direction-level margins and both fold-level margins are positive.

No secondary analysis can rescue a failed primary.

## Mandatory claim controls

Run after the primary regardless of its sign, using independent frozen seed
namespaces and the same complete product group.

1. **Population totals:** add log CA1, RSC, and CA3 spike totals to each area's
   temporal residual model. A statement beyond population vigor requires positive
   aggregate margin, `p <= 0.05`, and all three mouse margins positive.
2. **CA3 partial:** remove a cubic function of the aligned CA3 decoded residual from
   both CA1 and RSC and repeat the RSC trial-rotation test. A statement beyond the
   decoded CA3 deviation requires positive aggregate margin, `p <= 0.05`, and all
   three mouse margins positive.
3. **Tracking parity:** repeat the estimator using odd-indexed tracking frames. It
   must retain a positive aggregate margin and at least two of three positive mice.
4. **Literal RSC label:** repeat with `cell_type == "Pyramidal Cell"` and report
   evaluability and all outcomes. This is a low-power population-definition
   sensitivity and does not replace the source-convention primary.

The continuous-coordinate and fixed-bin-center residuals must agree within
`rtol=1e-10, atol=1e-10`; disagreement is an implementation failure because the
continuous within-bin coordinate is in the nuisance-design span.

## Interpretation matrix

- Primary failure: no support for CA1--RSC coherent position deviations in the
  held-out fixed sessions under this estimator.
- Primary pass, population-total failure: relation may be coupled to population
  vigor; withhold the beyond-count claim.
- Primary pass, CA3-partial failure: relation may reflect a broader hippocampal
  deviation shared with CA3; withhold CA1--RSC-specific language.
- Primary and both claim controls pass: support the narrow candidate contribution
  in `docs/novelty_boundary.md`.

Never use "subjective position," causal direction, or new-animal replication.

