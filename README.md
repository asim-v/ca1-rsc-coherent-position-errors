# Coherent position-decoding deviations in CA1 and RSC

This project asks whether independently decoded spatial estimates in hippocampal
CA1 and deep-layer retrosplenial cortex (RSC) fluctuate together during linear-track
navigation.

The intended empirical statement is deliberately narrow:

> After removing measured position, nonlinear speed effects, and trial-order
> trends, do within-traversal deviations of independently decoded position covary
> between CA1 and RSC?

This is not a test of causal communication, and a positive result would not by
itself establish a subjective position estimate. The design also asks whether any
CA1--RSC relation remains after accounting for the independently decoded CA3
deviation.

## Data separation

- **Discovery:** 12 single-maze sessions from four mice.
- **Sealed confirmation:** the fixed-condition sessions from M01, M02, and M03.
  They are different recordings, but not new animals.

The confirmation files must not be scored until the estimator, complete-group
trial-rotation null, simulation checks, and provenance manifest have been committed.

Data source: DANDI:001695, version `0.260319.2023`.

## Discovery result

The strict discovery estimator found a large separation from the complete-group
trial-rotation reference across 12 sessions:

- primary excess Fisher z: `0.2042`, Monte Carlo `p = 0.0001`, `4/4` mice positive;
- after CA1/RSC/CA3 population-total adjustment: `0.1757`, `p = 0.0001`, `4/4`;
- after cubic adjustment for decoded CA3 deviation: `0.1851`, `p = 0.0001`, `4/4`;
- odd tracking-frame reconstruction: primary excess `0.2024`, `p = 0.0001`, `4/4`.

Eleven of twelve sessions were positive. The exception was M02-20240312, which had
only seven active source-convention RSC units. These are exploratory results.

The estimator:

1. trains independent conditional-multinomial decoders in odd/even traversal folds;
2. retains zero-spike time bins instead of selecting simultaneous activity;
3. removes continuous position, cubic speed, and trial-order effects separately
   within each region;
4. tests the residual correlation against the full product group of within-fold
   cyclic traversal rotations, allowing zero shifts in individual strata.

The sealed fixed-condition sessions have not been opened for this outcome. The
frozen decision rule is in `docs/frozen_confirmation_spec.md`.

## Author

Javier Emilio Bazan Sanchez  
Facultad de Ciencias, Universidad Nacional Autonoma de Mexico  
<bazan@ciencias.unam.mx>
