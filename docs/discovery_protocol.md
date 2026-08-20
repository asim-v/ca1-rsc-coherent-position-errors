# Discovery protocol v0.1

## Question

Do moment-to-moment deviations of independently decoded position covary between
CA1 and RSC during locomotion on a linear track?

## Populations

- CA1: `cell_area == "CA1"` and `cell_type == "Pyramidal Cell"`.
- RSC: `cell_area == "RSC"` and `cell_type != "Narrow Interneuron"`.
- CA3 secondary: `cell_area == "CA3"` and `cell_type == "Pyramidal Cell"`.

Units are selected within a maze block using only their own locomotor spike count,
before any cross-area statistic is calculated. The primary decoder retains units
with at least 20 spikes across accepted traversals. A block requires at least ten
CA1 and five RSC units.

## Behavioral and neural bins

Traversal endpoints are 0.15 and 0.85 in a block-level physical coordinate. Each
accepted traversal covers at least 0.70 of the track and lasts at most 120 seconds.
Tracking is monotonized using PAVA fitted to even-indexed tracking frames within
each traversal. Neural observations are nonoverlapping 250-ms windows and require
speed greater than 2.5 cm/s. All such windows are retained, including windows with
zero CA1 or RSC spikes.

The decoder is trained on all 24 physical bins and scored only in bins 4--19.

## Decoder and cross-fitting

Separate conditional-multinomial decoders are fitted for CA1 and RSC, direction,
block, and traversal fold. Cell probabilities are estimated from the training
traversals with a fixed pseudocount and one-bin Gaussian smoothing. The spatial
prior is uniform. Odd traversals train even-traversal predictions and vice versa.

The signed error in every held-out 250-ms window is decoded posterior-mean position
minus the continuous PAVA coordinate. Measured covariates are removed at this
temporal resolution. Residuals are then averaged within traversal x true-position
bin, so within-bin speed mixtures cannot hide behind a single mean covariate.

## Removal of measured covariates

Within each block x direction x held-out fold, each area's signed errors are
residualized separately against the same fixed design:

- true-position-bin fixed effects;
- within-bin continuous position offset and its square and cube;
- standardized speed and its square and cube;
- normalized traversal order and its square.

Because continuous position is in the design span, the common subtraction of the
tracking coordinate is removed algebraically before the cross-area correlation.

## Statistic

The record-level statistic is Pearson correlation between CA1 and RSC residuals,
Fisher transformed before aggregation. Aggregation is equal fold, direction,
block, session within mouse, and finally equal mouse.

## Trial-rotation reference

RSC traversal identities are cyclically shifted within every block x direction x
held-out fold while true-position bins remain fixed. Each stratum shift is drawn
uniformly from `0, ..., n_trials - 1`; zero shifts in some or all strata are valid
members of the product group. The observed identity alignment is scored separately.
Monte Carlo p-values use the plus-one rule and include ties.

The primary effect size is observed hierarchical Fisher z minus the mean of the
rotation distribution. Circular traversal exchangeability is an observational
surrogate assumption, so linear and quadratic trial-order terms are removed and
temporal-drift diagnostics are mandatory.

## Prespecified discovery controls

- Refit the residual model with log-transformed CA1, RSC, and CA3 population spike
  totals from every 250-ms window.
- Remove a cubic function of the independently decoded CA3 residual from both CA1
  and RSC before their correlation; every rotation repeats that partial correlation
  with CA1 and CA3 kept aligned.

## Claim boundary

A positive discovery result supports only a candidate relation. Confirmation in
sealed sessions would support the statement that within-traversal signed deviations
of independently decoded position covary between CA1 and RSC in these recordings.
It would not establish CA1-to-RSC causality, a subjective position estimate, or
population-level generality beyond the recorded animals.
