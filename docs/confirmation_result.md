# Held-out confirmation result

## Outcome

The frozen held-out-session test confirmed coherent CA1--RSC position-decoding
deviations in the three fixed-condition recordings.

| Estimator | Observed z | Mean rotated z | Excess z | Monte Carlo p | Positive mice |
|---|---:|---:|---:|---:|---:|
| Primary | 0.1660 | 0.0195 | 0.1465 | 0.0001 | 3/3 |
| Population-total adjusted | 0.1495 | 0.0164 | 0.1331 | 0.0001 | 3/3 |
| CA3-deviation adjusted | 0.1586 | 0.0194 | 0.1392 | 0.0001 | 3/3 |
| Odd tracking frames | 0.1653 | 0.0129 | 0.1524 | 0.0001 | 3/3 |

All three leave-one-mouse-out estimates, both directions, both decoder folds, and
all six maze blocks had positive observed-minus-rotation margins. Nineteen of 24
block x direction x fold records were positive.

The discovery and confirmation effect sizes were similar: primary excess Fisher z
was 0.2042 in the 12 exploratory sessions and 0.1465 in the three sealed sessions.

## Plain-language interpretation

When the CA1 decoder placed the animal slightly ahead of or behind its measured
position, the independently trained RSC decoder tended to deviate in the same
direction within the same traversal and spatial bin. This coordination persisted
after removing measured position, nonlinear speed, trial order, population spike
totals, and a cubic function of the simultaneously decoded CA3 deviation.

The calculation does not arise from subtracting the same tracking value twice:
continuous-position and fixed-bin-center errors become identical after the frozen
position design is projected out, with maximum discrepancy `3.81e-12`. Rebuilding
the coordinate from the other half of tracking frames also preserved the result.

## Important limits

- The sealed sessions came from three mice already represented in discovery. This
  is held-out-session confirmation, not replication in new animals.
- The NWB files expose position and speed but not licking or behavioral choice
  errors. The result cannot be called a subjective position estimate.
- The analysis is correlational and does not establish CA1-to-RSC direction or
  causality.
- The literal `RSC Pyramidal Cell` sensitivity was not evaluable across all six
  blocks: M01 had four eligible literal-labeled cells per block, below the frozen
  minimum of five. The supported population is therefore the source-convention
  non-narrow RSC population, which includes wide-spiking cells.
- Circular traversal rotation is an observational surrogate reference, not a
  randomized intervention.

## Frozen provenance

- Analysis commit: `e35757987a7e7acca7d5000f75b088f6a4e36a6b`
- Execution commit: `a192a80a9785e048e92b90808bf91d77a8914339`
- Input manifest: `b594fa568f4ffdc174cdf040220fc468d8cb5fc447bb2397248f30a09ed71a44`
- Seed: `20260830`
- Trial-rotation draws: `9,999`

