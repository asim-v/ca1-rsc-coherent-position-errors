# Coherent position-decoding deviations in CA1 and RSC

[![Release](https://img.shields.io/github/v/release/asim-v/ca1-rsc-coherent-position-errors)](https://github.com/asim-v/ca1-rsc-coherent-position-errors/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the complete public analysis package for a reanalysis of
simultaneous CA3, CA1, and deep-layer retrosplenial cortex (RSC) recordings during
linear-track navigation.

> **Status:** public version 1.0.0. The manuscript is a research preprint and has
> not yet undergone external peer review.

## Main result

Independent CA1 and RSC position decoders tended to deviate in the same direction
within the same traversal and spatial bin. The relation remained after accounting
for measured position, nonlinear speed, trial order, population spike totals, and
the independently decoded CA3 deviation.

- **Discovery:** excess Fisher z `0.2042`, Monte Carlo `p = 0.0001`, `4/4` mice
  positive across 12 sessions.
- **Prespecified held-out sessions:** excess Fisher z `0.1465`, `p = 0.0001`,
  `3/3` mice positive across three reserved fixed-condition sessions.
- Both directions, both decoder folds, all six held-out maze blocks, and all
  leave-one-mouse-out estimates were positive.

The claim is deliberately narrow. This is evidence for coordinated decoded-position
deviations in these recordings, not subjective position, causal CA1-to-RSC
communication, or replication in new animals. See
[the confirmation report](docs/confirmation_result.md) and
[the novelty boundary](docs/novelty_boundary.md).

## Manuscript

- [Rendered manuscript (PDF)](output/pdf/ca1_rsc_coherent_position_deviations.pdf)
- [LaTeX source](manuscript/main.tex)
- [Main figure](manuscript/figure1.png)

## Verify the archived result

Python 3.10 is recommended.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pytest -q
python scripts/verify_confirmation.py
```

The verification script independently rebuilds the mouse-level hierarchy and all
four held-out statistics from the archived record tables and null arrays. It does
not require the raw NWB files.

For full raw-data reproduction, including the protected one-shot confirmation
workflow, see [Reproducibility](docs/reproducibility.md).

## Data

Raw recordings are not duplicated in this repository. They are publicly available
from [DANDI:001695, version 0.260319.2023](https://dandiarchive.org/dandiset/001695/0.260319.2023).
The exact input filenames and SHA-256 hashes used for confirmation are recorded in
[`outputs/provenance/frozen_manifest.json`](outputs/provenance/frozen_manifest.json).

## Repository contents

```text
docs/           frozen protocol, result interpretation, and reproduction guide
manuscript/     LaTeX source and publication figure
output/pdf/     rendered manuscript
outputs/        archived discovery, calibration, and confirmation results
scripts/        analysis, calibration, figure, provenance, and verification code
tests/          unit tests for the core estimator and null hierarchy
```

## Citation

Citation metadata are available in [`CITATION.cff`](CITATION.cff). GitHub's
**Cite this repository** button can export the citation directly.

## Author

Javier Emilio Bazan Sanchez  
Facultad de Ciencias, Universidad Nacional Autonoma de Mexico  
<bazan@ciencias.unam.mx>

## License

Code and repository-authored documentation are released under the [MIT License](LICENSE).
The source NWB recordings remain subject to the terms supplied by their original
authors and DANDI.
