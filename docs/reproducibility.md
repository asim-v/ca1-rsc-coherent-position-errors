# Reproducibility guide

This project separates three levels of reproducibility: checking the archived
result, regenerating exploratory outputs from raw data, and replaying the frozen
one-shot confirmation from its pre-result commit.

## 1. Environment

Python 3.10 and the exact package versions in `requirements.txt` are recommended.

```bash
git clone https://github.com/asim-v/ca1-rsc-coherent-position-errors.git
cd ca1-rsc-coherent-position-errors
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## 2. Verify the published result without raw data

```bash
python -m pytest -q
python scripts/verify_confirmation.py
```

The verifier reads the archived record-level CSV files and null arrays, rebuilds
the equal-fold, direction, block, session-within-mouse, and mouse hierarchy, and
checks every reported held-out statistic against `outputs/confirmation/summary.json`.

Expected primary values:

- observed Fisher z: `0.16600574744773336`;
- mean trial-rotation Fisher z: `0.01948046835787935`;
- excess Fisher z: `0.14652527908985402`;
- one-sided Monte Carlo p: `0.0001`;
- positive mice: `3/3`.

## 3. Obtain the raw recordings

The recordings are available from
[DANDI:001695, version 0.260319.2023](https://dandiarchive.org/dandiset/001695/0.260319.2023).
Place the required NWB files in one directory, such as `data/raw/`. That directory
is ignored by Git and must not be committed.

Discovery uses the 12 sessions listed in `scripts/run_discovery.py`. Confirmation
uses exactly these three files:

```text
sub-M01_ses-20240318T100000_behavior+ecephys.nwb
sub-M02_ses-20240318T100000_behavior+ecephys.nwb
sub-M03_ses-20240624T100000_behavior+ecephys.nwb
```

Their exact sizes and SHA-256 hashes are frozen in
`outputs/provenance/frozen_manifest.json`.

## 4. Reproduce discovery

Use a new output directory so the archived result remains unchanged.

```bash
python scripts/run_discovery.py \
  --raw-dir data/raw \
  --output-dir tmp/reproduced_discovery \
  --null-draws 9999 \
  --seed 20260820

python scripts/run_discovery.py \
  --raw-dir data/raw \
  --output-dir tmp/reproduced_discovery_tracking_odd \
  --null-draws 9999 \
  --seed 20260820 \
  --tracking-parity 1
```

## 5. Replay the frozen confirmation

The published branch already contains the write-once sentinel and confirmation
outputs, so `run_confirmation.py` intentionally refuses to run there. Replay the
analysis in a separate worktree at the frozen execution commit:

```bash
git worktree add ../ca1-rsc-confirmation-replay a192a80a9785e048e92b90808bf91d77a8914339
cd ../ca1-rsc-confirmation-replay
python -m venv .venv
# Activate the environment as above.
python -m pip install -r requirements.txt
python scripts/run_confirmation.py --raw-dir /path/to/nwb/files
python scripts/verify_confirmation.py
```

The script verifies the frozen manifest before opening the aligned held-out score,
creates a write-once sentinel, and refuses overwriting or a second run.

## 6. Rebuild the manuscript figure

From the current release checkout:

```bash
python scripts/make_manuscript_figure.py
```

The figure is generated only from archived CSV and NPZ outputs. The manuscript is
compiled from `manuscript/main.tex` with a standard `pdflatex` installation.

## Integrity anchors

- Analysis commit: `e35757987a7e7acca7d5000f75b088f6a4e36a6b`
- Confirmation execution commit: `a192a80a9785e048e92b90808bf91d77a8914339`
- Frozen input-manifest self-hash:
  `b594fa568f4ffdc174cdf040220fc468d8cb5fc447bb2397248f30a09ed71a44`
- Confirmation result-manifest self-hash:
  `8db8e516db42f4876fa9c2cd44874c9ae130ee2d71a668373e037c5ab34e5ec5`
- Primary seed: `20260830`
- Primary trial-rotation draws: `9,999`
