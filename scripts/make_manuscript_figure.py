#!/usr/bin/env python3
"""Build the main manuscript figure from archived discovery/confirmation tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    discovery = pd.read_csv(repo / "outputs/discovery/primary_mouse_scores.csv")
    confirmation = pd.read_csv(repo / "outputs/confirmation/primary_mouse_scores.csv")
    records = pd.read_csv(repo / "outputs/confirmation/primary_record_scores.csv")
    controls = []
    for cohort, root in (
        ("Discovery", repo / "outputs/discovery"),
        ("Held-out", repo / "outputs/confirmation"),
    ):
        for variant, label in (
            ("primary", "Primary"),
            ("population_totals", "Population\ntotals"),
            ("ca3_partial", "CA3\npartial"),
        ):
            table = pd.read_csv(root / f"{variant}_mouse_scores.csv")
            for row in table.itertuples():
                controls.append(
                    {
                        "cohort": cohort,
                        "variant": label,
                        "subject": row.subject,
                        "margin": row.excess_fisher_z,
                    }
                )
        if cohort == "Held-out":
            table = pd.read_csv(root / "tracking_odd_mouse_scores.csv")
        else:
            table = pd.read_csv(repo / "outputs/discovery_tracking_odd/primary_mouse_scores.csv")
        for row in table.itertuples():
            controls.append(
                {
                    "cohort": cohort,
                    "variant": "Tracking\nodd",
                    "subject": row.subject,
                    "margin": row.excess_fisher_z,
                }
            )
    controls = pd.DataFrame(controls)

    figure = plt.figure(figsize=(13.2, 6.8), constrained_layout=True)
    grid = figure.add_gridspec(2, 4, height_ratios=[0.88, 1.12])
    ax_a = figure.add_subplot(grid[0, :2])
    ax_b = figure.add_subplot(grid[0, 2])
    ax_c = figure.add_subplot(grid[0, 3])
    ax_d = figure.add_subplot(grid[1, :2])
    ax_e = figure.add_subplot(grid[1, 2:])

    # A: conceptual definition.
    ax_a.plot([0.05, 0.95], [0, 0], color="#4c566a", linewidth=5, solid_capstyle="round")
    actual, ca1, rsc = 0.52, 0.61, 0.59
    ax_a.scatter([actual], [0], s=170, color="#1f2937", zorder=3)
    ax_a.scatter([ca1], [0.20], s=145, color="#c45a4a", zorder=3)
    ax_a.scatter([rsc], [-0.20], s=145, color="#3978a8", zorder=3)
    ax_a.annotate("measured position", (actual, 0), (actual - 0.02, 0.37), ha="center",
                  arrowprops=dict(arrowstyle="->", color="#1f2937"), fontsize=10)
    ax_a.annotate("CA1 decoded", (ca1, 0.20), (0.73, 0.30), ha="center",
                  arrowprops=dict(arrowstyle="->", color="#c45a4a"), color="#9d3f33", fontsize=10)
    ax_a.annotate("RSC decoded", (rsc, -0.20), (0.73, -0.33), ha="center",
                  arrowprops=dict(arrowstyle="->", color="#3978a8"), color="#2d628c", fontsize=10)
    ax_a.text(0.57, 0.04, "same-signed deviation", ha="left", va="bottom", fontsize=10)
    ax_a.set(xlim=(0, 1), ylim=(-0.48, 0.48))
    ax_a.axis("off")
    ax_a.set_title("A  Within-traversal decoded-position deviations", loc="left", fontweight="bold")

    # B: mouse margins in discovery and held-out sessions.
    colors = {"Discovery": "#8da6c4", "Held-out": "#c45a4a"}
    for x, (label, table) in enumerate((("Discovery", discovery), ("Held-out", confirmation))):
        jitter = np.linspace(-0.09, 0.09, len(table))
        ax_b.scatter(x + jitter, table.excess_fisher_z, s=48, color=colors[label], edgecolor="white", zorder=3)
        ax_b.hlines(table.excess_fisher_z.mean(), x - 0.22, x + 0.22, color="#222", linewidth=2)
    ax_b.axhline(0, color="#777", linewidth=1, linestyle="--")
    ax_b.set_xticks([0, 1], ["Discovery\n12 sessions", "Held-out\n3 sessions"])
    ax_b.set_ylabel("Observed - rotated Fisher z")
    ax_b.set_title("B  Mouse-level replication", loc="left", fontweight="bold")

    # C: every held-out block.
    block = records.groupby(["subject", "block"]).agg(
        observed=("observed_fisher_z", "mean"), null=("null_mean_fisher_z", "mean")
    )
    block["margin"] = block.observed - block.null
    labels = [f"{mouse}\nB{number}" for mouse, number in block.index]
    ax_c.bar(np.arange(len(block)), block.margin, color="#5f86ad")
    ax_c.axhline(0, color="#777", linewidth=1)
    ax_c.set_xticks(np.arange(len(block)), labels, fontsize=8)
    ax_c.set_ylabel("Observed - rotated Fisher z")
    ax_c.set_title("C  Six held-out blocks", loc="left", fontweight="bold")

    # D: claim-control margins with individual mice.
    variants = ["Primary", "Population\ntotals", "CA3\npartial", "Tracking\nodd"]
    offsets = {"Discovery": -0.16, "Held-out": 0.16}
    for cohort in ("Discovery", "Held-out"):
        selected = controls[controls.cohort == cohort]
        means = []
        for x, variant in enumerate(variants):
            values = selected[selected.variant == variant].margin.to_numpy()
            means.append(values.mean())
            jitter = np.linspace(-0.05, 0.05, len(values))
            ax_d.scatter(x + offsets[cohort] + jitter, values, s=25, color=colors[cohort], alpha=0.88)
        ax_d.plot(np.arange(4) + offsets[cohort], means, "o-", color=colors[cohort],
                  label=cohort, linewidth=2, markersize=6)
    ax_d.axhline(0, color="#777", linewidth=1, linestyle="--")
    ax_d.set_xticks(np.arange(4), variants)
    ax_d.set_ylabel("Mouse observed - rotated Fisher z")
    ax_d.legend(frameon=False, ncol=2, loc="upper right")
    ax_d.set_title("D  Prespecified claim controls", loc="left", fontweight="bold")

    # E: confirmation primary randomization distribution.
    arrays = np.load(repo / "outputs/confirmation/primary_nulls.npz")
    summary = pd.read_json(repo / "outputs/confirmation/summary.json", typ="series")
    observed = float(summary["results"]["primary"]["observed_fisher_z"])
    ax_e.hist(arrays["overall"], bins=50, color="#b9c8d8", edgecolor="none")
    ax_e.axvline(observed, color="#b33e36", linewidth=2.5, label=f"Observed z = {observed:.3f}")
    ax_e.set_xlabel("Hierarchical Fisher z")
    ax_e.set_ylabel("Trial-rotation draws")
    ax_e.legend(frameon=False)
    ax_e.set_title("E  Held-out complete-group reference (p = 0.0001)", loc="left", fontweight="bold")

    for axis in (ax_b, ax_c, ax_d, ax_e):
        axis.spines[["top", "right"]].set_visible(False)
    output = repo / "manuscript/figure1.png"
    figure.savefig(output, dpi=320, facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()

