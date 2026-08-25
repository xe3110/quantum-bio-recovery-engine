"""Publication figures for the extended MS combination screen.

Run after ``experiments.ms.run_extended_screen``:
    python -m experiments.ms.plot_screen
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "experiments/ms/results"
FIGURES = ROOT / "figures"

# Validated categorical slots 1-3 (light mode, all-pairs safe).
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#dcdcd8"
SURFACE = "#fcfcfb"


def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.figure.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=8, length=3)
    ax.grid(axis="x", color=GRID, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)


def fig_strata(report):
    """Headline: mechanism strata ranked by median pair score, with bootstrap CIs."""
    strata = list(reversed(report["mechanism_strata"]["strata"]))
    labels = [s["axis_pair"].replace("_", " ") for s in strata]
    medians = [s["median_score"] for s in strata]
    los = [s["median_score"] - s["median_bootstrap_95ci"][0] for s in strata]
    his = [s["median_bootstrap_95ci"][1] - s["median_score"] for s in strata]
    significant = [s["q_greater"] is not None and s["q_greater"] < 0.05 for s in strata]

    fig, ax = plt.subplots(figsize=(9.6, 5.6))
    y = range(len(strata))
    ax.errorbar(medians, y, xerr=[los, his], fmt="none", ecolor=GRID, elinewidth=2, capsize=0, zorder=1)
    ax.scatter(
        medians, y, s=58, zorder=3,
        color=[BLUE if sig else "#b9c6d6" for sig in significant],
        edgecolor=SURFACE, linewidth=2,
    )
    # Anchor labels past the CI whisker so they never sit on top of the bars.
    for i, (stratum, high) in enumerate(zip(strata, [s["median_bootstrap_95ci"][1] for s in strata])):
        ax.annotate(f"{stratum['median_score']:.3f}  (n={stratum['n_pairs']})", (high, i),
                    xytext=(10, 0), textcoords="offset points", va="center",
                    fontsize=7.5, color=INK_2)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8.5, color=INK)
    ax.set_xlabel("median pair priority score (bars: bootstrap 95% CI)", fontsize=9, color=INK_2)
    ax.set_title(
        "CNS-innate + remyelination pairings outscore\n"
        "two peripheral immunomodulators",
        fontsize=12.5, color=INK, loc="left", pad=12, fontweight="semibold",
    )
    ax.annotate("filled = enriched at FDR q < 0.05", xy=(0.0, -0.155), xycoords="axes fraction",
                fontsize=7.5, color=INK_2)
    _style(ax)
    ax.margins(x=0.10)
    ax.set_xlim(right=ax.get_xlim()[1] + 0.11)
    fig.tight_layout()
    out = FIGURES / "fig_mechanism_strata.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def fig_rank_stability(report):
    """The central caveat: individual pair ranks are not stable."""
    items = list(report["stability"]["intervals"].items())[:20][::-1]
    labels = [k for k, _ in items]
    med = [v["rank_median"] for _, v in items]
    lo = [v["rank_median"] - v["rank_95ci"][0] for _, v in items]
    hi = [v["rank_95ci"][1] - v["rank_median"] for _, v in items]

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    y = range(len(items))
    ax.errorbar(med, y, xerr=[lo, hi], fmt="none", ecolor=ORANGE, elinewidth=2.2, alpha=0.55, capsize=0)
    ax.scatter(med, y, s=48, color=ORANGE, edgecolor=SURFACE, linewidth=2, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8, color=INK)
    ax.set_xscale("log")
    ax.set_xlabel("rank among primary pairs (log scale; bars: bootstrap 95% CI)", fontsize=9, color=INK_2)
    jac = report["stability"][f"top{report['stability']['top_k']}_jaccard_mean"]
    ax.set_title(
        "Individual pair ranks are not stable under target-effect uncertainty\n"
        f"top-{report['stability']['top_k']} bootstrap Jaccard = {jac:.2f}; "
        "read the strata, not the leaderboard",
        fontsize=12, color=INK, loc="left", pad=14, fontweight="semibold",
    )
    _style(ax)
    fig.tight_layout()
    out = FIGURES / "fig_rank_stability.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def fig_tradeoff(rows):
    """Efficacy against combined safety burden, with the Pareto front marked."""
    primary = [r for r in rows if not r["excluded_from_primary_ranking"]]
    front = [r for r in primary if r["pareto_optimal"]]
    rest = [r for r in primary if not r["pareto_optimal"]]

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    ax.scatter([r["safety_union"] for r in rest], [r["reversal_efficiency"] for r in rest],
               s=13, color="#c9cdc9", alpha=0.55, linewidth=0, label=f"dominated (n={len(rest)})")
    ax.scatter([r["safety_union"] for r in front], [r["reversal_efficiency"] for r in front],
               s=30, color=AQUA, edgecolor=SURFACE, linewidth=0.8,
               label=f"Pareto-optimal (n={len(front)})", zorder=3)

    # Label the corners of the front rather than the highest-scoring pairs, which
    # cluster together and overplot: most efficacious, lowest burden, and the knee.
    highlights = {
        "highest efficacy": max(front, key=lambda r: r["reversal_efficiency"]),
        "lowest burden": min(front, key=lambda r: (r["safety_union"], -r["reversal_efficiency"])),
        "best trade-off": max(front, key=lambda r: r["reversal_efficiency"] - r["safety_union"]),
    }
    offsets = {"highest efficacy": (8, 6), "lowest burden": (10, -4), "best trade-off": (8, -12)}
    for tag, row in highlights.items():
        ax.annotate(
            f"{row['drug_a']} + {row['drug_b']}\n({tag})",
            (row["safety_union"], row["reversal_efficiency"]),
            xytext=offsets[tag], textcoords="offset points", fontsize=7,
            color=INK, linespacing=1.35,
        )

    ax.set_xlabel("combined safety burden (weighted union across 7 risk domains)", fontsize=9, color=INK_2)
    ax.set_ylabel("reversal efficiency", fontsize=9, color=INK_2)
    ax.set_title(
        "Efficacy against safety burden across all primary pairs\n"
        "the Pareto front survives any choice of weights",
        fontsize=12, color=INK, loc="left", pad=14, fontweight="semibold",
    )
    legend = ax.legend(frameon=False, fontsize=8, loc="lower right")
    for text in legend.get_texts():
        text.set_color(INK_2)
    _style(ax)
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.9)
    fig.tight_layout()
    out = FIGURES / "fig_efficacy_safety_tradeoff.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def main() -> None:
    import csv

    report = json.loads((RESULTS / "ms_extended_screen.json").read_text())
    with (RESULTS / "ms_pairs_full.csv").open(newline="") as handle:
        rows = []
        for row in csv.DictReader(handle):
            for key in ("safety_union", "reversal_efficiency", "priority_score"):
                row[key] = float(row[key])
            for key in ("excluded_from_primary_ranking", "pareto_optimal"):
                row[key] = row[key] == "True"
            rows.append(row)

    FIGURES.mkdir(exist_ok=True)
    for out in (fig_strata(report), fig_rank_stability(report), fig_tradeoff(rows)):
        print(f"Wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
