#!/usr/bin/env python3
"""
Goldbach z/rho cellular automata  — version 2.

WHAT THIS SCRIPT DOES (plain English)
======================================
For every even number N from 4 to max_n the script computes two quantities:

    r_G(N)  — exact number of unordered prime pairs (p, q) with p + q = N
    h(N)    — Hardy-Littlewood heuristic prediction for r_G(N)

From these it derives the calibrated normalized residual

    z_h = (r_G(N) − h_cal(N)) / sqrt(c · h_cal(N))

which is analogous to a z-score: values near 0 mean "N behaves as the
heuristic predicts", large positive values mean "more prime pairs than
expected", and large negative values mean "fewer than expected".

ENCODING MAP  (how z_h becomes a discrete state label)
=======================================================
z_h is a continuous value. The automaton needs discrete states, so z_h is
rounded into fixed-width bins of size `step` and clipped at ±`clip`.

Run with --glossary to print the full encoding table for the default
parameters, or inspect the legend CSV written alongside the other outputs.

    z_h < −clip        → "<=−clip"   (well below expected)
    −clip ≤ z_h < 0    → e.g. "−1.0", "−0.5"  (below expected)
    z_h ≈ 0            → "+0.0"   (on target)
    0 < z_h ≤ clip     → e.g. "+0.5", "+1.0"  (above expected)
    z_h ≥ clip         → ">=+clip"   (well above expected)

SPATIAL STRUCTURE  (the rho30 dimension)
=========================================
rho30 = N mod 30. There are exactly 8 residue classes of even N mod 30
that can participate in Goldbach sums. Each rho30 class is one "cell" in
the 1-D automaton.

AUTOMATON RULES  (neighborhood notation)
=========================================
A rule is written as:

    left_state | center_state | right_state  →  predicted_next_center

where left/center/right are the state labels of adjacent rho30 cells at
time t, and the rule predicts the center cell's state at time t+1.

CANDIDATE SCORE  (what makes a good discretization)
=====================================================
    score = accuracy − 0.15 × entropy − 0.10 × dominant_share

    accuracy        — fraction of transitions correctly predicted by the mode rule
    entropy         — normalized Shannon entropy of rule outputs
                      (0 = fully deterministic, 1 = fully random)
    dominant_share  — fraction of cells in the most common state;
                      penalizes discretizations that trivially collapse
                      everything into one bucket

This is a visualization and operator-diagnostics tool, not a proof engine.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from shattering_mirrors import (
    build_dataset,
    ensure_parent_dir,
    residual_label_sort_key,
    z_label,
)


# ─────────────────────────────────────────────────────────────────────────────
# Symbol glossary
# ─────────────────────────────────────────────────────────────────────────────

SYMBOL_GLOSSARY: list[tuple[str, str]] = [
    ("N",              "Even integer being analyzed (Goldbach target)"),
    ("r_G(N)",         "Exact unordered Goldbach pair count: #{(p,q): p+q=N, p≤q, p,q prime}"),
    ("h(N)",           "Hardy-Littlewood heuristic for r_G(N): C₂·G(N)·∑ 1/(log k · log(N−k))"),
    ("C₂",             "Twin-prime constant ≈ 0.6602:  ∏_{p>2} (1 − 1/(p−1)²)"),
    ("G(N)",           "Local Goldbach boost:  ∏_{p|N, p>2} (p−1)/(p−2)"),
    ("h_cal(N)",       "Calibrated heuristic: h(N)·exp(α + β/log N), fitted by least squares"),
    ("α, β",           "Calibration coefficients fitted from the dataset"),
    ("c",              "Empirical variance scale: fitted so Var(r_G) ≈ c·h_cal"),
    ("z_h",            "Normalized residual: (r_G − h_cal) / √(c·h_cal)"),
    ("state_label",    "Discrete z_h bucket, e.g. '+0.5', '<=−2.0' (see encoding legend)"),
    ("rho30",          "N mod 30 — the spatial cell index in the automaton (8 possible values)"),
    ("step",           "Bin width for z_h discretization (default 0.5)"),
    ("clip",           "Absolute z_h value at which edge labels are applied (default 2.0)"),
    ("neighborhood",   "Compact key:  left_state|center_state|right_state"),
    ("left_state",     "State label of the rho30 cell immediately left of center"),
    ("center_state",   "State label of the cell whose next state is being predicted"),
    ("right_state",    "State label of the rho30 cell immediately right of center"),
    ("predicted_state","Most frequent next state observed for this neighborhood"),
    ("support",        "Number of observed transitions with this neighborhood pattern"),
    ("rule_accuracy",  "Fraction of support observations that produced predicted_state"),
    ("entropy",        "Normalized Shannon entropy of rule outputs (0=deterministic, 1=random)"),
    ("dominant_share", "Fraction of all automaton cells in the single most common state"),
    ("score",          "Quality metric: accuracy − 0.15·entropy − 0.10·dominant_share"),
    ("time_index",     "Window index: 0 is the earliest N-window, higher = later"),
    ("N_start/N_end",  "First and last even N in the automaton window"),
    ("cell_total",     "Number of even N values in this (window, rho30) cell"),
    ("state_count",    "Rows in the cell that have the dominant state label"),
    ("state_share",    "state_count / cell_total — how dominant the winning label is"),
    ("mean_z_h",       "Mean continuous z_h across all rows in this (window, rho30) cell"),
]


def print_glossary() -> None:
    """Print a human-readable symbol glossary to stdout."""
    col_w = max(len(sym) for sym, _ in SYMBOL_GLOSSARY)
    print()
    print("=" * 72)
    print("SYMBOL GLOSSARY")
    print("=" * 72)
    for symbol, definition in SYMBOL_GLOSSARY:
        print(f"  {symbol:<{col_w}}  {definition}")
    print("=" * 72)
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Encoding legend  (state label → z_h range → plain-English meaning)
# ─────────────────────────────────────────────────────────────────────────────

def build_encoding_legend(step: float, clip: float) -> pd.DataFrame:
    """
    Return a DataFrame that explicitly maps every possible state label to:
      - the z_h interval it covers  (discovered by probing z_label directly)
      - a plain-English description of what that interval means

    This is the authoritative encoding reference for end users who are not
    familiar with the z_label bucketing logic.
    """
    # Probe a fine grid so we discover the actual bin boundaries from z_label
    probe = np.linspace(-clip - 0.5, clip + 0.5, 20_000)
    label_to_values: dict[str, list[float]] = {}
    for z in probe:
        lab = z_label(z, step=step, clip=clip)
        label_to_values.setdefault(lab, []).append(float(z))

    half = step / 2.0
    rows: list[dict] = []

    for lab in sorted(label_to_values.keys(), key=residual_label_sort_key):
        values = label_to_values[lab]
        lo, hi = min(values), max(values)

        if lab.startswith("<="):
            z_range = f"z_h ≤ {-clip:.1f}"
            meaning = (
                f"Severely below expected: {clip:.0f}+ standard deviations "
                f"fewer Goldbach pairs than the heuristic predicts"
            )
        elif lab.startswith(">="):
            z_range = f"z_h ≥ {clip:.1f}"
            meaning = (
                f"Severely above expected: {clip:.0f}+ standard deviations "
                f"more Goldbach pairs than the heuristic predicts"
            )
        else:
            center = float(lab)
            z_range = f"{lo:.3f}  ≤  z_h  <  {hi:.3f}"
            if center < -half:
                meaning = (
                    f"Below expected (z_h ≈ {center:+.1f}): actual Goldbach count is "
                    f"slightly fewer than the calibrated heuristic"
                )
            elif center <= half:
                meaning = (
                    "On target (z_h ≈ 0): actual Goldbach count closely matches "
                    "the calibrated heuristic"
                )
            else:
                meaning = (
                    f"Above expected (z_h ≈ {center:+.1f}): actual Goldbach count is "
                    f"slightly more than the calibrated heuristic"
                )

        rows.append({
            "state_label":  lab,
            "z_h_range":    z_range,
            "approx_lower": round(lo, 4),
            "approx_upper": round(hi, 4),
            "meaning":      meaning,
        })

    rows.append({
        "state_label":  "missing",
        "z_h_range":    "N/A",
        "approx_lower": float("nan"),
        "approx_upper": float("nan"),
        "meaning":      "No even N values fall in this (window, rho30) cell",
    })

    return pd.DataFrame(rows)


def print_encoding_legend(legend_df: pd.DataFrame, step: float, clip: float) -> None:
    """Print the encoding legend table to stdout."""
    print("=" * 72)
    print(f"ENCODING LEGEND  (step={step:.2f}, clip={clip:.2f})")
    print("  Each state label is a discretized z_h bucket.  z_h is the")
    print("  normalized deviation of the exact Goldbach count from the")
    print("  calibrated Hardy-Littlewood heuristic.")
    print("=" * 72)
    print(f"  {'state_label':>12}  {'z_h_range':<30}  meaning")
    print(f"  {'-'*12}  {'-'*30}  {'-'*35}")
    for _, row in legend_df.iterrows():
        print(f"  {row['state_label']:>12}  {row['z_h_range']:<30}  {row['meaning']}")
    print("=" * 72)
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Core automaton dataclass
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class AutomatonSearchResult:
    best_step: float
    best_clip: float
    best_score: float
    best_accuracy: float
    best_entropy: float
    best_dominant_share: float
    states: pd.DataFrame
    rules: pd.DataFrame
    search: pd.DataFrame
    encoding_legend: pd.DataFrame = field(default_factory=pd.DataFrame)


# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

def parse_float_list(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def ordered_state_labels(labels: list[str]) -> list[str]:
    core = [label for label in labels if label != "missing"]
    core = sorted(core, key=residual_label_sort_key)
    if "missing" in labels:
        core.append("missing")
    return core


def dominant_label(series: pd.Series) -> str:
    counts = series.value_counts()
    top_count = int(counts.iloc[0])
    top_labels = sorted(
        [str(label) for label, count in counts.items() if int(count) == top_count],
        key=residual_label_sort_key,
    )
    return top_labels[0]


# ─────────────────────────────────────────────────────────────────────────────
# Automaton state builder
# ─────────────────────────────────────────────────────────────────────────────

def build_automaton_states(
    df: pd.DataFrame,
    step: float,
    clip: float,
    window_size: int,
    stride: int,
    rho_column: str = "rho30",
) -> pd.DataFrame:
    """
    Assign each (time-window, rho30-cell) pair its dominant z_h state label.

    Parameters
    ----------
    step        : bin width for z_h discretization
    clip        : absolute z_h at which edge labels are applied
    window_size : number of consecutive even N values per automaton time step
    stride      : N values to advance between successive windows
    rho_column  : column name for the residue-class cell index (default rho30)

    Returns
    -------
    DataFrame columns:
        time_index  — window index (0 = earliest)
        N_start     — first N in the window
        N_end       — last N in the window
        rho30       — N mod 30 residue class (spatial cell)
        state_label — dominant z_h bucket in this cell (see encoding legend)
        cell_total  — total even N values in this cell
        state_count — rows carrying the dominant state label
        state_share — state_count / cell_total
        mean_z_h    — mean continuous z_h across all rows in this cell
    """
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    if stride <= 0:
        raise ValueError("stride must be positive")

    work = df[["N", rho_column, "z_h"]].copy()
    work["z_state"] = [z_label(value, step=step, clip=clip) for value in work["z_h"]]

    rho_classes = sorted(work[rho_column].astype(int).unique())
    min_n = int(work["N"].min())
    max_n = int(work["N"].max())

    starts = list(range(min_n, max_n - window_size + 2, stride))
    if not starts:
        starts = [min_n]
    last_start = max(min_n, max_n - window_size + 1)
    if starts[-1] != last_start:
        starts.append(last_start)

    rows: list[dict] = []

    for time_index, start in enumerate(starts):
        end = min(start + window_size - 1, max_n)
        window = work[(work["N"] >= start) & (work["N"] <= end)].copy()

        for rho in rho_classes:
            cell = window[window[rho_column] == rho].copy()
            cell_total = int(len(cell))

            if cell_total == 0:
                rows.append({
                    "time_index":  time_index,
                    "N_start":     start,
                    "N_end":       end,
                    rho_column:    rho,
                    "state_label": "missing",
                    "cell_total":  0,
                    "state_count": 0,
                    "state_share": 0.0,
                    "mean_z_h":    float("nan"),
                })
                continue

            state_label = dominant_label(cell["z_state"])
            state_count = int((cell["z_state"] == state_label).sum())

            rows.append({
                "time_index":  time_index,
                "N_start":     start,
                "N_end":       end,
                rho_column:    rho,
                "state_label": state_label,
                "cell_total":  cell_total,
                "state_count": state_count,
                "state_share": state_count / cell_total,
                "mean_z_h":    float(cell["z_h"].mean()),
            })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Rule table
# ─────────────────────────────────────────────────────────────────────────────

def neighborhood_text(left: str, center: str, right: str) -> str:
    return f"{left}|{center}|{right}"


def entropy_from_counts(counts: dict[str, int]) -> float:
    total = float(sum(counts.values()))
    if total <= 0.0:
        return 0.0
    probs = np.array([count / total for count in counts.values()], dtype=np.float64)
    if len(probs) <= 1:
        return 0.0
    return float(-(probs * np.log(probs)).sum() / math.log(len(probs)))


def derive_rule_table(
    states_df: pd.DataFrame,
    rho_column: str = "rho30",
) -> tuple[pd.DataFrame, float, float, float]:
    """
    Build the local transition rule table.

    For each observed neighborhood pattern (left, center, right) at time t,
    predict the center cell's state at time t+1 as the most frequent outcome.

    Returns
    -------
    rule_df : DataFrame  — one row per unique neighborhood pattern, columns:
        neighborhood      compact string key "left|center|right"
        left_state        state label of the left neighbor at time t
        center_state      state label of the center cell at time t
        right_state       state label of the right neighbor at time t
        predicted_state   most frequent observed state for the center at t+1
        support           total transitions observed with this neighborhood
        rule_accuracy     fraction of support that produced predicted_state
        entropy           normalized Shannon entropy of next-state distribution
        next_state_counts raw counts "state:n;..." for full transparency
    global_accuracy : float   — support-weighted average accuracy across rules
    mean_entropy : float      — support-weighted average entropy across rules
    dominant_share : float    — fraction of cells in the most common state
    """
    rho_classes = sorted(states_df[rho_column].astype(int).unique())
    grid = (
        states_df.pivot(index="time_index", columns=rho_column, values="state_label")
        .reindex(columns=rho_classes)
        .sort_index()
    )

    transition_counts: dict[str, Counter[str]] = defaultdict(Counter)
    neighborhood_parts: dict[str, tuple[str, str, str]] = {}
    total_observations = 0
    total_correct = 0
    weighted_entropy = 0.0

    for t in range(len(grid) - 1):
        current_row = [str(v) for v in grid.iloc[t].tolist()]
        next_row = [str(v) for v in grid.iloc[t + 1].tolist()]

        for idx in range(len(rho_classes)):
            left   = current_row[(idx - 1) % len(rho_classes)]
            center = current_row[idx]
            right  = current_row[(idx + 1) % len(rho_classes)]
            key    = neighborhood_text(left, center, right)
            transition_counts[key][next_row[idx]] += 1
            neighborhood_parts[key] = (left, center, right)

    rows: list[dict] = []
    for neighborhood, counts in transition_counts.items():
        total = int(sum(counts.values()))
        predicted_state, best_count = sorted(
            counts.items(), key=lambda item: (-item[1], item[0])
        )[0]
        rule_accuracy = best_count / total
        entropy = entropy_from_counts(counts)
        total_observations += total
        total_correct += best_count
        weighted_entropy += entropy * total

        counts_text = ";".join(
            f"{state}:{count}"
            for state, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        )
        left, center, right = neighborhood_parts[neighborhood]

        rows.append({
            "neighborhood":      neighborhood,
            "left_state":        left,
            "center_state":      center,
            "right_state":       right,
            "predicted_state":   predicted_state,
            "support":           total,
            "rule_accuracy":     rule_accuracy,
            "entropy":           entropy,
            "next_state_counts": counts_text,
        })

    rule_df = pd.DataFrame(rows).sort_values(
        ["support", "rule_accuracy", "neighborhood"],
        ascending=[False, False, True],
        ignore_index=True,
    )

    global_accuracy  = total_correct / total_observations if total_observations else 0.0
    mean_entropy     = weighted_entropy / total_observations if total_observations else 0.0
    state_counts     = states_df["state_label"].value_counts()
    dominant_share   = float(state_counts.iloc[0] / state_counts.sum()) if len(state_counts) else 0.0
    return rule_df, global_accuracy, mean_entropy, dominant_share


# ─────────────────────────────────────────────────────────────────────────────
# Search over discretization parameters
# ─────────────────────────────────────────────────────────────────────────────

def candidate_score(accuracy: float, entropy: float, dominant_share: float) -> float:
    """
    Single-number quality metric for a candidate (step, clip) discretization.

    score = accuracy − 0.15·entropy − 0.10·dominant_share

    Rewards high accuracy (transitions are predictable).
    Penalizes high entropy (rule outcomes are non-deterministic) and high
    dominant_share (one bucket monopolizes all cells, making accuracy trivially
    high by always predicting the same state).
    """
    return float(accuracy - 0.15 * entropy - 0.10 * dominant_share)


def search_automaton(
    df: pd.DataFrame,
    steps: list[float],
    clips: list[float],
    window_size: int,
    stride: int,
    rho_column: str = "rho30",
) -> AutomatonSearchResult:
    """
    Grid-search over (step, clip) pairs and return the best-scoring automaton.

    The returned AutomatonSearchResult contains:
        states          — cell-level state table for the best discretization
        rules           — local transition rule table (with left/center/right columns)
        search          — full grid-search results sorted by score
        encoding_legend — human-readable mapping of every state label to its z_h range
    """
    candidate_rows: list[dict] = []
    best_states: pd.DataFrame | None = None
    best_rules: pd.DataFrame | None = None
    best_tuple: tuple[float, float, float, float, float, float] | None = None
    best_step = steps[0]
    best_clip = clips[0]
    best_score = best_accuracy = best_entropy = best_dominant_share = 0.0

    for step in steps:
        for clip in clips:
            states = build_automaton_states(
                df=df,
                step=step,
                clip=clip,
                window_size=window_size,
                stride=stride,
                rho_column=rho_column,
            )
            rules, accuracy, entropy, dominant_share = derive_rule_table(
                states, rho_column=rho_column
            )
            score       = candidate_score(accuracy, entropy, dominant_share)
            state_count = int(states["state_label"].nunique())

            candidate_rows.append({
                "step":           step,
                "clip":           clip,
                "score":          score,
                "accuracy":       accuracy,
                "entropy":        entropy,
                "dominant_share": dominant_share,
                "state_count":    state_count,
                "rule_count":     int(len(rules)),
            })

            current_tuple = (score, accuracy, -entropy, -dominant_share, -state_count, -clip)
            if best_tuple is None or current_tuple > best_tuple:
                best_tuple          = current_tuple
                best_states         = states
                best_rules          = rules
                best_step           = step
                best_clip           = clip
                best_score          = score
                best_accuracy       = accuracy
                best_entropy        = entropy
                best_dominant_share = dominant_share

    search_df = pd.DataFrame(candidate_rows).sort_values(
        ["score", "accuracy", "entropy"],
        ascending=[False, False, True],
        ignore_index=True,
    )

    assert best_states is not None
    assert best_rules is not None

    return AutomatonSearchResult(
        best_step=best_step,
        best_clip=best_clip,
        best_score=best_score,
        best_accuracy=best_accuracy,
        best_entropy=best_entropy,
        best_dominant_share=best_dominant_share,
        states=best_states,
        rules=best_rules,
        search=search_df,
        encoding_legend=build_encoding_legend(best_step, best_clip),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Dashboard plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_automaton_dashboard(
    states_df: pd.DataFrame,
    rules_df: pd.DataFrame,
    search_df: pd.DataFrame,
    out_path: str,
    rho_column: str = "rho30",
    top_rules: int = 12,
) -> None:
    ensure_parent_dir(out_path)

    mean_z_col = "mean_z_h" if "mean_z_h" in states_df.columns else "mean_z"

    rho_classes  = sorted(states_df[rho_column].astype(int).unique())
    state_labels = ordered_state_labels(states_df["state_label"].astype(str).unique().tolist())
    state_to_idx = {label: idx for idx, label in enumerate(state_labels)}

    state_grid = (
        states_df.assign(state_idx=states_df["state_label"].map(state_to_idx))
        .pivot(index="time_index", columns=rho_column, values="state_idx")
        .reindex(columns=rho_classes)
        .sort_index()
    )
    mean_z_grid = (
        states_df.pivot(index="time_index", columns=rho_column, values=mean_z_col)
        .reindex(columns=rho_classes)
        .sort_index()
    )

    cmap_base = plt.get_cmap("tab20", max(len(state_labels), 3))
    cmap = ListedColormap([cmap_base(i) for i in range(len(state_labels))])

    fig, axes = plt.subplots(2, 2, figsize=(16, 11), constrained_layout=True)
    ax_state, ax_mean, ax_occ, ax_rules = axes.ravel()

    # Panel 1 — state grid
    im_state = ax_state.imshow(
        state_grid.values, aspect="auto", interpolation="nearest", cmap=cmap
    )
    ax_state.set_title(
        "Automaton state grid\n"
        "Each cell = dominant z_h bucket in that N-window × rho30 class"
    )
    ax_state.set_xlabel(f"{rho_column}  (= N mod 30,  spatial cell index)")
    ax_state.set_ylabel("Window index  (0 = earliest N range)")
    ax_state.set_xticks(range(len(rho_classes)))
    ax_state.set_xticklabels([str(r) for r in rho_classes])
    cbar_state = fig.colorbar(im_state, ax=ax_state, fraction=0.046, pad=0.04)
    cbar_state.set_ticks(range(len(state_labels)))
    cbar_state.set_ticklabels(state_labels)
    cbar_state.set_label("State label  (z_h bucket;  negative = fewer pairs than predicted)")

    # Panel 2 — mean z_h heatmap
    im_mean = ax_mean.imshow(
        mean_z_grid.values, aspect="auto", interpolation="nearest", cmap="coolwarm"
    )
    ax_mean.set_title(
        "Mean z_h (continuous) per window × rho30 cell\n"
        "Blue = fewer Goldbach pairs than predicted,  Red = more"
    )
    ax_mean.set_xlabel(f"{rho_column}  (= N mod 30)")
    ax_mean.set_ylabel("Window index")
    ax_mean.set_xticks(range(len(rho_classes)))
    ax_mean.set_xticklabels([str(r) for r in rho_classes])
    fig.colorbar(im_mean, ax=ax_mean, fraction=0.046, pad=0.04, label="mean z_h")

    # Panel 3 — occupancy bar chart
    occupancy = (
        states_df["state_label"]
        .value_counts(normalize=True)
        .reindex(state_labels, fill_value=0.0)
        .reset_index()
    )
    occupancy.columns = ["state_label", "share"]
    ax_occ.bar(
        occupancy["state_label"],
        occupancy["share"],
        color=[cmap(state_to_idx[label]) for label in occupancy["state_label"]],
    )
    ax_occ.set_title(
        "State occupancy share\n"
        "How often each z_h bucket is the dominant label across all cells"
    )
    ax_occ.set_xlabel("State label (z_h bucket)")
    ax_occ.set_ylabel("Fraction of all automaton cells")
    ax_occ.tick_params(axis="x", rotation=35)

    # Panel 4 — top rules horizontal bar
    rules_plot = rules_df.head(top_rules).iloc[::-1]
    ax_rules.barh(rules_plot["neighborhood"], rules_plot["support"], color="#5b8e7d")
    for row in rules_plot.itertuples(index=False):
        ax_rules.text(
            row.support,
            row.neighborhood,
            f"  → {row.predicted_state}  acc={row.rule_accuracy:.2f}",
            va="center",
            ha="left",
            fontsize=8,
        )
    ax_rules.set_title(
        f"Top {top_rules} local transition rules\n"
        "Format: left_state | center_state | right_state  →  predicted next center state"
    )
    ax_rules.set_xlabel("Support  (number of observed transitions with this pattern)")
    ax_rules.set_ylabel("Neighborhood pattern  (left | center | right)")

    best = search_df.iloc[0]
    fig.suptitle(
        "Goldbach z/rho Cellular Automata  |  "
        f"Best discretization:  step={best['step']:.2f}  clip={best['clip']:.2f}\n"
        f"score={best['score']:.3f}  "
        f"accuracy={best['accuracy']:.3f}  "
        f"entropy={best['entropy']:.3f}  "
        f"dominant_share={best['dominant_share']:.3f}",
        fontsize=13,
    )
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Human-readable console summary
# ─────────────────────────────────────────────────────────────────────────────

def print_result_summary(result: AutomatonSearchResult) -> None:
    """Print a plain-English summary of the best automaton found."""
    print()
    print("=" * 72)
    print("BEST DISCRETIZATION FOUND")
    print("=" * 72)
    print(f"  step           = {result.best_step:.3f}"
          f"  (z_h bin width: values within ±{result.best_step/2:.3f} share a state label)")
    print(f"  clip           = {result.best_clip:.3f}"
          f"  (z_h values beyond ±{result.best_clip:.3f} are merged into edge labels)")
    print()
    print(f"  score          = {result.best_score:.6f}")
    print(f"                   = accuracy − 0.15·entropy − 0.10·dominant_share")
    print(f"                   = {result.best_accuracy:.6f}"
          f" − 0.15·{result.best_entropy:.6f}"
          f" − 0.10·{result.best_dominant_share:.6f}")
    print()
    print(f"  accuracy       = {result.best_accuracy:.6f}")
    print(f"                   fraction of cell transitions correctly predicted")
    print(f"  entropy        = {result.best_entropy:.6f}")
    print(f"                   rule determinism (0 = always correct, 1 = fully random)")
    print(f"  dominant_share = {result.best_dominant_share:.6f}")
    print(f"                   fraction of cells in the most common state")
    print("=" * 72)
    print()

    print_encoding_legend(result.encoding_legend, result.best_step, result.best_clip)

    top10 = result.rules.head(10)
    print("TOP 10 TRANSITION RULES")
    print("  Rule format:  left_state | center_state | right_state  →  predicted_next_center")
    print("-" * 72)
    print(f"  {'neighborhood':<24}  {'→ next':>8}  {'accuracy':>8}  {'support':>8}")
    print(f"  {'-'*24}  {'-'*8}  {'-'*8}  {'-'*8}")
    for _, row in top10.iterrows():
        print(
            f"  {row['neighborhood']:<24}  {row['predicted_state']:>8}  "
            f"{row['rule_accuracy']:>8.2f}  {row['support']:>8}"
        )
    print("-" * 72)
    print()


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--max-n", type=int, default=100_000,
        help="Analyze even integers from 4 to this value (default: 100000)",
    )
    parser.add_argument(
        "--window-size", type=int, default=1200,
        help="Even N values per automaton time step (default: 1200)",
    )
    parser.add_argument(
        "--stride", type=int, default=240,
        help="N values to advance between successive time steps (default: 240)",
    )
    parser.add_argument(
        "--search-steps", default="0.25,0.5,0.75,1.0",
        help="Comma-separated z_h bin widths to try (default: 0.25,0.5,0.75,1.0)",
    )
    parser.add_argument(
        "--search-clips", default="1.5,2.0,2.5,3.0",
        help="Comma-separated z_h saturation thresholds to try (default: 1.5,2.0,2.5,3.0)",
    )
    parser.add_argument(
        "--states-out", default="outputs/csv/goldbach_automata_states.csv",
        help="Output CSV: one row per (window, rho30) cell with state_label and mean_z_h",
    )
    parser.add_argument(
        "--rules-out", default="outputs/csv/goldbach_automata_rules.csv",
        help="Output CSV: local transition rules with left/center/right_state columns",
    )
    parser.add_argument(
        "--legend-out", default="outputs/csv/goldbach_automata_legend.csv",
        help="Output CSV: encoding legend — maps each state label to its z_h range and meaning",
    )
    parser.add_argument(
        "--summary-out", default="outputs/csv/goldbach_automata_search.csv",
        help="Output CSV: full grid-search results sorted by score",
    )
    parser.add_argument(
        "--plot-out", default="outputs/plots/goldbach_automata_dashboard.png",
    )
    parser.add_argument("--top-rules", type=int, default=12)
    parser.add_argument(
        "--glossary", action="store_true",
        help="Print the symbol glossary and exit",
    )
    parser.add_argument(
        "--no-glossary", action="store_true",
        help="Suppress the symbol glossary printed at startup",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.glossary:
        print_glossary()
        return

    if not args.no_glossary:
        print_glossary()

    steps = parse_float_list(args.search_steps)
    clips = parse_float_list(args.search_clips)

    print("[1/4] Build Goldbach residual dataset...")
    df = build_dataset(args.max_n, include_strings=True)

    print("[2/4] Search automata discretizations...")
    result = search_automaton(
        df=df,
        steps=steps,
        clips=clips,
        window_size=args.window_size,
        stride=args.stride,
        rho_column="rho30",
    )

    print("[3/4] Write CSV outputs...")
    for path in (
        args.states_out, args.rules_out, args.legend_out,
        args.summary_out, args.plot_out,
    ):
        ensure_parent_dir(path)
    result.states.to_csv(args.states_out, index=False)
    result.rules.to_csv(args.rules_out, index=False)
    result.encoding_legend.to_csv(args.legend_out, index=False)
    result.search.to_csv(args.summary_out, index=False)

    print("[4/4] Render dashboard...")
    plot_automaton_dashboard(
        states_df=result.states,
        rules_df=result.rules,
        search_df=result.search,
        out_path=args.plot_out,
        rho_column="rho30",
        top_rules=args.top_rules,
    )

    print_result_summary(result)

    print(f"states_csv  = {args.states_out}")
    print(f"rules_csv   = {args.rules_out}  ← now includes left_state / center_state / right_state columns")
    print(f"legend_csv  = {args.legend_out}  ← maps every state label to its z_h range and plain-English meaning")
    print(f"search_csv  = {args.summary_out}")
    print(f"plot_png    = {args.plot_out}")


if __name__ == "__main__":
    main()
