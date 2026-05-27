#!/usr/bin/env python3
"""
Goldbach z/rho cellular automata search.

This script builds a 1D cellular automaton on the rho lattice (`rho30` by
default). Each cell stores the dominant residual label `z_bucket` inside an
`N`-window, and the script searches over z-bucket discretizations to find the
most locally predictable automaton.

The "optimal z" in this file is operational, not theorem-level: it is the
discretization that maximizes a simple local-rule score built from

- one-step neighborhood prediction accuracy
- transition entropy
- occupancy balance across states

This is a visualization and operator-diagnostics tool for the Goldbach
framework. It is not a proof engine.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from dataclasses import dataclass

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
    top_labels = sorted([str(label) for label, count in counts.items() if int(count) == top_count], key=residual_label_sort_key)
    return top_labels[0]


def build_automaton_states(
    df: pd.DataFrame,
    step: float,
    clip: float,
    window_size: int,
    stride: int,
    rho_column: str = "rho30",
) -> pd.DataFrame:
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
                rows.append(
                    {
                        "time_index": time_index,
                        "N_start": start,
                        "N_end": end,
                        rho_column: rho,
                        "state_label": "missing",
                        "cell_total": 0,
                        "state_count": 0,
                        "state_share": 0.0,
                        "mean_z": np.nan,
                    }
                )
                continue

            state_label = dominant_label(cell["z_state"])
            state_count = int((cell["z_state"] == state_label).sum())

            rows.append(
                {
                    "time_index": time_index,
                    "N_start": start,
                    "N_end": end,
                    rho_column: rho,
                    "state_label": state_label,
                    "cell_total": cell_total,
                    "state_count": state_count,
                    "state_share": state_count / cell_total,
                    "mean_z": float(cell["z_h"].mean()),
                }
            )

    return pd.DataFrame(rows)


def neighborhood_text(left: str, center: str, right: str) -> str:
    return f"{left}|{center}|{right}"


def entropy_from_counts(counts: dict[str, int]) -> float:
    total = float(sum(counts.values()))
    if total <= 0.0:
        return 0.0
    probs = np.array([count / total for count in counts.values()], dtype=np.float64)
    if len(probs) <= 1:
        return 0.0
    entropy = float(-(probs * np.log(probs)).sum() / math.log(len(probs)))
    return entropy


def derive_rule_table(
    states_df: pd.DataFrame,
    rho_column: str = "rho30",
) -> tuple[pd.DataFrame, float, float, float]:
    rho_classes = sorted(states_df[rho_column].astype(int).unique())
    grid = (
        states_df.pivot(index="time_index", columns=rho_column, values="state_label")
        .reindex(columns=rho_classes)
        .sort_index()
    )

    transition_counts: dict[str, Counter[str]] = defaultdict(Counter)
    total_observations = 0
    total_correct = 0
    weighted_entropy = 0.0

    for t in range(len(grid) - 1):
        current_row = [str(value) for value in grid.iloc[t].tolist()]
        next_row = [str(value) for value in grid.iloc[t + 1].tolist()]

        for idx in range(len(rho_classes)):
            left = current_row[(idx - 1) % len(rho_classes)]
            center = current_row[idx]
            right = current_row[(idx + 1) % len(rho_classes)]
            next_state = next_row[idx]
            transition_counts[neighborhood_text(left, center, right)][next_state] += 1

    rows: list[dict] = []
    for neighborhood, counts in transition_counts.items():
        total = int(sum(counts.values()))
        predicted_state, best_count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
        rule_accuracy = best_count / total
        entropy = entropy_from_counts(counts)
        total_observations += total
        total_correct += best_count
        weighted_entropy += entropy * total
        counts_text = ";".join(f"{state}:{count}" for state, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])))

        rows.append(
            {
                "neighborhood": neighborhood,
                "predicted_state": predicted_state,
                "support": total,
                "rule_accuracy": rule_accuracy,
                "entropy": entropy,
                "next_state_counts": counts_text,
            }
        )

    rule_df = pd.DataFrame(rows).sort_values(
        ["support", "rule_accuracy", "neighborhood"],
        ascending=[False, False, True],
        ignore_index=True,
    )

    global_accuracy = total_correct / total_observations if total_observations else 0.0
    mean_entropy = weighted_entropy / total_observations if total_observations else 0.0
    state_counts = states_df["state_label"].value_counts()
    dominant_share = float(state_counts.iloc[0] / state_counts.sum()) if len(state_counts) else 0.0
    return rule_df, global_accuracy, mean_entropy, dominant_share


def candidate_score(accuracy: float, entropy: float, dominant_share: float) -> float:
    return float(accuracy - 0.15 * entropy - 0.10 * dominant_share)


def search_automaton(
    df: pd.DataFrame,
    steps: list[float],
    clips: list[float],
    window_size: int,
    stride: int,
    rho_column: str = "rho30",
) -> AutomatonSearchResult:
    candidate_rows: list[dict] = []
    best_states: pd.DataFrame | None = None
    best_rules: pd.DataFrame | None = None
    best_tuple: tuple[float, float, float, float, float, float] | None = None

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
            rules, accuracy, entropy, dominant_share = derive_rule_table(states, rho_column=rho_column)
            score = candidate_score(accuracy, entropy, dominant_share)
            state_count = int(states["state_label"].nunique())

            candidate_rows.append(
                {
                    "step": step,
                    "clip": clip,
                    "score": score,
                    "accuracy": accuracy,
                    "entropy": entropy,
                    "dominant_share": dominant_share,
                    "state_count": state_count,
                    "rule_count": int(len(rules)),
                }
            )

            current_tuple = (score, accuracy, -entropy, -dominant_share, -state_count, -clip)
            if best_tuple is None or current_tuple > best_tuple:
                best_tuple = current_tuple
                best_states = states
                best_rules = rules
                best_step = step
                best_clip = clip
                best_score = score
                best_accuracy = accuracy
                best_entropy = entropy
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
    )


def plot_automaton_dashboard(
    states_df: pd.DataFrame,
    rules_df: pd.DataFrame,
    search_df: pd.DataFrame,
    out_path: str,
    rho_column: str = "rho30",
    top_rules: int = 12,
) -> None:
    ensure_parent_dir(out_path)

    rho_classes = sorted(states_df[rho_column].astype(int).unique())
    state_labels = ordered_state_labels(states_df["state_label"].astype(str).unique().tolist())
    state_to_idx = {label: idx for idx, label in enumerate(state_labels)}

    state_grid = (
        states_df.assign(state_idx=states_df["state_label"].map(state_to_idx))
        .pivot(index="time_index", columns=rho_column, values="state_idx")
        .reindex(columns=rho_classes)
        .sort_index()
    )
    mean_z_grid = (
        states_df.pivot(index="time_index", columns=rho_column, values="mean_z")
        .reindex(columns=rho_classes)
        .sort_index()
    )

    cmap_base = plt.get_cmap("tab20", max(len(state_labels), 3))
    cmap = ListedColormap([cmap_base(i) for i in range(len(state_labels))])

    fig, axes = plt.subplots(2, 2, figsize=(16, 11), constrained_layout=True)
    ax_state, ax_mean, ax_occ, ax_rules = axes.ravel()

    im_state = ax_state.imshow(state_grid.values, aspect="auto", interpolation="nearest", cmap=cmap)
    ax_state.set_title("Automaton state zoo")
    ax_state.set_xlabel(rho_column)
    ax_state.set_ylabel("window index")
    ax_state.set_xticks(range(len(rho_classes)))
    ax_state.set_xticklabels([str(rho) for rho in rho_classes], rotation=0)
    cbar_state = fig.colorbar(im_state, ax=ax_state, fraction=0.046, pad=0.04)
    cbar_state.set_ticks(range(len(state_labels)))
    cbar_state.set_ticklabels(state_labels)

    im_mean = ax_mean.imshow(mean_z_grid.values, aspect="auto", interpolation="nearest", cmap="coolwarm")
    ax_mean.set_title("Mean z_h per window / rho cell")
    ax_mean.set_xlabel(rho_column)
    ax_mean.set_ylabel("window index")
    ax_mean.set_xticks(range(len(rho_classes)))
    ax_mean.set_xticklabels([str(rho) for rho in rho_classes], rotation=0)
    fig.colorbar(im_mean, ax=ax_mean, fraction=0.046, pad=0.04)

    occupancy = (
        states_df["state_label"]
        .value_counts(normalize=True)
        .reindex(state_labels, fill_value=0.0)
        .reset_index()
    )
    occupancy.columns = ["state_label", "share"]
    ax_occ.bar(occupancy["state_label"], occupancy["share"], color=[cmap(state_to_idx[label]) for label in occupancy["state_label"]])
    ax_occ.set_title("State occupancy share")
    ax_occ.set_ylabel("share of automaton cells")
    ax_occ.tick_params(axis="x", rotation=35)

    rules_plot = rules_df.head(top_rules).iloc[::-1]
    ax_rules.barh(rules_plot["neighborhood"], rules_plot["support"], color="#5b8e7d")
    for row in rules_plot.itertuples(index=False):
        ax_rules.text(
            row.support,
            row.neighborhood,
            f"  acc={row.rule_accuracy:.2f}",
            va="center",
            ha="left",
            fontsize=8,
        )
    ax_rules.set_title("Top local rules")
    ax_rules.set_xlabel("support")

    best = search_df.iloc[0]
    fig.suptitle(
        (
            "Goldbach z/rho cellular automata"
            f"  step={best['step']:.2f}"
            f"  clip={best['clip']:.2f}"
            f"  score={best['score']:.3f}"
            f"  acc={best['accuracy']:.3f}"
        ),
        fontsize=16,
    )
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=100_000)
    parser.add_argument("--window-size", type=int, default=1200)
    parser.add_argument("--stride", type=int, default=240)
    parser.add_argument("--search-steps", default="0.25,0.5,0.75,1.0")
    parser.add_argument("--search-clips", default="1.5,2.0,2.5,3.0")
    parser.add_argument("--states-out", default="outputs/csv/goldbach_automata_states.csv")
    parser.add_argument("--rules-out", default="outputs/csv/goldbach_automata_rules.csv")
    parser.add_argument("--summary-out", default="outputs/csv/goldbach_automata_search.csv")
    parser.add_argument("--plot-out", default="outputs/plots/goldbach_automata_dashboard.png")
    parser.add_argument("--top-rules", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
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
    for path in (args.states_out, args.rules_out, args.summary_out, args.plot_out):
        ensure_parent_dir(path)
    result.states.to_csv(args.states_out, index=False)
    result.rules.to_csv(args.rules_out, index=False)
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

    print("best_step={:.3f}".format(result.best_step))
    print("best_clip={:.3f}".format(result.best_clip))
    print("best_score={:.6f}".format(result.best_score))
    print("best_accuracy={:.6f}".format(result.best_accuracy))
    print("best_entropy={:.6f}".format(result.best_entropy))
    print("best_dominant_share={:.6f}".format(result.best_dominant_share))
    print(f"states_csv={args.states_out}")
    print(f"rules_csv={args.rules_out}")
    print(f"search_csv={args.summary_out}")
    print(f"plot_png={args.plot_out}")


if __name__ == "__main__":
    main()
