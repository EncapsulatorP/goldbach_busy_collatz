#!/usr/bin/env python3
"""
Goldbach super-compressed chamber analysis with cluster-filtered plots.

This script keeps the `q mod 30 + i p mod 30` pair-fiber compression from
`shattering_compressed.py`, but uses the native h-residual cluster filter from
`shattering_mirrors.py` to decide which families to visualize.
"""

import argparse
import math
from collections import Counter
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import patches

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shattering_mirrors import (
    build_dataset,
    make_cluster_palette,
    ordered_residual_buckets,
    resolve_cluster_filter,
    sieve_bool,
)


PRIME_RESIDUES_30 = [1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29]


def gi_label(real: int, imag: int) -> str:
    real = int(real)
    imag = int(imag)

    if imag == 0:
        return str(real)

    if real == 0:
        if imag == 1:
            return "i"
        if imag == -1:
            return "-i"
        return f"{imag}i"

    if imag == 1:
        return f"{real}+i"
    if imag == -1:
        return f"{real}-i"

    if imag > 0:
        return f"{real}+{imag}i"

    return f"{real}{imag}i"


def goldbach_pairs(N: int, is_prime) -> list[tuple[int, int]]:
    """
    Unordered Goldbach pairs p <= q with p + q = N.
    """
    pairs = []
    for p in range(2, N // 2 + 1):
        q = N - p
        if is_prime[p] and is_prime[q]:
            pairs.append((p, q))
    return pairs


def pair_chambers_i(N: int, is_prime) -> list[dict]:
    """
    Build a_q|m|b_p and i-compressed q30 + i*p30 forms.

    For each Goldbach pair p + q = N:

        a_q = q mod 30
        m   = N mod 30
        b_p = p mod 30

    chamber:
        a_q | m | b_p

    i form:
        a_q + i b_p [m]
    """
    m = N % 30
    out = []

    for p, q in goldbach_pairs(N, is_prime):
        a_q = q % 30
        b_p = p % 30

        chamber = f"{a_q}|{m}|{b_p}"
        z_pair = gi_label(a_q, b_p)

        check = (a_q + b_p) % 30

        out.append(
            {
                "N": N,
                "p": p,
                "q": q,
                "m": m,
                "a_q": a_q,
                "b_p": b_p,
                "pair_chamber": chamber,
                "pair_i": f"{z_pair}[{m}]",
                "chamber_check": check == m,
            }
        )

    return out


def super_compressed_N(N: int, r: int, h_floor: int, is_prime) -> dict:
    """
    Legacy helper for raw number-level + pair-fiber compression.
    """
    eps_h = r - h_floor
    rho30 = N % 30

    z_N = gi_label(eps_h, rho30)
    fiber = pair_chambers_i(N, is_prime)

    fiber_i = [x["pair_i"] for x in fiber]

    return {
        "N": N,
        "r": r,
        "h_floor": h_floor,
        "eps_h": eps_h,
        "rho30": rho30,
        "z_N": z_N,
        "pair_fiber_i": fiber_i,
        "pair_fiber_size": len(fiber_i),
    }


def build_selected_number_frame(
    max_n: int,
    cluster_filter: str,
    top_clusters: int,
) -> tuple[pd.DataFrame, list[str], pd.DataFrame, pd.DataFrame]:
    df = build_dataset(max_n, include_strings=True)
    selected_labels, filtered_counts = resolve_cluster_filter(df, cluster_filter, top_clusters)
    selected_numbers = df[df["native_cluster"].isin(selected_labels)].copy()
    selected_numbers["z_N"] = [
        f"{z_bucket}+{rho30}i"
        for z_bucket, rho30 in zip(selected_numbers["z_bucket"], selected_numbers["rho30"])
    ]
    selected_numbers["pair_fiber_size"] = selected_numbers["r"]
    return df, selected_labels, filtered_counts, selected_numbers


def ensure_parent_dir(path: str) -> None:
    parent = Path(path).expanduser().resolve().parent
    parent.mkdir(parents=True, exist_ok=True)


def build_pair_aggregates(
    selected_numbers: pd.DataFrame,
    selected_labels: list[str],
    is_prime,
) -> dict:
    residue_to_idx = {value: idx for idx, value in enumerate(PRIME_RESIDUES_30)}
    global_heat = np.zeros((len(PRIME_RESIDUES_30), len(PRIME_RESIDUES_30)), dtype=np.int64)
    per_cluster_heat = {
        label: np.zeros((len(PRIME_RESIDUES_30), len(PRIME_RESIDUES_30)), dtype=np.int64)
        for label in selected_labels
    }
    global_signatures = Counter()
    per_cluster_signatures = {label: Counter() for label in selected_labels}
    invalid_checks = 0

    for row in selected_numbers.itertuples(index=False):
        label = row.native_cluster
        pairs = pair_chambers_i(int(row.N), is_prime)

        for pair in pairs:
            if not pair["chamber_check"]:
                invalid_checks += 1

            a_q = int(pair["a_q"])
            b_p = int(pair["b_p"])
            x_idx = residue_to_idx.get(a_q)
            y_idx = residue_to_idx.get(b_p)

            if x_idx is not None and y_idx is not None:
                global_heat[y_idx, x_idx] += 1
                per_cluster_heat[label][y_idx, x_idx] += 1

            signature = pair["pair_i"]
            global_signatures[signature] += 1
            per_cluster_signatures[label][signature] += 1

    return {
        "global_heat": global_heat,
        "per_cluster_heat": per_cluster_heat,
        "global_signatures": global_signatures,
        "per_cluster_signatures": per_cluster_signatures,
        "invalid_checks": invalid_checks,
    }


def top_signature_table(
    per_cluster_signatures: dict[str, Counter],
    selected_labels: list[str],
    top_k: int,
) -> pd.DataFrame:
    total_counter = Counter()
    for counter in per_cluster_signatures.values():
        total_counter.update(counter)

    top_signatures = [signature for signature, _ in total_counter.most_common(top_k)]
    rows = []

    for signature in top_signatures:
        row = {"pair_i": signature, "total": total_counter[signature]}
        for label in selected_labels:
            row[label] = per_cluster_signatures[label][signature]
        rows.append(row)

    return pd.DataFrame(rows)


def format_signature_block(counter: Counter, lines: int = 5) -> str:
    top = counter.most_common(lines)
    if not top:
        return "top pair_i: none"
    return "\n".join(f"{signature}: {count}" for signature, count in top)


def plot_compressed_dashboard(
    df: pd.DataFrame,
    selected_labels: list[str],
    filtered_counts: pd.DataFrame,
    pair_agg: dict,
    out_path: str,
    top_k: int,
) -> pd.DataFrame:
    palette = make_cluster_palette(selected_labels)
    z_labels = ordered_residual_buckets(df)
    rho_labels = sorted(df["rho30"].astype(int).unique())

    heat = (
        df.groupby(["z_bucket", "rho30"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=z_labels, columns=rho_labels, fill_value=0)
    )
    signature_df = top_signature_table(pair_agg["per_cluster_signatures"], selected_labels, top_k=top_k)

    fig, axes = plt.subplots(3, 2, figsize=(16, 17))
    ax_heat, ax_bar, ax_eps, ax_r, ax_pair, ax_sig = axes.ravel()

    im = ax_heat.imshow(heat.values, aspect="auto", cmap="magma")
    ax_heat.set_title("Number-level normalized residual chambers")
    ax_heat.set_xlabel("rho30")
    ax_heat.set_ylabel("z_bucket")
    ax_heat.set_xticks(range(len(rho_labels)))
    ax_heat.set_xticklabels(rho_labels)
    ax_heat.set_yticks(range(len(z_labels)))
    ax_heat.set_yticklabels(z_labels)
    for i, label in enumerate(selected_labels):
        z_bucket, rho_text = label.split(";rho30=")
        jj = rho_labels.index(int(rho_text))
        ii = z_labels.index(z_bucket.replace("z=", ""))
        rect = patches.Rectangle((jj - 0.5, ii - 0.5), 1, 1, fill=False, lw=2.5, ec=palette[label])
        ax_heat.add_patch(rect)
        ax_heat.text(jj, ii, str(i + 1), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04, label="count")

    ax_bar.set_title("Selected number-cluster sizes")
    bar_positions = np.arange(len(filtered_counts))
    ax_bar.bar(
        bar_positions,
        filtered_counts["count"],
        color=[palette[label] for label in filtered_counts["native_cluster"]],
        alpha=0.92,
    )
    ax_bar.set_xticks(bar_positions)
    ax_bar.set_xticklabels([str(i + 1) for i in range(len(filtered_counts))])
    ax_bar.set_xlabel("selected cluster id")
    ax_bar.set_ylabel("count of N values")
    for pos, (_, row) in enumerate(filtered_counts.iterrows()):
        ax_bar.text(pos, row["count"], f"{row['rho30']}\n{row['z_bucket']}", ha="center", va="bottom", fontsize=8)

    ax_eps.set_title("N vs calibrated z_h with selected clusters highlighted")
    ax_eps.scatter(df["N"], df["z_h"], s=7, color="#b0b0b0", alpha=0.18, edgecolors="none")
    for label in selected_labels:
        subset = df[df["native_cluster"] == label]
        ax_eps.scatter(subset["N"], subset["z_h"], s=11, color=palette[label], alpha=0.82, edgecolors="none")
    ax_eps.axhline(0.0, color="#444444", lw=1.0)
    ax_eps.set_xlabel("N")
    ax_eps.set_ylabel("z_h = (r - alpha*h)/sqrt(alpha*h)")

    ax_r.set_title("N vs pair-fiber size r_G(N)")
    ax_r.scatter(df["N"], df["r"], s=7, color="#b0b0b0", alpha=0.18, edgecolors="none")
    for label in selected_labels:
        subset = df[df["native_cluster"] == label]
        ax_r.scatter(subset["N"], subset["r"], s=11, color=palette[label], alpha=0.82, edgecolors="none")
    ax_r.set_xlabel("N")
    ax_r.set_ylabel("pair_fiber_size = r_G(N)")

    pair_im = ax_pair.imshow(pair_agg["global_heat"], aspect="auto", cmap="viridis")
    ax_pair.set_title("Selected-pair chamber heatmap")
    ax_pair.set_xlabel("a_q = q mod 30")
    ax_pair.set_ylabel("b_p = p mod 30")
    ax_pair.set_xticks(range(len(PRIME_RESIDUES_30)))
    ax_pair.set_xticklabels(PRIME_RESIDUES_30)
    ax_pair.set_yticks(range(len(PRIME_RESIDUES_30)))
    ax_pair.set_yticklabels(PRIME_RESIDUES_30)
    pair_title = f"all selected pair fibers; chamber failures={pair_agg['invalid_checks']}"
    ax_pair.text(
        0.02,
        0.98,
        pair_title,
        transform=ax_pair.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.35, "edgecolor": "none"},
    )
    fig.colorbar(pair_im, ax=ax_pair, fraction=0.046, pad=0.04, label="pair count")

    ax_sig.set_title("Top pair_i signatures across selected clusters")
    if not signature_df.empty:
        x = np.arange(len(signature_df))
        bottom = np.zeros(len(signature_df))
        for label in selected_labels:
            values = signature_df[label].to_numpy()
            ax_sig.bar(x, values, bottom=bottom, color=palette[label], alpha=0.92)
            bottom += values
        ax_sig.set_xticks(x)
        ax_sig.set_xticklabels(signature_df["pair_i"], rotation=35, ha="right")
        ax_sig.set_ylabel("pair count")
    else:
        ax_sig.text(0.5, 0.5, "no pair_i signatures found", ha="center", va="center")
        ax_sig.set_xticks([])
        ax_sig.set_yticks([])

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=palette[label], label=f"{i + 1}: {label}")
        for i, label in enumerate(selected_labels)
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle("Shattering Compressed Cluster Dashboard", fontsize=17, y=0.995)
    fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return signature_df


def plot_cluster_gallery(
    selected_labels: list[str],
    filtered_counts: pd.DataFrame,
    pair_agg: dict,
    out_path: str,
) -> None:
    n_panels = len(selected_labels)
    ncols = 2
    nrows = int(math.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4.8 * nrows), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    vmax = max(
        int(pair_agg["per_cluster_heat"][label].max())
        for label in selected_labels
    )
    vmax = max(vmax, 1)

    for idx, label in enumerate(selected_labels):
        ax = axes[idx]
        meta = filtered_counts[filtered_counts["native_cluster"] == label].iloc[0]
        heat = pair_agg["per_cluster_heat"][label]
        counter = pair_agg["per_cluster_signatures"][label]

        im = ax.imshow(heat, aspect="auto", cmap="viridis", vmin=0, vmax=vmax)
        ax.set_title(f"{idx + 1}. {label}")
        ax.set_xlabel("a_q = q mod 30")
        ax.set_ylabel("b_p = p mod 30")
        ax.set_xticks(range(len(PRIME_RESIDUES_30)))
        ax.set_xticklabels(PRIME_RESIDUES_30)
        ax.set_yticks(range(len(PRIME_RESIDUES_30)))
        ax.set_yticklabels(PRIME_RESIDUES_30)
        ax.text(
            0.02,
            0.98,
            (
                f"N count={int(meta['count'])}\n"
                f"N range=[{int(meta['min_N'])}, {int(meta['max_N'])}]\n"
                f"mean_r={meta['mean_r']:.2f}, mean_h={meta['mean_h']:.2f}\n"
                f"mean_z={meta['mean_z']:.3f}, mean_boost={meta['mean_boost']:.3f}\n"
                f"total pairs={int(heat.sum())}\n"
                f"{format_signature_block(counter, lines=4)}"
            ),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8.8,
            color="white",
            bbox={"facecolor": "black", "alpha": 0.42, "edgecolor": "none"},
        )

    for ax in axes[n_panels:]:
        ax.axis("off")

    fig.colorbar(im, ax=axes[:n_panels].tolist(), fraction=0.025, pad=0.02, label="pair count")
    fig.suptitle("Shattering Compressed Cluster Gallery", fontsize=17, y=0.995)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=100_000, help="Largest even N to include.")
    parser.add_argument("--plot", action="store_true", help="Generate the compressed dashboard and gallery PNGs.")
    parser.add_argument(
        "--plot-prefix",
        type=str,
        default="outputs/plots/shattering_compressed",
        help="Prefix for PNG outputs.",
    )
    parser.add_argument(
        "--numbers-out",
        type=str,
        default="outputs/csv/shattering_compressed_numbers.csv",
        help="CSV path for the selected number rows.",
    )
    parser.add_argument(
        "--summary-out",
        type=str,
        default="outputs/csv/shattering_compressed_cluster_summary.csv",
        help="CSV path for the selected cluster summary.",
    )
    parser.add_argument(
        "--signatures-out",
        type=str,
        default="outputs/csv/shattering_compressed_pair_signatures.csv",
        help="CSV path for the pair_i signature table.",
    )
    parser.add_argument(
        "--plot-top-clusters",
        type=int,
        default=8,
        help="If no explicit cluster filter is given, plot the top K native clusters.",
    )
    parser.add_argument(
        "--cluster-filter",
        type=str,
        default="",
        help="Comma-separated native_cluster labels such as 'z=+0.5;rho30=6,z=-0.5;rho30=0'.",
    )
    parser.add_argument(
        "--top-signatures",
        type=int,
        default=10,
        help="How many pair_i signatures to show in the dashboard bar chart.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("[1/7] Build number-level h-residual dataset...")
    df, selected_labels, filtered_counts, selected_numbers = build_selected_number_frame(
        max_n=args.max_n,
        cluster_filter=args.cluster_filter,
        top_clusters=args.plot_top_clusters,
    )

    counts_path = args.summary_out
    numbers_path = args.numbers_out
    signature_path = args.signatures_out
    ensure_parent_dir(counts_path)
    ensure_parent_dir(numbers_path)
    ensure_parent_dir(signature_path)
    filtered_counts.to_csv(counts_path, index=False)
    selected_numbers.to_csv(numbers_path, index=False)
    print(f"[2/7] Saved selected number rows to {numbers_path}")
    print(f"[3/7] Saved cluster summary to {counts_path}")
    print("[4/7] Selected native clusters:")
    for idx, label in enumerate(selected_labels, start=1):
        print(f"  {idx}. {label}")

    if not args.plot:
        print("[5/7] Plot generation skipped. Use --plot to render the dashboard and gallery.")
        return

    print(f"[5/7] Sieve primes up to {args.max_n} for pair-fiber aggregation...")
    is_prime = sieve_bool(args.max_n)

    print(f"[6/7] Aggregate pair fibers for {len(selected_numbers)} selected N values...")
    pair_agg = build_pair_aggregates(selected_numbers, selected_labels, is_prime)

    dashboard_path = f"{args.plot_prefix}_dashboard.png"
    gallery_path = f"{args.plot_prefix}_gallery.png"
    ensure_parent_dir(dashboard_path)
    ensure_parent_dir(gallery_path)

    print("[7/7] Render compressed plots...")
    signature_df = plot_compressed_dashboard(
        df=df,
        selected_labels=selected_labels,
        filtered_counts=filtered_counts,
        pair_agg=pair_agg,
        out_path=dashboard_path,
        top_k=args.top_signatures,
    )
    signature_df.to_csv(signature_path, index=False)
    plot_cluster_gallery(
        selected_labels=selected_labels,
        filtered_counts=filtered_counts,
        pair_agg=pair_agg,
        out_path=gallery_path,
    )

    print(f"Saved dashboard: {dashboard_path}")
    print(f"Saved gallery:   {gallery_path}")
    print(f"Saved pair_i table: {signature_path}")


if __name__ == "__main__":
    main()
