#!/usr/bin/env python3
"""
Goldbach volume splash generator.

The original file in this workspace was truncated to a one-liner. This rebuild
reconstructs the checked-in artifact shape from the surrounding repo outputs:

- exact Goldbach counts `r_G(N)` via FFT
- raw Hardy-Littlewood-style heuristic `h(N)`
- raw normalized residual `z_h = (r - h) / sqrt(h)`
- residue/family labels tied to the odd prime support of `N / 2`
- pair-signature counts of the form `q mod m + i p mod m [N mod m]`

The resulting CSV layout matches the current `goldbach_volume_*` artifacts, and
the script can also generate a static PNG dashboard plus an interactive Plotly
HTML splash when requested.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shattering_mirrors import (
    exact_goldbach_counts_fft,
    h_goldbach,
    local_goldbach_boost,
    reciprocal_log_density_convolution,
    spf_sieve,
    sieve_bool,
)


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


def eps_label(eps: int, clip: int = 5) -> str:
    if eps < -clip:
        return f"<-{clip}"
    if eps > clip:
        return f">{clip}"
    return str(eps)


def ensure_parent_dir(path: str) -> None:
    parent = Path(path).expanduser().resolve().parent
    parent.mkdir(parents=True, exist_ok=True)


def human_k(max_n: int) -> str:
    if max_n % 1000 == 0:
        return f"{max_n // 1000}k"
    return str(max_n)


def default_stem(max_n: int, modulus: int) -> str:
    base = f"goldbach_volume_{human_k(max_n)}"
    if modulus == 30:
        return f"{base}_v2"
    return f"{base}_mod{modulus}"


def odd_prime_support(n: int, spf: np.ndarray) -> list[int]:
    support: list[int] = []
    last_p = None
    x = n

    while x > 1:
        p = int(spf[x])
        if p == 0:
            break
        if p != last_p and p > 2:
            support.append(p)
        last_p = p
        while x % p == 0:
            x //= p

    return support


def family_label(support: list[int]) -> str:
    if not support:
        return "spine_none"
    if support == [3]:
        return "cone_3"
    if support == [5]:
        return "ring_5"
    if support == [3, 5]:
        return "mobius_3x5"
    return "shell_other"


def goldbach_pairs(N: int, prime_list: np.ndarray, is_prime: np.ndarray) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    half = N // 2
    for p in prime_list:
        p = int(p)
        if p > half:
            break
        q = N - p
        if is_prime[q]:
            pairs.append((p, q))
    return pairs


def volume_coordinates(df: pd.DataFrame, modulus: int, residue_column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    residue = df[residue_column].to_numpy(dtype=float)
    angle = (2.0 * math.pi * residue) / float(modulus)
    radius = np.log1p(df["h"].to_numpy(dtype=float)) * (
        1.0 + 0.25 * np.log1p(df["boost_G"].to_numpy(dtype=float))
    )
    shear = df["eps_h"].to_numpy(dtype=float) / np.sqrt(np.maximum(df["N"].to_numpy(dtype=float), 4.0))
    z_h = df["z_h"].to_numpy(dtype=float)

    x = -radius * np.cos(angle) - 4.0 * shear
    y = radius * np.sin(angle) + 0.60 * z_h
    z = z_h * np.sqrt(np.log1p(df["N"].to_numpy(dtype=float))) + 0.20 * df["eps_h"].to_numpy(dtype=float)
    return x, y, z


def build_volume_dataset(max_n: int, modulus: int) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    if max_n < 4:
        raise ValueError("max_n must be at least 4")
    if max_n % 2:
        raise ValueError("max_n must be even")
    if modulus <= 1:
        raise ValueError("modulus must be greater than 1")

    print(f"[1/6] Sieve primes up to {max_n}...")
    is_prime = sieve_bool(max_n)
    prime_list = np.flatnonzero(is_prime)

    print(f"[2/6] Build SPF up to {max_n}...")
    spf = spf_sieve(max_n)

    print("[3/6] Build reciprocal-log density convolution...")
    density_conv = reciprocal_log_density_convolution(max_n)

    print(f"[4/6] Compute exact Goldbach counts by FFT...")
    r_counts = exact_goldbach_counts_fft(max_n, is_prime)

    print(f"[5/6] Build raw volume rows modulo {modulus}...")
    rows: list[dict] = []
    residue_column = "rho30" if modulus == 30 else "rho_mod"
    residue_label = f"rho{modulus}"

    for N in range(4, max_n + 1, 2):
        n = N // 2
        r = int(r_counts[N])
        boost = local_goldbach_boost(N, spf)
        h = float(h_goldbach(N, boost, density_conv=density_conv))
        h_floor = int(math.floor(h))
        eps_h = int(r - h_floor)
        z_h = float((r - h) / math.sqrt(h)) if h > 0 else 0.0
        residue = int(N % modulus)
        support = odd_prime_support(n, spf)
        boost_primes = ",".join(str(p) for p in support) or None
        family = family_label(support)
        eps_bucket = eps_label(eps_h)

        row = {
            "N": N,
            "n": n,
            "r": r,
            "h": h,
            "h_floor": h_floor,
            "eps_h": eps_h,
            "eps_bucket": eps_bucket,
            "z_h": z_h,
            "boost_G": float(boost),
            "boost_primes": boost_primes,
            "family": family,
            "native_cluster": f"eps={eps_bucket};{residue_label}={residue}",
            "z_native_i": gi_label(eps_h, residue),
            residue_column: residue,
        }
        if modulus != 30:
            row["modulus"] = modulus

        rows.append(row)

    df = pd.DataFrame(rows)
    x, y, z = volume_coordinates(df, modulus, residue_column)
    df["x"] = x
    df["y"] = y
    df["z"] = z

    if modulus == 30:
        ordered_columns = [
            "N",
            "n",
            "r",
            "h",
            "h_floor",
            "eps_h",
            "eps_bucket",
            "z_h",
            "rho30",
            "boost_G",
            "boost_primes",
            "family",
            "native_cluster",
            "z_native_i",
            "x",
            "y",
            "z",
        ]
    else:
        ordered_columns = [
            "N",
            "n",
            "r",
            "h",
            "h_floor",
            "eps_h",
            "eps_bucket",
            "z_h",
            "modulus",
            "rho_mod",
            "boost_G",
            "boost_primes",
            "family",
            "native_cluster",
            "z_native_i",
            "x",
            "y",
            "z",
        ]

    print("[6/6] Build summaries...")
    df = df[ordered_columns]
    return df, is_prime, prime_list


def summarize_families(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby("family", dropna=False)
        .agg(
            count=("N", "size"),
            min_N=("N", "min"),
            max_N=("N", "max"),
            mean_r=("r", "mean"),
            mean_h=("h", "mean"),
            mean_z=("z_h", "mean"),
            mean_boost=("boost_G", "mean"),
        )
        .reset_index()
        .sort_values(["count", "family"], ascending=[False, True], ignore_index=True)
    )
    return summary


def summarize_clusters(df: pd.DataFrame, residue_column: str) -> pd.DataFrame:
    summary = (
        df.groupby(["eps_bucket", residue_column, "native_cluster"], dropna=False)
        .agg(
            count=("N", "size"),
            min_N=("N", "min"),
            max_N=("N", "max"),
            mean_r=("r", "mean"),
            mean_h=("h", "mean"),
            mean_z=("z_h", "mean"),
            mean_boost=("boost_G", "mean"),
        )
        .reset_index()
        .sort_values(["count", residue_column], ascending=[False, True], ignore_index=True)
    )
    return summary


def pair_signature_counts(
    max_n: int,
    modulus: int,
    prime_list: np.ndarray,
    is_prime: np.ndarray,
) -> pd.DataFrame:
    counter: Counter[str] = Counter()

    for N in range(4, max_n + 1, 2):
        residue = N % modulus
        for p, q in goldbach_pairs(N, prime_list, is_prime):
            signature = f"{gi_label(q % modulus, p % modulus)}[{residue}]"
            counter[signature] += 1

    rows = [{"pair_i": pair_i, "count": count} for pair_i, count in counter.most_common()]
    return pd.DataFrame(rows)


def plot_dashboard(
    df: pd.DataFrame,
    family_summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    pair_df: pd.DataFrame,
    residue_column: str,
    plot_out: str,
    top_k: int,
) -> None:
    ensure_parent_dir(plot_out)

    family_colors = {
        "spine_none": "#264653",
        "cone_3": "#d1495b",
        "ring_5": "#edae49",
        "mobius_3x5": "#00798c",
        "shell_other": "#5c4d7d",
    }

    fig, axes = plt.subplots(2, 2, figsize=(15, 11), constrained_layout=True)
    ax_xy, ax_heat, ax_family, ax_pair = axes.ravel()

    for family, group in df.groupby("family", sort=False):
        ax_xy.scatter(
            group["x"],
            group["y"],
            s=16,
            alpha=0.65,
            label=family,
            color=family_colors.get(family, "#444444"),
        )
    ax_xy.set_title("Goldbach volume splash")
    ax_xy.set_xlabel("volume x")
    ax_xy.set_ylabel("volume y")
    ax_xy.grid(alpha=0.20, linewidth=0.6)
    ax_xy.legend(loc="best", fontsize=8)

    heat = (
        cluster_summary.pivot(index="eps_bucket", columns=residue_column, values="count")
        .fillna(0)
        .sort_index(axis=1)
    )
    im = ax_heat.imshow(heat.values, aspect="auto", cmap="magma")
    ax_heat.set_title("Native cluster occupancy")
    ax_heat.set_xlabel(residue_column)
    ax_heat.set_ylabel("eps_bucket")
    ax_heat.set_xticks(range(len(heat.columns)))
    ax_heat.set_xticklabels([str(int(x)) for x in heat.columns], rotation=90)
    ax_heat.set_yticks(range(len(heat.index)))
    ax_heat.set_yticklabels(heat.index.tolist())
    fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04)

    fam_plot = family_summary.head(5).copy()
    ax_family.bar(fam_plot["family"], fam_plot["count"], color="#4c78a8")
    ax_family.set_title("Family counts")
    ax_family.set_ylabel("count")
    ax_family.tick_params(axis="x", rotation=20)

    pair_plot = pair_df.head(top_k).iloc[::-1]
    ax_pair.barh(pair_plot["pair_i"], pair_plot["count"], color="#59a14f")
    ax_pair.set_title(f"Top {min(top_k, len(pair_df))} pair signatures")
    ax_pair.set_xlabel("count")

    fig.savefig(plot_out, dpi=180)
    plt.close(fig)


def write_html_splash(
    df: pd.DataFrame,
    family_summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    pair_df: pd.DataFrame,
    residue_column: str,
    modulus: int,
    html_out: str,
    top_k: int,
) -> None:
    try:
        from plotly import graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as exc:
        raise RuntimeError(
            "Plotly is required for --html output but is not installed in this environment."
        ) from exc

    ensure_parent_dir(html_out)

    fig = make_subplots(
        rows=2,
        cols=2,
        specs=[[{"type": "scene"}, {"type": "xy"}], [{"type": "heatmap"}, {"type": "xy"}]],
        column_widths=[0.63, 0.37],
        row_heights=[0.68, 0.32],
        horizontal_spacing=0.06,
        vertical_spacing=0.10,
        subplot_titles=(
            "Volume splash",
            "Family counts",
            "Native cluster occupancy",
            "Top pair signatures",
        ),
    )

    family_colors = {
        "spine_none": "#264653",
        "cone_3": "#d1495b",
        "ring_5": "#edae49",
        "mobius_3x5": "#00798c",
        "shell_other": "#5c4d7d",
    }

    for family, group in df.groupby("family", sort=False):
        hover = np.stack(
            [
                group["N"].to_numpy(),
                group["r"].to_numpy(),
                group["h"].round(3).to_numpy(),
                group["eps_h"].to_numpy(),
                group[residue_column].to_numpy(),
            ],
            axis=-1,
        )
        fig.add_trace(
            go.Scatter3d(
                x=group["x"],
                y=group["y"],
                z=group["z"],
                mode="markers",
                name=family,
                legendgroup=family,
                marker={
                    "size": 3.4,
                    "opacity": 0.68,
                    "color": family_colors.get(family, "#444444"),
                },
                customdata=hover,
                hovertemplate=(
                    "family=%{fullData.name}<br>"
                    "N=%{customdata[0]}<br>"
                    "r=%{customdata[1]}<br>"
                    "h=%{customdata[2]}<br>"
                    "eps_h=%{customdata[3]}<br>"
                    f"{residue_column}=%{{customdata[4]}}"
                    "<extra></extra>"
                ),
            ),
            row=1,
            col=1,
        )

    fam_plot = family_summary.head(5)
    fig.add_trace(
        go.Bar(
            x=fam_plot["family"],
            y=fam_plot["count"],
            marker_color="#4c78a8",
            showlegend=False,
            hovertemplate="family=%{x}<br>count=%{y}<extra></extra>",
        ),
        row=1,
        col=2,
    )

    heat = (
        cluster_summary.pivot(index="eps_bucket", columns=residue_column, values="count")
        .fillna(0)
        .sort_index(axis=1)
    )
    fig.add_trace(
        go.Heatmap(
            z=heat.values,
            x=[str(int(x)) for x in heat.columns],
            y=heat.index.tolist(),
            colorscale="Magma",
            showscale=True,
            colorbar={"title": "count"},
            hovertemplate="eps=%{y}<br>rho=%{x}<br>count=%{z}<extra></extra>",
        ),
        row=2,
        col=1,
    )

    pair_plot = pair_df.head(top_k).iloc[::-1]
    fig.add_trace(
        go.Bar(
            x=pair_plot["count"],
            y=pair_plot["pair_i"],
            orientation="h",
            marker_color="#59a14f",
            showlegend=False,
            hovertemplate="pair_i=%{y}<br>count=%{x}<extra></extra>",
        ),
        row=2,
        col=2,
    )

    fig.update_layout(
        title=(
            "Goldbach h-Residual Volume Splash"
            f"<br><sup>mod {modulus} native chambers, boost volumes, and pair signatures</sup>"
        ),
        scene={
            "xaxis_title": "volume x",
            "yaxis_title": "volume y",
            "zaxis_title": "volume z",
        },
        legend={"itemsizing": "constant"},
        margin={"l": 10, "r": 10, "t": 80, "b": 10},
        paper_bgcolor="white",
    )

    fig.write_html(html_out, include_plotlyjs=True, full_html=True)


def resolve_outputs(args: argparse.Namespace) -> dict[str, str]:
    stem = args.stem or default_stem(args.max_n, args.modulus)
    return {
        "csv": args.csv_out or f"outputs/csv/{stem}.csv",
        "family": args.family_summary_out or f"outputs/csv/{stem}_family_summary.csv",
        "cluster": args.cluster_summary_out or f"outputs/csv/{stem}_cluster_summary.csv",
        "pairs": args.signatures_out or f"outputs/csv/{stem}_pair_signatures.csv",
        "plot": args.plot_out or f"outputs/plots/{stem}_dashboard.png",
        "html": args.html_out or f"outputs/html/{stem}.html",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=5000)
    parser.add_argument("--modulus", type=int, default=30)
    parser.add_argument("--stem", default=None)
    parser.add_argument("--csv-out", default=None)
    parser.add_argument("--family-summary-out", default=None)
    parser.add_argument("--cluster-summary-out", default=None)
    parser.add_argument("--signatures-out", default=None)
    parser.add_argument("--plot-out", default=None)
    parser.add_argument("--html-out", default=None)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--html", action="store_true")
    parser.add_argument("--pair-top-k", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = resolve_outputs(args)
    residue_column = "rho30" if args.modulus == 30 else "rho_mod"

    df, is_prime, prime_list = build_volume_dataset(args.max_n, args.modulus)
    family_summary = summarize_families(df)
    cluster_summary = summarize_clusters(df, residue_column)

    print("[6/6] Enumerate pair signatures...")
    pair_df = pair_signature_counts(args.max_n, args.modulus, prime_list, is_prime)

    for path in (outputs["csv"], outputs["family"], outputs["cluster"], outputs["pairs"]):
        ensure_parent_dir(path)

    df.to_csv(outputs["csv"], index=False)
    family_summary.to_csv(outputs["family"], index=False)
    cluster_summary.to_csv(outputs["cluster"], index=False)
    pair_df.to_csv(outputs["pairs"], index=False)

    if args.plot or args.plot_out:
        plot_dashboard(
            df=df,
            family_summary=family_summary,
            cluster_summary=cluster_summary,
            pair_df=pair_df,
            residue_column=residue_column,
            plot_out=outputs["plot"],
            top_k=args.pair_top_k,
        )

    if args.html or args.html_out:
        write_html_splash(
            df=df,
            family_summary=family_summary,
            cluster_summary=cluster_summary,
            pair_df=pair_df,
            residue_column=residue_column,
            modulus=args.modulus,
            html_out=outputs["html"],
            top_k=args.pair_top_k,
        )

    print("Done.")
    print(f"rows_csv={outputs['csv']}")
    print(f"family_summary_csv={outputs['family']}")
    print(f"cluster_summary_csv={outputs['cluster']}")
    print(f"pair_signatures_csv={outputs['pairs']}")
    if args.plot or args.plot_out:
        print(f"plot_png={outputs['plot']}")
    if args.html or args.html_out:
        print(f"html_out={outputs['html']}")


if __name__ == "__main__":
    main()
