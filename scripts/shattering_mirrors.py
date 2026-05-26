#!/usr/bin/env python3
"""
Goldbach h-residual mirror clustering.

This script builds the hierarchy:

    fake mirror:
        r | delta10 | r

    h-residue mirror:
        r | eps_h | r

    root decimal:
        eps_h | delta10 | eps_h

    root native:
        eps_h | rho30 | eps_h

where:

    r       = exact unordered Goldbach count r_G(N)
    h       = Hardy-Littlewood expected Goldbach count
    eps_h   = r_G(N) - floor(h(N))
    delta10 = N - 10*r_G(N)
    rho30   = N mod 30

Use this as a diagnostic / clustering tool, not as a proof engine.
"""

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches


GOLDBACH_C2 = 0.6601618158468696
# Product over p > 2 of (1 - 1/(p-1)^2)

Z_BUCKET_STEP = 0.5
Z_BUCKET_CLIP = 2.0


def sieve_bool(limit: int) -> np.ndarray:
    """Boolean prime sieve up to limit inclusive."""
    is_prime = np.ones(limit + 1, dtype=bool)
    if limit >= 0:
        is_prime[0] = False
    if limit >= 1:
        is_prime[1] = False

    max_p = int(limit ** 0.5)
    for p in range(2, max_p + 1):
        if is_prime[p]:
            is_prime[p * p : limit + 1 : p] = False

    return is_prime


def spf_sieve(limit: int) -> np.ndarray:
    """Smallest-prime-factor sieve up to limit inclusive."""
    spf = np.arange(limit + 1, dtype=np.int64)

    if limit >= 0:
        spf[0] = 0
    if limit >= 1:
        spf[1] = 1

    for p in range(2, int(limit ** 0.5) + 1):
        if spf[p] == p:
            start = p * p
            mask = spf[start : limit + 1 : p] == np.arange(
                start, limit + 1, p, dtype=np.int64
            )
            spf[start : limit + 1 : p][mask] = p

    return spf


def exact_goldbach_counts_fft(max_n: int, is_prime: np.ndarray) -> np.ndarray:
    """
    Compute unordered Goldbach counts r_G(N) for all N <= max_n.

    Uses FFT convolution of the prime indicator with itself.

    Ordered count:
        ordered[N] = # ordered prime pairs (p, q) with p + q = N

    Unordered count:
        r_G(N) = (ordered[N] + self_pair[N]) // 2

    where self_pair[N] = 1 if N/2 is prime, else 0.
    """
    arr = is_prime.astype(np.float64)

    size = 1 << ((2 * max_n + 1) - 1).bit_length()

    fft_arr = np.fft.rfft(arr, size)
    conv = np.fft.irfft(fft_arr * fft_arr, size)
    ordered = np.rint(conv[: max_n + 1]).astype(np.int64)

    r = np.zeros(max_n + 1, dtype=np.int64)

    for n in range(4, max_n + 1, 2):
        self_pair = 1 if is_prime[n // 2] else 0
        r[n] = (ordered[n] + self_pair) // 2

    return r


def local_goldbach_boost(n: int, spf: np.ndarray) -> float:
    """
    Compute:

        G(n) = product_{p | n, p > 2} (p - 1)/(p - 2)

    using unique odd prime divisors of n.
    """
    if n <= 1:
        return 1.0

    boost = 1.0
    last_p = None
    x = n

    while x > 1:
        p = int(spf[x])
        if p == 0:
            break

        if p != last_p and p > 2:
            boost *= (p - 1) / (p - 2)

        last_p = p

        while x % p == 0:
            x //= p

    return boost


def h_goldbach(N: int, boost: float) -> float:
    """
    Hardy-Littlewood-style expected unordered Goldbach count:

        h(N) ≈ N / log(N)^2 * C2 * G(N/2)
    """
    if N < 4:
        return 0.0

    return (N / (math.log(N) ** 2)) * GOLDBACH_C2 * boost


def eps_label(eps: int, clip: int = 5) -> str:
    """
    Bucket residuals so large outliers do not explode the cluster tree.
    """
    if eps < -clip:
        return f"<-{clip}"
    if eps > clip:
        return f">{clip}"
    return str(eps)


def z_label(z: float, step: float = Z_BUCKET_STEP, clip: float = Z_BUCKET_CLIP) -> str:
    """
    Bucket normalized residuals on a fixed scale.

    The clipped edge buckets exist only to keep the chamber table finite; the
    meaningful comparison is whether occupancy concentrates away from zero.
    """
    if z <= -clip:
        return f"<={-clip:.1f}"
    if z >= clip:
        return f">={clip:.1f}"

    bucket = round(z / step) * step
    if abs(bucket) < step / 2:
        bucket = 0.0
    return f"{bucket:+.1f}"


@dataclass
class CompressionRow:
    N: int
    n: int
    r: int
    h: float
    h_floor: int
    eps_h: int
    z_h_raw: float
    h_scale: float
    h_cal: float
    z_h: float
    boost_G: float
    delta10: int
    rho30: int
    fake_mirror: str
    h_mirror: str
    root_decimal: str
    root_native: str
    native_core: str
    eps_bucket: str
    z_bucket: str
    native_cluster: str


def build_row(N: int, r: int, spf: np.ndarray) -> CompressionRow:
    n = N // 2
    boost = local_goldbach_boost(n, spf)
    h = h_goldbach(N, boost)

    h_floor = math.floor(h)
    eps_h = int(r - h_floor)

    z_h_raw = 0.0
    if h > 0:
        z_h_raw = (r - h) / math.sqrt(h)

    delta10 = int(N - 10 * r)
    rho30 = int(N % 30)

    fake_mirror = f"{r}|{delta10}|{r}"
    h_mirror = f"{r}|{eps_h}|{r}"
    root_decimal = f"{eps_h}|{delta10}|{eps_h}"
    root_native = f"{eps_h}|{rho30}|{eps_h}"
    native_core = f"{r}|({rho30},{eps_h})|{r}"

    bucket = eps_label(eps_h)
    return CompressionRow(
        N=N,
        n=n,
        r=int(r),
        h=float(h),
        h_floor=int(h_floor),
        eps_h=int(eps_h),
        z_h_raw=float(z_h_raw),
        h_scale=1.0,
        h_cal=float(h),
        z_h=float(z_h_raw),
        boost_G=float(boost),
        delta10=delta10,
        rho30=rho30,
        fake_mirror=fake_mirror,
        h_mirror=h_mirror,
        root_decimal=root_decimal,
        root_native=root_native,
        native_core=native_core,
        eps_bucket=bucket,
        z_bucket="",
        native_cluster="",
    )


def calibrated_h_scale(df: pd.DataFrame) -> float:
    """
    Fit one global multiplicative scale for h(N) so mean normalized residual is 0.

    This does not upgrade the heuristic into a theorem; it only removes the
    single global bias before bucketing residuals.
    """
    positive = df["h"] > 0
    sqrt_h = np.sqrt(df.loc[positive, "h"].to_numpy())
    r = df.loc[positive, "r"].to_numpy()

    if len(sqrt_h) == 0:
        return 1.0

    return float((r / sqrt_h).sum() / sqrt_h.sum())


def apply_residual_calibration(df: pd.DataFrame) -> pd.DataFrame:
    """Attach calibrated h and normalized residual buckets to the dataset."""
    out = df.copy()
    scale = calibrated_h_scale(out)
    out["h_scale"] = scale
    out["h_cal"] = out["h"] * scale

    z = np.zeros(len(out), dtype=np.float64)
    positive = out["h_cal"] > 0
    z[positive] = (
        out.loc[positive, "r"].to_numpy() - out.loc[positive, "h_cal"].to_numpy()
    ) / np.sqrt(out.loc[positive, "h_cal"].to_numpy())
    out["z_h"] = z
    out["z_bucket"] = [z_label(value) for value in out["z_h"]]

    if "native_cluster" in out.columns:
        out["native_cluster"] = [
            f"z={z_bucket};rho30={rho30}"
            for z_bucket, rho30 in zip(out["z_bucket"], out["rho30"])
        ]

    return out


def ensure_parent_dir(path: str) -> None:
    parent = Path(path).expanduser().resolve().parent
    parent.mkdir(parents=True, exist_ok=True)


def build_dataset(max_n: int, include_strings: bool = True) -> pd.DataFrame:
    """
    Build full dataset for even N from 4 to max_n.

    For very large max_n, set include_strings=False to reduce memory.
    """
    if max_n < 4:
        raise ValueError("max_n must be >= 4")

    print(f"[1/5] Sieve primes up to {max_n}...")
    is_prime = sieve_bool(max_n)

    print(f"[2/5] Build SPF up to {max_n // 2}...")
    spf = spf_sieve(max_n // 2)

    print("[3/5] Compute exact Goldbach counts by FFT...")
    r_counts = exact_goldbach_counts_fft(max_n, is_prime)

    print("[4/5] Build h-residual compression rows...")
    rows = []

    for N in range(4, max_n + 1, 2):
        row = build_row(N, int(r_counts[N]), spf)

        if include_strings:
            rows.append(row.__dict__)
        else:
            rows.append(
                {
                    "N": row.N,
                    "n": row.n,
                    "r": row.r,
                    "h": row.h,
                    "h_floor": row.h_floor,
                    "eps_h": row.eps_h,
                    "z_h": row.z_h,
                    "boost_G": row.boost_G,
                    "delta10": row.delta10,
                    "rho30": row.rho30,
                    "eps_bucket": row.eps_bucket,
                }
            )

    df = pd.DataFrame(rows)
    df = apply_residual_calibration(df)

    print("[5/5] Done.")
    return df


def summarize_clusters(df: pd.DataFrame, top: int = 30) -> pd.DataFrame:
    """
    Summary by normalized residual bucket and native mod-30 chamber.
    """
    summary = (
        df.groupby(["z_bucket", "rho30"])
        .agg(
            count=("N", "count"),
            min_N=("N", "min"),
            max_N=("N", "max"),
            mean_r=("r", "mean"),
            mean_h=("h", "mean"),
            mean_h_cal=("h_cal", "mean"),
            mean_eps=("eps_h", "mean"),
            mean_z=("z_h", "mean"),
            mean_boost=("boost_G", "mean"),
        )
        .reset_index()
        .sort_values(["count"], ascending=False)
    )

    return summary.head(top)


def full_cluster_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Full native-cluster counts, not truncated to the top rows."""
    return (
        df.groupby(["z_bucket", "rho30", "native_cluster"])
        .agg(
            count=("N", "count"),
            min_N=("N", "min"),
            max_N=("N", "max"),
            mean_r=("r", "mean"),
            mean_h=("h", "mean"),
            mean_h_cal=("h_cal", "mean"),
            mean_eps=("eps_h", "mean"),
            mean_z=("z_h", "mean"),
            mean_boost=("boost_G", "mean"),
        )
        .reset_index()
        .sort_values(["count", "z_bucket", "rho30"], ascending=[False, True, True])
    )


def residual_bucket_sort_key(label: str) -> tuple[int, float]:
    """Stable ordering for clipped normalized-residual buckets."""
    if label.startswith("<="):
        return (-10_000, float(label[2:]))
    if label.startswith(">="):
        return (10_000, float(label[2:]))
    return (0, float(label))


def ordered_residual_buckets(df: pd.DataFrame) -> list[str]:
    labels = sorted(df["z_bucket"].astype(str).unique(), key=residual_bucket_sort_key)
    return labels


def resolve_cluster_filter(
    df: pd.DataFrame,
    requested: str,
    top_clusters: int,
) -> tuple[list[str], pd.DataFrame]:
    """Resolve either an explicit cluster filter or the top native clusters."""
    counts = full_cluster_counts(df)
    known = set(counts["native_cluster"])

    if requested.strip():
        labels = [item.strip() for item in requested.split(",") if item.strip()]
        missing = [label for label in labels if label not in known]
        if missing:
            raise ValueError(
                "unknown native_cluster label(s): "
                + ", ".join(missing)
                + ". Use labels like 'z=+0.5;rho30=6'."
            )
    else:
        labels = counts["native_cluster"].head(top_clusters).tolist()

    filtered_counts = counts[counts["native_cluster"].isin(labels)].copy()
    label_order = {label: i for i, label in enumerate(labels)}
    filtered_counts["_order"] = filtered_counts["native_cluster"].map(label_order)
    filtered_counts = filtered_counts.sort_values("_order").drop(columns="_order")
    return labels, filtered_counts


def filtered_cluster_frame(df: pd.DataFrame, selected_labels: list[str]) -> pd.DataFrame:
    """Annotate the dataframe with a selection mask for plotting."""
    out = df.copy()
    selected = set(selected_labels)
    out["cluster_view"] = np.where(out["native_cluster"].isin(selected), out["native_cluster"], "other")
    return out


def make_cluster_palette(labels: list[str]) -> dict[str, tuple[float, float, float, float]]:
    cmap = plt.get_cmap("tab10" if len(labels) <= 10 else "tab20")
    return {label: cmap(i % cmap.N) for i, label in enumerate(labels)}


def build_cluster_share_table(
    df: pd.DataFrame,
    labels: list[str],
    bins: int,
) -> pd.DataFrame:
    """Share of each selected cluster in linear N-bins."""
    work = df[["N", "native_cluster"]].copy()
    work["N_bin"] = pd.cut(work["N"], bins=bins, include_lowest=True)
    grouped = (
        work.groupby(["N_bin", "native_cluster"], observed=False)
        .size()
        .rename("count")
        .reset_index()
    )
    totals = grouped.groupby("N_bin", observed=False)["count"].sum().rename("bin_total").reset_index()
    grouped = grouped.merge(totals, on="N_bin", how="left")
    grouped["share"] = grouped["count"] / grouped["bin_total"]
    grouped = grouped[grouped["native_cluster"].isin(labels)].copy()
    grouped["bin_mid"] = grouped["N_bin"].apply(lambda interval: 0.5 * (interval.left + interval.right))
    return grouped


def build_cluster_z_table(
    df: pd.DataFrame,
    labels: list[str],
    bins: int,
) -> pd.DataFrame:
    """Mean z_h by selected cluster in linear N-bins."""
    work = df[df["native_cluster"].isin(labels)][["N", "native_cluster", "z_h"]].copy()
    work["N_bin"] = pd.cut(work["N"], bins=bins, include_lowest=True)
    grouped = (
        work.groupby(["N_bin", "native_cluster"], observed=False)
        .agg(mean_z=("z_h", "mean"), count=("z_h", "size"))
        .reset_index()
    )
    grouped["bin_mid"] = grouped["N_bin"].apply(lambda interval: 0.5 * (interval.left + interval.right))
    return grouped


def plot_cluster_dashboard(
    df: pd.DataFrame,
    selected_labels: list[str],
    filtered_counts: pd.DataFrame,
    out_path: str,
    bins: int = 50,
) -> None:
    """Create a multi-panel explanatory dashboard for the native cluster filter."""
    plot_df = filtered_cluster_frame(df, selected_labels)
    palette = make_cluster_palette(selected_labels)
    z_labels = ordered_residual_buckets(df)
    rho_labels = sorted(df["rho30"].astype(int).unique())

    heat = (
        df.groupby(["z_bucket", "rho30"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=z_labels, columns=rho_labels, fill_value=0)
    )

    share_table = build_cluster_share_table(df, selected_labels, bins=bins)
    z_table = build_cluster_z_table(df, selected_labels, bins=bins)

    fig, axes = plt.subplots(3, 2, figsize=(15, 16))
    ax_heat, ax_bar, ax_z, ax_delta, ax_share, ax_centroid = axes.ravel()

    im = ax_heat.imshow(heat.values, aspect="auto", cmap="magma")
    ax_heat.set_title("Normalized residual chambers: count heatmap")
    ax_heat.set_xlabel("rho30")
    ax_heat.set_ylabel("z_bucket")
    ax_heat.set_xticks(range(len(rho_labels)))
    ax_heat.set_xticklabels(rho_labels)
    ax_heat.set_yticks(range(len(z_labels)))
    ax_heat.set_yticklabels(z_labels)
    for i, label in enumerate(selected_labels):
        z_bucket, rho_text = label.split(";rho30=")
        j = rho_labels.index(int(rho_text))
        ii = z_labels.index(z_bucket.replace("z=", ""))
        rect = patches.Rectangle((j - 0.5, ii - 0.5), 1, 1, fill=False, lw=2.5, ec=palette[label])
        ax_heat.add_patch(rect)
        ax_heat.text(j, ii, str(i + 1), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04, label="count")

    ax_bar.set_title("Selected cluster sizes")
    bar_positions = np.arange(len(filtered_counts))
    bar_colors = [palette[label] for label in filtered_counts["native_cluster"]]
    ax_bar.bar(bar_positions, filtered_counts["count"], color=bar_colors, alpha=0.9)
    ax_bar.set_xticks(bar_positions)
    ax_bar.set_xticklabels([f"{i + 1}" for i in range(len(filtered_counts))])
    ax_bar.set_xlabel("selected cluster id")
    ax_bar.set_ylabel("count")
    for pos, (_, row) in enumerate(filtered_counts.iterrows()):
        ax_bar.text(pos, row["count"], f"{row['rho30']}\n{row['z_bucket']}", ha="center", va="bottom", fontsize=8)

    ax_z.set_title("N vs calibrated z_h with cluster filter highlight")
    ax_z.scatter(df["N"], df["z_h"], s=7, color="#b0b0b0", alpha=0.18, edgecolors="none", label="all points")
    for label in selected_labels:
        subset = plot_df[plot_df["native_cluster"] == label]
        ax_z.scatter(subset["N"], subset["z_h"], s=10, color=palette[label], alpha=0.8, edgecolors="none", label=label)
    ax_z.axhline(0.0, color="#444444", lw=1.0)
    ax_z.set_xlabel("N")
    ax_z.set_ylabel("z_h = (r - alpha*h)/sqrt(alpha*h)")

    ax_delta.set_title("N vs delta10 with cluster filter highlight")
    ax_delta.scatter(df["N"], df["delta10"], s=7, color="#b0b0b0", alpha=0.18, edgecolors="none")
    for label in selected_labels:
        subset = plot_df[plot_df["native_cluster"] == label]
        ax_delta.scatter(subset["N"], subset["delta10"], s=10, color=palette[label], alpha=0.8, edgecolors="none")
    ax_delta.axhline(0.0, color="#444444", lw=1.0)
    ax_delta.set_xlabel("N")
    ax_delta.set_ylabel("delta10 = N - 10 r")

    ax_share.set_title("Selected cluster share across N-bins")
    for label in selected_labels:
        subset = share_table[share_table["native_cluster"] == label]
        ax_share.plot(subset["bin_mid"], subset["share"], color=palette[label], lw=2, label=label)
    ax_share.set_xlabel("N-bin midpoint")
    ax_share.set_ylabel("share of all rows in bin")

    ax_centroid.set_title("Selected cluster centroids")
    ax_centroid.scatter(
        filtered_counts["mean_boost"],
        filtered_counts["mean_z"],
        s=np.clip(filtered_counts["count"] / 2, 30, 500),
        color=[palette[label] for label in filtered_counts["native_cluster"]],
        alpha=0.85,
        edgecolors="black",
        linewidths=0.4,
    )
    for i, (_, row) in enumerate(filtered_counts.iterrows()):
        ax_centroid.text(
            row["mean_boost"],
            row["mean_z"],
            f"{i + 1}",
            ha="center",
            va="center",
            color="white",
            fontsize=9,
            fontweight="bold",
        )
    ax_centroid.set_xlabel("mean Hardy-Littlewood boost G")
    ax_centroid.set_ylabel("mean z_h")

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=palette[label], label=f"{i + 1}: {label}")
        for i, label in enumerate(selected_labels)
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle("Shattering Mirrors Cluster Dashboard", fontsize=17, y=0.995)
    fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_cluster_gallery(
    df: pd.DataFrame,
    selected_labels: list[str],
    filtered_counts: pd.DataFrame,
    out_path: str,
    bins: int = 50,
) -> None:
    """Create one panel per selected cluster with local trend summaries."""
    palette = make_cluster_palette(selected_labels)
    z_table = build_cluster_z_table(df, selected_labels, bins=bins)

    n_panels = len(selected_labels)
    ncols = 2
    nrows = int(math.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4.2 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for idx, label in enumerate(selected_labels):
        ax = axes[idx]
        subset = df[df["native_cluster"] == label]
        meta = filtered_counts[filtered_counts["native_cluster"] == label].iloc[0]
        z_subset = z_table[z_table["native_cluster"] == label]

        ax.scatter(subset["N"], subset["z_h"], s=10, color=palette[label], alpha=0.55, edgecolors="none")
        ax.plot(z_subset["bin_mid"], z_subset["mean_z"], color="black", lw=2.0, alpha=0.85)
        ax.axhline(0.0, color="#444444", lw=1.0)
        ax.set_title(f"{idx + 1}. {label}")
        ax.set_xlabel("N")
        ax.set_ylabel("z_h")
        ax.text(
            0.02,
            0.98,
            (
                f"count={int(meta['count'])}\n"
                f"N range=[{int(meta['min_N'])}, {int(meta['max_N'])}]\n"
                f"mean_r={meta['mean_r']:.2f}, mean_h={meta['mean_h']:.2f}\n"
                f"mean_h_cal={meta['mean_h_cal']:.2f}, mean_eps={meta['mean_eps']:.2f}\n"
                f"mean_boost={meta['mean_boost']:.3f}"
            ),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
        )

    for ax in axes[n_panels:]:
        ax.axis("off")

    fig.suptitle("Shattering Mirrors Selected Cluster Gallery", fontsize=17, y=0.995)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def find_decimal_mirror_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Strict fake mirror hits:

        r_G(N) = floor(N/10)

    These are base-10 artifacts / negative controls.
    """
    return df[df["r"] == (df["N"] // 10)].copy()


def add_kmeans_labels(df: pd.DataFrame, k: int = 8) -> pd.DataFrame:
    """
    Optional k-means on numeric features.

    Symbolic hierarchy remains:
        z_bucket + rho30

    K-means is only a secondary numeric grouping.
    """
    try:
        from sklearn.cluster import KMeans
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for k-means. Install with: pip install scikit-learn"
        ) from exc

    out = df.copy()

    features = pd.DataFrame(
        {
            "log_N": np.log(out["N"]),
            "r": out["r"],
            "h": out["h"],
            "z_h": out["z_h"],
            "boost_G": out["boost_G"],
            "rho30": out["rho30"],
            "delta10": out["delta10"],
        }
    )

    X = StandardScaler().fit_transform(features)

    km = KMeans(n_clusters=k, random_state=42, n_init="auto")
    out["kmeans_cluster"] = km.fit_predict(X)

    out["hybrid_cluster"] = (
        out["z_bucket"].astype(str)
        + "|rho30="
        + out["rho30"].astype(str)
        + "|k="
        + out["kmeans_cluster"].astype(str)
    )

    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=100_000)
    parser.add_argument("--out", type=str, default="outputs/csv/goldbach_h_clusters.csv")
    parser.add_argument("--summary-out", type=str, default="outputs/csv/goldbach_h_cluster_summary.csv")
    parser.add_argument("--mirror-out", type=str, default="outputs/csv/goldbach_decimal_mirror_hits.csv")
    parser.add_argument("--no-strings", action="store_true")
    parser.add_argument("--kmeans", action="store_true")
    parser.add_argument("--k", type=int, default=8)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument(
        "--plot-prefix",
        type=str,
        default="outputs/plots/goldbach_h_cluster_plots",
        help="Prefix for dashboard/gallery plot files.",
    )
    parser.add_argument(
        "--plot-top-clusters",
        type=int,
        default=8,
        help="If no explicit cluster filter is provided, highlight this many native clusters by count.",
    )
    parser.add_argument(
        "--cluster-filter",
        type=str,
        default="",
        help="Comma-separated native_cluster labels to highlight, e.g. 'z=+0.5;rho30=6,z=-0.5;rho30=0'.",
    )
    parser.add_argument(
        "--plot-bins",
        type=int,
        default=48,
        help="Number of linear N-bins used in the occupancy and mean-z plots.",
    )

    args = parser.parse_args()

    df = build_dataset(args.max_n, include_strings=not args.no_strings)

    if args.kmeans:
        print(f"Running k-means with k={args.k}...")
        df = add_kmeans_labels(df, k=args.k)

    summary = summarize_clusters(df)
    mirror_hits = find_decimal_mirror_hits(df)

    ensure_parent_dir(args.out)
    ensure_parent_dir(args.summary_out)
    ensure_parent_dir(args.mirror_out)
    df.to_csv(args.out, index=False)
    summary.to_csv(args.summary_out, index=False)
    mirror_hits.to_csv(args.mirror_out, index=False)

    print(f"\nCalibrated h-scale alpha={df['h_scale'].iloc[0]:.6f}")
    print("\nTop normalized residual clusters:")
    print(summary.to_string(index=False))

    print("\nStrict decimal mirror hits:")
    cols = ["N", "r", "h", "h_cal", "eps_h", "z_h", "z_bucket", "delta10", "rho30"]
    extra_cols = [
        "fake_mirror",
        "h_mirror",
        "root_decimal",
        "root_native",
        "native_core",
    ]
    display_cols = [c for c in cols + extra_cols if c in mirror_hits.columns]

    if len(mirror_hits) == 0:
        print("No strict mirror hits found.")
    else:
        print(mirror_hits[display_cols].to_string(index=False))

    print(f"\nSaved full dataset to: {args.out}")
    print(f"Saved summary to: {args.summary_out}")
    print(f"Saved mirror hits to: {args.mirror_out}")

    if args.plot:
        if args.no_strings:
            raise ValueError("--plot requires the string cluster columns; rerun without --no-strings")

        selected_labels, filtered_counts = resolve_cluster_filter(
            df,
            requested=args.cluster_filter,
            top_clusters=args.plot_top_clusters,
        )
        dashboard_path = f"{args.plot_prefix}_dashboard.png"
        gallery_path = f"{args.plot_prefix}_gallery.png"
        ensure_parent_dir(dashboard_path)
        ensure_parent_dir(gallery_path)

        print("\nGenerating explanatory plots with native cluster filter:")
        for i, label in enumerate(selected_labels, start=1):
            print(f"  {i}. {label}")

        plot_cluster_dashboard(
            df,
            selected_labels=selected_labels,
            filtered_counts=filtered_counts,
            out_path=dashboard_path,
            bins=args.plot_bins,
        )
        plot_cluster_gallery(
            df,
            selected_labels=selected_labels,
            filtered_counts=filtered_counts,
            out_path=gallery_path,
            bins=args.plot_bins,
        )
        print(f"Saved dashboard plot to: {dashboard_path}")
        print(f"Saved gallery plot to: {gallery_path}")


if __name__ == "__main__":
    main()
