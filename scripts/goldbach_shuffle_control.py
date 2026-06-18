#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


C2 = 0.6601618158468695  # twin-prime constant


def ensure_parent_dir(path: str) -> None:
    Path(path).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)


def build_prime_tables(max_n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sieve = np.ones(max_n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(max_n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    primes = np.nonzero(sieve)[0]
    small_primes = primes[primes <= int(max_n**0.5) + 2]
    return sieve, primes, small_primes


def odd_prime_factors(n: int, small_primes: np.ndarray) -> list[int]:
    facs = []
    m = n
    while m % 2 == 0:
        m //= 2
    for p in small_primes:
        if p == 2:
            continue
        p = int(p)
        if p * p > m:
            break
        if m % p == 0:
            facs.append(p)
            while m % p == 0:
                m //= p
    if m > 1 and m % 2 == 1:
        facs.append(int(m))
    return facs


def singular_series(n: int, small_primes: np.ndarray) -> float:
    series = 2.0 * C2
    for p in odd_prime_factors(n, small_primes):
        series *= (p - 1) / (p - 2)
    return series


def run_experiment(max_n: int, sample_size: int, seed: int) -> dict[str, np.ndarray | float]:
    sieve, primes, small_primes = build_prime_tables(max_n)
    rng = np.random.default_rng(seed)

    candidates = np.arange(2000, max_n, 2)
    if sample_size > len(candidates):
        raise ValueError("sample_size exceeds the number of available even targets")
    sample = np.sort(rng.choice(candidates, size=sample_size, replace=False))

    reps = np.empty(len(sample))
    series = np.empty(len(sample))
    for k, n in enumerate(sample):
        jmax = np.searchsorted(primes, n - 2, side="right")
        pl = primes[:jmax]
        reps[k] = np.count_nonzero(sieve[n - pl])
        series[k] = singular_series(int(n), small_primes)

    ln_n = np.log(sample.astype(float))
    detrended = reps * ln_n**2 / sample
    detrended_shuffled = rng.permutation(detrended)

    return {
        "sample": sample,
        "series": series,
        "detrended": detrended,
        "detrended_shuffled": detrended_shuffled,
        "r_real": float(np.corrcoef(series, detrended)[0, 1]),
        "r_shuf": float(np.corrcoef(series, detrended_shuffled)[0, 1]),
    }


def plot_results(result: dict[str, np.ndarray | float], out: str) -> None:
    sample = result["sample"]
    series = result["series"]
    detrended = result["detrended"]
    detrended_shuffled = result["detrended_shuffled"]
    r_real = float(result["r_real"])
    r_shuf = float(result["r_shuf"])

    div3 = sample % 3 == 0

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 6.0))
    lo = min(series.min(), detrended.min()) - 0.1
    hi = max(series.max(), detrended.max()) + 0.1

    for axis in ax:
        axis.plot([lo, hi], [lo, hi], "k--", lw=1.3, zorder=5, label="y = x  (= S(N))")
        axis.set_xlim(lo, hi)
        axis.set_ylim(lo, hi)
        axis.set_xlabel("singular series  S(N)   (Möbius / inclusion-exclusion built)")

    c_no3 = "#c0392b"
    c_y3 = "#2471a3"

    ax[0].scatter(series[~div3], detrended[~div3], s=6, alpha=0.45, c=c_no3, label="3 ∤ N")
    ax[0].scatter(series[div3], detrended[div3], s=6, alpha=0.45, c=c_y3, label="3 | N")
    ax[0].set_ylabel("detrended r(N) = r(N)·(ln N)^2 / N")
    ax[0].set_title(f"REAL Goldbach counts — residual tracks S(N)\nr = {r_real:.4f}", fontsize=12)
    ax[0].legend(loc="upper left", fontsize=9, framealpha=0.9)

    ax[1].scatter(series[~div3], detrended_shuffled[~div3], s=6, alpha=0.45, c=c_no3)
    ax[1].scatter(series[div3], detrended_shuffled[div3], s=6, alpha=0.45, c=c_y3)
    ax[1].set_ylabel("shuffled detrended r(N)  (same marginal, pairing destroyed)")
    ax[1].set_title(f"SHUFFLED control — diagonal evaporates\nr = {r_shuf:.4f}", fontsize=12)

    fig.suptitle(
        "Shuffle control: is the S(N) alignment data, or an artifact of the transform?",
        fontsize=13,
        y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    ensure_parent_dir(out)
    fig.savefig(out, dpi=130)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=500_000)
    parser.add_argument("--sample-size", type=int, default=4_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="outputs/plots/goldbach_shuffle_control.png")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = run_experiment(max_n=args.max_n, sample_size=args.sample_size, seed=args.seed)
    print(f"corr(S, detrended r(N))      = {float(result['r_real']):.4f}")
    print(f"corr(S, shuffled detrended)  = {float(result['r_shuf']):.4f}")
    plot_results(result, args.out)
    print("saved", args.out)


if __name__ == "__main__":
    main()
