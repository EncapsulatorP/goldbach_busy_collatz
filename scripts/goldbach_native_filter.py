#!/usr/bin/env python3
"""
Goldbach-native filter and tier stratification.

This script provides a lightweight triage layer for even targets `2n`:

- a local singular-series boost `G`
- a small-prime tier label
- a heuristic ordering score `h`
- a simple keep/drop filter based on `h` and `G`

The implementation is aligned with the repo's current raw Goldbach heuristic:

- unordered Goldbach-count convention
- singular-series boost over the odd prime support
- reciprocal-log density convolution when a finite search bound is known

This is still a heuristic prioritization tool, not a sieve with guarantees.
"""

from __future__ import annotations

import argparse
from math import log, prod
from pathlib import Path

import pandas as pd
from sympy import factorint

from shattering_mirrors import (
    GOLDBACH_C2,
    ensure_parent_dir,
    reciprocal_log_density_convolution,
)


def odd_prime_factors(n: int) -> list[int]:
    """Distinct odd prime divisors of n."""
    return sorted(int(p) for p in factorint(n) if int(p) > 2)


def local_boost(n: int) -> float:
    """
    Goldbach singular-series boost from the odd prime support of n.

    For an even target `2n`, the odd-prime support of `n` and `2n` is the same,
    so the native `n`-based view and the even-target view agree here.
    """
    factors = odd_prime_factors(n)
    if not factors:
        return 1.0
    return float(prod((p - 1) / (p - 2) for p in factors))


def pointwise_goldbach_heuristic(n: int) -> float:
    """
    Legacy pointwise approximation for the unordered Goldbach count at 2n.

    Kept as a fallback when no precomputed density convolution is available.
    """
    m = 2 * n
    if m < 4:
        return 0.0
    return float((m / log(m) ** 2) * GOLDBACH_C2 * local_boost(n))


def goldbach_heuristic(n: int, density_conv=None) -> float:
    """
    Raw heuristic Goldbach count for the even target 2n.

    If `density_conv` is available, use the repo's current finite-range raw
    heuristic:

        h(2n) = C2 * G(n) * sum_k 1/log(k)1/log(2n-k)

    Otherwise fall back to the older pointwise approximation.
    """
    m = 2 * n
    if m < 4:
        return 0.0
    if density_conv is not None and m < len(density_conv):
        return float(GOLDBACH_C2 * local_boost(n) * density_conv[m])
    return pointwise_goldbach_heuristic(n)


def tier(n: int) -> int:
    """
    Classify the even target 2n by the small-prime structure of n.

      Tier 0 : n = 2^k                  -- no odd-prime boost, very sparse
      Tier 1 : exactly one odd prime    -- single odd-prime boost
      Tier 2 : >= 2 odd primes, 3 | n   -- includes the maximal p=3 factor
      Tier 3 : 15 | n  (3 and 5 | n)    -- two strong low-prime boosts
      Tier 4 : >= 3 distinct odd primes -- high singular-series boost

    A number lands in exactly one tier; checked most-specific first.
    """
    opf = set(odd_prime_factors(n))

    if not opf:
        return 0
    if len(opf) >= 3:
        return 4
    if {3, 5} <= opf:
        return 3
    if 3 in opf and len(opf) >= 2:
        return 2
    if len(opf) == 1:
        return 1
    return 2


def resolve_threshold(max_even_target: int, threshold_rule) -> float:
    """
    Resolve a scalar threshold from either a float or a callable of the bound.
    """
    if callable(threshold_rule):
        return float(threshold_rule(max_even_target))
    return float(threshold_rule)


def native_filter(
    n: int,
    threshold: float,
    gamma: float = 1.0,
    density_conv=None,
) -> bool:
    """
    Keep / prioritize the even target 2n if both conditions hold:

        h(2n) > threshold
        G(n)  > gamma
    """
    return goldbach_heuristic(n, density_conv=density_conv) > threshold and local_boost(n) > gamma


def build_native_filter_frame(
    max_even_target: int,
    threshold_rule=None,
    gamma: float = 1.0,
) -> tuple[pd.DataFrame, float]:
    """
    Build a triage frame for all even targets `2n <= max_even_target`.
    """
    if max_even_target < 4:
        raise ValueError("max_even_target must be at least 4")
    if max_even_target % 2:
        raise ValueError("max_even_target must be even")

    if threshold_rule is None:
        threshold_rule = lambda bound: log(bound)

    threshold = resolve_threshold(max_even_target, threshold_rule)
    density_conv = reciprocal_log_density_convolution(max_even_target)

    rows: list[dict] = []
    for n in range(2, max_even_target // 2 + 1):
        factors = odd_prime_factors(n)
        rows.append(
            {
                "n": n,
                "N": 2 * n,
                "boost_G": local_boost(n),
                "boost_primes": ",".join(str(p) for p in factors) or None,
                "h": goldbach_heuristic(n, density_conv=density_conv),
                "tier": tier(n),
                "keep": native_filter(n, threshold, gamma=gamma, density_conv=density_conv),
            }
        )

    df = pd.DataFrame(rows)
    return df, threshold


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-even-target", type=int, default=2000)
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Scalar h-threshold. Defaults to log(max-even-target).",
    )
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument(
        "--out",
        type=str,
        default="outputs/csv/goldbach_native_filter.csv",
    )
    parser.add_argument(
        "--sample-step",
        type=int,
        default=250,
        help="Print every kth row after the first 30 n values.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    threshold_rule = args.threshold if args.threshold is not None else None
    df, threshold = build_native_filter_frame(
        max_even_target=args.max_even_target,
        threshold_rule=threshold_rule,
        gamma=args.gamma,
    )

    ensure_parent_dir(args.out)
    df.to_csv(args.out, index=False)

    tier_counts = df["tier"].value_counts().sort_index()
    kept = int(df["keep"].sum())
    total = int(len(df))

    print(f"threshold={threshold:.6f}")
    print(f"gamma={args.gamma:.6f}")
    print(f"{'n':>5} {'2n':>6} {'G(n)':>8} {'h(2n)':>10} {'tier':>5} {'keep':>5}")
    print("-" * 46)

    for row in df.itertuples(index=False):
        if row.n <= 30 or row.n % args.sample_step == 0:
            print(
                f"{row.n:>5} {row.N:>6} {row.boost_G:>8.4f} "
                f"{row.h:>10.2f} {row.tier:>5} {str(bool(row.keep)):>5}"
            )

    print("-" * 46)
    print(f"\nTotals over 2n <= {args.max_even_target}  ({total} even targets):")
    for t in range(5):
        count = int(tier_counts.get(t, 0))
        print(f"  Tier {t}: {count:>5}  ({count / total:6.1%})")
    print(f"  Kept by native_filter: {kept}  ({kept / total:.1%})")
    print(f"  Saved CSV to: {args.out}")


if __name__ == "__main__":
    main()
