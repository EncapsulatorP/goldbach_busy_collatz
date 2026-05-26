"""
Goldbach-native filter and tier stratification.

Implements the *heuristic* Hardy-Littlewood predictor for the number of
Goldbach representations of an even target 2n, plus the local-boost factor
G(n) and a stratification of even targets by small-prime structure.

NOTE ON STATUS
--------------
Nothing here proves anything about the Goldbach conjecture. h(2n) is a
heuristic *expected count*, useful for ranking/prioritising a search, not a
sieve with guarantees. Treat "filter" as "ordering / triage", never "exclusion
with certainty".
"""

from __future__ import annotations
from math import log, prod
from sympy import primerange, factorint

# Hardy-Littlewood twin-prime constant C_2 = prod_{p>2} (1 - 1/(p-1)^2)
TWIN_PRIME_C2 = 0.6601618158468695739278121100145


# ---------------------------------------------------------------------------
# Local boost G(n)
# ---------------------------------------------------------------------------
def odd_prime_factors(n: int) -> list[int]:
    """Distinct odd prime divisors of n."""
    return [p for p in factorint(n) if p > 2]


def local_boost(n: int) -> float:
    """
    G(n) = prod_{p | n, p > 2} (p - 1) / (p - 2).

    Each odd prime divisor p inflates the expected Goldbach count.
    p = 3 contributes the maximal single factor: (3-1)/(3-2) = 2.
    A pure power of 2 has no odd prime divisors -> G(n) = 1.
    """
    return prod((p - 1) / (p - 2) for p in odd_prime_factors(n))


# ---------------------------------------------------------------------------
# Heuristic Goldbach representation count h(2n)
# ---------------------------------------------------------------------------
def goldbach_heuristic(n: int) -> float:
    """
    h(2n) ~ (2n / log^2(2n)) * C_2 * G(n).

    Hardy-Littlewood heuristic for the number of ways 2n can be written
    as an (ordered-ish) sum of two primes. Undefined / meaningless for
    2n <= 2; guard the caller.
    """
    m = 2 * n
    if m < 4:
        return 0.0
    return (m / log(m) ** 2) * 2.0 * TWIN_PRIME_C2 * local_boost(n)


# ---------------------------------------------------------------------------
# Tier stratification
# ---------------------------------------------------------------------------
def tier(n: int) -> int:
    """
    Classify the even target 2n by the small-prime structure of n.

      Tier 0 : n = 2^k                  -- no odd-prime boost, very sparse
      Tier 1 : exactly one odd prime    -- single odd-prime boost
      Tier 2 : >= 2 odd primes, 3 | n   -- includes the maximal p=3 factor
      Tier 3 : 15 | n  (3 and 5 | n)    -- two strong low-prime boosts
      Tier 4 : >= 3 distinct odd primes -- high singular-series boost

    A number lands in exactly one tier; checked most-specific first.
    Tiers 2-4 overlap conceptually, so order matters: 4 before 3 before 2.
    """
    opf = set(odd_prime_factors(n))

    if not opf:                               # n is a power of 2 (incl. n = 1)
        return 0
    if len(opf) >= 3:                         # rich singular series
        return 4
    if {3, 5} <= opf:                         # exactly the 3,5 pair
        return 3
    if 3 in opf and len(opf) >= 2:            # 3 plus one other odd prime
        return 2
    if len(opf) == 1:                         # one distinct odd prime (any power)
        return 1
    return 2                                  # two odd primes, neither is 3


def native_filter(n: int, tau, gamma: float = 1.0) -> bool:
    """
    F_Goldbach(n) = [ h(2n) > tau(X) ]  AND  [ G(n) > gamma ].

    `tau` is a callable threshold (e.g. lambda X: log(X)) so it can scale
    with the search bound. Returns True = keep / prioritise.
    """
    return goldbach_heuristic(n) > tau and local_boost(n) > gamma


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    X = 2000                       # search even targets 2n <= X, i.e. n <= X//2
    tau_X = log(X)                 # threshold scales with the bound

    tier_counts = {t: 0 for t in range(5)}
    kept = 0

    print(f"{'n':>5} {'2n':>6} {'G(n)':>8} {'h(2n)':>10} {'tier':>5} {'keep':>5}")
    print("-" * 46)
    for n in range(2, X // 2 + 1):
        g = local_boost(n)
        h = goldbach_heuristic(n)
        t = tier(n)
        keep = native_filter(n, tau_X)
        tier_counts[t] += 1
        kept += keep
        if n <= 30 or n % 250 == 0:          # sample rows only
            print(f"{n:>5} {2*n:>6} {g:>8.4f} {h:>10.2f} {t:>5} {str(keep):>5}")

    total = X // 2 - 1
    print("-" * 46)
    print(f"\nTotals over 2n <= {X}  ({total} even targets):")
    for t, c in tier_counts.items():
        print(f"  Tier {t}: {c:>5}  ({c/total:6.1%})")
    print(f"  Kept by native_filter (h > log X, G > 1): {kept}  ({kept/total:.1%})")
