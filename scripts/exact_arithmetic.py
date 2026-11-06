"""Integer Goldbach convolution and shared modular/Mobius research kernels.

NTT coefficients are reconstructed uniquely in a single prime field because
the input is binary and its length is strictly below the modulus. No rounding.
This module proves neither Goldbach nor a global Mertens bound.
"""
from __future__ import annotations

from math import isqrt
from operator import index

import numpy as np

NTT_MODULUS = 998244353  # prime, 119 * 2**23 + 1
NTT_ROOT = 3
MAX_TRANSFORM = 1 << 23
MAX_LIMIT = (MAX_TRANSFORM - 1) // 2


def integer(value, name, minimum=0, maximum=MAX_LIMIT):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        value = index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if not minimum <= value <= maximum:
        raise ValueError(f"{name} must be in [{minimum}, {maximum}]")
    return value


def prime_indicator(limit):
    limit = integer(limit, "limit")
    flags = np.ones(limit + 1, dtype=np.bool_)
    flags[:2] = False
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p::p] = False
    return flags


def _ntt(values, inverse=False):
    """Radix-2 field transform; int64 products stay below 2**63."""
    n = len(values)
    if n < 1 or n & (n - 1) or n > MAX_TRANSFORM:
        raise ValueError("unsupported transform length")
    indices = np.arange(n, dtype=np.uint32)
    rev = np.zeros(n, dtype=np.uint32)
    work = indices.copy()
    for _ in range(n.bit_length() - 1):
        rev = (rev << 1) | (work & 1)
        work >>= 1
    a = np.asarray(values, dtype=np.int64)[rev].copy()
    width = 2
    while width <= n:
        half = width // 2
        root = pow(NTT_ROOT, (NTT_MODULUS - 1) // width, NTT_MODULUS)
        if inverse:
            root = pow(root, NTT_MODULUS - 2, NTT_MODULUS)
        roots = np.ones(half, dtype=np.int64)
        step = 1
        while step < half:
            roots[step:2 * step] = roots[:step] * pow(root, step, NTT_MODULUS) % NTT_MODULUS
            step *= 2
        blocks = a.reshape(-1, width)
        left = blocks[:, :half].copy()
        right = blocks[:, half:] * roots % NTT_MODULUS
        blocks[:, :half] = (left + right) % NTT_MODULUS
        blocks[:, half:] = (left - right) % NTT_MODULUS
        width *= 2
    if inverse:
        a = a * pow(n, NTT_MODULUS - 2, NTT_MODULUS) % NTT_MODULUS
    return a


def binary_self_convolution(flags):
    """Exact ordinary self-convolution of a bounded binary integer vector."""
    arr = np.asarray(flags)
    if arr.ndim != 1 or not 1 <= len(arr) <= MAX_LIMIT + 1:
        raise ValueError("binary vector has unsupported shape or length")
    if arr.dtype.kind not in "biu" or not np.all((arr == 0) | (arr == 1)):
        raise ValueError("expected a binary integer or boolean vector")
    size = 1 << (2 * len(arr) - 2).bit_length()
    padded = np.zeros(size, dtype=np.int64)
    padded[:len(arr)] = arr
    spectrum = _ntt(padded)
    result = _ntt(spectrum * spectrum % NTT_MODULUS, inverse=True)
    # Each integer coefficient <= len(arr) < modulus, so no CRT ambiguity.
    return result[:2 * len(arr) - 1]


def exact_goldbach_counts(max_n, is_prime=None):
    """Unordered prime pairs for even N <= max_n; other entries are zero.

    Optional is_prime is a caller-supplied prime indicator, not itself a
    primality certificate. Shape and binary values are checked. O(N log N)
    time and O(N) memory, with max_n <= MAX_LIMIT.
    """
    max_n = integer(max_n, "max_n")
    flags = prime_indicator(max_n) if is_prime is None else np.asarray(is_prime)
    if flags.ndim != 1 or len(flags) != max_n + 1:
        raise ValueError("prime indicator must have length max_n + 1")
    ordered = binary_self_convolution(flags)[:max_n + 1]
    if any(flags[:2]):
        raise ValueError("zero and one must not be marked prime")
    result = np.zeros(max_n + 1, dtype=np.int64)
    targets = np.arange(4, max_n + 1, 2)
    diagonal = flags[targets // 2].astype(np.int64)
    result[targets] = (ordered[targets] + diagonal) // 2
    return result


def unit_pair_count(target, modulus):
    """Ordered (a,b) units mod m with a+b=target, via prime-power CRT.

    For p**e: p**(e-1)*(p-1 if p divides target else p-2).
    These are local admissibility counts, not counts of actual prime pairs.
    """
    target = integer(target, "target")
    modulus = integer(modulus, "modulus", 2)
    remaining, p, count = modulus, 2, 1
    while p * p <= remaining:
        exponent = 0
        while remaining % p == 0:
            remaining //= p
            exponent += 1
        if exponent:
            count *= p ** (exponent - 1) * (p - 1 if target % p == 0 else p - 2)
        p = 3 if p == 2 else p + 2
    if remaining > 1:
        count *= remaining - 1 if target % remaining == 0 else remaining - 2
    return count


def mobius_prefix(limit):
    """Return exact mu[0..limit] and M[n]=sum(mu[1..n]); mu[0]=0."""
    limit = integer(limit, "limit")
    mu = np.ones(limit + 1, dtype=np.int64)
    mu[0] = 0
    flags = prime_indicator(limit)
    for p in np.flatnonzero(flags):
        mu[p::p] *= -1
        if p * p <= limit:
            mu[p * p::p * p] = 0
    return mu, np.cumsum(mu, dtype=np.int64)


def chi4(n):
    """Primitive real Dirichlet character modulo four."""
    n = integer(n, "n")
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)


def signed_prime_pair_sum(target, flags):
    """Direct oracle for sum chi4(p)*chi4(q) over unordered prime pairs."""
    target = integer(target, "target", 4)
    if target % 2 or len(flags) <= target:
        raise ValueError("expected an even target covered by the indicator")
    return sum(chi4(p) * chi4(target - p)
               for p in range(2, target // 2 + 1) if flags[p] and flags[target - p])
