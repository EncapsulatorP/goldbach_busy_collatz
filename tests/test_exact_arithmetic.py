import sys
import unittest
from math import gcd, isqrt
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from exact_arithmetic import (MAX_LIMIT, NTT_MODULUS, NTT_ROOT,
    binary_self_convolution, chi4, exact_goldbach_counts, mobius_prefix,
    prime_indicator, signed_prime_pair_sum, unit_pair_count)


def trial_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


class ExactArithmeticTests(unittest.TestCase):
    def test_ntt_field_and_primitive_root(self):
        self.assertTrue(trial_prime(NTT_MODULUS))
        for factor in (2, 7, 17):
            self.assertNotEqual(pow(NTT_ROOT, (NTT_MODULUS - 1) // factor, NTT_MODULUS), 1)
        self.assertLess((NTT_MODULUS - 1)**2, 2**63)

    def test_ntt_against_independent_integer_convolution(self):
        rng = np.random.default_rng(518)
        for length in (1, 2, 3, 7, 8, 9, 31, 32, 33, 127, 129):
            for arr in (np.ones(length, dtype=np.int64),
                        np.zeros(length, dtype=np.int64),
                        rng.integers(0, 2, length)):
                expected = [sum(int(arr[i])*int(arr[k-i])
                                for i in range(length) if 0 <= k-i < length)
                            for k in range(2*length-1)]
                self.assertEqual(binary_self_convolution(arr).tolist(), expected)

    def test_counts_against_trial_division_oracle(self):
        for limit in (0, 1, 2, 3, 4, 5, 8, 30, 31, 127, 257, 1024):
            actual = exact_goldbach_counts(limit)
            for n in range(limit + 1):
                expected = sum(trial_prime(p) and trial_prime(n-p)
                               for p in range(2, n//2+1)) if n % 2 == 0 else 0
                self.assertEqual(int(actual[n]), expected, (limit, n))

    def test_reject_invalid_or_unbounded_inputs(self):
        for value in (-1, 1.5, True, MAX_LIMIT+1):
            with self.assertRaises(ValueError):
                exact_goldbach_counts(value)
        for flags in ([], [0., 1.], [0, 2], [[0, 1]], [np.nan]):
            with self.assertRaises(ValueError):
                binary_self_convolution(flags)
        with self.assertRaises(ValueError):
            exact_goldbach_counts(4, [0, 0, 1])
        with self.assertRaises(ValueError):
            exact_goldbach_counts(4, [0, 1, 1, 1, 0])

    def test_prime_power_crt_against_exhaustive_units(self):
        for modulus in range(2, 121):
            for target in range(modulus):
                expected = sum(gcd(a, modulus) == 1 and gcd(target-a, modulus) == 1
                               for a in range(modulus))
                self.assertEqual(unit_pair_count(target, modulus), expected)

    def test_mobius_dirichlet_inverse_identity(self):
        mu, mertens = mobius_prefix(500)
        for n in range(1, 501):
            self.assertEqual(sum(int(mu[d]) for d in range(1, n+1) if n % d == 0), int(n == 1))
            self.assertEqual(int(mertens[n]), sum(mu[1:n+1]))
        self.assertEqual(mobius_prefix(0)[1].tolist(), [0])

    def test_character_bridge_including_exceptional_prime_two(self):
        flags = np.array([trial_prime(n) for n in range(1001)])
        counts = exact_goldbach_counts(1000, flags)
        for n in range(4, 1001, 2):
            expected = (-1)**(n//2+1)*(int(counts[n])-int(n == 4))
            self.assertEqual(signed_prime_pair_sum(n, flags), expected)
        self.assertEqual([chi4(n) for n in range(8)], [0,1,0,-1,0,1,0,-1])


if __name__ == '__main__':
    unittest.main()
