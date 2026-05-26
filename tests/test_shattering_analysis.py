import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from shattering_compressed import pair_chambers_i
from shattering_mirrors import (
    build_dataset,
    exact_goldbach_counts_fft,
    sieve_bool,
    summarize_clusters,
)


def brute_goldbach_count(limit: int, is_prime):
    counts = [0] * (limit + 1)
    for n in range(4, limit + 1, 2):
        total = 0
        for p in range(2, n // 2 + 1):
            q = n - p
            if is_prime[p] and is_prime[q]:
                total += 1
        counts[n] = total
    return counts


class ShatteringAnalysisTests(unittest.TestCase):
    def test_exact_goldbach_counts_fft_matches_bruteforce(self):
        limit = 200
        is_prime = sieve_bool(limit)
        fft_counts = exact_goldbach_counts_fft(limit, is_prime)
        brute_counts = brute_goldbach_count(limit, is_prime)

        for n in range(4, limit + 1, 2):
            self.assertEqual(fft_counts[n], brute_counts[n], msg=f"Mismatch at N={n}")

    def test_dataset_uses_normalized_cluster_axis(self):
        df = build_dataset(2_000, include_strings=True)

        self.assertIn("z_h_raw", df.columns)
        self.assertIn("h_cal", df.columns)
        self.assertIn("z_bucket", df.columns)
        self.assertTrue(df["native_cluster"].str.startswith("z=").all())
        self.assertAlmostEqual(float(df["z_h"].mean()), 0.0, places=10)

    def test_summary_groups_by_z_bucket(self):
        df = build_dataset(2_000, include_strings=True)
        summary = summarize_clusters(df, top=10)

        self.assertIn("z_bucket", summary.columns)
        self.assertNotIn("eps_bucket", summary.columns)
        self.assertTrue(summary["z_bucket"].astype(str).str.len().gt(0).all())

    def test_pair_chambers_preserve_mod_30_constraint(self):
        is_prime = sieve_bool(500)
        for n in (10, 28, 100, 250):
            chambers = pair_chambers_i(n, is_prime)
            self.assertTrue(chambers, msg=f"No Goldbach pairs found for N={n}")
            self.assertTrue(all(item["chamber_check"] for item in chambers), msg=f"Bad chamber for N={n}")


if __name__ == "__main__":
    unittest.main()
