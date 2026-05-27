import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from goldbach_volume import (
    build_volume_dataset,
    family_label,
    pair_signature_counts,
)


class GoldbachVolumeTests(unittest.TestCase):
    def test_family_label_matches_expected_support_patterns(self):
        self.assertEqual(family_label([]), "spine_none")
        self.assertEqual(family_label([3]), "cone_3")
        self.assertEqual(family_label([5]), "ring_5")
        self.assertEqual(family_label([3, 5]), "mobius_3x5")
        self.assertEqual(family_label([7]), "shell_other")

    def test_mod30_dataset_matches_current_volume_schema(self):
        df, _, _ = build_volume_dataset(100, 30)
        self.assertEqual(
            df.columns.tolist(),
            [
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
            ],
        )
        self.assertEqual(df.iloc[0]["family"], "spine_none")
        self.assertEqual(int(df.iloc[1]["rho30"]), 6)

    def test_pair_signature_counts_exist_for_small_range(self):
        df, is_prime, prime_list = build_volume_dataset(60, 30)
        pair_df = pair_signature_counts(60, 30, prime_list, is_prime)
        self.assertFalse(pair_df.empty)
        self.assertIn("pair_i", pair_df.columns)
        self.assertIn("count", pair_df.columns)
        self.assertEqual(int(pair_df["count"].sum()), int(df["r"].sum()))


if __name__ == "__main__":
    unittest.main()
