import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from goldbach_native_filter import (
    build_native_filter_frame,
    goldbach_heuristic,
    local_boost,
    native_filter,
    resolve_threshold,
    tier,
)


class GoldbachNativeFilterTests(unittest.TestCase):
    def test_local_boost_matches_expected_small_cases(self):
        self.assertEqual(local_boost(4), 1.0)
        self.assertEqual(local_boost(3), 2.0)
        self.assertAlmostEqual(local_boost(15), 8.0 / 3.0)

    def test_tier_small_cases(self):
        self.assertEqual(tier(4), 0)
        self.assertEqual(tier(3), 1)
        self.assertEqual(tier(21), 2)
        self.assertEqual(tier(15), 3)
        self.assertEqual(tier(105), 4)

    def test_threshold_resolution_accepts_callable(self):
        threshold = resolve_threshold(2000, lambda bound: 2.0)
        self.assertEqual(threshold, 2.0)

    def test_filter_frame_has_expected_columns(self):
        df, threshold = build_native_filter_frame(200)
        self.assertGreater(threshold, 0.0)
        self.assertEqual(
            df.columns.tolist(),
            ["n", "N", "boost_G", "boost_primes", "h", "tier", "keep"],
        )
        self.assertTrue(df["keep"].isin([True, False]).all())

    def test_native_filter_uses_scalar_threshold(self):
        h = goldbach_heuristic(15)
        self.assertTrue(native_filter(15, threshold=h - 1e-9))
        self.assertFalse(native_filter(15, threshold=h + 1e-9))


if __name__ == "__main__":
    unittest.main()
