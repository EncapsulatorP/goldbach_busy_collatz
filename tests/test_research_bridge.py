import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from research_bridge import beta_interval, character_identity, formal_readiness, run


class BridgeTests(unittest.TestCase):
    def test_beta_bounds_nest_and_include_independent_value(self):
        # Reference decimal is an independent low-precision check, not the bound.
        anchor = Fraction('0.915965594177219015')
        previous = (Fraction(0), Fraction(2))
        for terms in (1, 2, 4, 16, 64, 128):
            low, high = beta_interval(terms)
            self.assertLess(low, anchor)
            self.assertGreater(high, anchor)
            self.assertLessEqual(previous[0], low)
            self.assertLessEqual(high, previous[1])
            self.assertEqual(high-low, Fraction(1, (2*terms+1)**2))
            previous = low, high

    def test_finite_channel_split_varied_parities(self):
        for terms in (1, 2, 3, 8, 31):
            for weight in (1, 2, 3, 4, 8):
                result = character_identity(terms, weight)
                self.assertTrue(result['real_matches_minus_eta_over_two_power_s'])
                self.assertTrue(result['imag_matches_beta'])

    def test_fixture_cannot_promote_a_theorem(self):
        manifests = [None, {}, {'production_ready': True},
                     {'corpus_version': 'fixture:1', 'production_ready': True},
                     {'backends': {'kernel_runner': 'fixture_only'}}]
        for manifest in manifests:
            self.assertFalse(formal_readiness(manifest)['formal_proof_verified'])
        self.assertEqual(formal_readiness(manifests[-1])['status'], 'fixture_only')
        for invalid in ([1], {'backends': None}, {'backends': []}):
            with self.assertRaises(ValueError): formal_readiness(invalid)

    def test_local_capsule_is_reproducible_and_explicitly_finite(self):
        a, b = run(200), run(200)
        self.assertEqual(a['goldbach']['counts_sha256'], b['goldbach']['counts_sha256'])
        self.assertTrue(a['computation_checks_passed'])
        self.assertEqual(a['parameters']['oracle_target_count'], 99)
        self.assertEqual(a['tower']['status'], 'not_run')
        self.assertFalse(a['formal_proof_verified'])
        self.assertFalse(a['character_bridge']['independent_predictive_signal'])
        json.dumps(a, allow_nan=False)

    def test_budget_and_domain_are_enforced(self):
        for limit in (0, 5, 100002, True):
            with self.assertRaises(ValueError): run(limit)


if __name__ == '__main__':
    unittest.main()
