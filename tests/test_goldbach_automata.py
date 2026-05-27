import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from goldbach_automata import build_automaton_states, derive_rule_table, search_automaton
from shattering_mirrors import build_dataset


class GoldbachAutomataTests(unittest.TestCase):
    def test_build_automaton_states_has_one_cell_per_window_and_rho(self):
        df = build_dataset(400, include_strings=True)
        states = build_automaton_states(df, step=0.5, clip=2.0, window_size=80, stride=40)

        rho_count = df["rho30"].nunique()
        window_count = states["time_index"].nunique()
        self.assertEqual(len(states), rho_count * window_count)
        self.assertIn("state_label", states.columns)
        self.assertIn("state_share", states.columns)

    def test_rule_table_returns_probability_like_metrics(self):
        df = build_dataset(400, include_strings=True)
        states = build_automaton_states(df, step=0.5, clip=2.0, window_size=80, stride=40)
        rules, accuracy, entropy, dominant_share = derive_rule_table(states)

        self.assertFalse(rules.empty)
        self.assertTrue(0.0 <= accuracy <= 1.0)
        self.assertTrue(0.0 <= entropy <= 1.0)
        self.assertTrue(0.0 <= dominant_share <= 1.0)

    def test_search_automaton_selects_a_best_candidate(self):
        df = build_dataset(400, include_strings=True)
        result = search_automaton(
            df=df,
            steps=[0.25, 0.5],
            clips=[1.5, 2.0],
            window_size=80,
            stride=40,
        )

        self.assertFalse(result.states.empty)
        self.assertFalse(result.rules.empty)
        self.assertFalse(result.search.empty)
        self.assertIn(result.best_step, [0.25, 0.5])
        self.assertIn(result.best_clip, [1.5, 2.0])


if __name__ == "__main__":
    unittest.main()
