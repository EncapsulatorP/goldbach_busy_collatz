import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from busy_beaver_waring_goldbach import (
    build_cell_records,
    build_row_records,
    build_transition_plaquettes,
    builtin_machines,
    history_span,
    simulate_machine,
)


class BusyBeaverWaringGoldbachTests(unittest.TestCase):
    def test_bb2_reaches_known_small_busy_beaver_totals(self):
        machine = builtin_machines()["bb2"]
        history = simulate_machine(machine, max_steps=32)

        self.assertTrue(history[-1].halted)
        self.assertEqual(len(history) - 1, 6)
        self.assertEqual(sum(history[-1].tape.values()), 4)

    def test_encoding_uses_unique_primes_for_each_cell(self):
        machine = builtin_machines()["bb3"]
        history = simulate_machine(machine, max_steps=64)
        span = history_span(history, padding=2)
        cell_records, _ = build_cell_records(history, machine, span)

        primes = [record["prime"] for record in cell_records]
        self.assertEqual(len(primes), len(set(primes)))
        self.assertTrue(all(record["exponent"] >= 2 for record in cell_records))

    def test_transition_plaquettes_match_executed_steps(self):
        machine = builtin_machines()["bb3"]
        history = simulate_machine(machine, max_steps=64)
        span = history_span(history, padding=2)
        _, lookup = build_cell_records(history, machine, span)
        row_records = build_row_records(history, span, lookup)
        tiles = build_transition_plaquettes(history, machine, lookup)

        self.assertEqual(len(tiles), len(history) - 1)
        self.assertEqual(len(row_records), len(history))


if __name__ == "__main__":
    unittest.main()
