import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


class OutputPathConventionTests(unittest.TestCase):
    def test_scripts_do_not_hardcode_external_output_roots(self):
        banned = re.compile(r"['\"](?:/mnt/|[A-Za-z]:\\\\)[^'\"]*outputs[^'\"]*['\"]")

        for path in SCRIPTS.glob("*.py"):
            source = path.read_text(encoding="utf-8")
            self.assertIsNone(
                banned.search(source),
                msg=f"{path.name} hardcodes an external output path",
            )


if __name__ == "__main__":
    unittest.main()
