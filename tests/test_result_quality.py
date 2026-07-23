from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "contourWorkflow" / "result_quality.py"
SPEC = importlib.util.spec_from_file_location("result_quality_under_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
QUALITY_MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = QUALITY_MODULE
SPEC.loader.exec_module(QUALITY_MODULE)


class ResultQualityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.case = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_passes_finite_evidence_within_conservative_gates(self) -> None:
        (self.case / "log").write_text(
            "i dt t ke Xc Vcm\n"
            "0 0.001 0 12.5 0 0\n"
            "1 0.001 0.001 999 0 0\n"
        )
        (self.case / "interface-latest.dat").write_text("0 1\n\n2 3\n")

        self.assertEqual(
            QUALITY_MODULE.case_quality(self.case),
            {
                "max_ke": "999",
                "facet_lines": "2",
                "quality_state": "pass",
                "quality_reason": "within_gates",
            },
        )

    def test_fails_large_ke_and_facet_count(self) -> None:
        (self.case / "log").write_text("1 0.001 0.1 1000.01 0 0\n")
        (self.case / "interface-latest.dat").write_text("0 1\n" * 8001)

        quality = QUALITY_MODULE.case_quality(self.case)
        self.assertEqual(quality["quality_state"], "fail")
        self.assertEqual(quality["max_ke"], "1000.01")
        self.assertEqual(quality["facet_lines"], "8001")
        self.assertEqual(
            quality["quality_reason"],
            "max_ke_exceeds_1000;facet_lines_exceeds_8000",
        )

    def test_fails_missing_or_nonfinite_evidence(self) -> None:
        missing = QUALITY_MODULE.case_quality(self.case)
        self.assertEqual(missing["quality_state"], "fail")
        self.assertEqual(missing["quality_reason"], "missing_log;missing_facets")

        (self.case / "log").write_text("1 0.001 0.1 nan 0 0\n")
        (self.case / "interface-latest.dat").write_text("nan 1\n")
        nonfinite = QUALITY_MODULE.case_quality(self.case)
        self.assertEqual(nonfinite["quality_state"], "fail")
        self.assertEqual(
            nonfinite["quality_reason"],
            "nonfinite_ke;missing_ke;nonfinite_facets",
        )


if __name__ == "__main__":
    unittest.main()
