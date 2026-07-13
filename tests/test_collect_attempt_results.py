from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "contourWorkflow" / "collect_attempt_results.py"
SPEC = importlib.util.spec_from_file_location("collect_attempt_results_under_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class CollectAttemptResultsTests(unittest.TestCase):
    def test_collects_retained_complete_and_interrupted_cases(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = []
            for index in range(2):
                case_dir = root / "cases" / f"case-{7000 + index}"
                case_dir.mkdir(parents=True)
                (case_dir / "classification.status").write_text(
                    "id,reason,t,drop_volume,drop_radius,drop_axial_position,"
                    "drop_axial_velocity,threshold,consecutive\n"
                    + (
                        "0,observation_horizon_without_drop,1,0,0,0,0,0.015625,0\n"
                        if index == 0
                        else "-1,running,0.6,0,0,0,0,0.015625,0\n"
                    )
                )
                (case_dir / "runner.status").write_text(
                    "state=complete\nexit_code=0\n"
                    if index == 0
                    else "state=running\n"
                )
                (case_dir / "log").write_text("1 0.001 0.1 1 0 0\n")
                (case_dir / "interface-latest.dat").write_text("0 1\n")
                manifest.append(
                    {
                        "CaseNo": str(7000 + index),
                        "Rr": "6",
                        "OhOut": f"{0.04 + index * 0.001:g}",
                        "case_dir": str(case_dir),
                    }
                )
            (root / "case-manifest.json").write_text(
                json.dumps(manifest) + "\n"
            )

            destination = MODULE.collect(root)
            with destination.open(newline="") as stream:
                rows = list(csv.DictReader(stream))

            self.assertEqual([row["caseId"] for row in rows], ["7000", "7001"])
            self.assertEqual(rows[0]["runner_state"], "complete")
            self.assertEqual(rows[0]["id"], "0")
            self.assertEqual(rows[1]["runner_state"], "running")
            self.assertEqual(rows[1]["id"], "-1")


if __name__ == "__main__":
    unittest.main()
