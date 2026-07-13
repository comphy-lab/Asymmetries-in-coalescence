from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "contourWorkflow" / "materialize_cases.py"


class MaterializeCasesTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.data_dir = self.root / "DataFiles"
        self.data_dir.mkdir()
        (self.data_dir / "InitialConditionRr-4.00.dat").write_text("shape\n")
        self.source = self.root / "coalescenceBubble.c"
        self.source.write_text("int main(void) { return 0; }\n")
        self.proposal = self.root / "cases.csv"
        self.case_root = self.root / "iteration" / "cases"
        self.rows = [
            {
                "caseId": str(5000 + index),
                "x": "4",
                "y": f"{0.01 + 0.004 * index:.3f}",
                "id": "-1",
            }
            for index in range(16)
        ]

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_materializer(
        self, rows: list[dict[str, str]], *, expected: int = 16
    ) -> subprocess.CompletedProcess[str]:
        with self.proposal.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("caseId", "x", "y", "id"))
            writer.writeheader()
            writer.writerows(rows)
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                str(self.proposal),
                str(self.case_root),
                "--data-dir",
                str(self.data_dir),
                "--source",
                str(self.source),
                "--expected",
                str(expected),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_validates_entire_batch_before_writing(self) -> None:
        rows = deepcopy(self.rows)
        rows[-1]["x"] = "4.004"
        result = self.run_materializer(rows)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not exactly match", result.stderr)
        self.assertFalse(self.case_root.exists())

    def test_rejects_invalid_oh_and_duplicate_coordinates(self) -> None:
        for change, expected in (
            ((15, "y", "nan"), "must be finite"),
            ((15, "y", "0.009"), "must be finite"),
            ((15, "caseId", "5000"), "duplicate caseId"),
            ((15, "y", self.rows[0]["y"]), "duplicate proposed point"),
        ):
            with self.subTest(change=change):
                rows = deepcopy(self.rows)
                index, key, value = change
                rows[index][key] = value
                result = self.run_materializer(rows)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected, result.stderr)
                self.assertFalse(self.case_root.exists())

    def test_reuses_byte_identical_cases_with_runtime_outputs(self) -> None:
        first = self.run_materializer(self.rows)
        self.assertEqual(first.returncode, 0, first.stderr)
        case = self.case_root / "case-5000"
        source_before = (case / "coalescenceBubble.c").read_bytes()
        params_before = (case / "case.params").read_bytes()
        (case / "restart").write_bytes(b"checkpoint")
        (case / "runner.status").write_text("state=running\n")
        (case / "intermediate" / "snapshot-0.5").write_bytes(b"snapshot")

        second = self.run_materializer(self.rows)
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("materialised=0 reused=16", second.stdout)
        self.assertEqual((case / "coalescenceBubble.c").read_bytes(), source_before)
        self.assertEqual((case / "case.params").read_bytes(), params_before)
        self.assertEqual((case / "restart").read_bytes(), b"checkpoint")

    def test_rejects_changed_immutable_inputs(self) -> None:
        first = self.run_materializer(self.rows)
        self.assertEqual(first.returncode, 0, first.stderr)
        changed = self.case_root / "case-5007" / "case.params"
        changed.write_text(changed.read_text() + "tampered=1\n")

        second = self.run_materializer(self.rows)
        self.assertNotEqual(second.returncode, 0)
        self.assertIn("different parameters", second.stderr)

    def test_materialises_selective_retry_subset(self) -> None:
        rows = [self.rows[index] for index in (0, 1, 12)]
        result = self.run_materializer(rows, expected=3)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("materialised=3 reused=0", result.stdout)
        manifest = self.case_root.parent / "case-manifest.json"
        self.assertTrue(manifest.exists())
        self.assertEqual(
            sorted(path.name for path in self.case_root.iterdir()),
            ["case-5000", "case-5001", "case-5012"],
        )


if __name__ == "__main__":
    unittest.main()
