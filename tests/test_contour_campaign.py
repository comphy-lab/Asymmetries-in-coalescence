from __future__ import annotations

import copy
import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = Path(__file__).parents[1] / "contourWorkflow" / "contour_campaign.py"
SPEC = importlib.util.spec_from_file_location("contour_campaign_under_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
CAMPAIGN_MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CAMPAIGN_MODULE
SPEC.loader.exec_module(CAMPAIGN_MODULE)


class ContourCampaignTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.workspace = Path(self.temporary.name)
        self.project_root = self.workspace / "project"
        self.predictor_root = self.workspace / "predictor"
        data_dir = self.project_root / "simulationCases" / "DataFiles"
        data_dir.mkdir(parents=True)
        self.predictor_root.mkdir()
        for value in (1.0, 2.0, 4.0):
            (data_dir / f"InitialConditionRr-{value:.2f}.dat").write_text("shape\n")
        self.campaign = CAMPAIGN_MODULE.Campaign(
            root=self.workspace / "campaign",
            project_root=self.project_root,
            predictor_root=self.predictor_root,
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def create_layout(self) -> None:
        for directory in (
            self.campaign.root,
            self.campaign.completed,
            self.campaign.proposals,
            self.campaign.contours,
            self.campaign.iterations,
        ):
            directory.mkdir(parents=True, exist_ok=True)
        (self.campaign.root / "campaign-config.json").write_text(
            json.dumps({"allowed_x_values": [1.0, 2.0, 4.0]}) + "\n"
        )

    def valid_proposal(self) -> list[dict[str, str]]:
        return [
            {
                "caseId": str(5000 + index),
                "x": str((1.0, 2.0, 4.0)[index % 3]),
                "y": f"{0.011 + 0.003 * index:.3f}",
                "id": "-1",
            }
            for index in range(CAMPAIGN_MODULE.BATCH_SIZE)
        ]

    def mark_completed_through(self, iteration: int) -> None:
        self.campaign.completed.mkdir(parents=True, exist_ok=True)
        for value in range(iteration + 1):
            (self.campaign.completed / f"Sweep-{value}_completed.csv").write_text(
                "caseId,x,y,id\n"
            )

    def test_initialise_is_idempotent_and_rejects_config_conflict(self) -> None:
        seed = self.workspace / "seed.csv"
        self.write_rows(
            seed,
            [
                {"caseId": "10", "Rr": "1", "Oh": "0.02", "id": "0"},
                {"caseId": "11", "Rr": "4", "Oh": "0.03", "id": "1"},
            ],
        )
        predictor_help = subprocess.CompletedProcess(
            [], 0, stdout="usage: predictor --x-candidates VALUES\n", stderr=""
        )
        predictor_commit = subprocess.CompletedProcess(
            [], 0, stdout="abc123\n", stderr=""
        )

        with (
            mock.patch.object(CAMPAIGN_MODULE, "run", return_value=predictor_help),
            mock.patch.object(
                CAMPAIGN_MODULE.subprocess, "run", return_value=predictor_commit
            ),
        ):
            CAMPAIGN_MODULE.initialise(self.campaign, seed, [])
            marker = self.campaign.iterations / "runtime-evidence.txt"
            marker.write_text("preserve me\n")
            original = {
                path.relative_to(self.campaign.root): path.read_bytes()
                for path in (
                    self.campaign.root / "campaign-config.json",
                    self.campaign.root / "campaign-state.json",
                    self.campaign.completed / "Sweep-0_completed.csv",
                )
            }

            CAMPAIGN_MODULE.initialise(self.campaign, seed, [])
            self.assertEqual(marker.read_text(), "preserve me\n")
            for relative, contents in original.items():
                self.assertEqual((self.campaign.root / relative).read_bytes(), contents)

            with self.assertRaisesRegex(FileExistsError, "config differs"):
                CAMPAIGN_MODULE.initialise(self.campaign, seed, [4.0])
            self.assertEqual(marker.read_text(), "preserve me\n")
            for relative, contents in original.items():
                self.assertEqual((self.campaign.root / relative).read_bytes(), contents)

    def test_validate_proposal_rejects_duplicate_and_out_of_domain_rows(self) -> None:
        self.create_layout()
        rows = self.valid_proposal()
        CAMPAIGN_MODULE.validate_proposal(self.campaign, rows)

        failures = []
        duplicate_id = copy.deepcopy(rows)
        duplicate_id[-1]["caseId"] = duplicate_id[0]["caseId"]
        failures.append((duplicate_id, "duplicate caseId"))
        duplicate_point = copy.deepcopy(rows)
        duplicate_point[-1]["x"] = duplicate_point[0]["x"]
        duplicate_point[-1]["y"] = duplicate_point[0]["y"]
        failures.append((duplicate_point, "duplicate proposal point"))
        unavailable_x = copy.deepcopy(rows)
        unavailable_x[-1]["x"] = "3"
        failures.append((unavailable_x, "unavailable Rr"))
        too_small_y = copy.deepcopy(rows)
        too_small_y[-1]["y"] = "0.009"
        failures.append((too_small_y, "lies outside"))
        too_large_y = copy.deepcopy(rows)
        too_large_y[-1]["y"] = "0.076"
        failures.append((too_large_y, "lies outside"))

        for invalid, message in failures:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    CAMPAIGN_MODULE.validate_proposal(self.campaign, invalid)

    def test_manual_candidate_and_advance_require_iteration_eight(self) -> None:
        self.create_layout()
        self.mark_completed_through(7)
        with self.assertRaisesRegex(ValueError, "only after iteration 8"):
            CAMPAIGN_MODULE.propose_manual_batch(self.campaign)

        self.mark_completed_through(8)
        candidate = self.campaign.proposals / "Sweep-9_manual-candidate.csv"
        with mock.patch.object(
            CAMPAIGN_MODULE, "proposal_for", return_value=candidate
        ) as propose:
            self.assertEqual(
                CAMPAIGN_MODULE.propose_manual_batch(self.campaign), candidate
            )
            propose.assert_called_once_with(self.campaign, 9, output=candidate)

        with (
            mock.patch.object(CAMPAIGN_MODULE, "proposal_for") as propose,
            mock.patch.object(CAMPAIGN_MODULE, "submit") as submit,
        ):
            CAMPAIGN_MODULE.advance(self.campaign, submit_jobs=True)
            propose.assert_not_called()
            submit.assert_not_called()
        state = json.loads((self.campaign.root / "campaign-state.json").read_text())
        self.assertEqual(state["state"], "manual_checkpoint")
        self.assertEqual(state["last_completed_iteration"], 8)

    def test_advance_finishing_iteration_sixteen_stops_without_new_work(self) -> None:
        self.create_layout()
        self.mark_completed_through(15)
        iteration_root = self.campaign.iterations / "iteration-16"
        iteration_root.mkdir()
        (iteration_root / "current-attempt.json").write_text(
            json.dumps({"job_id": "1600", "iteration": 16, "attempt": 1}) + "\n"
        )

        completed = self.campaign.completed / "Sweep-16_completed.csv"

        def collect_iteration(*_args: object, **_kwargs: object) -> Path:
            completed.write_text("caseId,x,y,id\n")
            return completed

        with (
            mock.patch.object(CAMPAIGN_MODULE, "slurm_state", return_value="COMPLETED"),
            mock.patch.object(
                CAMPAIGN_MODULE, "collect", side_effect=collect_iteration
            ) as collect,
            mock.patch.object(CAMPAIGN_MODULE, "proposal_for") as propose,
            mock.patch.object(CAMPAIGN_MODULE, "submit") as submit,
        ):
            CAMPAIGN_MODULE.advance(self.campaign, submit_jobs=True)

        collect.assert_called_once_with(self.campaign, 16, 1)
        propose.assert_not_called()
        submit.assert_not_called()
        state = json.loads((self.campaign.root / "campaign-state.json").read_text())
        self.assertEqual(state, {"state": "complete", "last_completed_iteration": 16})

    def test_manual_approval_is_gated_canonical_and_non_destructive(self) -> None:
        self.create_layout()
        source = self.workspace / "reviewed.csv"
        rows = [
            {
                "Rr": str((1.0, 2.0, 4.0)[index % 3]),
                "Oh": f"{0.012 + 0.003 * index:.3f}",
            }
            for index in range(CAMPAIGN_MODULE.BATCH_SIZE)
        ]
        self.write_rows(source, rows)

        with self.assertRaisesRegex(ValueError, "only after iteration 8"):
            CAMPAIGN_MODULE.approve_manual_batch(self.campaign, source)

        self.mark_completed_through(8)
        CAMPAIGN_MODULE.approve_manual_batch(self.campaign, source)
        approved = self.campaign.proposals / "Sweep-9_manual-approved.csv"
        approved_rows = CAMPAIGN_MODULE.read_rows(approved)
        self.assertEqual(
            [row["caseId"] for row in approved_rows],
            [str(5128 + index) for index in range(CAMPAIGN_MODULE.BATCH_SIZE)],
        )
        self.assertTrue(all(row["id"] == "-1" for row in approved_rows))

        original = approved.read_bytes()
        CAMPAIGN_MODULE.approve_manual_batch(self.campaign, source)
        self.assertEqual(approved.read_bytes(), original)

        changed = copy.deepcopy(rows)
        changed[-1]["Oh"] = "0.074"
        self.write_rows(source, changed)
        with self.assertRaises(FileExistsError):
            CAMPAIGN_MODULE.approve_manual_batch(self.campaign, source)
        self.assertEqual(approved.read_bytes(), original)

    def test_submit_reuses_current_job_and_numbers_isolated_attempts(self) -> None:
        self.create_layout()
        proposal = self.campaign.proposals / "Sweep-1_proposed.csv"
        self.write_rows(proposal, self.valid_proposal())
        submitted = [
            subprocess.CompletedProcess([], 0, "Submitted batch job 101\n", ""),
            subprocess.CompletedProcess([], 0, "Submitted batch job 202\n", ""),
        ]

        with mock.patch.object(
            CAMPAIGN_MODULE, "run", side_effect=submitted
        ) as submit_command:
            self.assertEqual(CAMPAIGN_MODULE.submit(self.campaign, 1, proposal), "101")
            attempt_one = self.campaign.iterations / "iteration-01" / "attempt-01"
            (attempt_one / "runtime-evidence.txt").write_text("first attempt\n")

            self.assertEqual(CAMPAIGN_MODULE.submit(self.campaign, 1, proposal), "101")
            self.assertEqual(submit_command.call_count, 1)
            self.assertFalse(
                (self.campaign.iterations / "iteration-01" / "attempt-02").exists()
            )

            self.assertEqual(
                CAMPAIGN_MODULE.submit(
                    self.campaign, 1, proposal, new_attempt=True
                ),
                "202",
            )

        attempt_two = self.campaign.iterations / "iteration-01" / "attempt-02"
        self.assertEqual((attempt_one / "runtime-evidence.txt").read_text(), "first attempt\n")
        self.assertEqual(
            (attempt_one / "cases.csv").read_bytes(),
            (attempt_two / "cases.csv").read_bytes(),
        )
        current = json.loads(
            (
                self.campaign.iterations
                / "iteration-01"
                / "current-attempt.json"
            ).read_text()
        )
        self.assertEqual(current, {"job_id": "202", "iteration": 1, "attempt": 2})


if __name__ == "__main__":
    unittest.main()
