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

    @staticmethod
    def result_rows(
        cases: list[dict[str, str]], unresolved: set[str] | None = None
    ) -> list[dict[str, str]]:
        unresolved = unresolved or set()
        return [
            {
                "caseId": row["caseId"],
                "x": row["x"],
                "y": row["y"],
                "id": (
                    "-1"
                    if row["caseId"] in unresolved
                    else str((int(row["caseId"]) - 5000) % 2)
                ),
                "runner_state": "failed" if row["caseId"] in unresolved else "complete",
                "exit_code": "1" if row["caseId"] in unresolved else "0",
            }
            for row in cases
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

    def test_proposal_uses_configured_mix_and_noncolliding_case_ids(self) -> None:
        self.create_layout()
        (self.campaign.root / "campaign-config.json").write_text(
            json.dumps(
                {
                    "allowed_x_values": [1.0, 2.0, 4.0],
                    "case_id_start": 6000,
                    "n_new": 14,
                    "n_repeats": 2,
                    "posterior_samples": 64,
                }
            )
            + "\n"
        )
        self.mark_completed_through(0)
        output = self.campaign.proposals / "Sweep-1_proposed.csv"

        def predictor(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
            building = Path(command[command.index("--outfile") + 1])
            self.write_rows(building, self.valid_proposal())
            return subprocess.CompletedProcess(command, 0, "proposed\n", "")

        with mock.patch.object(CAMPAIGN_MODULE, "run", side_effect=predictor) as run:
            CAMPAIGN_MODULE.proposal_for(self.campaign, 1, output=output)

        command = run.call_args.args[0]
        self.assertEqual(command[command.index("--n-simulations") + 1], "16")
        self.assertEqual(command[command.index("--n-new") + 1], "14")
        self.assertEqual(command[command.index("--n-repeats") + 1], "2")
        self.assertEqual(command[command.index("--posterior-samples") + 1], "64")
        rows = CAMPAIGN_MODULE.read_rows(output)
        self.assertEqual(rows[0]["caseId"], "6000")
        self.assertEqual(rows[-1]["caseId"], "6015")

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

    def test_unattended_config_skips_checkpoint_and_runs_to_twenty(self) -> None:
        self.create_layout()
        (self.campaign.root / "campaign-config.json").write_text(
            json.dumps(
                {
                    "allowed_x_values": [1.0, 2.0, 4.0],
                    "manual_checkpoint_after": None,
                    "final_iteration": 20,
                }
            )
            + "\n"
        )
        self.mark_completed_through(8)
        proposal = self.campaign.proposals / "Sweep-9_proposed.csv"
        with (
            mock.patch.object(
                CAMPAIGN_MODULE, "proposal_for", return_value=proposal
            ) as propose,
            mock.patch.object(CAMPAIGN_MODULE, "submit") as submit,
        ):
            CAMPAIGN_MODULE.advance(self.campaign, submit_jobs=True)
        propose.assert_called_once_with(self.campaign, 9, output=proposal)
        submit.assert_called_once_with(self.campaign, 9, proposal)

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
        self.assertEqual(
            current,
            {
                "job_id": "202",
                "iteration": 1,
                "attempt": 2,
                "backend": "slurm",
            },
        )

    def test_local_submit_uses_bounded_systemd_runner(self) -> None:
        self.create_layout()
        (self.campaign.root / "campaign-config.json").write_text(
            json.dumps(
                {
                    "allowed_x_values": [1.0, 2.0, 4.0],
                    "unit_prefix": "dropinj-l11",
                    "max_level": 11,
                    "drop_radius_min": 0.015625,
                    "workers": 3,
                    "threads_per_case": 8,
                    "max_threads": 48,
                }
            )
            + "\n"
        )
        campaign = CAMPAIGN_MODULE.Campaign(
            root=self.campaign.root,
            project_root=self.campaign.project_root,
            predictor_root=self.campaign.predictor_root,
            backend="local",
        )
        proposal = campaign.proposals / "Sweep-1_proposed.csv"
        self.write_rows(proposal, self.valid_proposal())
        launched = subprocess.CompletedProcess([], 0, "Running as unit\n", "")

        with mock.patch.object(CAMPAIGN_MODULE, "run", return_value=launched) as run:
            job_id = CAMPAIGN_MODULE.submit(campaign, 1, proposal)

        self.assertEqual(job_id, "local:dropinj-l11-i01-a01.service")
        command = run.call_args.args[0]
        self.assertEqual(command[:4], ["systemd-run", "--user", "--unit", "dropinj-l11-i01-a01"])
        self.assertIn("Nice=10", command)
        self.assertIn("CONTOUR_MAXLEVEL=11", command)
        self.assertIn("CONTOUR_DROP_RADIUS_MIN=0.015625", command)
        self.assertIn("CONTOUR_WORKERS=3", command)
        self.assertIn("CONTOUR_THREADS_PER_CASE=8", command)
        self.assertIn("CONTOUR_MAX_THREADS=48", command)
        self.assertIn(str(campaign.project_root / "runContourLocal.sh"), command)
        current = json.loads(
            (campaign.iterations / "iteration-01" / "current-attempt.json").read_text()
        )
        self.assertEqual(current["backend"], "local")

    def test_local_state_maps_systemd_results(self) -> None:
        states = (
            ("LoadState=loaded\nActiveState=active\nResult=success\n", "RUNNING"),
            ("LoadState=loaded\nActiveState=inactive\nResult=success\n", "COMPLETED"),
            ("LoadState=loaded\nActiveState=failed\nResult=exit-code\n", "FAILED"),
            ("LoadState=not-found\nActiveState=inactive\nResult=success\n", "UNKNOWN"),
        )
        for stdout, expected in states:
            with self.subTest(expected=expected), mock.patch.object(
                CAMPAIGN_MODULE.subprocess,
                "run",
                return_value=subprocess.CompletedProcess([], 0, stdout, ""),
            ):
                self.assertEqual(
                    CAMPAIGN_MODULE.local_state("local:dropinj-i01-a01.service"),
                    expected,
                )

    def test_inline_submit_runs_with_scoped_numerical_environment(self) -> None:
        self.create_layout()
        (self.campaign.root / "campaign-config.json").write_text(
            json.dumps(
                {
                    "allowed_x_values": [1.0, 2.0, 4.0],
                    "unit_prefix": "dropinj-l11",
                    "max_level": 11,
                    "drop_radius_min": 0.015625,
                    "workers": 3,
                    "threads_per_case": 8,
                    "max_threads": 48,
                }
            )
            + "\n"
        )
        campaign = CAMPAIGN_MODULE.Campaign(
            root=self.campaign.root,
            project_root=self.campaign.project_root,
            predictor_root=self.campaign.predictor_root,
            backend="inline",
        )
        proposal = campaign.proposals / "Sweep-1_proposed.csv"
        self.write_rows(proposal, self.valid_proposal())
        completed = subprocess.CompletedProcess([], 0, "resolved=16/16\n", "")

        with mock.patch.object(
            CAMPAIGN_MODULE.subprocess, "run", return_value=completed
        ) as run:
            job_id = CAMPAIGN_MODULE.submit(campaign, 1, proposal)

        self.assertEqual(job_id, "inline:dropinj-l11-i01-a01")
        environment = run.call_args.kwargs["env"]
        self.assertEqual(environment["CONTOUR_MAXLEVEL"], "11")
        self.assertEqual(environment["CONTOUR_DROP_RADIUS_MIN"], "0.015625")
        self.assertEqual(environment["CONTOUR_WORKERS"], "3")
        self.assertEqual(environment["CONTOUR_THREADS_PER_CASE"], "8")
        self.assertEqual(environment["CONTOUR_MAX_THREADS"], "48")
        self.assertEqual(CAMPAIGN_MODULE.execution_state(job_id), "COMPLETED")

    def test_retry_submits_only_cases_unresolved_across_prior_attempts(self) -> None:
        self.create_layout()
        proposal = self.campaign.proposals / "Sweep-1_proposed.csv"
        cases = self.valid_proposal()
        self.write_rows(proposal, cases)
        with mock.patch.object(
            CAMPAIGN_MODULE,
            "run",
            return_value=subprocess.CompletedProcess(
                [], 0, "Submitted batch job 101\n", ""
            ),
        ):
            CAMPAIGN_MODULE.submit(self.campaign, 1, proposal)

        unresolved = {"5000", "5001", "5012"}
        attempt_one = self.campaign.iterations / "iteration-01" / "attempt-01"
        self.write_rows(attempt_one / "results.csv", self.result_rows(cases, unresolved))

        with (
            mock.patch.object(CAMPAIGN_MODULE, "slurm_state", return_value="COMPLETED"),
            mock.patch.object(
                CAMPAIGN_MODULE,
                "run",
                return_value=subprocess.CompletedProcess(
                    [], 0, "Submitted batch job 202\n", ""
                ),
            ),
        ):
            self.assertEqual(
                CAMPAIGN_MODULE.retry_current_iteration(
                    self.campaign, submit_job=True
                ),
                "202",
            )

        retry_rows = CAMPAIGN_MODULE.read_rows(
            self.campaign.iterations
            / "iteration-01"
            / "attempt-02"
            / "cases.csv"
        )
        self.assertEqual(
            [row["caseId"] for row in retry_rows],
            ["5000", "5001", "5012"],
        )

    def test_collect_merges_attempts_in_canonical_order(self) -> None:
        self.create_layout()
        proposal = self.campaign.proposals / "Sweep-1_proposed.csv"
        cases = self.valid_proposal()
        self.write_rows(proposal, cases)
        iteration = CAMPAIGN_MODULE.stage_proposal(self.campaign, 1, proposal)
        unresolved = {"5000", "5001", "5012"}

        attempt_one = iteration / "attempt-01"
        attempt_one.mkdir()
        self.write_rows(attempt_one / "cases.csv", cases)
        self.write_rows(attempt_one / "results.csv", self.result_rows(cases, unresolved))

        retry_cases = [row for row in cases if row["caseId"] in unresolved]
        attempt_two = iteration / "attempt-02"
        attempt_two.mkdir()
        self.write_rows(attempt_two / "cases.csv", retry_cases)
        self.write_rows(attempt_two / "results.csv", self.result_rows(retry_cases))

        with mock.patch.object(CAMPAIGN_MODULE, "assess") as assess:
            destination = CAMPAIGN_MODULE.collect(self.campaign, 1, 2)

        completed = CAMPAIGN_MODULE.read_rows(destination)
        self.assertEqual(
            [row["caseId"] for row in completed],
            [row["caseId"] for row in cases],
        )
        self.assertEqual(
            [row["id"] for row in completed],
            [str(index % 2) for index in range(CAMPAIGN_MODULE.BATCH_SIZE)],
        )
        assess.assert_called_once()
        measurements = CAMPAIGN_MODULE.read_rows(
            self.campaign.measurements / "Sweep-1_measurements.csv"
        )
        self.assertEqual(len(measurements), CAMPAIGN_MODULE.BATCH_SIZE)
        self.assertTrue(all(row["iteration"] == "1" for row in measurements))
        self.assertEqual(measurements[0]["source_attempt"], "2")
        self.assertEqual(measurements[3]["source_attempt"], "1")

    def test_quality_fail_and_exact_quarantine_are_excluded_from_merge(self) -> None:
        self.create_layout()
        proposal = self.campaign.proposals / "Sweep-1_proposed.csv"
        cases = self.valid_proposal()
        self.write_rows(proposal, cases)
        iteration = CAMPAIGN_MODULE.stage_proposal(self.campaign, 1, proposal)

        attempt_one = iteration / "attempt-01"
        attempt_one.mkdir()
        self.write_rows(attempt_one / "cases.csv", cases)
        results = self.result_rows(cases)
        results[0]["quality_state"] = "fail"
        results[0]["quality_reason"] = "max_ke_exceeds_1000"
        self.write_rows(attempt_one / "results.csv", results)
        self.write_rows(
            self.campaign.root / "quality-quarantine.csv",
            [
                {
                    "iteration": "1",
                    "attempt": "1",
                    "caseId": "5001",
                    "reason": "manual facet review",
                }
            ],
        )

        self.assertEqual(
            [
                row["caseId"]
                for row in CAMPAIGN_MODULE.unresolved_cases(self.campaign, 1, 1)
            ],
            ["5000", "5001"],
        )

        attempt_two = iteration / "attempt-02"
        attempt_two.mkdir()
        retry_case = [cases[1]]
        self.write_rows(attempt_two / "cases.csv", retry_case)
        self.write_rows(attempt_two / "results.csv", self.result_rows(retry_case))
        self.assertEqual(
            [
                row["caseId"]
                for row in CAMPAIGN_MODULE.unresolved_cases(self.campaign, 1, 2)
            ],
            ["5000"],
        )

    def test_qualityless_attempt_is_backfilled_from_retained_case_evidence(self) -> None:
        self.create_layout()
        proposal = self.campaign.proposals / "Sweep-1_proposed.csv"
        cases = self.valid_proposal()
        self.write_rows(proposal, cases)
        iteration = CAMPAIGN_MODULE.stage_proposal(self.campaign, 1, proposal)

        attempt = iteration / "attempt-01"
        attempt.mkdir()
        self.write_rows(attempt / "cases.csv", cases)
        self.write_rows(attempt / "results.csv", self.result_rows(cases))
        for row in cases:
            case_dir = attempt / "cases" / f"case-{row['caseId']}"
            case_dir.mkdir(parents=True)
            maximum = "1000.01" if row["caseId"] == "5000" else "1"
            (case_dir / "log").write_text(f"1 0.001 0.1 {maximum} 0 0\n")
            (case_dir / "interface-latest.dat").write_text("0 1\n")

        resolved = CAMPAIGN_MODULE.merged_resolved_results(self.campaign, 1, 1)
        self.assertNotIn("5000", resolved)
        self.assertEqual(len(resolved), CAMPAIGN_MODULE.BATCH_SIZE - 1)


if __name__ == "__main__":
    unittest.main()
