"""Invariants of the two R_r -> infinity launchers.

These are the properties whose violation has already cost this project
measurements, so they are asserted rather than reviewed:

* the drop detector is pinned to a physical radius, never left at 0 (which
  makes the solver derive `2*Ldomain/2^MAXlevel` and ties the criterion to the
  mesh, so runs at different resolutions stop being comparable);
* the gas-to-liquid viscosity ratio is passed explicitly, because the
  drill-solver drop-map ladder ran 2e-2 against the manuscript's 1e-2 and
  nothing in its argv or run directory recorded the difference;
* the discriminator's parallel arrays stay aligned, since a silent
  off-by-one there would attribute a divergence to the wrong solver stack.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


ROOT = Path(__file__).parents[1]
LADDER = ROOT / "runBurstingBubbleInfiniteRr.sbatch"
PROBE = ROOT / "runInfiniteRrStackProbe.sbatch"
SOLVER_ARGS = ROOT / "src-local" / "solver_args.sh"

# The one value AGENTS.md forbids changing without explicit instruction.
PINNED_DETECTOR = "0.021005127"

# argv positions, 0-based, in the 26-argument contract.
ARGV_OH = 0
ARGV_MAXLEVEL = 3
ARGV_DROP_RADIUS_MIN = 6
ARGV_DRILL_AMR = 9
ARGV_GEOMETRY = 22
ARGV_MURIN = 25

BASILISK_REF = "v2026-07-20"


def bash_array(source: str, name: str) -> list[str]:
    """Return every element of every assignment to `name=( ... )`."""
    elements: list[str] = []
    for body in re.findall(rf"\b{name}=\(([^)]*)\)", source):
        elements.extend(body.split())
    return elements


class LauncherHarness:
    """Run a launcher for real, with the scheduler and compiler stubbed.

    Reading the launcher as text cannot tell a value that reaches the solver
    from one that is shadowed later, spelt only in a comment, or dropped by a
    quoting slip -- and a ladder point measured with the wrong argument looks
    exactly like a ladder point measured with the right one. So the script is
    executed and the solver's argv is captured from a stub.
    """

    STUBS = {
        # `module` is a shell function on the cluster and simply absent here.
        "module": "#!/bin/bash\nexit 0\n",
        # qcc must produce the file named by -o, since the launcher copies it.
        "qcc": (
            "#!/bin/bash\n"
            'out=""\n'
            'while [[ $# -gt 0 ]]; do\n'
            '  if [[ "$1" == "-o" ]]; then out="$2"; shift 2; continue; fi\n'
            '  shift\n'
            'done\n'
            '[[ -n "$out" ]] || exit 1\n'
            'printf "#!/bin/bash\\ntrue\\n" > "$out"\n'
            'chmod +x "$out"\n'
        ),
        # srun records the solver argv and writes the outputs the launcher's
        # own completeness check looks for.
        "srun": (
            "#!/bin/bash\n"
            "while [[ $# -gt 0 ]]; do\n"
            '  case "$1" in\n'
            "    --exclusive|--exact) shift ;;\n"
            "    -N1|-n1|--cpu-bind=cores) shift ;;\n"
            "    -c) shift 2 ;;\n"
            "    *) break ;;\n"
            "  esac\n"
            "done\n"
            'binary="$1"; shift\n'
            'printf "%s\\n" "$binary" > solver_binary.txt\n'
            'printf "%s\\n" "$#" "$@" > argv.txt\n'
            'printf "i dt t\\n0 1e-6 0\\n" > log\n'
            'printf "t,z_tip\\n%s,0\\n" "${LAUNCHER_TEST_TMAX:-1.5}" > jet.log\n'
            'printf "t,component\\n0,1\\n" > components.log\n'
            'printf "0 0\\n" > foot.dat\n'
            "exit 0\n"
        ),
    }

    def __init__(self, launcher: Path, environment: dict[str, str]) -> None:
        self.launcher = launcher
        self.environment = environment

    def __enter__(self) -> "LauncherHarness":
        self._directory = tempfile.TemporaryDirectory()
        self.root = Path(self._directory.name)

        (self.root / "simulationCases" / "DataFiles").mkdir(parents=True)
        (self.root / "src-local").mkdir()
        shutil.copy(SOLVER_ARGS, self.root / "src-local" / "solver_args.sh")
        # The stub compiler ignores the source, but the launcher expects it.
        (self.root / "simulationCases" / "burstingBubbleInfiniteRr.c").write_text("")
        shutil.copy(self.launcher, self.root / self.launcher.name)

        stub_dir = self.root / "stubs"
        stub_dir.mkdir()
        for name, body in self.STUBS.items():
            path = stub_dir / name
            path.write_text(body)
            path.chmod(0o755)
        self.stub_dir = stub_dir

        # A Basilisk tree that satisfies the launcher's ref lock.
        basilisk = self.root / "basilisk"
        (basilisk / "src").mkdir(parents=True)
        qcc = basilisk / "src" / "qcc"
        shutil.copy(stub_dir / "qcc", qcc)
        qcc.chmod(0o755)
        (basilisk / ".comphy-lock").write_text(f"ref={BASILISK_REF}\n")
        self.basilisk = basilisk
        return self

    def __exit__(self, *exception: object) -> None:
        self._directory.cleanup()

    def run(self, **overrides: str) -> subprocess.CompletedProcess[str]:
        environment = dict(os.environ)
        environment.update(
            {
                "PATH": f"{self.stub_dir}:{environment['PATH']}",
                "SLURM_SUBMIT_DIR": str(self.root),
                "SLURM_CPUS_PER_TASK": "1",
                "BASILISK_ROOT": str(self.basilisk),
                "EXPECTED_BASILISK_REF": BASILISK_REF,
            }
        )
        environment.update(self.environment)
        environment.update(overrides)
        environment["LAUNCHER_TEST_TMAX"] = environment.get("TMAX", "1.5")
        return subprocess.run(
            ["bash", str(self.root / self.launcher.name)],
            capture_output=True,
            text=True,
            check=False,
            env=environment,
            timeout=300,
        )

    def case_argv(self, case_id: str) -> list[str]:
        recorded = (self.root / "simulationCases" / case_id / "argv.txt").read_text()
        fields = recorded.split()
        count, argv = int(fields[0]), fields[1:]
        assert count == len(argv) == 26, f"case {case_id}: {count} arguments, {argv}"
        return argv

    def case_params(self, case_id: str) -> dict[str, str]:
        text = (self.root / "simulationCases" / case_id / "case.params").read_text()
        return dict(
            line.split("=", 1) for line in text.splitlines() if "=" in line
        )


class LauncherContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.ladder = LADDER.read_text()
        self.probe = PROBE.read_text()

    def test_both_launchers_parse(self) -> None:
        for script in (LADDER, PROBE):
            result = subprocess.run(
                ["bash", "-n", str(script)], capture_output=True, text=True, check=False
            )
            self.assertEqual(result.returncode, 0, f"{script.name}: {result.stderr}")

    def test_detector_radius_is_pinned_not_mesh_derived(self) -> None:
        for name, source in (("ladder", self.ladder), ("probe", self.probe)):
            with self.subTest(launcher=name):
                self.assertIn(f'"dropRadiusMin={PINNED_DETECTOR}"', source)
                self.assertNotIn('"dropRadiusMin=0"', source)

    def test_gas_viscosity_ratio_is_passed_explicitly(self) -> None:
        self.assertIn('"MuRin=$MURIN"', self.ladder)
        self.assertIn('MURIN="${MURIN:-1e-2}"', self.ladder)
        self.assertIn('"MuRin=$murin"', self.probe)

    def test_launchers_record_the_resolved_argv_beside_the_outputs(self) -> None:
        for name, source in (("ladder", self.ladder), ("probe", self.probe)):
            with self.subTest(launcher=name):
                self.assertIn("case.params", source)
                self.assertIn("argv=%s", source)

    def test_ladder_case_ids_are_unique_across_groups(self) -> None:
        ids = bash_array(self.ladder, "CASE_IDS")
        self.assertEqual(sorted(ids), sorted(set(ids)), "case id reused between groups")

    def test_ladder_covers_the_forward_injection_boundary(self) -> None:
        """Group F exists to re-measure Oh_c^(U) on the corrected stack."""
        values = {v for v in bash_array(self.ladder, "OH_VALUES") if v[0].isdigit()}
        for oh in ("0.0340", "0.0345", "0.0350", "0.0355"):
            self.assertIn(oh, values)

    # The design the job is meant to run. A reordering of any one array would
    # attribute a divergence to the wrong solver stack or the wrong viscosity
    # ratio, which is the single way this experiment can mislead, so the whole
    # mapping is pinned rather than its cardinality.
    EXPECTED_PROBE_DESIGN = [
        # case,   Oh,       MuRin,  qcc flags
        ("6551", "0.0325", "2e-2", ""),  # decisive: 6330's physics on Stack 1
        ("6552", "0.0280", "2e-2", ""),  # decisive: 6326's physics on Stack 1
        ("6553", "0.0280", "1e-2", "-DSINGLE_PROJECTION=1"),
        ("6554", "0.0280", "1e-2", "-DUSE_CONSERVING=1"),
    ]

    def probe_columns(self) -> dict[str, list[str]]:
        # STACK_FLAGS holds an empty string for Stack 1, which `split()` drops,
        # so it is read with a quote-aware scan instead.
        body = re.search(r"STACK_FLAGS=\(([^)]*)\)", self.probe)
        self.assertIsNotNone(body, "STACK_FLAGS array not found")
        return {
            "CASE_IDS": bash_array(self.probe, "CASE_IDS"),
            "OH_VALUES": bash_array(self.probe, "OH_VALUES"),
            "MURIN_VALUES": bash_array(self.probe, "MURIN_VALUES"),
            "STACK_FLAGS": re.findall(r'"([^"]*)"', body.group(1)),
        }

    def probe_design(self) -> list[tuple[str, str, str, str]]:
        columns = self.probe_columns()
        # `zip` truncates to the shortest input, so a surplus entry in any one
        # column would vanish here instead of being caught. Check first.
        lengths = {name: len(values) for name, values in columns.items()}
        self.assertEqual(len(set(lengths.values())), 1, f"misaligned arrays: {lengths}")
        return list(zip(*columns.values()))

    def test_probe_design_is_exactly_the_intended_experiment(self) -> None:
        self.assertEqual(self.probe_design(), self.EXPECTED_PROBE_DESIGN)

    def test_probe_varies_one_factor_at_a_time(self) -> None:
        """A discriminator that changed two things at once would prove nothing."""
        design = self.probe_design()
        by_case = {row[0]: row for row in design}

        # 6552 and 6551 differ from the two dead drill-solver cases in the
        # solver stack alone: same Oh, same gas viscosity, same everything the
        # launcher pins.
        self.assertEqual(by_case["6552"][1:3], ("0.0280", "2e-2"))
        self.assertEqual(by_case["6551"][1:3], ("0.0325", "2e-2"))
        self.assertEqual(by_case["6552"][3], "", "6552 must run the default Stack 1")

        # 6553 and 6554 differ from the running crossdriver case 6506
        # (Oh=0.0280, MuRin=1e-2, Stack 1) in the stack alone.
        for case in ("6553", "6554"):
            self.assertEqual(by_case[case][1:3], ("0.0280", "1e-2"))

        # And Oh=0.0280 is present at both viscosity ratios, so the gas-viscosity
        # sensitivity is readable at fixed stack.
        pairs = {row[1:3] for row in design}
        self.assertIn(("0.0280", "2e-2"), pairs)
        self.assertIn(("0.0280", "1e-2"), pairs)

    def test_probe_arrays_are_aligned(self) -> None:
        lengths = {name: len(values) for name, values in self.probe_columns().items()}
        self.assertEqual(len(set(lengths.values())), 1, lengths)
        self.assertEqual(next(iter(lengths.values())), 4)
        # STACK_NAMES only labels the summary table, but a mislabelled summary
        # is how a correct experiment gets read backwards.
        self.assertEqual(len(re.findall(r'"([^"]*)"', re.search(
            r"STACK_NAMES=\((.*?)\)", self.probe, re.S).group(1))), 4)


class LauncherExecutionTests(unittest.TestCase):
    """What the solver is actually handed, not what the script appears to say."""

    def test_ladder_group_f_reaches_the_solver_intact(self) -> None:
        with LauncherHarness(LADDER, {"GROUP": "F", "TMAX": "1.5"}) as harness:
            result = harness.run()
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            expected = {
                "6513": "0.0340",
                "6514": "0.0345",
                "6515": "0.0350",
                "6516": "0.0355",
            }
            for case_id, oh in expected.items():
                with self.subTest(case=case_id):
                    argv = harness.case_argv(case_id)
                    self.assertEqual(argv[ARGV_OH], oh)
                    self.assertEqual(argv[ARGV_MAXLEVEL], "13")
                    self.assertEqual(argv[ARGV_DROP_RADIUS_MIN], PINNED_DETECTOR)
                    self.assertEqual(argv[ARGV_DRILL_AMR], "0")
                    self.assertEqual(argv[ARGV_GEOMETRY], "halfspace")
                    self.assertEqual(argv[ARGV_MURIN], "1e-2")

                    params = harness.case_params(case_id)
                    self.assertEqual(params["Oh"], oh)
                    self.assertEqual(params["MuRin"], "1e-2")
                    self.assertEqual(params["solverStack"],
                                     "filtered+double-projection")
                    self.assertEqual(params["argv"].split(), argv)

    def test_ladder_honours_an_overridden_gas_viscosity_ratio(self) -> None:
        with LauncherHarness(LADDER, {"GROUP": "B", "MURIN": "2e-2"}) as harness:
            result = harness.run()
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(harness.case_argv("6506")[ARGV_MURIN], "2e-2")
            self.assertEqual(harness.case_params("6506")["MuRin"], "2e-2")

    def test_ladder_rejects_an_unknown_group(self) -> None:
        with LauncherHarness(LADDER, {"GROUP": "Z"}) as harness:
            result = harness.run()
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("GROUP must be", result.stderr)

    def test_probe_runs_its_declared_design(self) -> None:
        with LauncherHarness(PROBE, {"TMAX": "0.75"}) as harness:
            result = harness.run()
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            for case_id, oh, murin, flags in \
                    LauncherContractTests.EXPECTED_PROBE_DESIGN:
                with self.subTest(case=case_id):
                    argv = harness.case_argv(case_id)
                    self.assertEqual(argv[ARGV_OH], oh)
                    self.assertEqual(argv[ARGV_MURIN], murin)
                    self.assertEqual(argv[ARGV_DROP_RADIUS_MIN], PINNED_DETECTOR)

                    params = harness.case_params(case_id)
                    self.assertEqual(params["MuRin"], murin)
                    self.assertEqual(params["qccFlags"], flags)
                    self.assertEqual(params["argv"].split(), argv)

            # Each case must run its own binary: the stacks are mutually
            # exclusive at compile time and cannot share an executable.
            binaries = {
                (harness.root / "simulationCases" / case / "solver_binary.txt")
                .read_text().strip()
                for case, *_ in LauncherContractTests.EXPECTED_PROBE_DESIGN
            }
            self.assertEqual(len(binaries), 4, binaries)

    def test_concurrent_groups_do_not_share_a_build_output(self) -> None:
        """Groups share one staging root, so a fixed binary name races.

        Submitted together on 2026-08-30, group C died 32 s in with
        `cp: cannot open 'burstingBubbleInfiniteRr' for reading` because
        another group's `qcc` was rewriting the file it was copying.
        """
        groups = (("A", "6501", "111"), ("C", "6509", "222"))
        with LauncherHarness(LADDER, {"TMAX": "1.5"}) as harness:
            # Both launchers run at once, against one staging root, which is
            # the situation that actually broke. With a per-job build output
            # the two share no mutable file, so this stays deterministic.
            with ThreadPoolExecutor(max_workers=len(groups)) as pool:
                futures = {
                    group: pool.submit(harness.run, GROUP=group, SLURM_JOB_ID=job)
                    for group, _, job in groups
                }
                results = {group: future.result() for group, future in futures.items()}

            names: dict[str, str] = {}
            for group, case, _ in groups:
                result = results[group]
                self.assertEqual(
                    result.returncode, 0,
                    f"group {group}:\n{result.stdout}{result.stderr}",
                )
                names[group] = (
                    harness.root / "simulationCases" / case / "solver_binary.txt"
                ).read_text().strip()

            self.assertNotEqual(
                names["A"], names["C"], f"both groups built into {names['A']}"
            )
            for group, name in names.items():
                self.assertIn(group, name)

    def test_probe_summary_reports_a_divergence_rather_than_failing(self) -> None:
        """A case dying the way 6326 and 6330 died is this job's result."""
        with LauncherHarness(PROBE, {"TMAX": "0.75"}) as harness:
            # Make one case stop early with the solver's own ke-gate message.
            srun = harness.stub_dir / "srun"
            srun.write_text(
                srun.read_text().replace(
                    'printf "t,z_tip\\n%s,0\\n" "${LAUNCHER_TEST_TMAX:-1.5}" > jet.log',
                    'if [[ "$PWD" == *6554 ]]; then\n'
                    '  printf "t,z_tip\\n0.50,0\\n" > jet.log\n'
                    '  echo "The kinetic energy blew up" > runner.err\n'
                    "else\n"
                    '  printf "t,z_tip\\n%s,0\\n" "${LAUNCHER_TEST_TMAX:-1.5}" > jet.log\n'
                    "fi",
                )
            )
            result = harness.run()
            self.assertEqual(result.returncode, 0, "a divergence is a result, not a failure")
            self.assertRegex(result.stdout, r"6554.*DIVERGED")
            self.assertRegex(result.stdout, r"6551.*reached horizon")


if __name__ == "__main__":
    unittest.main()
