from __future__ import annotations

import re
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
LIBRARY = ROOT / "src-local" / "solver_args.sh"
PARSE_PARAMS = ROOT / "src-local" / "parse_params.sh"
SOLVER = ROOT / "simulationCases" / "coalescenceBubble.c"
RUNNER = ROOT / "BayesianWorkflow" / "run_one_contour_case.sh"

# `OhOut = atof(argv[1]);` and `snprintf (geometryMode, sizeof(geometryMode), "%s", argv[23]);`
DIRECT_ARGV = re.compile(r"^\s*(\w+)\s*=\s*ato[fi]\s*\(argv\[(\d+)\]\)")
SNPRINTF_ARGV = re.compile(r"^\s*snprintf\s*\(\s*(\w+),.*argv\[(\d+)\]")


def solver_argv_order() -> list[str]:
    """Read the solver's own argv assignments as the authoritative order."""
    found: dict[int, str] = {}
    for line in SOLVER.read_text().splitlines():
        match = DIRECT_ARGV.match(line) or SNPRINTF_ARGV.match(line)
        if match:
            found[int(match.group(2))] = match.group(1)
    return [found[index] for index in sorted(found)]


def run_bash(script: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", "-c", script], text=True, capture_output=True, check=False, cwd=ROOT
    )


class SolverArgumentContractTests(unittest.TestCase):
    def library_keys(self) -> list[str]:
        result = run_bash(
            f'source "{LIBRARY}"; printf "%s\\n" "${{COALESCENCE_SOLVER_ARG_KEYS[@]}}"'
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.split()

    def test_key_order_matches_the_solver(self) -> None:
        """The shared contract is the reason four launchers drifted out of date."""
        self.assertEqual(self.library_keys(), solver_argv_order())

    def test_builds_every_argument_from_a_parameter_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            params = Path(temporary) / "case.params"
            params.write_text("OhOut=0.03\nRr=4.00\nMAXlevel=13\ntmax=1.0\n")
            result = run_bash(
                f'source "{PARSE_PARAMS}"; source "{LIBRARY}"; '
                f'parse_param_file "{params}" >/dev/null; '
                "coalescence_solver_args_from_params; "
                'printf "%s\\n" "${#COALESCENCE_SOLVER_ARGS[@]}" '
                '"${COALESCENCE_SOLVER_ARGS[0]}" "${COALESCENCE_SOLVER_ARGS[4]}" '
                '"${COALESCENCE_SOLVER_ARGS[22]}" "${COALESCENCE_SOLVER_ARGS[24]}"'
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        count, oh, tmax, geometry, floor = result.stdout.split()
        self.assertEqual(count, "25")
        self.assertEqual(oh, "0.03")
        self.assertEqual(tmax, "1.0")
        # Unset keys fall back to the solver's own compiled defaults.
        self.assertEqual(geometry, "finite")
        self.assertEqual(floor, "1")

    def test_overrides_replace_the_parameter_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            params = Path(temporary) / "case.params"
            params.write_text("OhOut=0.03\ntmax=1.0\ninterfaceFloor=0\n")
            result = run_bash(
                f'source "{PARSE_PARAMS}"; source "{LIBRARY}"; '
                f'parse_param_file "{params}" >/dev/null; '
                'coalescence_solver_args_from_params "tmax=0.05"; '
                'printf "%s\\n" "${COALESCENCE_SOLVER_ARGS[4]}" '
                '"${COALESCENCE_SOLVER_ARGS[24]}"'
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.split(), ["0.05", "0"])

    def test_omitted_keys_do_not_leak_from_the_previous_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first.params"
            second = Path(temporary) / "second.params"
            first.write_text("OhOut=0.03\nRhoIn=9e-2\nRr=2.0\nMAXlevel=13\ntmax=1.0\nzWall=0.05\n")
            second.write_text("OhOut=0.04\nRr=3.0\nMAXlevel=13\ntmax=1.0\nzWall=0.05\n")
            result = run_bash(
                f'source "{PARSE_PARAMS}"; source "{LIBRARY}"; '
                f'parse_param_file "{first}" >/dev/null; '
                f'parse_param_file "{second}" >/dev/null; '
                "coalescence_solver_args_from_params; "
                'printf "%s\\n" "${COALESCENCE_SOLVER_ARGS[0]}" '
                '"${COALESCENCE_SOLVER_ARGS[1]}" "${COALESCENCE_SOLVER_ARGS[2]}"'
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        oh, rho, rr = result.stdout.split()
        self.assertEqual(oh, "0.04")
        self.assertEqual(rho, "1e-3")
        self.assertEqual(rr, "3.0")

    def test_parse_leaves_unrelated_param_file_variables_alone(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            params = Path(temporary) / "case.params"
            params.write_text("OhOut=0.03\n")
            result = run_bash(
                f'PARAM_FILE=keep-me; PARAM_FILES=also-keep; '
                f'source "{PARSE_PARAMS}"; '
                f'parse_param_file "{params}" >/dev/null; '
                'printf "%s\\n" "$PARAM_FILE" "$PARAM_FILES"'
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.split(), ["keep-me", "also-keep"])

    def test_parse_ignores_reserved_file_keys(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            params = Path(temporary) / "case.params"
            params.write_text("OhOut=0.03\nFILE=clobber\nFILES=also-clobber\n")
            result = run_bash(
                f'PARAM_FILE=keep-me; PARAM_FILES=also-keep; '
                f'source "{PARSE_PARAMS}"; '
                f'parse_param_file "{params}" >/dev/null; '
                'printf "%s\\n" "$PARAM_FILE" "$PARAM_FILES" "$PARAM_OhOut"'
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.split(), ["keep-me", "also-keep", "0.03"])
        self.assertIn("reserved parameter key FILE", result.stderr)
        self.assertIn("reserved parameter key FILES", result.stderr)

    def test_failed_parse_leaves_the_previous_file_loaded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first.params"
            bad = Path(temporary) / "bad.params"
            first.write_text("OhOut=0.03\nRhoIn=9e-2\n")
            bad.write_text("OhOut=0.04\nbad-key=1\n")
            result = run_bash(
                f'source "{PARSE_PARAMS}"; '
                f'parse_param_file "{first}" >/dev/null; '
                f'parse_param_file "{bad}" >/dev/null; '
                'printf "%s\\n" "$?"; '
                'printf "%s\\n" "${PARAM_OhOut}" "${PARAM_RhoIn}"'
            )
        lines = result.stdout.split()
        self.assertEqual(lines[0], "1")
        self.assertEqual(lines[1:], ["0.03", "9e-2"])
        self.assertIn("invalid parameter key", result.stderr)

    def test_reports_a_missing_required_argument(self) -> None:
        result = run_bash(
            f'source "{LIBRARY}"; nothing() {{ :; }}; '
            "coalescence_solver_args nothing"
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing required solver argument OhOut", result.stderr)


class ContourRunnerTests(unittest.TestCase):
    PARAMS = {
        "OhOut": "0.03",
        "RhoIn": "1e-3",
        "Rr": "4.00",
        "MAXlevel": "13",
        "tmax": "1.0",
        "zWall": "0.05",
        "dropRadiusMin": "0.021005127",
        "dropPersistence": "3",
        "snapshotInterval": "0.01",
        "drillAMR": "0",
        "drillMaxlevelStart": "6",
        "drillMaxlevelFocus": "9",
        "drillNcells": "5",
        "drillRegionMinX": "-3",
        "drillArmSteps": "5",
        "drillArmTime": "0.45",
        "drillCoarsenTime": "0",
        "drillRegionMaxX": "3",
        "drillRegionRadius": "1.5",
        "drillFireX": "0.25",
        "drillTipRadius": "0.25",
        "drillRegionalOnly": "0",
        "geometryMode": "finite",
        "wallClearance": "0.027",
    }

    def run_case(self, params: dict[str, str]) -> tuple[subprocess.CompletedProcess[str], Path]:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        case_dir = root / "case"
        case_dir.mkdir()
        (case_dir / "case.params").write_text(
            "# case\n" + "".join(f"{key}={value}\n" for key, value in params.items())
        )
        stub = root / "stub.sh"
        stub.write_text(
            '#!/bin/bash\nprintf "%s\\n" "$#" "$@" > argv.txt\n'
            'printf "id,reason\\n1,ok\\n" > classification.status\n'
        )
        stub.chmod(0o755)
        result = subprocess.run(
            ["bash", str(RUNNER), str(case_dir), str(stub)],
            text=True,
            capture_output=True,
            check=False,
        )
        return result, case_dir

    def test_forwards_all_twenty_five_arguments(self) -> None:
        result, case_dir = self.run_case(dict(self.PARAMS))
        self.assertEqual(result.returncode, 0, result.stderr)
        recorded = (case_dir / "argv.txt").read_text().split()
        self.assertEqual(recorded[0], "25")
        argv = recorded[1:]
        self.assertEqual(argv[22], "finite")
        self.assertEqual(argv[23], "0.027")
        # argv 25 is absent from older case tables and defaults to the
        # solver's compiled value.
        self.assertEqual(argv[24], "1")

    def test_honours_an_explicit_interface_floor(self) -> None:
        params = dict(self.PARAMS)
        params["interfaceFloor"] = "0"
        result, case_dir = self.run_case(params)
        self.assertEqual(result.returncode, 0, result.stderr)
        argv = (case_dir / "argv.txt").read_text().split()[1:]
        self.assertEqual(argv[24], "0")

    def test_rejects_a_truncated_case_table(self) -> None:
        params = dict(self.PARAMS)
        del params["drillNcells"]
        result, _ = self.run_case(params)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing drillNcells", result.stderr)


if __name__ == "__main__":
    unittest.main()
