from __future__ import annotations

import re
import shutil
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
# `if (argc > 26 && !parse_positive_double (argv[26], "MuRin", &MuRin))` --
# the strictly validated form used where a bad value would otherwise be silent.
VALIDATED_ARGV = re.compile(r"parse_\w+\s*\(\s*argv\[(\d+)\]\s*,\s*\"(\w+)\"")


def solver_argv_order() -> list[str]:
    """Read the solver's own argv assignments as the authoritative order."""
    found: dict[int, str] = {}
    for line in SOLVER.read_text().splitlines():
        validated = VALIDATED_ARGV.search(line)
        if validated:
            found[int(validated.group(1))] = validated.group(2)
            continue
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

    def test_the_viscosity_ratio_uses_the_strict_parser(self) -> None:
        """`StrictNumericParserTests` shows the parser rejects bad input; this
        shows argv 26 actually goes through it and not through `atof`."""
        code = re.sub(r"/\*.*?\*/|//[^\n]*", "", SOLVER.read_text(), flags=re.DOTALL)
        self.assertRegex(
            code,
            r'\bparse_positive_double\s*\(\s*argv\[26\]\s*,\s*"MuRin"\s*,\s*&MuRin\s*\)',
        )
        self.assertNotRegex(code, r"\bMuRin\s*=\s*atof\b")

    def test_the_bond_number_uses_the_nonnegative_parser(self) -> None:
        """Bond = 0 is the existing ladder; a negative value is not a Bond."""
        code = re.sub(r"/\*.*?\*/|//[^\n]*", "", SOLVER.read_text(), flags=re.DOTALL)
        self.assertRegex(
            code,
            r'\bparse_nonnegative_double\s*\(\s*argv\[27\]\s*,\s*"Bond"\s*,\s*&Bond\s*\)',
        )
        self.assertNotRegex(code, r"\bBond\s*=\s*atof\b")


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
        self.assertEqual(count, "27")
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

    def test_forwards_all_twenty_seven_arguments(self) -> None:
        result, case_dir = self.run_case(dict(self.PARAMS))
        self.assertEqual(result.returncode, 0, result.stderr)
        recorded = (case_dir / "argv.txt").read_text().split()
        self.assertEqual(recorded[0], "27")
        argv = recorded[1:]
        self.assertEqual(argv[22], "finite")
        self.assertEqual(argv[23], "0.027")
        # argv 25 is absent from older case tables and defaults to the
        # solver's compiled value.
        self.assertEqual(argv[24], "1")
        # argv 26 likewise. Its default is the gas-to-liquid viscosity ratio
        # the manuscript states; the drill-solver ladder ran 2e-2 with nothing
        # in argv to show it, which is why the value is passed at all.
        self.assertEqual(argv[25], "1e-2")
        # argv 27 is Bond. Zero keeps Bo0.0000.dat and G.x = 0 in reduced.h.
        self.assertEqual(argv[26], "0")

    def test_honours_an_explicit_interface_floor(self) -> None:
        params = dict(self.PARAMS)
        params["interfaceFloor"] = "0"
        result, case_dir = self.run_case(params)
        self.assertEqual(result.returncode, 0, result.stderr)
        argv = (case_dir / "argv.txt").read_text().split()[1:]
        self.assertEqual(argv[24], "0")

    def test_honours_an_explicit_gas_viscosity_ratio(self) -> None:
        """The default alone would pass even if the key were never wired up."""
        params = dict(self.PARAMS)
        params["MuRin"] = "2e-2"
        result, case_dir = self.run_case(params)
        self.assertEqual(result.returncode, 0, result.stderr)
        recorded = (case_dir / "argv.txt").read_text().split()
        self.assertEqual(recorded[0], "27")
        self.assertEqual(recorded[1:][25], "2e-2")

    def test_rejects_a_truncated_case_table(self) -> None:
        params = dict(self.PARAMS)
        del params["drillNcells"]
        result, _ = self.run_case(params)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing drillNcells", result.stderr)


if __name__ == "__main__":
    unittest.main()


class StrictNumericParserTests(unittest.TestCase):
    """Exercise `parse_positive_double` itself.

    `atof` would accept `1e-2x` as 0.01 and `1e9999` as an infinity that passes
    an ordinary positivity test. Every other solver argument shows a wrong
    value immediately -- in the banner, the domain size or the first snapshot.
    A wrong property ratio does not: the run looks entirely normal and simply
    answers a different question, which is how the drop-map ladder came to sit
    at a gas viscosity the manuscript never states. So the parser is compiled
    and run rather than pattern-matched.

    The function is plain C, so it is lifted out of the Basilisk source and
    built with the system compiler; `ferr` is Basilisk's stderr alias and is
    stubbed.
    """

    HARNESS = """
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
static FILE * ferr;
%s
int main (int argc, char ** argv) {
  ferr = stderr;
  double value = -12345.;
  bool ok = parse_positive_double (argc > 1 ? argv[1] : NULL, "MuRin", &value);
  printf ("%%d %%.17g\\n", (int) ok, value);
  return ok ? 0 : 1;
}
"""

    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which("cc") or shutil.which("gcc")
        if compiler is None:
            raise unittest.SkipTest("no C compiler available")

        source = SOLVER.read_text()
        start = source.index("static bool parse_positive_double")
        end = source.index("\n}\n", start) + len("\n}\n")
        function = source[start:end]

        cls._directory = tempfile.TemporaryDirectory()
        root = Path(cls._directory.name)
        (root / "harness.c").write_text(cls.HARNESS % function)
        cls.binary = root / "harness"
        build = subprocess.run(
            [compiler, "-std=c11", "-Wall", "-Werror", str(root / "harness.c"),
             "-o", str(cls.binary), "-lm"],
            capture_output=True, text=True, check=False,
        )
        if build.returncode != 0:
            raise AssertionError(f"harness did not build:\n{build.stderr}")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._directory.cleanup()

    def parse(self, text: str | None) -> subprocess.CompletedProcess[str]:
        command = [str(self.binary)] + ([] if text is None else [text])
        return subprocess.run(command, capture_output=True, text=True, check=False)

    def test_accepts_well_formed_positive_values(self) -> None:
        for text, expected in (("2e-2", 0.02), ("1e-2", 0.01), ("0.018", 0.018)):
            with self.subTest(text=text):
                result = self.parse(text)
                self.assertEqual(result.returncode, 0, result.stderr)
                accepted, value = result.stdout.split()
                self.assertEqual(accepted, "1")
                self.assertAlmostEqual(float(value), expected, places=15)

    def test_rejects_values_atof_would_have_swallowed(self) -> None:
        cases = {
            "1e-2x": "is not a number",     # atof returns 0.01
            "abc": "is not a number",       # atof returns 0
            "": "empty argument",           # atof returns 0
            "1e9999": "is not a finite number",  # atof returns +inf, which is > 0
            "nan": "is not a finite number",      # and NaN fails every comparison
            "1e-9999": "is out of range",         # underflows to a denormal or 0
            "-1": "must be strictly positive",
            "0": "must be strictly positive",
        }
        for text, message in cases.items():
            with self.subTest(text=text):
                result = self.parse(text)
                self.assertNotEqual(result.returncode, 0, f"{text!r} was accepted")
                self.assertIn(message, result.stderr)
                # A rejected parse must leave the caller's value untouched, so
                # a compiled default survives a bad argument rather than being
                # half-overwritten.
                self.assertEqual(result.stdout.split(), ["0", "-12345"])

    def test_rejects_a_missing_argument(self) -> None:
        result = self.parse(None)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("empty argument", result.stderr)
