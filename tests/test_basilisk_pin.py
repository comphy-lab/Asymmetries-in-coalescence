"""Fail-closed version binding for current and legacy launch paths."""
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "src-local" / "basilisk_pin.sh"
LEGACY = ROOT / "src-local" / "basilisk_version.sh"
CURRENT = "v2026-07-20"
OLD = "v2026-01-13"
LAUNCHERS = (
    "runInfiniteRrCampaign.sbatch", "runCoalescenceCampaign.sbatch",
    "runBurstingBubbleInfiniteRr.sbatch", "runInfiniteRrRestart.sbatch",
    "runInfiniteRrDeltaAnchor.sbatch", "runInfiniteRrStackProbe.sbatch",
    "runInfiniteRrBondProbe.sbatch",
)


class PinTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "src-local").mkdir()
        shutil.copy(HELPER, self.root / "src-local" / HELPER.name)
        shutil.copy(LEGACY, self.root / "src-local" / LEGACY.name)
        shutil.copy(ROOT / "src-local" / "solver_args.sh", self.root / "src-local" / "solver_args.sh")
        (self.root / "basilisk.ref").write_text(CURRENT + "\n")
        shutil.copy(ROOT / "basilisk-legacy.ref", self.root / "basilisk-legacy.ref")
        self.install = self.root / "install"
        (self.install / "src").mkdir(parents=True)
        self.lock = self.install / ".comphy-lock"
        self.lock.write_text("ref=" + CURRENT + "\n")
        self.stub = self.root / "stub"
        self.stub.mkdir()
        module = self.stub / "module"
        module.write_text("#!/bin/sh\nexit 0\n")
        module.chmod(0o755)
        self.curl_marker = self.root / "curl-was-called"
        curl = self.stub / "curl"
        curl.write_text(f"#!/bin/sh\ntouch '{self.curl_marker}'\n")
        curl.chmod(0o755)
        self.env = dict(os.environ, PATH=f"{self.stub}:{os.environ['PATH']}",
                        SLURM_SUBMIT_DIR=str(self.root), BASILISK_ROOT=str(self.install),
                        SLURM_JOB_ID="999", CASE="6564", TMAX="1.5",
                        MANIFEST="campaign.tsv")
        self.env.pop("EXPECTED_BASILISK_REF", None)
        (self.root / "campaign.tsv").write_text("# guard test\n")

    def helper(self, command, **extra):
        return subprocess.run(["bash", "-eu", "-c", command, "bash", str(self.root)],
                              env=dict(self.env, **extra), capture_output=True, text=True)

    def run_launcher(self, name, **extra):
        shutil.copy(ROOT / name, self.root / name)
        return subprocess.run(["bash", str(self.root / name)], env=dict(self.env, **extra),
                              cwd=self.root, capture_output=True, text=True, timeout=10)

    def test_reader_accepts_one_commented_pin_and_matching_override(self):
        (self.root / "basilisk.ref").write_text(" # comment\n  v2026-07-20   # pinned\n\n")
        result = self.helper('source "$1/src-local/basilisk_pin.sh"; '
                             'basilisk_require_campaign_pin "$1"; printf "%s" "$EXPECTED_BASILISK_REF"',
                             EXPECTED_BASILISK_REF=CURRENT)
        self.assertEqual((result.returncode, result.stdout), (0, CURRENT), result.stderr)

    def test_reader_rejects_missing_malformed_and_multiple_pins(self):
        for content in (None, "", "v1\nv2\n", "$(touch injected)\n", "../bad\n", "v2026-07- 20\n"):
            with self.subTest(content=content):
                path = self.root / "basilisk.ref"
                path.unlink(missing_ok=True)
                if content is not None:
                    path.write_text(content)
                result = self.helper('source "$1/src-local/basilisk_pin.sh"; basilisk_read_pin "$1/basilisk.ref"')
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse((self.root / "injected").exists())

    def test_lock_requires_one_matching_ref(self):
        command = 'source "$1/src-local/basilisk_pin.sh"; basilisk_require_lock "$1/install/.comphy-lock" v2026-07-20'
        for content, ok in (("ref=v2026-07-20\nos=linux\n", True),
                            ("ref=v2026-01-13\n", False),
                            ("ref=v2026-07-20\nref=v2026-01-13\n", False),
                            ("os=linux\n", False)):
            with self.subTest(content=content):
                self.lock.write_text(content)
                self.assertEqual(self.helper(command).returncode == 0, ok)
        self.lock.unlink()
        self.assertNotEqual(self.helper(command).returncode, 0)

    def test_all_seven_launchers_reject_conflicting_override_before_qcc(self):
        for name in LAUNCHERS:
            with self.subTest(name=name):
                result = self.run_launcher(name, EXPECTED_BASILISK_REF=OLD)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("conflicts with the component pin", result.stderr)

    def test_all_seven_launchers_reject_missing_or_invalid_pin(self):
        for value in (None, "v1\nv2\n"):
            for name in LAUNCHERS:
                with self.subTest(name=name, value=value):
                    p = self.root / "basilisk.ref"
                    p.unlink(missing_ok=True)
                    if value:
                        p.write_text(value)
                    result = self.run_launcher(name)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("Basilisk pin", result.stderr)

    def test_all_seven_launchers_read_valid_component_pin(self):
        # No qcc is installed. Reaching the existing qcc check proves the
        # actual launcher passed the pin reader, before any compile can run.
        for name in LAUNCHERS:
            with self.subTest(name=name):
                result = self.run_launcher(name)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("no qcc", result.stderr)
                self.assertNotIn("Basilisk pin", result.stderr)

    def test_all_seven_launchers_reject_incorrect_or_ambiguous_lock(self):
        qcc = self.install / "src" / "qcc"
        marker = self.root / "qcc-was-called"
        qcc.write_text(f"#!/bin/sh\ntouch '{marker}'\n")
        qcc.chmod(0o755)
        for contents in ("ref=" + OLD + "\n", "ref=" + CURRENT + "\nref=" + OLD + "\n"):
            for name in LAUNCHERS:
                with self.subTest(name=name, contents=contents):
                    self.lock.write_text(contents)
                    result = self.run_launcher(name)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("expected Basilisk ref=", result.stderr)
                    self.assertFalse(marker.exists())

    def test_legacy_mismatch_refuses_before_installer(self):
        marker = self.root / "curl-was-called"
        curl = self.stub / "curl"
        curl.write_text(f"#!/bin/sh\ntouch '{marker}'\n")
        curl.chmod(0o755)
        (self.root / "basilisk-legacy.ref").write_text(CURRENT + "\n")
        result = self.helper('cd "$1"; source "$1/src-local/basilisk_version.sh"; install_basilisk --hard')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("legacy workflows require", result.stderr)
        self.assertFalse(marker.exists())

    def test_legacy_matching_pin_keeps_january_and_rechecks_before_install(self):
        result = self.helper('source "$1/src-local/basilisk_version.sh"; printf "%s" "$BASILISK_REF"')
        self.assertEqual((result.returncode, result.stdout), (0, OLD), result.stderr)
        result = self.helper('cd "$1"; source "$1/src-local/basilisk_version.sh"; '
                             'printf "v2026-07-20\\n" > "$1/basilisk-legacy.ref"; install_basilisk --hard')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("installation refused", result.stderr)
        self.assertFalse(self.curl_marker.exists())

    def test_real_single_run_help_works_with_both_version_pins(self):
        result = subprocess.run(["bash", str(ROOT / "runSimulation.sh"), "--help"],
                                cwd=ROOT, env=self.env, capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Usage:", result.stdout)

    def test_legacy_pin_missing_or_invalid_never_calls_curl(self):
        for value in (None, "", "v1\nv2\n", "bad!\n", CURRENT + "\n"):
            with self.subTest(value=value):
                pin = self.root / "basilisk-legacy.ref"
                pin.unlink(missing_ok=True)
                if value is not None:
                    pin.write_text(value)
                result = self.helper('cd "$1"; source "$1/src-local/basilisk_version.sh"; install_basilisk --hard')
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.curl_marker.exists())

    def test_legacy_installer_rejects_another_working_directory(self):
        marker = self.root / "curl-was-called"
        curl = self.stub / "curl"
        curl.write_text(f"#!/bin/sh\ntouch '{marker}'\n")
        curl.chmod(0o755)
        (self.root / "other").mkdir()
        result = self.helper('source "$1/src-local/basilisk_version.sh"; cd "$1/other"; install_basilisk --hard')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("own component root", result.stderr)
        self.assertFalse(marker.exists())

    def test_legacy_installer_never_overwrites_unlocked_or_other_ref(self):
        install = self.root / "basilisk"
        install.mkdir()
        marker = self.root / "curl-was-called"
        curl = self.stub / "curl"
        curl.write_text(f"#!/bin/sh\ntouch '{marker}'\n")
        curl.chmod(0o755)
        command = 'cd "$1"; source "$1/src-local/basilisk_version.sh"; install_basilisk --hard'
        for contents in (None, "ref=" + CURRENT + "\n", "ref=" + OLD + "\nref=" + CURRENT + "\n"):
            with self.subTest(contents=contents):
                lock = install / ".comphy-lock"
                lock.unlink(missing_ok=True)
                if contents is not None:
                    lock.write_text(contents)
                self.assertNotEqual(self.helper(command).returncode, 0)
                self.assertFalse(marker.exists())
        (install / ".comphy-lock").write_text("ref=" + OLD + "\n")
        result = self.helper(command)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(marker.exists())

    def test_local_contour_refuses_wrong_or_duplicate_lock_before_qcc(self):
        shutil.copy(ROOT / "runContourLocal.sh", self.root / "runContourLocal.sh")
        (self.root / "basilisk" / "src").mkdir(parents=True)
        qcc = self.root / "basilisk" / "src" / "qcc"
        marker = self.root / "qcc-was-called"
        qcc.write_text(f"#!/bin/sh\ntouch '{marker}'\n")
        qcc.chmod(0o755)
        run = self.root / "attempt"
        run.mkdir()
        (run / "cases.csv").write_text("id\n1\n")
        for contents in ("ref=" + CURRENT + "\n", "ref=" + CURRENT + "\nref=" + OLD + "\n"):
            with self.subTest(contents=contents):
                (self.root / "basilisk" / ".comphy-lock").write_text(contents)
                env = dict(self.env, CONTOUR_QCC=str(qcc), CONTOUR_THREADS_PER_CASE="1",
                           CONTOUR_WORKERS="1", CONTOUR_CPU_COUNT="1", CONTOUR_MAX_THREADS="1")
                result = subprocess.run(["bash", str(self.root / "runContourLocal.sh"), str(run)],
                                        env=env, capture_output=True, text=True, timeout=10)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("expected Basilisk ref=", result.stderr)
                self.assertFalse(marker.exists())

    def test_local_contour_rejects_external_qcc_even_when_root_lock_matches(self):
        shutil.copy(ROOT / "runContourLocal.sh", self.root / "runContourLocal.sh")
        install = self.root / "basilisk"
        (install / "src").mkdir(parents=True)
        (install / ".comphy-lock").write_text("ref=" + OLD + "\n")
        expected = install / "src" / "qcc"
        local_marker = self.root / "local-qcc-called"
        expected.write_text(f"#!/bin/sh\n: <<'END'\n{install / 'src'}\nEND\ntouch '{local_marker}'\nexit 42\n")
        expected.chmod(0o755)
        other = self.root / "other" / "qcc"
        other.parent.mkdir()
        marker = self.root / "external-qcc-called"
        other.write_text(f"#!/bin/sh\ntouch '{marker}'\n")
        other.chmod(0o755)
        attempt = self.root / "attempt"
        attempt.mkdir()
        (attempt / "cases.csv").write_text("id\n1\n")
        for override in (str(other), "qcc"):
            with self.subTest(override=override):
                env = dict(self.env, PATH=f"{other.parent}:{self.env['PATH']}",
                           CONTOUR_QCC=override, CONTOUR_THREADS_PER_CASE="1",
                           CONTOUR_WORKERS="1", CONTOUR_CPU_COUNT="1", CONTOUR_MAX_THREADS="1")
                result = subprocess.run(["bash", str(self.root / "runContourLocal.sh"), str(attempt)],
                                        env=env, capture_output=True, text=True, timeout=10)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(marker.exists())
        env = dict(self.env, PATH=f"{other.parent}:{self.env['PATH']}",
                   CONTOUR_THREADS_PER_CASE="1", CONTOUR_WORKERS="1",
                   CONTOUR_CPU_COUNT="1", CONTOUR_MAX_THREADS="1")
        env.pop("CONTOUR_QCC", None)
        result = subprocess.run(["bash", str(self.root / "runContourLocal.sh"), str(attempt)],
                                env=env, capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 42, result.stderr)
        self.assertTrue(local_marker.exists())
        self.assertFalse(marker.exists())
        expected.unlink()
        expected.symlink_to(other)
        local_marker.unlink()
        env["CONTOUR_QCC"] = str(expected)
        result = subprocess.run(["bash", str(self.root / "runContourLocal.sh"), str(attempt)],
                                env=env, capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("outside the verified", result.stderr)
        self.assertFalse(marker.exists())
        self.assertFalse(local_marker.exists())

    def test_local_contour_rejects_relocated_qcc_with_wrong_embedded_path(self):
        shutil.copy(ROOT / "runContourLocal.sh", self.root / "runContourLocal.sh")
        install = self.root / "basilisk"
        (install / "src").mkdir(parents=True)
        (install / ".comphy-lock").write_text("ref=" + OLD + "\n")
        qcc = install / "src" / "qcc"
        marker = self.root / "wrong-qcc-called"
        qcc.write_text(f"#!/bin/sh\n: <<'END'\n/old/install/src\nEND\ntouch '{marker}'\n")
        qcc.chmod(0o755)
        attempt = self.root / "attempt"
        attempt.mkdir()
        (attempt / "cases.csv").write_text("id\n1\n")
        env = dict(self.env, CONTOUR_THREADS_PER_CASE="1", CONTOUR_WORKERS="1",
                   CONTOUR_CPU_COUNT="1", CONTOUR_MAX_THREADS="1")
        result = subprocess.run(["bash", str(self.root / "runContourLocal.sh"), str(attempt)],
                                env=env, capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not embed", result.stderr)
        self.assertFalse(marker.exists())

    def test_docs_installer_uses_component_pin_and_refuses_unsafe_install(self):
        spec = importlib.util.spec_from_file_location("generate_docs", ROOT / ".github/scripts/generate_docs.py")
        module = importlib.util.module_from_spec(spec)
        with mock.patch.object(sys, "argv", ["generate_docs"]):
            spec.loader.exec_module(module)
        module.REPO_ROOT = self.root
        with mock.patch.object(module.subprocess, "run") as run:
            for value in ("bad!\n", "", None):
                with self.subTest(value=value):
                    pin = self.root / "basilisk.ref"
                    pin.unlink(missing_ok=True)
                    if value is not None:
                        pin.write_text(value)
                    self.assertFalse(module.install_basilisk())
            run.assert_not_called()
            (self.root / "basilisk.ref").write_text(CURRENT + "\n")
            existing = self.root / "basilisk"
            existing.mkdir()
            for contents in (None, "ref=" + OLD + "\n", "ref=" + OLD + "\nref=" + CURRENT + "\n"):
                lock = existing / ".comphy-lock"
                lock.unlink(missing_ok=True)
                if contents is not None:
                    lock.write_text(contents)
                self.assertFalse(module.install_basilisk())
                run.assert_not_called()
            (existing / ".comphy-lock").write_text("ref=" + CURRENT + "\n")
            run.return_value = SimpleNamespace(returncode=0)
            self.assertTrue(module.install_basilisk())
            self.assertEqual(run.call_count, 3)
            self.assertIn(f"--ref={CURRENT}", run.call_args_list[-1].args[0])

    def test_source_commit_marks_dirty_and_failed_status_unverified(self):
        spec = importlib.util.spec_from_file_location("source_provenance", ROOT / "src-local/source_provenance.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        head = SimpleNamespace(returncode=0, stdout="abcdef\n")
        for status, expected in ((SimpleNamespace(returncode=0, stdout=""), "abcdef"),
                                 (SimpleNamespace(returncode=0, stdout=" M file\n"), "abcdef-dirty"),
                                 (SimpleNamespace(returncode=1, stdout=""), "abcdef-unverified")):
            with self.subTest(expected=expected), mock.patch.object(module.subprocess, "run", side_effect=[head, status]):
                self.assertEqual(module.source_commit(self.root), expected)
        with mock.patch.object(module.subprocess, "run", return_value=SimpleNamespace(returncode=1, stdout="")):
            self.assertEqual(module.source_commit(self.root), "unavailable")

    def test_contour_metadata_records_dirty_commit_and_selected_pin(self):
        import json
        source = ROOT / "src-local/source_provenance.py"
        shutil.copy(source, self.root / "src-local" / source.name)
        git = self.stub / "git"
        git.write_text("#!/bin/sh\ncase $1 in rev-parse) echo abcdef ;; status) echo ' M changed' ;; *) exit 2 ;; esac\n")
        git.chmod(0o755)
        attempt = self.root / "attempt"
        attempt.mkdir()
        for name in ("runContourHamilton.sbatch", "runContourSnellius.sbatch", "runContourLocal.sh"):
            with self.subTest(name=name):
                script = (ROOT / name).read_text()
                start = script.index("python3 - \"$iteration_dir\" \"$project_root\"")
                start = script.index("<<'PY'\n", start) + len("<<'PY'\n")
                end = script.index("\nPY\n", start)
                body = script[start:end]
                if name == "runContourLocal.sh":
                    args = [str(attempt), str(self.root), OLD, "1", "1", "1", "1", "", "0.021", "/qcc",
                            "0", "9", "10", "5", "-2.1", "5", "0", "0", "3", "1.5", ".25", ".25",
                            "0", "finite", "-1", "0", "1"]
                else:
                    args = [str(attempt), str(self.root), OLD, "1"]
                env = dict(self.env, SLURM_CPUS_PER_TASK="1")
                result = subprocess.run([sys.executable, "-c", body, *args], env=env,
                                        capture_output=True, text=True, timeout=10)
                self.assertEqual(result.returncode, 0, result.stderr)
                data = json.loads((attempt / "run-metadata.json").read_text())
                self.assertEqual(data["source_commit"], "abcdef-dirty")
                if name == "runContourLocal.sh":
                    self.assertEqual(data["qcc_lock_ref"], OLD)
                    self.assertEqual(data["basilisk_pin_file"], "basilisk-legacy.ref")
                else:
                    self.assertEqual(data["basilisk_pin_file"], "basilisk-legacy.ref")
                    self.assertEqual(data["basilisk_ref"], OLD)

    def test_docs_main_and_build_propagate_generation_failure(self):
        spec = importlib.util.spec_from_file_location("generate_docs", ROOT / ".github/scripts/generate_docs.py")
        module = importlib.util.module_from_spec(spec)
        with mock.patch.object(sys, "argv", ["generate_docs"]):
            spec.loader.exec_module(module)
        with mock.patch.object(module, "validate_config", return_value=False):
            with self.assertRaises(SystemExit) as raised:
                module.main()
            self.assertEqual(raised.exception.code, 1)
        module.DOCS_DIR = self.root / "docs-output"
        module.TEMPLATE_PATH = self.root / "template.html"
        with mock.patch.object(module, "validate_config", return_value=True), \
                mock.patch.object(module, "copy_assets", return_value=False):
            with self.assertRaises(SystemExit) as raised:
                module.main()
            self.assertEqual(raised.exception.code, 1)
        scripts = self.root / ".github" / "scripts"
        scripts.mkdir(parents=True)
        shutil.copy(ROOT / ".github/scripts/build.sh", scripts / "build.sh")
        (scripts / "generate_docs.py").write_text("raise SystemExit(7)\n")
        marker = self.root / "anchors-were-run"
        (scripts / "fix_empty_anchors.py").write_text(f"from pathlib import Path\nPath({str(marker)!r}).touch()\n")
        git = self.stub / "git"
        git.write_text("#!/bin/sh\nexit 1\n")
        git.chmod(0o755)
        result = subprocess.run(["bash", str(scripts / "build.sh")], cwd=self.root,
                                env=self.env, capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Documentation generation failed", result.stdout)
        self.assertFalse(marker.exists())
        self.assertFalse(self.curl_marker.exists())


if __name__ == "__main__":
    unittest.main()
