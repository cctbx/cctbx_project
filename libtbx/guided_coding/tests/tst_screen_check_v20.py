"""Controls for the claims that screen_check.py actually makes."""

import compileall
import contextlib
import importlib.util
import io
import os
import py_compile
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "payload" / "tools"))
import screen_check as checker


TREE = "b" * 40
BASE = "a" * 40
COMMIT = "c" * 40
REMOTE_URL = "ssh://git@example.invalid/phenix.git"
OUTGOING_KEYS = ("repository", "remote", "remote_url", "base", "commit", "tree", "refspec")


def outgoing(omit=None, **fields):
    """One OUTGOING.txt block (SPEC_A 1.2); `omit` drops a key, fields override."""
    values = dict(repository="phenix", remote="origin", remote_url=REMOTE_URL,
                  base=BASE, commit=COMMIT, tree=TREE)
    values.update(fields)
    values.setdefault("refspec", f"{values['commit']}:refs/heads/master")
    return "".join(f"{key}: {values[key]}\n" for key in OUTGOING_KEYS if key != omit)


class SourceInventoryChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.source = Path(self.temporary.name) / "source"
        tools = self.source / "payload" / "tools"
        tools.mkdir(parents=True)
        shutil.copy2(Path(checker.__file__), tools / "screen_check.py")
        shutil.copy2(Path(checker.__file__).with_name("review_bundle.py"),
                     tools / "review_bundle.py")
        (self.source / "SKILL.md").write_text("reviewed entry\n")
        names = ("SKILL.md", "payload/tools/review_bundle.py",
                 "payload/tools/screen_check.py")
        (self.source / "SOURCE_MANIFEST.sha256").write_text(
            "".join(f"{checker.digest((self.source / name).read_bytes())}  ./{name}\n"
                    for name in names))

    def run_source_check(self):
        return subprocess.run(
            [sys.executable, "-I", "-B", str(self.source / "payload/tools/screen_check.py"),
             "verify-source", str(self.source)],
            cwd=self.source, text=True, capture_output=True)

    def test_clean_source_passes_and_changed_source_fails(self):
        self.assertEqual(self.run_source_check().returncode, 0)
        (self.source / "SKILL.md").write_text("altered\n")
        result = self.run_source_check()
        self.assertEqual(result.returncode, 2)
        self.assertIn("changed source file", result.stderr)

    def test_unlisted_module_refused_without_executing_it(self):
        marker = self.source.parent / "executed"
        (self.source / "payload/tools/hashlib.py").write_text(
            f"from pathlib import Path\nPath({str(marker)!r}).write_text('executed')\n")
        result = self.run_source_check()
        self.assertEqual(result.returncode, 2)
        self.assertIn("extra source file", result.stderr)
        self.assertFalse(marker.exists())
        direct = subprocess.run(
            [sys.executable, "-B", str(self.source / "payload/tools/screen_check.py"),
             "--help"], cwd=self.source, text=True, capture_output=True)
        self.assertEqual(direct.returncode, 0)
        self.assertFalse(marker.exists())
        transport = subprocess.run(
            [sys.executable, "-B", str(self.source / "payload/tools/review_bundle.py"),
             "--help"], cwd=self.source, text=True, capture_output=True)
        self.assertEqual(transport.returncode, 0)
        self.assertFalse(marker.exists())

    def test_unlisted_regular_file_and_link_refused(self):
        extra = self.source / "unlisted.txt"
        extra.write_text("unlisted\n")
        self.assertIn("extra source file", self.run_source_check().stderr)
        extra.unlink()
        extra.symlink_to(self.source / "SKILL.md")
        self.assertIn("symbolic link", self.run_source_check().stderr)

    def precompile(self):
        # The real installer pattern: libtbx.py_compile_all -i calls
        # compileall.compile_dir(dir, 100, ...), which writes
        # __pycache__/<stem>.<cache tag>.pyc beside every module.
        self.assertTrue(compileall.compile_dir(str(self.source), 100, quiet=2))
        tag = sys.implementation.cache_tag
        self.assertEqual(
            sorted(p.relative_to(self.source).as_posix()
                   for p in self.source.rglob("*.pyc")),
            [f"payload/tools/__pycache__/review_bundle.{tag}.pyc",
             f"payload/tools/__pycache__/screen_check.{tag}.pyc"])
        return tag

    # A configured bytecode cache prefix (Apple's python3, PYTHONPYCACHEPREFIX)
    # moves bytecode out of __pycache__; these checks need the default layout.
    @unittest.skipIf(getattr(sys, "pycache_prefix", None), "bytecode cache prefix is set")
    def test_real_precompilation_accepted_and_reported(self):
        self.precompile()
        result = self.run_source_check()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("VERIFIED complete source", result.stdout)
        self.assertIn("ignored 2 precompiled bytecode file(s)", result.stdout)

    @unittest.skipIf(getattr(sys, "pycache_prefix", None), "bytecode cache prefix is set")
    def test_unlisted_code_still_refused_after_precompilation(self):
        tag = self.precompile()
        tools = self.source / "payload" / "tools"
        cached = tools / "__pycache__" / f"screen_check.{tag}.pyc"
        cases = {
            "unlisted module": (tools / "extra.py", b"print('extra')\n"),
            "bytecode for an unlisted module":
                (tools / "__pycache__" / f"extra.{tag}.pyc", cached.read_bytes()),
            "sourceless bytecode": (tools / "screen_helper.pyc", cached.read_bytes()),
            "non-bytecode file in __pycache__":
                (tools / "__pycache__" / "notes.txt", b"note\n"),
            "bytecode beside an unlisted stem":
                (self.source / "__pycache__" / f"SKILL.{tag}.pyc", cached.read_bytes()),
            "listed stem in another directory":
                (self.source / "__pycache__" / f"screen_check.{tag}.pyc", cached.read_bytes()),
        }
        for label, (path, data) in cases.items():
            with self.subTest(label):
                path.parent.mkdir(exist_ok=True)
                path.write_bytes(data)
                result = self.run_source_check()
                self.assertEqual(result.returncode, 2, label)
                self.assertIn("extra source file", result.stderr)
                path.unlink()
        with self.subTest("unchecked-hash bytecode for a listed module"):
            # Python runs unchecked-hash bytecode without consulting the source.
            py_compile.compile(str(tools / "screen_check.py"), cfile=str(cached),
                               invalidation_mode=py_compile.PycInvalidationMode.UNCHECKED_HASH)
            result = self.run_source_check()
            self.assertEqual(result.returncode, 2)
            self.assertIn("extra source file", result.stderr)
            cached.unlink()
        self.assertEqual(self.run_source_check().returncode, 0)

    @unittest.skipIf(getattr(sys, "pycache_prefix", None), "bytecode cache prefix is set")
    def test_tolerated_bytecode_is_never_loaded_by_the_tools(self):
        marker = self.source.parent / "cached-checker-ran"
        source = self.source / "payload" / "tools" / "screen_check.py"
        stat_result = source.stat()
        code = compile(f"from pathlib import Path\nPath({str(marker)!r}).write_text('ran')\n",
                       str(source), "exec")
        cache = Path(importlib.util.cache_from_source(str(source)))
        cache.parent.mkdir(exist_ok=True)
        # A header matching the source's mtime and size, as Python checks it.
        cache.write_bytes(importlib._bootstrap_external._code_to_timestamp_pyc(
            code, stat_result.st_mtime, stat_result.st_size))
        self.assertEqual(self.run_source_check().returncode, 0)
        for tool in ("screen_check.py", "review_bundle.py"):
            result = subprocess.run(
                [sys.executable, "-I", "-B", str(self.source / "payload/tools" / tool), "--help"],
                cwd=self.source, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(marker.exists())
        # Positive control: an ordinary cached import of the same file does
        # load this bytecode, so the check above is meaningful.
        subprocess.run(
            [sys.executable, "-B", "-c",
             "import importlib.util as u, sys; s = u.spec_from_file_location("
             "'screen_check', sys.argv[1]); s.loader.exec_module(u.module_from_spec(s))",
             str(source)], cwd=self.source, capture_output=True)
        self.assertTrue(marker.exists())


class SkillRegistrationChecks(unittest.TestCase):
    def setUp(self):
        SourceInventoryChecks.setUp(self)
        self.config = self.source.parent / "config"
        self.bin = self.source.parent / "bin"
        self.bin.mkdir()
        self.claude = self.bin / "claude"
        self.claude.write_text("#!/bin/sh\necho '2.1.284 (Claude Code)'\n")
        self.claude.chmod(0o755)

    def register(self, config=None):
        return subprocess.run(
            [sys.executable, "-I", "-B", str(self.source / "payload/tools/screen_check.py"),
             "register-skill", str(self.source)],
            env=dict(os.environ, PATH=str(self.bin),
                     CLAUDE_CONFIG_DIR=str(config or self.config)),
            text=True, capture_output=True)

    def test_fresh_registration_and_second_run_refused_without_source_change(self):
        first = self.register()
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertIn("/gc setup", first.stdout)
        self.assertIn("Registration has not adopted a project", first.stdout)
        self.assertFalse((self.source / ".claude").exists())
        link = self.config / "skills" / "guided_coding"
        self.assertTrue(link.is_symlink())
        self.assertEqual(link.resolve(), self.source.resolve())
        second = self.register()
        self.assertEqual(second.returncode, 2)
        self.assertIn("already exists", second.stderr)
        self.assertTrue(link.is_symlink())
        self.assertFalse((self.source / "guided_coding").exists())
        self.assertEqual(SourceInventoryChecks.run_source_check(self).returncode, 0)

    def test_other_link_directory_and_dangling_link_refused_without_nested_link(self):
        link = self.config / "skills" / "guided_coding"
        link.parent.mkdir(parents=True)
        other = self.source.parent / "other_package"
        other.mkdir()
        for kind in ("other_link", "directory", "dangling"):
            with self.subTest(kind=kind):
                if kind == "directory":
                    link.mkdir()
                else:
                    link.symlink_to(other if kind == "other_link" else
                                    self.source.parent / "not_here", target_is_directory=True)
                result = self.register()
                self.assertEqual(result.returncode, 2)
                self.assertIn("already exists", result.stderr)
                self.assertFalse((link / "guided_coding").exists())
                self.assertFalse((other / "guided_coding").exists())
                self.assertEqual(SourceInventoryChecks.run_source_check(self).returncode, 0)
                if kind == "directory":
                    link.rmdir()
                else:
                    link.unlink()

    def test_old_version_and_symlinked_skills_parent_refused(self):
        self.claude.write_text("#!/bin/sh\necho '2.1.268 (Claude Code)'\n")
        self.assertEqual(self.register().returncode, 2)
        self.assertFalse(self.config.exists())
        self.claude.write_text("#!/bin/sh\necho '2.1.284 (Claude Code)'\n")
        self.config.mkdir()
        other = self.source.parent / "other_skills"
        other.mkdir()
        (self.config / "skills").symlink_to(other, target_is_directory=True)
        result = self.register()
        self.assertEqual(result.returncode, 2)
        self.assertIn("symbolic link", result.stderr)
        self.assertFalse((other / "guided_coding").exists())

    def test_missing_prefix_dotdot_refused_before_any_source_or_parent_write(self):
        def snapshot(root):
            return {str(p.relative_to(root)):
                    (p.lstat().st_mode,
                     os.readlink(p) if p.is_symlink() else
                     p.read_bytes() if p.is_file() else None)
                    for p in root.rglob("*")}

        other = self.source.parent / "other_package"
        other.mkdir()
        for kind in ("fresh", "occupied", "config_link", "skills_link"):
            with self.subTest(kind=kind):
                endpoint = self.source.parent / f"outside_{kind}"
                prefix = self.source / f"created_{kind}"
                if kind == "occupied":
                    (endpoint / "skills").mkdir(parents=True)
                    (endpoint / "skills" / "guided_coding").symlink_to(self.source)
                elif kind == "config_link":
                    endpoint.symlink_to(other, target_is_directory=True)
                elif kind == "skills_link":
                    endpoint.mkdir()
                    (endpoint / "skills").symlink_to(other, target_is_directory=True)
                config = prefix / ".." / ".." / endpoint.name
                source_before = snapshot(self.source)
                other_before = snapshot(other)
                result = self.register(config)
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn("must not contain '..'", result.stderr)
                self.assertFalse(prefix.exists())
                self.assertEqual(source_before, snapshot(self.source))
                self.assertEqual(other_before, snapshot(other))
                if kind == "fresh":
                    self.assertFalse(endpoint.exists())
                if kind == "occupied":
                    self.assertTrue((endpoint / "skills" / "guided_coding").is_symlink())
                if kind == "config_link":
                    self.assertTrue(endpoint.is_symlink())
                if kind == "skills_link":
                    self.assertTrue((endpoint / "skills").is_symlink())


class ClaudeVersionChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.bin = Path(self.temporary.name)
        self.executable = self.bin / "claude"

    def run_gate(self, version=None, exit_code=0):
        if version is not None:
            self.executable.write_text(
                "#!/bin/sh\n" +
                f"printf '%s\\n' '{version}'\n" +
                f"exit {exit_code}\n")
            self.executable.chmod(0o755)
        environment = dict(os.environ, PATH=str(self.bin))
        return subprocess.run(
            [sys.executable, "-I", "-B", str(Path(checker.__file__)),
             "check-claude-version"],
            env=environment, text=True, capture_output=True)

    def test_minimum_and_current_mac_version(self):
        for version, allowed in (("2.1.268 (Claude Code)", False),
                                 ("2.1.280 (Claude Code)", False),
                                 ("2.1.281 (Claude Code)", True),
                                 ("2.1.284 (Claude Code)", True),
                                 ("2.2.0 (Claude Code)", True)):
            with self.subTest(version=version):
                result = self.run_gate(version)
                self.assertEqual(result.returncode, 0 if allowed else 2)
                self.assertIn("VERIFIED Claude Code CLI" if allowed else
                              "below 2.1.281", result.stdout if allowed else result.stderr)

    def test_missing_unreadable_or_failing_executable_refused(self):
        missing = self.run_gate()
        self.assertEqual(missing.returncode, 2)
        self.assertIn("not found on PATH", missing.stderr)
        malformed = self.run_gate("Claude Code unknown")
        self.assertEqual(malformed.returncode, 2)
        self.assertIn("unrecognized", malformed.stderr)
        failed = self.run_gate("2.1.284 (Claude Code)", exit_code=1)
        self.assertEqual(failed.returncode, 2)
        self.assertIn("exit 1", failed.stderr)


class ScreenChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.evidence = self.root / "evidence"
        self.evidence.mkdir()
        (self.evidence / "CODE_IDENTITY.txt").write_text(
            f"base_commit: {BASE}\ntested_tree: {TREE}\ncheck: python3 -m py_compile sample.py\n"
            "outcome: PASS\ncriterion_by: Developer\ntest_written_by: Worker\nrun_by: Worker\n"
        )
        (self.evidence / "RUN_LOG.txt").write_text("check exited 0\n")
        (self.evidence / "CHANGE.diff").write_text("+ # sample comment\n")
        (self.evidence / "APPROVAL_REPORT.md").write_text(
            "# Approval report\n## Exact tested changes\n```diff\n"
            "+ # sample comment\n```\n## New tests\nNone added\n"
        )
        (self.evidence / "SERVER_SUITE.txt").write_text("synthetic suite log\n")
        (self.evidence / "ROSTER_COMPARISON.txt").write_text("synthetic roster comparison\n")
        # SPEC_A 1.2: one valid block; its commit is named in the publication BATCH.
        (self.evidence / "OUTGOING.txt").write_text(outgoing())
        with contextlib.redirect_stdout(io.StringIO()):
            checker.freeze(self.evidence)
        self.identity = checker.verify(self.evidence)
        self.reading = self.root / "reading.txt"
        self.write_reading()

    def write_reading(self, verdict="PROCEED", scope="integration and publication",
                      identity=None, extra=""):
        """SPEC_A 1.1: a reading carries exactly one Scope line (None omits it)."""
        self.reading.write_text(
            f"Packet identity:  {identity or self.identity}\nVerdict:          {verdict}\n"
            + (f"Scope: {scope}\n" if scope is not None else "") + extra)

    def refreeze(self):
        """Freeze the evidence again after a fixture change and re-point the reading."""
        manifest = self.evidence / "MANIFEST.sha256"
        if manifest.exists():
            manifest.unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            checker.freeze(self.evidence)
        self.identity = checker.verify(self.evidence)
        self.write_reading()

    def put(self, body):
        path = self.root / "screen.md"
        path.write_text(body)
        return path

    def plan(self):
        return self.put("PLAN | comment-1 | light\nGOAL\nExplain the confusing comment.\n"
                        "Why it matters: Readers may misunderstand the explanation.\n"
                        "SCOPE\nOne comment in sample.py; approval authorizes a local draft.\n"
                        "CHECK\nInspect diff and compile; criterion by Developer.\n"
                        "LIMITS\nThe largest effect if wrong is a misleading explanation.\n"
                        "DECISION\nApprove draft, revise scope, or stop.\n"
                        "Recommendation: APPROVE because the change is limited to the comment.\n"
                        "Next: If approved, edit the comment and run the stated check.\n"
                        "ACTION: APPROVE / REVISE / STOP\n")

    def result(self, path="light", evidence_id=None, reading_id=None):
        evidence_id = evidence_id or self.identity
        if path == "light":
            gate = "Reading: deferred to publication\nVerdict: deferred to publication"
        else:
            reading_id = reading_id or checker.digest(self.reading.read_bytes())
            gate = f"Reading SHA256: {reading_id}\nVerdict: PROCEED"
        return self.put(f"RESULT | comment-1 | {path} | {evidence_id}\n"
                        "CHANGED\nClarified a comment in sample.py.\n"
                        "Why it matters: Readers can understand the comment more easily.\n"
                        "Full diff and exact new tests: evidence/APPROVAL_REPORT.md\n"
                        f"CHECKED\nTested tree: {TREE}\nThe log says compile passed.\n"
                        "LIMITS\nIt did not test runtime behavior.\n"
                        f"GATE\n{gate}\nDECISION\n"
                        "Integrate applies this tree; revise updates it; discard abandons it.\n"
                        "Recommendation: INTEGRATE because the check passed and scope is contained.\n"
                        "Next: If integrated, record the local change and close this task.\n"
                        "ACTION: INTEGRATE / REVISE / DISCARD\n")

    def stop(self):
        return self.put("STOP | comment-1\nWHY\nThe intended meaning is uncertain.\n"
                        "PRESERVED\nThe isolated draft is saved; no integration.\n"
                        "NEEDED\nDecide intended meaning, or cancel this change.\n"
                        "ACTION: DECIDE / CANCEL\n")

    def publication(self, batch=None, commit=COMMIT):
        reading_id = checker.digest(self.reading.read_bytes())
        if batch is None:
            batch = f"One local change to shared master.\nphenix commit {commit} to master."
        return self.put(f"PUBLICATION | batch-1 | {self.identity}\n"
                        f"BATCH\n{batch}\n"
                        "SERVER CHECK\nSuite and roster are in evidence.\n"
                        "LIMITS\nNo material gap claimed.\n"
                        f"GATE\nReading SHA256: {reading_id}\nVerdict: PROCEED\n"
                        "DECISION\nPublish this batch, or hold locally.\n"
                        "ACTION: PUBLISH / HOLD\n")

    def ok(self, kind, path, reading=None, evidence=None, disposition=None):
        with contextlib.redirect_stdout(io.StringIO()):
            checker.check(kind, path, evidence, reading, disposition)

    def test_four_valid_shapes(self):
        self.ok("plan", self.plan())
        self.ok("result", self.result(), evidence=self.evidence)
        self.ok("result", self.result("full"), self.reading, self.evidence)
        self.ok("stop", self.stop())
        self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_present_checks_screen_and_hides_internal_identifiers(self):
        for kind, make_screen, reading, evidence in (
            ("plan", self.plan, None, None),
            ("result", self.result, None, self.evidence),
            ("result", lambda: self.result("full"), self.reading, self.evidence),
            ("stop", self.stop, None, None),
            ("publication", self.publication, self.reading, self.evidence),
        ):
            with self.subTest(kind=kind, reading=bool(reading)):
                path = make_screen()
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    checker.present(kind, path, evidence, reading)
                display = output.getvalue()
                self.assertTrue(display.startswith(
                    "**PLEASE READ — YOUR DECISION IS NEEDED**\n\n```\n"))
                self.assertTrue(display.endswith("\n```\n"))
                self.assertIn("\n" + kind.upper() + " | ", display)
                self.assertNotIn(self.identity, display)
                self.assertNotIn(TREE, display)
                self.assertNotIn(COMMIT, display)  # SPEC_A 1.2: shown as [recorded]
                self.assertNotIn(checker.digest(self.reading.read_bytes()), display)
                self.assertNotIn("SHA256:", display)
                self.assertIn("ACTION:", display)
        screen = self.result()
        screen.write_text(screen.read_text().replace(
            "The log says compile passed.",
            "The log says compile passed; commit 1faab028cc, tree fe5a57cb…614a."))
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            checker.present("result", screen, self.evidence)
        self.assertNotIn("1faab028cc", output.getvalue())
        self.assertNotIn("fe5a57cb", output.getvalue())
        with self.assertRaisesRegex(ValueError, "does not identify"):
            checker.present("result", self.result("full", reading_id="0" * 64),
                            self.evidence, self.reading)

    def test_plan_display_has_blank_line_before_each_section(self):
        screen = self.plan()
        saved = screen.read_bytes()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            checker.present("plan", screen)
        display = output.getvalue()
        for heading in checker.HEADINGS["plan"]:
            self.assertIn("\n\n" + heading + "\n", display)
        self.assertEqual(screen.read_bytes(), saved)
        self.ok("plan", screen)

    def test_present_keeps_code_fences_inside_the_screen_body(self):
        screen = self.plan()
        screen.write_text(screen.read_text().replace(
            "One comment in sample.py; approval authorizes a local draft.",
            "One comment in sample.py; approval authorizes a local draft.\n"
            "Example: ```python should remain in the screen body."))
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            checker.present("plan", screen)
        display = output.getvalue()
        self.assertIn("**PLEASE READ — YOUR DECISION IS NEEDED**\n\n````\n", display)
        self.assertTrue(display.endswith("\n````\n"))
        self.assertIn("Example: ```python should remain in the screen body.", display)

    def test_headings_action_and_placeholders(self):
        source = self.plan().read_text()
        for altered, reason in (
            (source.replace("CHECK\n", "TEST\n"), "CHECK heading"),
            (source.replace("GOAL\n", "SCOPE\nGOAL\n"), "SCOPE heading"),
            (source + "unchecked postscript\n", "last line"),
            (source.replace("One comment", "<one comment>"), "placeholder"),
            (source.replace("Explain the confusing comment.", "many " * 230), "one-page"),
        ):
            with self.subTest(reason=reason):
                with self.assertRaisesRegex(ValueError, reason):
                    self.ok("plan", self.put(altered))

    def test_plan_and_result_require_useful_decision_fields(self):
        for kind, original in (("plan", self.plan().read_text()),
                               ("result", self.result().read_text())):
            for prefix, reason in (("Why it matters:", "Why it matters"),
                                   ("Recommendation:", "recommendation and reason"),
                                   ("Next:", "next action")):
                with self.subTest(kind=kind, missing=prefix):
                    lines = original.splitlines(keepends=True)
                    changed = "".join(line for line in lines if not line.startswith(prefix))
                    with self.assertRaisesRegex(ValueError, reason):
                        self.ok(kind, self.put(changed), evidence=self.evidence if kind == "result" else None)

    def test_manifest_rejects_change_or_extra(self):
        (self.evidence / "RUN_LOG.txt").write_text("check failed\n")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            checker.verify(self.evidence)
        (self.evidence / "RUN_LOG.txt").write_text("check exited 0\n")
        (self.evidence / "unlisted.txt").write_text("extra")
        with self.assertRaisesRegex(ValueError, "missing or extra"):
            checker.verify(self.evidence)

    def test_no_self_hash_overwrite_or_missing_evidence(self):
        with self.assertRaisesRegex(ValueError, "already exists"):
            checker.freeze(self.evidence)
        (self.evidence / "MANIFEST.sha256").write_text(
            self.evidence.joinpath("MANIFEST.sha256").read_text() + "0" * 64 + "  MANIFEST.sha256\n"
        )
        with self.assertRaisesRegex(ValueError, "self-listed"):
            checker.verify(self.evidence)
        (self.evidence / "MANIFEST.sha256").unlink()
        with self.assertRaisesRegex(ValueError, "missing manifest"):
            checker.verify(self.evidence)

    def test_rejects_symlink_and_hardlink(self):
        link = self.evidence / "link"
        link.symlink_to(self.evidence / "RUN_LOG.txt")
        with self.assertRaisesRegex(ValueError, "symbolic link"):
            checker.verify(self.evidence)
        link.unlink()
        os.link(self.evidence / "RUN_LOG.txt", link)
        with self.assertRaisesRegex(ValueError, "linked file"):
            checker.verify(self.evidence)

    def test_freeze_refuses_symlinked_parent_before_writing(self):
        outside = self.root / "outside"
        target = outside / "fresh"
        target.mkdir(parents=True)
        (target / "RUN_LOG.txt").write_text("some evidence\n")
        link = self.root / "parent_link"
        link.symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "symbolic link in evidence path"):
            checker.freeze(link / "fresh")
        self.assertFalse((target / "MANIFEST.sha256").exists())

    def test_wrong_packet_or_tested_tree(self):
        with self.assertRaisesRegex(ValueError, "current evidence"):
            self.ok("result", self.result(evidence_id="0" * 64), evidence=self.evidence)
        screen = self.result().read_text().replace(TREE, "0" * 40)
        with self.assertRaisesRegex(ValueError, "tested tree"):
            self.ok("result", self.put(screen), evidence=self.evidence)

    def test_reading_identity_and_verdict_are_required(self):
        self.write_reading(identity="0" * 64)
        with self.assertRaisesRegex(ValueError, "does not name"):
            self.ok("result", self.result("full"), self.reading, self.evidence)
        self.write_reading("PROCEED IF repaired")
        with self.assertRaisesRegex(ValueError, "conditional verdict"):
            self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_conditional_reading_links_exact_disposition(self):
        self.write_reading("PROCEED IF suite matches")
        disposition = self.root / "disposition.txt"
        disposition.write_text("Condition: suite matches\nStatus: SATISFIED\nEvidence: roster review\n")
        gate = (f"Verdict: PROCEED IF suite matches\nDisposition SHA256: "
                f"{checker.digest(disposition.read_bytes())}")
        screen = self.put(self.result("full").read_text().replace("Verdict: PROCEED", gate))
        self.ok("result", screen, self.reading, self.evidence, disposition)
        disposition.write_text("Condition: a different condition\nStatus: SATISFIED\nEvidence: roster\n")
        with self.assertRaisesRegex(ValueError, "exact condition"):
            self.ok("result", screen, self.reading, self.evidence, disposition)

    def test_waiver_and_pending_final_decision_are_truthful_and_visible(self):
        self.write_reading("PROCEED IF a fresh reading is done")
        note = self.root / "disposition.txt"
        note.write_text(
            'Condition: a fresh reading is done\nStatus: WAIVED\n'
            'Developer authorization: "I waive a second reading"\n'
            'Evidence: user decision recorded in change record\n'
        )
        gate = ("Verdict: PROCEED IF a fresh reading is done\nDisposition SHA256: "
                f"{checker.digest(note.read_bytes())}")
        result = self.result("full")
        result.write_text(result.read_text().replace("Verdict: PROCEED", gate)
                          .replace("It did not test runtime behavior.",
                                   "It did not test runtime behavior. WAIVED second reading by Developer."))
        self.ok("result", result, self.reading, self.evidence, note)
        note.write_text(note.read_text().replace(
            'Developer authorization: "I waive a second reading"\n', ""))
        with self.assertRaisesRegex(ValueError, "quoted authorization"):
            self.ok("result", result, self.reading, self.evidence, note)

        self.write_reading("PROCEED IF Developer chooses PUBLISH")
        note.write_text("Condition: Developer chooses PUBLISH\nStatus: PENDING\n"
                        "Pending: Developer PUBLISH or HOLD\nEvidence: other findings resolved\n")
        gate = ("Verdict: PROCEED IF Developer chooses PUBLISH\nDisposition SHA256: "
                f"{checker.digest(note.read_bytes())}")
        publication = self.publication()
        publication.write_text(publication.read_text().replace("Verdict: PROCEED", gate)
                               .replace("No material gap claimed.",
                                        "PENDING Developer choice; nothing pushed."))
        self.ok("publication", publication, self.reading, self.evidence, note)
        publication.write_text(publication.read_text().replace(
            "PENDING Developer choice; nothing pushed.", "Nothing pushed."))
        with self.assertRaisesRegex(ValueError, "visible in Developer View"):
            self.ok("publication", publication, self.reading, self.evidence, note)

    def test_result_needs_exact_report_diff_and_usable_path(self):
        report = self.evidence / "APPROVAL_REPORT.md"
        original = report.read_text()
        report.write_text(original.replace("+ # sample comment", "+ # other comment"))
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            self.ok("result", self.result(), evidence=self.evidence)
        (self.evidence / "MANIFEST.sha256").unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            checker.freeze(self.evidence)
        self.identity = checker.verify(self.evidence)
        with self.assertRaisesRegex(ValueError, "complete exact CHANGE.diff"):
            self.ok("result", self.result(), evidence=self.evidence)
        report.write_text(original)
        (self.evidence / "MANIFEST.sha256").unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            checker.freeze(self.evidence)
        self.identity = checker.verify(self.evidence)
        screen = self.result()
        screen.write_text(screen.read_text().replace(
            "Full diff and exact new tests: evidence/APPROVAL_REPORT.md\n", ""))
        with self.assertRaisesRegex(ValueError, "point to the complete approval report"):
            self.ok("result", screen, evidence=self.evidence)

    def test_publication_requires_named_logs(self):
        (self.evidence / "SERVER_SUITE.txt").unlink()
        self.refreeze()
        with self.assertRaisesRegex(ValueError, "missing or empty publication evidence"):
            self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_stale_reading_sha_and_light_reading(self):
        with self.assertRaisesRegex(ValueError, "does not identify"):
            self.ok("result", self.result("full", reading_id="0" * 64), self.reading, self.evidence)
        with self.assertRaisesRegex(ValueError, "belongs to publication"):
            self.ok("result", self.result(), self.reading, self.evidence)

    def test_failed_record_cannot_be_a_result(self):
        identity_file = self.evidence / "CODE_IDENTITY.txt"
        identity_file.write_text(identity_file.read_text().replace("outcome: PASS", "outcome: FAIL"))
        # A newly frozen packet still cannot present the failed run as a passing result.
        (self.evidence / "MANIFEST.sha256").unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            checker.freeze(self.evidence)
        with self.assertRaisesRegex(ValueError, "requires a passing"):
            self.ok("result", self.result(evidence_id=checker.verify(self.evidence)), evidence=self.evidence)

    def test_passing_record_can_explain_the_result(self):
        identity_file = self.evidence / "CODE_IDENTITY.txt"
        identity_file.write_text(identity_file.read_text().replace(
            "outcome: PASS\n", "outcome: PASS (details in RUN_LOG.txt)\n"))
        (self.evidence / "MANIFEST.sha256").unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            checker.freeze(self.evidence)
        self.ok("result", self.result(evidence_id=checker.verify(self.evidence)), evidence=self.evidence)

    def test_counted_success_does_not_require_pass_prefix(self):
        identity_file = self.evidence / "CODE_IDENTITY.txt"
        original = identity_file.read_text()
        for outcome, passes in (
                ("before 4/4 exit 0; after 4/4 exit 0", True),
                ("4/4 exit 0", True),
                ("before 3/4 exit 0; after 4/4 exit 0", False),
                ("before 4/4 exit 0; after 4/4 exit 1", False),
                ("before 4/4 exit 0; FAILED", False)):
            with self.subTest(outcome=outcome):
                identity_file.write_text(original.replace("outcome: PASS", f"outcome: {outcome}"))
                (self.evidence / "MANIFEST.sha256").unlink()
                with contextlib.redirect_stdout(io.StringIO()):
                    checker.freeze(self.evidence)
                self.identity = checker.verify(self.evidence)
                screen = self.result(evidence_id=self.identity)
                if passes:
                    self.ok("result", screen, evidence=self.evidence)
                else:
                    with self.assertRaisesRegex(ValueError, "requires a passing"):
                        self.ok("result", screen, evidence=self.evidence)

    # ---- SPEC_A 1.1: reading scope ------------------------------------------

    def test_reading_needs_exactly_one_scope_line(self):
        """SPEC_A 1.1, counted as in SPEC_A 10: zero or several Scope lines are
        refused on the full RESULT and on the PUBLICATION screen with "exactly
        one Scope line"; a single line with an unrecognized value is refused
        with "Scope value not recognized"."""
        for label, scope_lines in (
                ("none", ""),
                ("two different", "Scope: integration\nScope: publication\n"),
                ("two identical", "Scope: integration and publication\n" * 2)):
            with self.subTest(scope=label):
                self.write_reading(scope=None, extra=scope_lines)
                with self.assertRaisesRegex(ValueError, "exactly one Scope line"):
                    self.ok("result", self.result("full"), self.reading, self.evidence)
                with self.assertRaisesRegex(ValueError, "exactly one Scope line"):
                    self.ok("publication", self.publication(), self.reading, self.evidence)
        for label, scope_lines in (
                ("unrecognized value", "Scope: everything\n"),
                ("trailing words", "Scope: integration and publication and more\n")):
            with self.subTest(scope=label):
                self.write_reading(scope=None, extra=scope_lines)
                with self.assertRaisesRegex(ValueError, "Scope value not recognized"):
                    self.ok("result", self.result("full"), self.reading, self.evidence)
                with self.assertRaisesRegex(ValueError, "Scope value not recognized"):
                    self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_publication_only_reading_does_not_cover_a_result(self):
        """SPEC_A 1.1: kind result needs scope integration or integration and
        publication; a publication-only reading still serves PUBLICATION."""
        self.write_reading(scope="publication")
        with self.assertRaisesRegex(ValueError, "does not cover integration"):
            self.ok("result", self.result("full"), self.reading, self.evidence)
        self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_integration_only_reading_does_not_cover_a_publication(self):
        """SPEC_A 1.1: kind publication needs scope publication or integration
        and publication; an integration-only reading still serves RESULT."""
        self.write_reading(scope="integration")
        with self.assertRaisesRegex(ValueError, "does not cover publication"):
            self.ok("publication", self.publication(), self.reading, self.evidence)
        self.ok("result", self.result("full"), self.reading, self.evidence)

    def test_dual_scope_reading_serves_both_screens_for_the_same_packet(self):
        """SPEC_A 1.1: `integration and publication` is accepted by RESULT and
        PUBLICATION for the same packet; surrounding whitespace and a
        Supersedes reading SHA256 line are allowed."""
        for label, scope, extra in (
                ("plain", "integration and publication", ""),
                ("padded", "   integration and publication   ", ""),
                ("supersedes", "integration and publication",
                 f"Supersedes reading SHA256: {'1' * 64}\n")):
            with self.subTest(form=label):
                self.write_reading(scope=scope, extra=extra)
                self.ok("result", self.result("full"), self.reading, self.evidence)
                self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_light_result_is_unchanged_by_the_scope_rule(self):
        """SPEC_A 1.1: a light RESULT defers the reading; it takes no reading
        file and does not care what reading is on disk."""
        self.ok("result", self.result(), evidence=self.evidence)
        self.write_reading(scope=None)
        self.ok("result", self.result(), evidence=self.evidence)
        self.write_reading(scope="publication")
        self.ok("result", self.result(), evidence=self.evidence)
        self.write_reading()
        with self.assertRaisesRegex(ValueError, "belongs to publication"):
            self.ok("result", self.result(), self.reading, self.evidence)

    def test_altered_reading_after_screens_written_fails_both_gates(self):
        """SPEC_A 1.1 (GATE hashes unchanged): one added byte in an accepted
        dual-scope reading makes RESULT and PUBLICATION fail at the GATE."""
        result = self.result("full").rename(self.root / "result.md")
        publication = self.publication().rename(self.root / "publication.md")
        self.ok("result", result, self.reading, self.evidence)
        self.ok("publication", publication, self.reading, self.evidence)
        self.reading.write_bytes(self.reading.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "does not identify"):
            self.ok("result", result, self.reading, self.evidence)
        with self.assertRaisesRegex(ValueError, "does not identify"):
            self.ok("publication", publication, self.reading, self.evidence)

    # ---- SPEC_A 1.2: OUTGOING.txt (PUBLICATION only) ------------------------

    def publication_refused(self, reason, **screen):
        with self.assertRaisesRegex(ValueError, reason):
            self.ok("publication", self.publication(**screen), self.reading, self.evidence)

    def test_publication_requires_a_non_empty_outgoing_record(self):
        """SPEC_A 1.2: a missing or empty OUTGOING.txt refuses PUBLICATION;
        RESULT does not need it."""
        record = self.evidence / "OUTGOING.txt"
        record.unlink()
        self.refreeze()
        self.publication_refused("missing or empty publication evidence: OUTGOING.txt")
        self.ok("result", self.result("full"), self.reading, self.evidence)
        record.write_text("")
        self.refreeze()
        self.publication_refused("missing or empty publication evidence: OUTGOING.txt")

    def test_outgoing_record_must_be_utf8(self):
        """SPEC_A 1.2: a non-UTF-8 OUTGOING.txt is refused."""
        (self.evidence / "OUTGOING.txt").write_bytes(b"repository: phenix\n\xff\xfe\n")
        self.refreeze()
        self.publication_refused("outgoing record is not UTF-8")

    def test_outgoing_record_needs_every_key_in_every_block(self):
        """SPEC_A 1.2: a block missing any of the seven keys, or a record with
        no repository block at all, is refused."""
        needs = ("outgoing record needs repository, remote, remote_url, base, "
                 "commit, tree and refspec")
        for key in OUTGOING_KEYS:
            with self.subTest(missing=key):
                (self.evidence / "OUTGOING.txt").write_text(outgoing(omit=key))
                self.refreeze()
                self.publication_refused(needs)
        with self.subTest(missing="second block's tree"):
            (self.evidence / "OUTGOING.txt").write_text(
                outgoing() + outgoing(omit="tree", repository="cctbx"))
            self.refreeze()
            self.publication_refused(needs)
        with self.subTest(missing="any block"):
            (self.evidence / "OUTGOING.txt").write_text("# comments only\n\n")
            self.refreeze()
            self.publication_refused(needs)

    def test_outgoing_record_ignores_comments_and_blank_lines(self):
        """SPEC_A 1.2: blank lines and lines starting with # are not entries."""
        (self.evidence / "OUTGOING.txt").write_text(
            "# header comment\n\n" + outgoing().replace("base:", "# note\n\nbase:")
            + "\n# trailing comment\n")
        self.refreeze()
        self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_outgoing_identities_must_be_forty_lowercase_hex(self):
        """SPEC_A 1.2: base, commit and tree are exactly 40 lowercase hex
        characters (the refspec and BATCH are kept consistent so only this
        rule is violated)."""
        for field, value in (("base", "A" * 40), ("tree", "b" * 39),
                             ("commit", "c" * 41), ("base", "g" * 40),
                             ("tree", "B" * 40)):
            with self.subTest(field=field, value=value):
                (self.evidence / "OUTGOING.txt").write_text(outgoing(**{field: value}))
                self.refreeze()
                self.publication_refused("invalid Git identity in outgoing record",
                                         commit=value if field == "commit" else COMMIT)

    def test_outgoing_refspec_must_name_the_commit_and_a_branch(self):
        """SPEC_A 1.2: refspec is <commit>:refs/heads/<branch> with the block's
        own commit, a [A-Za-z0-9._/-]+ branch and no *."""
        for refspec in (f"{'d' * 40}:refs/heads/master", f"{COMMIT}:refs/heads/*",
                        f"{COMMIT}:refs/heads/feature*", f"{COMMIT}:refs/tags/v1",
                        f"{COMMIT}:master", f"{COMMIT}:refs/heads/",
                        "HEAD:refs/heads/master", f"{COMMIT}:refs/heads/a b",
                        f"{COMMIT}", "refs/heads/master"):
            with self.subTest(refspec=refspec):
                (self.evidence / "OUTGOING.txt").write_text(outgoing(refspec=refspec))
                self.refreeze()
                self.publication_refused("outgoing refspec must be")
        (self.evidence / "OUTGOING.txt").write_text(
            outgoing(refspec=f"{COMMIT}:refs/heads/release/1.0-rc_2"))
        self.refreeze()
        self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_outgoing_commit_must_be_a_whole_word_in_batch(self):
        """SPEC_A 1.2: the block's commit appears as a whole word on some BATCH
        line; a longer token or another section does not count."""
        for label, batch in (
                ("absent", "One local change to shared master.\nNo identifiers here."),
                ("inside a longer token",
                 f"One local change.\nphenix commit {COMMIT}0 to master."),
                ("prefixed", f"One local change.\nphenix commit x{COMMIT} to master.")):
            with self.subTest(batch=label):
                self.publication_refused("outgoing commit is not named in BATCH", batch=batch)
        with self.subTest(batch="named only in SERVER CHECK"):
            screen = self.publication(batch="One local change to shared master.")
            screen.write_text(screen.read_text().replace(
                "Suite and roster are in evidence.",
                f"Suite and roster are in evidence for {COMMIT}."))
            with self.assertRaisesRegex(ValueError, "outgoing commit is not named in BATCH"):
                self.ok("publication", screen, self.reading, self.evidence)
        with self.subTest(batch="followed by punctuation"):
            self.ok("publication", self.publication(
                batch=f"One local change to shared master.\nCommit {COMMIT}."),
                self.reading, self.evidence)

    def test_outgoing_record_accepts_several_blocks_and_refuses_duplicates(self):
        """SPEC_A 1.2: with two repositories both commits must be in BATCH; a
        repeated repository name is refused."""
        second = "d" * 40
        (self.evidence / "OUTGOING.txt").write_text(
            outgoing() + "\n" + outgoing(repository="cctbx", commit=second,
                                         remote_url="ssh://git@example.invalid/cctbx.git"))
        self.refreeze()
        self.ok("publication", self.publication(
            batch=f"Two changes.\nphenix commit {COMMIT} and cctbx commit {second}."),
            self.reading, self.evidence)
        self.publication_refused("outgoing commit is not named in BATCH")
        (self.evidence / "OUTGOING.txt").write_text(outgoing() + "\n" + outgoing())
        self.refreeze()
        self.publication_refused("duplicate repository in outgoing record")

    # ---- SPEC_A 1.3: suite waiver form (PUBLICATION only) -------------------

    def suite(self, text):
        (self.evidence / "SERVER_SUITE.txt").write_text(text)
        self.refreeze()

    def test_suite_not_run_requires_quoted_developer_waiver(self):
        """SPEC_A 1.3: a first line of exactly SERVER_SUITE: NOT RUN needs a
        `Waiver (Developer...):` line immediately followed by `> ` lines."""
        reason = "suite not run without the Developer's quoted waiver"
        for label, text in (
                ("no waiver", "SERVER_SUITE: NOT RUN\nReason: no server reachable.\n"),
                ("header then unquoted line",
                 "SERVER_SUITE: NOT RUN\nWaiver (Developer, 2026-10-04, batch x):\n"
                 "I waive the suite for this batch.\n"),
                ("header then blank then quoted",
                 "SERVER_SUITE: NOT RUN\nWaiver (Developer, 2026-10-04, batch x):\n\n"
                 "> quoted words\n"),
                ("header at end of file",
                 "SERVER_SUITE: NOT RUN\nWaiver (Developer, 2026-10-04, batch x):\n"),
                ("quote without header", "SERVER_SUITE: NOT RUN\n> quoted words\n"),
                ("header not by the Developer",
                 "SERVER_SUITE: NOT RUN\nWaiver (Worker, 2026-10-04):\n> quoted words\n"),
                ("header with trailing text",
                 "SERVER_SUITE: NOT RUN\nWaiver (Developer): see below\n> quoted words\n"),
                ("surrounding whitespace is stripped",
                 "  SERVER_SUITE: NOT RUN  \nReason: offline.\n")):
            with self.subTest(form=label):
                self.suite(text)
                self.publication_refused(reason)

    def test_suite_not_run_with_quoted_waiver_is_accepted(self):
        """SPEC_A 1.3: the Waiver (Developer...) header immediately followed by
        one or more `> ` lines with words passes the form check."""
        for label, text in (
                ("minimal", "SERVER_SUITE: NOT RUN\nWaiver (Developer):\n"
                            "> I waive the suite for this batch.\n"),
                ("dated with two quoted lines",
                 "SERVER_SUITE: NOT RUN\nWaiver (Developer, 2026-10-04, batch x):\n"
                 "> quoted words\n> more quoted words\n"),
                ("waiver later in the file",
                 "SERVER_SUITE: NOT RUN\nReason: machine offline.\n\n"
                 "Waiver (Developer, 2026-10-04, batch x):\n> quoted words\n"
                 "Note: recorded by the Guide.\n")):
            with self.subTest(form=label):
                self.suite(text)
                self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_suite_form_check_applies_only_to_not_run_and_only_to_publication(self):
        """SPEC_A 1.3 (prefix rule, section 5): a first line that does not start
        with SERVER_SUITE: NOT RUN needs no waiver; one that starts with it,
        whatever the suffix, needs the quoted waiver; RESULT never applies the
        check."""
        for text in ("synthetic suite log\n", "SERVER_SUITE: PASS 12/12\n",
                     "Summary\nSERVER_SUITE: NOT RUN\n"):
            with self.subTest(first_line=text.splitlines()[0]):
                self.suite(text)
                self.ok("publication", self.publication(), self.reading, self.evidence)
        with self.subTest(first_line="SERVER_SUITE: NOT RUN (see below)", waiver=False):
            self.suite("SERVER_SUITE: NOT RUN (see below)\nReason: offline.\n")
            self.publication_refused("suite not run without the Developer's quoted waiver")
        with self.subTest(first_line="SERVER_SUITE: NOT RUN (see below)", waiver=True):
            self.suite("SERVER_SUITE: NOT RUN (see below)\n"
                       "Waiver (Developer, 2026-10-04, batch x):\n> quoted words\n")
            self.ok("publication", self.publication(), self.reading, self.evidence)
        self.suite("SERVER_SUITE: NOT RUN\nNo waiver here.\n")
        self.ok("result", self.result("full"), self.reading, self.evidence)

    # ---- SPEC_A 6: byte-order mark before SERVER_SUITE: NOT RUN -------------

    def suite_bytes(self, data):
        (self.evidence / "SERVER_SUITE.txt").write_bytes(data)
        self.assertEqual((self.evidence / "SERVER_SUITE.txt").read_bytes()[:3], b"\xef\xbb\xbf")
        self.refreeze()

    def test_bom_before_not_run_without_waiver_is_refused(self):
        """SPEC_A 6 (1.3): a UTF-8 byte-order mark immediately before
        `SERVER_SUITE: NOT RUN` does not evade the waiver rule: with no waiver
        block the PUBLICATION is refused."""
        reason = "suite not run without the Developer's quoted waiver"
        for label, data in (
                ("BOM then NOT RUN only", b"\xef\xbb\xbfSERVER_SUITE: NOT RUN\n"),
                ("BOM then NOT RUN and a reason",
                 b"\xef\xbb\xbfSERVER_SUITE: NOT RUN\nReason: no server reachable.\n"),
                ("BOM then NOT RUN with a suffix",
                 b"\xef\xbb\xbfSERVER_SUITE: NOT RUN (see below)\nReason: offline.\n"),
                ("BOM then NOT RUN with an unquoted waiver",
                 b"\xef\xbb\xbfSERVER_SUITE: NOT RUN\nWaiver (Developer, 2026-10-04, batch x):\n"
                 b"I waive the suite for this batch.\n")):
            with self.subTest(form=label):
                self.suite_bytes(data)
                self.publication_refused(reason)

    def test_bom_before_not_run_with_quoted_waiver_is_accepted(self):
        """SPEC_A 6 (1.3): the same BOM-prefixed file with a `Waiver
        (Developer, ...):` line immediately followed by a `> ` quoted line
        passes the form check."""
        self.suite_bytes(b"\xef\xbb\xbfSERVER_SUITE: NOT RUN\n"
                         b"Waiver (Developer, 2026-10-04, batch x):\n> quoted words\n")
        self.ok("publication", self.publication(), self.reading, self.evidence)
        self.suite_bytes(b"\xef\xbb\xbfSERVER_SUITE: NOT RUN\nReason: offline.\n\n"
                         b"Waiver (Developer, 2026-10-04, batch x):\n> quoted words\n"
                         b"> more quoted words\n")
        self.ok("publication", self.publication(), self.reading, self.evidence)

    # ---- SPEC_A 10: Scope lines are counted first; no repeated OUTGOING key --

    def test_extra_scope_line_with_any_wording_is_refused(self):
        """SPEC_A 10 (1.1): all lines beginning `Scope:` are counted first, so
        a valid line followed by `Scope: integration only` is refused at RESULT
        and at PUBLICATION with "exactly one Scope line"; a single `Scope:
        integration only` is refused with "Scope value not recognized"; each
        single recognized value still serves its screen(s)."""
        self.write_reading(scope="integration and publication",
                           extra="Scope: integration only\n")
        with self.assertRaisesRegex(ValueError, "exactly one Scope line"):
            self.ok("result", self.result("full"), self.reading, self.evidence)
        with self.assertRaisesRegex(ValueError, "exactly one Scope line"):
            self.ok("publication", self.publication(), self.reading, self.evidence)
        self.write_reading(scope="integration only")
        with self.assertRaisesRegex(ValueError, "Scope value not recognized"):
            self.ok("result", self.result("full"), self.reading, self.evidence)
        with self.assertRaisesRegex(ValueError, "Scope value not recognized"):
            self.ok("publication", self.publication(), self.reading, self.evidence)
        for scope, kinds in (("integration", ("result",)),
                             ("publication", ("publication",)),
                             ("integration and publication", ("result", "publication"))):
            with self.subTest(scope=scope):
                self.write_reading(scope=scope)
                for kind in kinds:
                    screen = self.result("full") if kind == "result" else self.publication()
                    self.ok(kind, screen, self.reading, self.evidence)

    def test_duplicate_key_within_an_outgoing_block_is_refused(self):
        """SPEC_A 10 (1.2): a key repeated within one repository block, required
        (remote_url: first a different destination, then the correct one) or
        optional (companions twice), is refused at PUBLICATION with "duplicate
        key in outgoing record: <key>"; a single optional key and the plain
        single block still pass."""
        record = self.evidence / "OUTGOING.txt"
        with self.subTest(key="remote_url"):
            record.write_text(outgoing().replace(
                "remote_url:", "remote_url: ssh://git@elsewhere.invalid/other.git\nremote_url:", 1))
            self.refreeze()
            self.publication_refused("duplicate key in outgoing record: remote_url")
        with self.subTest(key="companions"):
            record.write_text(outgoing() + "companions: cctbx\ncompanions: dxtbx\n")
            self.refreeze()
            self.publication_refused("duplicate key in outgoing record: companions")
        with self.subTest(key="single optional key"):
            record.write_text(outgoing() + "companions: cctbx\n")
            self.refreeze()
            self.ok("publication", self.publication(), self.reading, self.evidence)
        with self.subTest(key="plain single block"):
            record.write_text(outgoing())
            self.refreeze()
            self.ok("publication", self.publication(), self.reading, self.evidence)


if __name__ == "__main__":
    unittest.main()
