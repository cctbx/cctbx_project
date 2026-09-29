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
        with contextlib.redirect_stdout(io.StringIO()):
            checker.freeze(self.evidence)
        self.identity = checker.verify(self.evidence)
        self.reading = self.root / "reading.txt"
        self.reading.write_text(
            f"Packet identity:  {self.identity}\nVerdict:          PROCEED\n"
        )

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

    def publication(self):
        reading_id = checker.digest(self.reading.read_bytes())
        return self.put(f"PUBLICATION | batch-1 | {self.identity}\n"
                        "BATCH\nOne local change to shared master.\n"
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
        self.reading.write_text(f"Packet identity: {'0' * 64}\nVerdict: PROCEED\n")
        with self.assertRaisesRegex(ValueError, "does not name"):
            self.ok("result", self.result("full"), self.reading, self.evidence)
        self.reading.write_text(f"Packet identity: {self.identity}\nVerdict: PROCEED IF repaired\n")
        with self.assertRaisesRegex(ValueError, "conditional verdict"):
            self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_conditional_reading_links_exact_disposition(self):
        self.reading.write_text(f"Packet identity: {self.identity}\nVerdict: PROCEED IF suite matches\n")
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
        self.reading.write_text(
            f"Packet identity: {self.identity}\nVerdict: PROCEED IF a fresh reading is done\n"
        )
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

        self.reading.write_text(
            f"Packet identity: {self.identity}\nVerdict: PROCEED IF Developer chooses PUBLISH\n"
        )
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
        (self.evidence / "MANIFEST.sha256").unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            checker.freeze(self.evidence)
        self.identity = checker.verify(self.evidence)
        self.reading.write_text(f"Packet identity: {self.identity}\nVerdict: PROCEED\n")
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


if __name__ == "__main__":
    unittest.main()
