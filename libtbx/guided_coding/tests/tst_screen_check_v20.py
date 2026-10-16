"""Controls for the claims that screen_check.py actually makes."""

import compileall
import contextlib
import importlib.util
import io
import json
import os
import py_compile
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "payload" / "tools"))
import screen_check as checker


TREE = "b" * 40
BASE = "a" * 40
COMMIT = "c" * 40
REMOTE_URL = "ssh://git@example.invalid/phenix.git"
OUTGOING_KEYS = ("repository", "remote", "remote_url", "base", "commit", "tree", "refspec")


def stand_in(lines, code=0):
    """Revision 15: an executable stand-in run by this test's own interpreter (no /bin/sh):
    prints `lines` and exits `code`."""
    return ("#!%s\nimport sys\n" % sys.executable
            + "".join("print(%r)\n" % line for line in lines) + "sys.exit(%d)\n" % code)


def outgoing(omit=None, **fields):
    """One OUTGOING.txt block (SPEC_A 1.2); `omit` drops a key, fields override."""
    values = dict(repository="phenix", remote="origin", remote_url=REMOTE_URL,
                  base=BASE, commit=COMMIT, tree=TREE)
    values.update(fields)
    values.setdefault("refspec", f"{values['commit']}:refs/heads/master")
    return "".join(f"{key}: {values[key]}\n" for key in OUTGOING_KEYS if key != omit)


APP_ENTRYPOINT = "claude-desktop"


def session_environment(path, **variables):
    """A subprocess environment built from scratch, not inherited: `path` is
    the whole PATH, the interpreter keeps the few variables it needs, and the
    Claude session kind (CLAUDE_CODE_ENTRYPOINT, CLAUDE_CODE_EXECPATH) comes
    only from `variables` (None omits one), never from the session this test
    process itself runs in, which may be a Claude app session."""
    environment = {"PATH": str(path), "HOME": os.environ["HOME"],
                   "LANG": os.environ.get("LANG", "C.UTF-8")}
    for name in ("LC_ALL", "LC_CTYPE", "TMPDIR", "TEMP", "TMP", "SYSTEMROOT"):
        if name in os.environ:
            environment[name] = os.environ[name]
    environment.update((name, value) for name, value in variables.items()
                       if value is not None)
    return environment


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
        self.claude.write_text(stand_in(["2.1.284 (Claude Code)"]))
        self.claude.chmod(0o755)

    def register(self, config=None, entrypoint=None):
        """Register from a Terminal session (the default), or from a Claude
        app session when `entrypoint` is APP_ENTRYPOINT; CLAUDE_CODE_EXECPATH
        is never set, and this test process's own session never leaks in."""
        return subprocess.run(
            [sys.executable, "-I", "-B", str(self.source / "payload/tools/screen_check.py"),
             "register-skill", str(self.source)],
            env=session_environment(self.bin, CLAUDE_CONFIG_DIR=str(config or self.config),
                                    CLAUDE_CODE_ENTRYPOINT=entrypoint),
            text=True, capture_output=True)

    def test_app_session_without_engine_path_registers_as_not_checked(self):
        self.claude.unlink()
        result = self.register(entrypoint=APP_ENTRYPOINT)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("NOT CHECKED", result.stdout)
        self.assertIn("REGISTERED", result.stdout)
        link = self.config / "skills" / "guided_coding"
        self.assertTrue(link.is_symlink())
        self.assertEqual(link.resolve(), self.source.resolve())

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
        self.claude.write_text(stand_in(["2.1.268 (Claude Code)"]))
        self.assertEqual(self.register().returncode, 2)
        self.assertFalse(self.config.exists())
        self.claude.write_text(stand_in(["2.1.284 (Claude Code)"]))
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
    """check-claude-version (minimum 2.1.281). The session kind is decided
    only by CLAUDE_CODE_ENTRYPOINT: exactly "claude-desktop" is a Claude app
    session, whose engine is the program named by CLAUDE_CODE_EXECPATH; any
    other or missing value is a Terminal session, which reads the `claude`
    command on PATH and ignores CLAUDE_CODE_EXECPATH."""

    NEW = "2.1.284 (Claude Code)"
    OLD = "2.1.268 (Claude Code)"

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.bin = root / "bin"                # the only directory on PATH
        self.engine_dir = root / "app-engine"  # never on PATH
        self.bin.mkdir()
        self.engine_dir.mkdir()
        self.command = self.bin / "claude"
        self.engine = self.engine_dir / "claude"

    @staticmethod
    def write_fake(path, version, exit_code):
        path.write_text(stand_in([version], exit_code))
        path.chmod(0o755)

    def run_gate(self, path_version=None, path_exit=0, entrypoint=None, engine=None):
        """Run check-claude-version in a built-from-scratch environment.

        path_version: what the fake `claude` on PATH prints (None: no such
        command); path_exit: its exit code. entrypoint: CLAUDE_CODE_ENTRYPOINT
        (None: unset). engine: CLAUDE_CODE_EXECPATH, None for unset, "missing"
        for a path with no file, or (version, exit code) for a fake engine in
        a directory that is not on PATH."""
        if path_version is None:
            self.command.unlink(missing_ok=True)
        else:
            self.write_fake(self.command, path_version, path_exit)
        self.engine.unlink(missing_ok=True)
        if isinstance(engine, tuple):
            self.write_fake(self.engine, *engine)
        return subprocess.run(
            [sys.executable, "-I", "-B", str(Path(checker.__file__)), "check-claude-version"],
            env=session_environment(
                self.bin, CLAUDE_CODE_ENTRYPOINT=entrypoint,
                CLAUDE_CODE_EXECPATH=None if engine is None else str(self.engine)),
            text=True, capture_output=True)

    def assert_refused(self, result, *mentions):
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertTrue(result.stderr.startswith("ERROR: "), result.stderr)
        self.assertNotIn("VERIFIED", result.stdout)
        for text in mentions:
            self.assertIn(text, result.stderr)

    def assert_not_checked(self, result, *mentions):
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("NOT CHECKED", result.stdout)
        self.assertNotIn("VERIFIED", result.stdout)
        for text in mentions:
            self.assertIn(text, result.stdout)

    # ---- Terminal session ---------------------------------------------------

    def test_terminal_minimum_and_newer_versions_verified(self):
        for version in ("2.1.281 (Claude Code)", "2.1.284 (Claude Code)", "2.2.0 (Claude Code)"):
            with self.subTest(version=version):
                result = self.run_gate(path_version=version)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("VERIFIED Claude Code CLI", result.stdout)
                self.assertIn(version, result.stdout)
                self.assertIn(str(self.command), result.stdout)
                self.assertIn("the claude command on this shell's PATH", result.stdout)

    def test_terminal_versions_below_minimum_refused(self):
        for version in ("2.1.268 (Claude Code)", "2.1.280 (Claude Code)"):
            with self.subTest(version=version):
                self.assert_refused(self.run_gate(path_version=version),
                                    "below 2.1.281", "claude update")

    def test_terminal_missing_malformed_or_failing_command_refused(self):
        with self.subTest(command="not on PATH"):
            self.assert_refused(self.run_gate(), "not found on PATH")
        with self.subTest(command="unrecognized answer"):
            self.assert_refused(self.run_gate(path_version="Claude Code unknown"), "unrecognized")
        with self.subTest(command="exit 1 with a plausible version"):
            self.assert_refused(self.run_gate(path_version=self.NEW, path_exit=1), "exit 1")

    def test_terminal_session_ignores_the_app_engine_variable(self):
        for entrypoint in (None, "cli", "vscode", "claude-desktop-3p"):
            with self.subTest(entrypoint=entrypoint):
                result = self.run_gate(entrypoint=entrypoint, engine=(self.NEW, 0))
                self.assert_refused(result, "not found on PATH")

    def test_this_sessions_claude_variables_do_not_reach_the_gate(self):
        """Control for the environment every case here builds: even when this
        test process runs inside a Claude app session, a Terminal case with
        no `claude` on PATH is still refused rather than NOT CHECKED."""
        with unittest.mock.patch.dict(os.environ, {
                "CLAUDE_CODE_ENTRYPOINT": APP_ENTRYPOINT,
                "CLAUDE_CODE_EXECPATH": str(self.engine)}):
            environment = session_environment(self.bin)
            self.assertNotIn("CLAUDE_CODE_ENTRYPOINT", environment)
            self.assertNotIn("CLAUDE_CODE_EXECPATH", environment)
            self.assert_refused(self.run_gate(), "not found on PATH")

    # ---- Claude app session -------------------------------------------------

    def test_app_engine_at_or_above_minimum_verified(self):
        for version in ("2.1.281 (Claude Code)", "2.1.284 (Claude Code)", "2.2.0 (Claude Code)"):
            with self.subTest(version=version):
                result = self.run_gate(entrypoint=APP_ENTRYPOINT, engine=(version, 0))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("VERIFIED Claude app's Claude Code engine", result.stdout)
                self.assertIn(version, result.stdout)
                self.assertIn(str(self.engine), result.stdout)
                self.assertNotIn("CLI", result.stdout.replace(str(self.engine), ""))

    def test_app_engine_below_minimum_refused(self):
        for version in ("2.1.268 (Claude Code)", "2.1.280 (Claude Code)"):
            with self.subTest(version=version):
                result = self.run_gate(entrypoint=APP_ENTRYPOINT, engine=(version, 0))
                self.assert_refused(result, "below 2.1.281", "update the Claude app")
                self.assertNotIn("claude update", result.stderr)

    def test_app_without_engine_path_is_not_checked_whatever_is_on_path(self):
        with self.subTest(path_command="absent"):
            result = self.run_gate(entrypoint=APP_ENTRYPOINT)
            self.assert_not_checked(result, "CLAUDE_CODE_EXECPATH")
            self.assertTrue(result.stdout.startswith("NOT CHECKED"), result.stdout)
            self.assertEqual(len(result.stdout.splitlines()), 1, result.stdout)
            self.assertEqual(result.stderr, "")
        with self.subTest(path_command="new"):
            result = self.run_gate(path_version="2.2.0 (Claude Code)", entrypoint=APP_ENTRYPOINT)
            self.assert_not_checked(result, "CLAUDE_CODE_EXECPATH")

    def test_app_engine_that_gives_no_version_is_not_checked(self):
        for label, engine, mention in (
                ("missing file", "missing", "could not be run"),
                ("exit 1 with a plausible version", ("2.1.290 (Claude Code)", 1), "exit 1"),
                ("unrecognized answer", ("Claude Code unknown", 0), "unrecognized")):
            with self.subTest(engine=label):
                result = self.run_gate(path_version="2.2.0 (Claude Code)",
                                       entrypoint=APP_ENTRYPOINT, engine=engine)
                self.assert_not_checked(result, mention)

    def test_app_engine_decides_over_the_path_command(self):
        with self.subTest(engine="old", path_command="new"):
            result = self.run_gate(path_version="2.2.0 (Claude Code)", entrypoint=APP_ENTRYPOINT,
                                   engine=("2.1.270 (Claude Code)", 0))
            self.assert_refused(result, "below 2.1.281", "update the Claude app")
        with self.subTest(engine="new", path_command="old"):
            result = self.run_gate(path_version=self.OLD, entrypoint=APP_ENTRYPOINT,
                                   engine=("2.1.288 (Claude Code)", 0))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("VERIFIED Claude app's Claude Code engine", result.stdout)


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
        (self.evidence / "SERVER_SUITE.txt").write_text("SERVER_SUITE: PASS 12/12\n"
                                                        "synthetic suite log\n")
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

    def result(self, path="light", evidence_id=None, reading_id=None,
               limits="It did not test runtime behavior."):
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
                        f"LIMITS\n{limits}\n"
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

    def publication(self, batch=None, commit=COMMIT, limits="No material gap claimed."):
        reading_id = checker.digest(self.reading.read_bytes())
        if batch is None:
            batch = f"One local change to shared master.\nphenix commit {commit} to master."
        return self.put(f"PUBLICATION | batch-1 | {self.identity}\n"
                        f"BATCH\n{batch}\n"
                        "SERVER CHECK\nSuite and roster are in evidence.\n"
                        f"LIMITS\n{limits}\n"
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
        """SPEC_A 1.1, counted as in SPEC_A 10: zero Scope lines, or several
        with different values, are refused on the full RESULT and on the
        PUBLICATION screen with "exactly one Scope line"; a single line with an
        unrecognized value is refused with "Scope value not recognized".
        Changed in checkpoint 1 (rev19): two identical Scope lines are accepted
        as one (this test refused them before)."""
        for label, scope_lines in (
                ("none", ""),
                ("two different", "Scope: integration\nScope: publication\n")):
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
        # Checkpoint 1 (rev19): two identical Scope lines now count as one.
        self.write_reading(scope=None, extra="Scope: integration and publication\n" * 2)
        self.ok("result", self.result("full"), self.reading, self.evidence)
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
        check. Checkpoint 3: the first line must be `SERVER_SUITE: <outcome>`
        (see test_suite_first_line_form_is_checked_in_draft_and_final)."""
        for text in ("SERVER_SUITE: PASS 12/12\n", "SERVER_SUITE: RUN, 3 failed\nSummary\n"):
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

    def test_suite_first_line_form_is_checked_in_draft_and_final(self):
        """Checkpoint 3 (review finding 4): SERVER_SUITE.txt's first line must be
        `SERVER_SUITE: <outcome>` (after a byte-order mark and spaces) at the
        final PUBLICATION check, as in check-draft; the forms that passed
        before because only a first line starting NOT RUN was read are refused."""
        reason = "first line must be 'SERVER_SUITE: <outcome>'"
        for label, text in (
                ("free text", "synthetic suite log\n"),
                ("heading first", "# server suite\nSERVER_SUITE: NOT RUN\n"),
                ("summary first", "Summary\nSERVER_SUITE: NOT RUN\n"),
                ("prose not run", "not run, he said it was fine\n"),
                ("no outcome", "SERVER_SUITE:\n"),
                ("lowercase key", "server_suite: PASS\n")):
            with self.subTest(form=label):
                self.suite(text)
                self.publication_refused(reason)
                packet = self.root / "draft"
                shutil.rmtree(packet, ignore_errors=True)
                shutil.copytree(self.evidence, packet)
                (packet / "MANIFEST.sha256").unlink()
                with self.assertRaisesRegex(ValueError, reason), \
                        contextlib.redirect_stdout(io.StringIO()):
                    checker.check_draft("publication", packet)
        self.suite("\ufeff  SERVER_SUITE: PASS 12/12\n")
        self.ok("publication", self.publication(), self.reading, self.evidence)

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
        """SPEC_A 10 (1.2): a key repeated within one repository block
        (remote_url: first a different destination, then the correct one) is
        refused at PUBLICATION with "duplicate key in outgoing record: <key>";
        the plain single block passes. Checkpoint 3 (review finding 3): a key
        outside the seven, such as the formerly tolerated `companions`, is
        refused once or twice."""
        record = self.evidence / "OUTGOING.txt"
        with self.subTest(key="remote_url"):
            record.write_text(outgoing().replace(
                "remote_url:", "remote_url: ssh://git@elsewhere.invalid/other.git\nremote_url:", 1))
            self.refreeze()
            self.publication_refused("duplicate key in outgoing record: remote_url")
        with self.subTest(key="companions"):
            record.write_text(outgoing() + "companions: cctbx\ncompanions: dxtbx\n")
            self.refreeze()
            self.publication_refused("not one of the seven outgoing keys")
        with self.subTest(key="single extra key"):
            record.write_text(outgoing() + "companions: cctbx\n")
            self.refreeze()
            self.publication_refused("not one of the seven outgoing keys")
        with self.subTest(key="plain single block"):
            record.write_text(outgoing())
            self.refreeze()
            self.ok("publication", self.publication(), self.reading, self.evidence)

    # ---- Checkpoint 1 (rev19): repeated reading fields ----------------------

    def both_gates_accept(self):
        self.ok("result", self.result("full"), self.reading, self.evidence)
        self.ok("publication", self.publication(), self.reading, self.evidence)

    def both_gates_refuse(self, reason):
        with self.assertRaisesRegex(ValueError, reason):
            self.ok("result", self.result("full"), self.reading, self.evidence)
        with self.assertRaisesRegex(ValueError, reason):
            self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_identical_repeated_reading_fields_count_once(self):
        """Identical repeated `Packet identity:`, `Verdict:` and `Scope:` lines
        are accepted as one at RESULT and PUBLICATION, including a whole
        closing block repeated in a summary."""
        block = (f"Packet identity:  {self.identity}\nVerdict:          PROCEED\n"
                 "Scope:            integration and publication\n")
        for label, extra in (
                ("identity", f"Packet identity: {self.identity}\n"),
                ("verdict", "Verdict: PROCEED\n"),
                ("scope", "Scope: integration and publication\n"),
                ("whole closing block in a summary",
                 "\nSummary of the reading\n" + block + "Findings: 0 blocking\n")):
            with self.subTest(repeated=label):
                self.write_reading(extra=extra)
                self.both_gates_accept()

    def test_conflicting_repeated_reading_fields_are_refused(self):
        """Repeats with a different value are still refused: identity,
        verdict and scope."""
        for label, extra, reason in (
                ("identity", f"Packet identity: {'0' * 64}\n", "does not name this packet"),
                ("identity not a hash", "Packet identity: see the message\n",
                 "does not name this packet"),
                ("verdict", "Verdict: DO NOT PROCEED\n", "does not name this packet"),
                ("verdict condition", "Verdict: PROCEED IF suite matches\n",
                 "does not name this packet"),
                ("scope", "Scope: integration\n", "exactly one Scope line")):
            with self.subTest(conflict=label):
                self.write_reading(extra=extra)
                self.both_gates_refuse(reason)

    def test_wrong_identity_refused_even_when_repeated_identically(self):
        """A reading naming another packet is refused, once or repeated."""
        other = "0" * 64
        for extra in ("", f"Packet identity: {other}\n"):
            with self.subTest(repeated=bool(extra)):
                self.write_reading(identity=other, extra=extra)
                self.both_gates_refuse("does not name this packet")

    def test_repeated_conditional_verdict_still_needs_its_disposition(self):
        """An identical repeated PROCEED IF is one condition: it still needs
        its disposition, and is accepted with it."""
        self.write_reading("PROCEED IF suite matches",
                           extra="Verdict: PROCEED IF suite matches\n")
        with self.assertRaisesRegex(ValueError, "conditional verdict"):
            self.ok("result", self.result("full"), self.reading, self.evidence)
        disposition = self.root / "disposition.txt"
        disposition.write_text("Condition: suite matches\nStatus: SATISFIED\nEvidence: roster review\n")
        gate = (f"Verdict: PROCEED IF suite matches\nDisposition SHA256: "
                f"{checker.digest(disposition.read_bytes())}")
        screen = self.put(self.result("full").read_text().replace("Verdict: PROCEED", gate))
        self.ok("result", screen, self.reading, self.evidence, disposition)

    def test_reading_with_repeats_is_never_modified_and_gate_names_its_exact_bytes(self):
        """The checker only reads the reading: its bytes are unchanged after
        both gates, and a GATE naming a de-duplicated copy is refused."""
        self.write_reading(extra=f"Packet identity: {self.identity}\nVerdict: PROCEED\n")
        before = self.reading.read_bytes()
        self.both_gates_accept()
        self.assertEqual(self.reading.read_bytes(), before)
        cleaned = before.replace(f"Packet identity: {self.identity}\nVerdict: PROCEED\n".encode(), b"")
        self.assertNotEqual(cleaned, before)
        with self.assertRaisesRegex(ValueError, "does not identify"):
            self.ok("result", self.result("full", reading_id=checker.digest(cleaned)),
                    self.reading, self.evidence)
        self.assertEqual(self.reading.read_bytes(), before)

    # ---- Checkpoint 3: final check agrees with the draft; decorated reading lines

    def test_final_publication_check_refuses_prose_and_extra_keys_like_the_draft(self):
        """Checkpoint 3 (item 3b, review finding 3): the final PUBLICATION check
        applies check-draft's OUTGOING.txt rule: prose (with or without ': '),
        an indented key and a key outside the seven are refused in both."""
        for label, extra, reason in (
                ("sentence", "This batch fixes the map sharpening bug.\n", "prose or non-key line"),
                ("capitalized prose key", "Public description: fixes the bug.\n",
                 "prose or non-key line"),
                ("indented key", " remote: origin\n", "prose or non-key line"),
                ("description key", "description: fixes the map sharpening bug for Miao\n",
                 "not one of the seven outgoing keys")):
            with self.subTest(case=label):
                (self.evidence / "OUTGOING.txt").write_text(outgoing() + extra)
                self.refreeze()
                self.publication_refused(reason)
                packet = self.root / "draft"
                shutil.rmtree(packet, ignore_errors=True)
                shutil.copytree(self.evidence, packet)
                (packet / "MANIFEST.sha256").unlink()
                with self.assertRaisesRegex(ValueError, reason), \
                        contextlib.redirect_stdout(io.StringIO()):
                    checker.check_draft("publication", packet)
        (self.evidence / "OUTGOING.txt").write_text("# comment\n\n" + outgoing())
        self.refreeze()
        self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_final_publication_check_refuses_companions_inside_the_packet(self):
        """Checkpoint 3 (review finding 14): a companion file inside the frozen
        packet is refused by the final check, as by check-draft."""
        for name in ("ERRATA.md", "notes/APPROVALS.md", "CHECKER.md"):
            with self.subTest(name=name):
                path = self.evidence / name
                path.parent.mkdir(exist_ok=True)
                path.write_text("companion\n")
                self.refreeze()
                self.publication_refused("companion file inside the packet")
                path.unlink()
        self.refreeze()
        self.ok("publication", self.publication(), self.reading, self.evidence)

    def test_decorated_conflicting_reading_fields_are_refused(self):
        """Checkpoint 3 (review finding 14): an indented, bulleted, quoted or
        bold Verdict, Scope or Packet identity line whose value differs from
        the column-0 value is refused at both gates; the same decorations
        with the same value are accepted. Checkpoint 5: Markdown headings and
        table rows too."""
        other = "0" * 64
        for label, extra in (
                ("indented verdict", "    Verdict: DO NOT PROCEED\n"),
                ("bulleted verdict", "- Verdict: DO NOT PROCEED\n"),
                ("bold bulleted verdict", "* **Verdict:** DO NOT PROCEED\n"),
                ("quoted verdict", "> Verdict: DO NOT PROCEED\n"),
                ("numbered verdict", "1. Verdict: PROCEED IF the suite passes\n"),
                ("bold scope", "**Scope:** integration\n"),
                ("code identity", f"`Packet identity: {other}`\n"),
                ("heading verdict", "## Verdict: DO NOT PROCEED\n"),
                ("table verdict", "| Verdict | DO NOT PROCEED |\n"),
                ("table bold scope", "| **Scope** | integration |\n")):
            with self.subTest(conflict=label):
                self.write_reading(extra=extra)
                self.both_gates_refuse("indented, bulleted, quoted or bold")
        for label, extra in (
                ("bold verdict", "**Verdict:** PROCEED\n"),
                ("bulleted block", f"- Packet identity: {self.identity}\n- Verdict: PROCEED\n"
                                   "- Scope: integration and publication\n"),
                ("unrelated prose", "  The verdict above stands; scope is as stated.\n"),
                ("heading same verdict", "### Verdict: PROCEED\n"),
                ("table same verdict", "| Verdict | PROCEED |\n|---|---|\n")):
            with self.subTest(same=label):
                self.write_reading(extra=extra)
                self.both_gates_accept()

    # ---- Checkpoint 3 (1b): misplaced verified evidence, correction companion

    def misplace(self, name="SERVER_SUITE.txt", where="previous_packet"):
        """Freeze the packet with `name` only under `where/` (the batch 4 case)."""
        (self.evidence / where).mkdir(exist_ok=True)
        (self.evidence / name).rename(self.evidence / where / name)
        self.refreeze()
        corrections = self.root / "corrections"
        corrections.mkdir(exist_ok=True)
        return corrections

    def correct(self, corrections, name, source, data):
        (corrections / "ERRATA.md").write_text(
            f"# Errata\nClerical: verified file frozen in the wrong place.\n"
            f"Relocated: {name} from {source}\n")
        (corrections / name).write_bytes(data)

    ERRATA_LIMITS = "SERVER_SUITE.txt was relocated after the reading; see ERRATA.md."

    def publication_with(self, corrections, limits=ERRATA_LIMITS):
        with contextlib.redirect_stdout(io.StringIO()) as output:
            checker.check("publication", self.publication(limits=limits), self.evidence,
                          self.reading, corrections=corrections)
        return output.getvalue()

    def test_misplaced_verified_file_accepted_from_identical_correction_copy(self):
        """1b: SERVER_SUITE.txt frozen only under previous_packet/ is refused
        without a companion; with ERRATA.md naming it and a byte-identical copy
        beside it, PUBLICATION passes (and says it used the copy) and the
        packet identity is unchanged. RESULT accepts a misplaced RUN_LOG.txt
        the same way."""
        data = (self.evidence / "SERVER_SUITE.txt").read_bytes()
        corrections = self.misplace()
        identity = self.identity
        self.publication_refused("missing or empty publication evidence: SERVER_SUITE.txt")
        self.correct(corrections, "SERVER_SUITE.txt", "previous_packet/SERVER_SUITE.txt", data)
        output = self.publication_with(corrections)
        self.assertIn("CHECK OK publication", output)
        self.assertIn("NOTE SERVER_SUITE.txt read from the correction companion", output)
        self.assertEqual(checker.verify(self.evidence), identity)
        log = (self.evidence / "RUN_LOG.txt").read_bytes()
        corrections = self.misplace("RUN_LOG.txt", "logs")
        with self.assertRaisesRegex(ValueError, "missing or empty change evidence: RUN_LOG.txt"):
            self.ok("result", self.result("full"), self.reading, self.evidence)
        self.correct(corrections, "RUN_LOG.txt", "logs/RUN_LOG.txt", log)
        with contextlib.redirect_stdout(io.StringIO()):
            checker.check("result", self.result("full", limits="RUN_LOG.txt: see ERRATA.md."),
                          self.evidence, self.reading, corrections=corrections)

    def test_correction_source_must_have_the_same_name_and_serve_one_name(self):
        """Checkpoint 5 (review CP3 finding 4): a frozen notes/draft.txt whose text
        is a passing SERVER_SUITE line cannot stand in for SERVER_SUITE.txt, even
        with a byte-identical copy; one frozen file cannot serve two names (the
        reviewer's previous_packet/SERVER_SUITE.txt for ROSTER_COMPARISON.txt too)."""
        corrections = self.misplace()
        (self.evidence / "notes").mkdir()
        fake = b"SERVER_SUITE: PASS 999/999\n"
        (self.evidence / "notes" / "draft.txt").write_bytes(fake)
        self.refreeze()
        self.correct(corrections, "SERVER_SUITE.txt", "notes/draft.txt", fake)
        with self.assertRaisesRegex(ValueError, "from a file of another name: notes/draft.txt"):
            self.publication_with(corrections)
        data = (self.evidence / "previous_packet" / "SERVER_SUITE.txt").read_bytes()
        (self.evidence / "ROSTER_COMPARISON.txt").rename(self.evidence / "roster.txt")
        self.refreeze()
        (corrections / "ERRATA.md").write_text(
            "Relocated: SERVER_SUITE.txt from previous_packet/SERVER_SUITE.txt\n"
            "Relocated: ROSTER_COMPARISON.txt from previous_packet/SERVER_SUITE.txt\n")
        (corrections / "SERVER_SUITE.txt").write_bytes(data)
        (corrections / "ROSTER_COMPARISON.txt").write_bytes(data)
        with self.assertRaisesRegex(ValueError, "relocates ROSTER_COMPARISON.txt from a file "
                                                "of another name"):
            self.publication_with(corrections)
        # Control: the same-name relocation alone is accepted once the roster is back.
        (self.evidence / "roster.txt").rename(self.evidence / "ROSTER_COMPARISON.txt")
        self.refreeze()
        self.correct(corrections, "SERVER_SUITE.txt", "previous_packet/SERVER_SUITE.txt", data)
        self.assertIn("CHECK OK publication", self.publication_with(corrections))

    def test_correction_companion_use_is_named_in_limits_and_in_present(self):
        """Checkpoint 5 (review CP3 finding 5): with a correction companion used,
        a screen whose LIMITS does not name ERRATA.md is refused, and `present`
        says which file came from the companion; without a companion, present
        says nothing of it."""
        data = (self.evidence / "SERVER_SUITE.txt").read_bytes()
        corrections = self.misplace()
        self.correct(corrections, "SERVER_SUITE.txt", "previous_packet/SERVER_SUITE.txt", data)
        with self.assertRaisesRegex(ValueError, "LIMITS must name ERRATA.md"):
            self.publication_with(corrections, limits="No material gap claimed.")
        with self.assertRaisesRegex(ValueError, "LIMITS must name ERRATA.md"):
            with contextlib.redirect_stdout(io.StringIO()):
                checker.present("publication", self.publication(), self.evidence, self.reading,
                                corrections=corrections)
        with contextlib.redirect_stdout(io.StringIO()) as output:
            checker.present("publication", self.publication(limits=self.ERRATA_LIMITS),
                            self.evidence, self.reading, corrections=corrections)
        shown = output.getvalue()
        self.assertIn("NOTE: correction companion used: SERVER_SUITE.txt was read from "
                      + str(corrections / "SERVER_SUITE.txt"), shown)
        self.assertIn("ERRATA.md", shown)
        self.assertNotIn("CHECK OK", shown)
        (self.evidence / "previous_packet" / "SERVER_SUITE.txt").rename(
            self.evidence / "SERVER_SUITE.txt")
        self.refreeze()
        with contextlib.redirect_stdout(io.StringIO()) as output:
            checker.present("publication", self.publication(), self.evidence, self.reading)
        self.assertNotIn("correction companion", output.getvalue())

    def test_correction_copy_one_byte_different_is_refused(self):
        data = (self.evidence / "SERVER_SUITE.txt").read_bytes()
        corrections = self.misplace()
        self.correct(corrections, "SERVER_SUITE.txt", "previous_packet/SERVER_SUITE.txt",
                     data[:-1] + b"!")
        with self.assertRaisesRegex(ValueError, "differs from the frozen"):
            self.publication_with(corrections)

    def test_correction_naming_a_file_not_in_the_manifest_is_refused(self):
        data = (self.evidence / "SERVER_SUITE.txt").read_bytes()
        corrections = self.misplace()
        for source in ("previous_packet/OTHER.txt", "SERVER_SUITE.txt"):
            with self.subTest(source=source):
                self.correct(corrections, "SERVER_SUITE.txt", source, data)
                with self.assertRaisesRegex(ValueError, "not listed in the frozen manifest"):
                    self.publication_with(corrections)

    def test_correction_naming_a_path_outside_the_packet_is_refused(self):
        data = (self.evidence / "SERVER_SUITE.txt").read_bytes()
        corrections = self.misplace()
        outside = self.root / "SERVER_SUITE.txt"
        outside.write_bytes(data)
        for source in ("../SERVER_SUITE.txt", str(outside), "previous_packet/../../SERVER_SUITE.txt",
                       "./previous_packet/SERVER_SUITE.txt"):
            with self.subTest(source=source):
                self.correct(corrections, "SERVER_SUITE.txt", source, data)
                with self.assertRaisesRegex(ValueError, "path outside the packet"):
                    self.publication_with(corrections)

    def test_correction_companion_limits(self):
        """1b limits: the companion directory must lie outside the packet; it may
        not stand in for a file the packet has at the top level or for a file the
        checker does not read; plan screens take none."""
        data = (self.evidence / "SERVER_SUITE.txt").read_bytes()
        corrections = self.misplace()
        self.correct(corrections, "SERVER_SUITE.txt", "previous_packet/SERVER_SUITE.txt", data)
        inside = self.evidence / "previous_packet"
        with self.assertRaisesRegex(ValueError, "outside the frozen packet"):
            self.publication_with(inside)
        roster = (self.evidence / "ROSTER_COMPARISON.txt").read_bytes()
        (corrections / "ERRATA.md").write_text(
            "Relocated: ROSTER_COMPARISON.txt from previous_packet/SERVER_SUITE.txt\n")
        (corrections / "ROSTER_COMPARISON.txt").write_bytes(roster)
        with self.assertRaisesRegex(ValueError, "has it at the top level"):
            self.publication_with(corrections)
        (corrections / "ERRATA.md").write_text(
            "Relocated: README_FIRST.md from previous_packet/SERVER_SUITE.txt\n")
        with self.assertRaisesRegex(ValueError, "does not read"):
            self.publication_with(corrections)
        with self.assertRaisesRegex(ValueError, "only result, publication and accept"):
            with contextlib.redirect_stdout(io.StringIO()):
                checker.check("plan", self.plan(), corrections=corrections)

    # ---- Checkpoint 1 (rev19): ACCEPT screen ---------------------------------

    def accept(self, evidence_id=None, candidate_extra="", pending="Pending: combined outside "
               "reading (integration and publication) of this packet.",
               decision="Accept records this exact candidate; master and origin stay unchanged.",
               name="screen.md"):
        path = self.root / name
        path.write_text(f"ACCEPT | comment-1 | {evidence_id or self.identity}\n"
                        "CANDIDATE\nClarified a comment in sample.py.\n"
                        "Why it matters: Readers can understand the comment more easily.\n"
                        "Full diff and exact new tests: evidence/APPROVAL_REPORT.md\n"
                        f"{candidate_extra}"
                        f"CHECKED\nTested tree: {TREE}\nThe log says compile passed.\n"
                        f"PENDING\n{pending}\n"
                        "LIMITS\nIt did not test runtime behavior.\n"
                        f"DECISION\n{decision}\n"
                        "Recommendation: ACCEPT because the check passed and scope is contained.\n"
                        "Next: After the reading, show the RESULT screen for this packet.\n"
                        "ACTION: ACCEPT / REVISE / DISCARD\n")
        return path

    def test_accept_screen_valid_and_presented_without_identities(self):
        """A valid ACCEPT names this exact frozen packet and is presented with
        its identities hidden."""
        screen = self.accept()
        self.ok("accept", screen, evidence=self.evidence)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            checker.present("accept", screen, self.evidence)
        display = output.getvalue()
        self.assertTrue(display.startswith("**PLEASE READ — YOUR DECISION IS NEEDED**\n\n```\n"))
        self.assertIn("\nACCEPT | comment-1\n", display)
        self.assertNotIn(self.identity, display)
        self.assertNotIn(TREE, display)
        self.assertIn("ACTION: ACCEPT / REVISE / DISCARD", display)

    def test_accept_screen_claiming_integration_or_publication_is_refused(self):
        """ACCEPT neither offers nor claims integration or publication.
        Checkpoint 3 (review findings 1, 2): each screen is written to its own
        file and checked in its own subtest, and the claim cases must fail for
        the claim (not for another rule): every verb form of integrate, merge,
        push and publish is refused. The action-line case is separate, with
        its own reason. The nouns integration and publication stay allowed."""
        claims = (
            ("integrated", "Integrated into master.", "candidate"),
            ("merged", "Merged on master.", "candidate"),
            ("published", "Published to origin.", "candidate"),
            ("pushed", "Pushed to origin.", "candidate"),
            ("integrate offered", "Accept and integrate; master and origin stay unchanged.",
             "decision"),
            ("will merge and push", "Next: On ACCEPT I will merge to master and push to origin.",
             "candidate"),
            ("integrating", "Now integrating into master.", "candidate"),
            ("will publish", "Accepting will publish this to origin.", "candidate"),
            ("merges", "It merges the fix.", "candidate"),
            ("push noun-verb", "Then push.", "candidate"),
            ("publishes", "It publishes nothing yet but publishes later.", "candidate"))
        for index, (label, text, where) in enumerate(claims):
            with self.subTest(claim=label):
                if where == "candidate":
                    screen = self.accept(candidate_extra=text + "\n", name=f"claim-{index}.md")
                else:
                    screen = self.accept(decision=text, name=f"claim-{index}.md")
                self.assertIn(text, screen.read_text())
                with self.assertRaisesRegex(ValueError, "must not offer or claim"):
                    self.ok("accept", screen, evidence=self.evidence)
        with self.subTest(case="integration action line"):
            screen = self.accept(name="action.md")
            screen.write_text(screen.read_text().replace(
                "ACTION: ACCEPT / REVISE / DISCARD", "ACTION: INTEGRATE / REVISE / DISCARD"))
            with self.assertRaisesRegex(ValueError, "last line must be the one exact action line"):
                self.ok("accept", screen, evidence=self.evidence)
        with self.subTest(case="nouns allowed"):
            self.ok("accept", self.accept(
                candidate_extra="Integration and publication decisions follow the reading.\n",
                name="nouns.md"), evidence=self.evidence)
        with self.assertRaisesRegex(ValueError, "master and origin stay unchanged"):
            self.ok("accept", self.accept(decision="Accept records this exact candidate.",
                                          name="no-unchanged.md"), evidence=self.evidence)
        with self.assertRaisesRegex(ValueError, "Pending:"):
            self.ok("accept", self.accept(pending="The reading.", name="no-pending.md"),
                    evidence=self.evidence)

    def test_accept_screen_without_its_frozen_evidence_is_refused(self):
        """ACCEPT needs the frozen packet it names: no evidence, an unfrozen
        directory, another identity and a reading are refused."""
        with self.assertRaisesRegex(ValueError, "evidence directory and manifest identity"):
            self.ok("accept", self.accept())
        with self.assertRaisesRegex(ValueError, "current evidence"):
            self.ok("accept", self.accept(evidence_id="0" * 64), evidence=self.evidence)
        with self.assertRaisesRegex(ValueError, "takes no reading"):
            self.ok("accept", self.accept(), self.reading, self.evidence)
        screen = self.accept()
        (self.evidence / "MANIFEST.sha256").unlink()
        with self.assertRaisesRegex(ValueError, "missing manifest"):
            self.ok("accept", screen, evidence=self.evidence)


class DraftPublicationChecks(unittest.TestCase):
    """Checkpoint 1 (rev19): `check-draft publication DIR` on an unfrozen packet."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.packet = Path(self.temp.name) / "packet"
        self.packet.mkdir()
        (self.packet / "OUTGOING.txt").write_text(outgoing())
        (self.packet / "SERVER_SUITE.txt").write_text("SERVER_SUITE: PASS 12/12\n")
        (self.packet / "ROSTER_COMPARISON.txt").write_text(
            "baseline: run 1\ncandidate: run 2\nnew failures: none\n")
        (self.packet / "PUBLICATION_NOTES.md").write_text("Public description: fix.\n")

    def draft(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            checker.check_draft("publication", self.packet)
        return output.getvalue()

    def refused(self, reason):
        with self.assertRaisesRegex(ValueError, reason):
            self.draft()

    def snapshot(self):
        """Every path under the packet, with the bytes of each file (checkpoint 3,
        review finding 13: subdirectories and contents, not only top-level names)."""
        return sorted((p.relative_to(self.packet).as_posix(),
                       p.read_bytes() if p.is_file() else None)
                      for p in self.packet.rglob("*"))

    def test_draft_passes_names_skipped_checks_and_writes_nothing(self):
        (self.packet / "logs").mkdir()
        (self.packet / "logs" / "suite.log").write_text("suite log\n")
        before = self.snapshot()
        output = self.draft()
        self.assertIn("DRAFT OK publication", output)
        self.assertIn("skipped: freeze and packet identity", output)
        self.assertIn("skipped: outgoing commits named in BATCH", output)
        self.assertIn("skipped: outside reading", output)
        self.assertEqual(self.snapshot(), before)
        result = subprocess.run(
            [sys.executable, "-I", "-B", checker.__file__, "check-draft", "publication",
             str(self.packet)], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("DRAFT OK publication", result.stdout)
        self.assertEqual(self.snapshot(), before)

    def test_draft_refuses_missing_outgoing(self):
        (self.packet / "OUTGOING.txt").unlink()
        self.refused("missing or empty publication evidence at the packet top level: OUTGOING.txt")

    def test_draft_refuses_prose_in_outgoing(self):
        for prose in ("This batch fixes the map sharpening bug.\n",
                      "Public description: fixes the map sharpening bug.\n",
                      " remote: origin\n"):
            with self.subTest(prose=prose):
                (self.packet / "OUTGOING.txt").write_text(outgoing() + prose)
                self.refused("prose or non-key line")

    def test_draft_refuses_prose_disguised_as_an_extra_key(self):
        """Checkpoint 3 (review finding 3): a lowercase key outside the seven,
        for example a description, is refused."""
        for extra in ("description: fixes the map sharpening bug for Miao\n",
                      "companions: cctbx\n"):
            with self.subTest(extra=extra):
                (self.packet / "OUTGOING.txt").write_text(outgoing() + extra)
                self.refused("not one of the seven outgoing keys")

    def test_draft_refuses_suite_or_roster_not_at_top_level(self):
        for name in ("SERVER_SUITE.txt", "ROSTER_COMPARISON.txt"):
            with self.subTest(name=name):
                (self.packet / "logs").mkdir(exist_ok=True)
                (self.packet / name).rename(self.packet / "logs" / name)
                self.refused(f"at the packet top level: {name}")
                (self.packet / "logs" / name).rename(self.packet / name)

    def test_draft_refuses_not_run_suite_without_quoted_waiver(self):
        (self.packet / "SERVER_SUITE.txt").write_text(
            "SERVER_SUITE: NOT RUN\nReason: no server reachable.\n")
        self.refused("suite not run without the Developer's quoted waiver")
        (self.packet / "SERVER_SUITE.txt").write_text(
            "SERVER_SUITE: NOT RUN\nWaiver (Developer, 2026-10-12, batch x):\n> quoted words\n")
        self.assertIn("DRAFT OK publication", self.draft())

    def test_draft_refuses_companions_inside_the_packet(self):
        (self.packet / "ERRATA.md").write_text("errata\n")
        self.refused("companion file inside the packet")


class VersionGateScopeChecks(unittest.TestCase):
    """Only check-claude-version (and register-skill, covered above) run the
    version gate: the source and screen commands succeed in a Terminal
    environment that has no `claude` at all."""
    put = ScreenChecks.put
    plan = ScreenChecks.plan

    def setUp(self):
        SourceInventoryChecks.setUp(self)
        self.root = Path(self.temporary.name)
        self.empty_bin = self.root / "empty-bin"
        self.empty_bin.mkdir()

    def run_tool(self, *arguments):
        return subprocess.run(
            [sys.executable, "-I", "-B", str(self.source / "payload/tools/screen_check.py"),
             *arguments],
            cwd=self.source, env=session_environment(self.empty_bin),
            text=True, capture_output=True)

    def test_only_the_version_command_needs_claude_on_path(self):
        gate = self.run_tool("check-claude-version")
        self.assertEqual(gate.returncode, 2, gate.stdout)
        self.assertIn("not found on PATH", gate.stderr)
        source = self.run_tool("verify-source", str(self.source))
        self.assertEqual(source.returncode, 0, source.stderr)
        self.assertIn("VERIFIED complete source", source.stdout)
        shown = self.run_tool("present", "plan", str(self.plan()))
        self.assertEqual(shown.returncode, 0, shown.stderr)
        self.assertIn("PLAN | comment-1 | light", shown.stdout)


AUTO_QUEUE = Path(__file__).resolve().parent.parent / "auto_queue.py"
# AUTO_QUEUE Revision 12: `run` needs macOS; off macOS (real, or simulated inside the
# test processes with GC_TEST_SIMULATE_NON_MACOS=1, read only by the tests) the queue
# is built through the test-only launcher, whose sandbox stand-in runs the local
# command `true` unrestricted.
AUTO_QUEUE_STANDIN = Path(__file__).resolve().parent / "auto_queue_standin.py"
# Revision 12 (P8): the launcher runs only queues in a temporary directory with this marker.
STANDIN_MARKER = "GC_KIT_TEST_QUEUE_ROOT"
STANDIN_MARKER_TEXT = "guided_coding kit tests: test-only queue root for this test run; never a real queue\n"
NON_MACOS = sys.platform != "darwin" or os.environ.get("GC_TEST_SIMULATE_NON_MACOS") == "1"
# Revision 2 (S5): the one STOP command the checker accepts names the
# auto_queue.py in the same procedure root as the screen_check.py being run.
CHECKER_AUTO_QUEUE = (Path(checker.__file__).parent.parent.parent / "auto_queue.py").resolve()
AUTO_START_HEADINGS = ("JOBS", "CHECKS", "INSTALLATION", "RESULTS", "STOP", "NOT AUTHORIZED")
AUTO_REPORT_HEADINGS = ("JOBS", "IN SIMPLE WORDS", "APPROVAL SUMMARIES", "NEXT")
AUTO_START_ACTION = "ACTION: NONE NEEDED (reply stop to stop the queue)"
AUTO_REPORT_ACTION = "ACTION: APPROVE / INSPECT / REVISE / DISCARD"
# AUTO_SCREENS Revision 9: a row is ready when its step begins with one of these.
READY_STEPS = ("approve for integration", "decide P")
SUMMARY_LABELS = ("Bug: ", "Fix: ", "Test: ", "Criterion: ", "Limits: ", "Approval: ")
APPROVAL_LINE = "Approval: accepts this exact packet and merges nothing."


# Revision 16: the auto-queue helper needs psutil, and building a queue needs git; tests that
# need either skip with these reasons when it is missing.
HAS_PSUTIL = importlib.util.find_spec("psutil") is not None
PSUTIL_SKIP = ("psutil not importable by %s: the auto-queue helper stops without it (exit 18); "
               "run with the cctbx conda Python (libtbx.python) or pip install psutil"
               % sys.executable)
GIT_SKIP = "git not on PATH: building a test queue needs a git repository (candidate)"
# Revision 16 (repair 1, finding 3): the start screen prints the interpreter's absolute path.
# Repair 2 (reading 2, finding 2): the checker requires an existing executable interpreter, so
# the tests make stand-in interpreter files (never run) in a temporary build tree.
import atexit
STOP_BUILD = Path(os.path.realpath(tempfile.mkdtemp(prefix="gc-stop-build-")))
atexit.register(shutil.rmtree, str(STOP_BUILD), True)
for _name in ("build/bin/libtbx.python", "conda_base/bin/python3", "conda_base/bin/python3.11",
              "noexec/build/bin/libtbx.python"):
    (STOP_BUILD / _name).parent.mkdir(parents=True, exist_ok=True)
    (STOP_BUILD / _name).write_text("#!/bin/sh\nexit 0\n")
    (STOP_BUILD / _name).chmod(0o644 if _name.startswith("noexec") else 0o755)
(STOP_BUILD / "dir" / "bin" / "libtbx.python").mkdir(parents=True)  # executable, not a file
(STOP_BUILD / "link" / "bin").mkdir(parents=True)
(STOP_BUILD / "link" / "bin" / "libtbx.python").symlink_to(STOP_BUILD / "build/bin/libtbx.python")
STOP_PYTHON = f"{STOP_BUILD}/build/bin/libtbx.python -I -B"
_SANDBOX_STATE = []


def sandbox_problem():
    """Revision 16: None when /usr/bin/sandbox-exec can apply a profile (the helper's network
    block for `run` on macOS); else why not (missing, or refused as inside another
    sandbox-exec). Checked once."""
    if not _SANDBOX_STATE:
        problem = None
        if not os.path.exists("/usr/bin/sandbox-exec"):
            problem = "/usr/bin/sandbox-exec missing: building a test queue runs it (`run`)"
        else:
            probe = subprocess.run(["/usr/bin/sandbox-exec", "-p", "(version 1)(allow default)",
                                    sys.executable, "-c", "pass"], capture_output=True, text=True,
                                   timeout=60)
            if probe.returncode != 0:
                problem = ("sandbox-exec cannot apply a profile here (exit %d: %s), as inside "
                           "another sandbox-exec; building a test queue runs it (`run`)"
                           % (probe.returncode, probe.stderr.strip()[:120] or "no message"))
        _SANDBOX_STATE.append(problem)
    return _SANDBOX_STATE[0]


def word_count(text):
    return len(text.split())


class AutoScreenChecks(unittest.TestCase):
    """AUTO_SCREENS_SPEC: the auto_start and auto_report kinds, checked only
    through the screen_check.py command line. auto_report screens are compared
    with the report of a real queue built with auto_queue.py in a temporary
    directory (temporary --home, GC_AUTO_HOME and --lock-dir)."""

    QUEUE_ID = "night-1"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        (self.root / STANDIN_MARKER).write_text(STANDIN_MARKER_TEXT)  # Revision 12 (P8)
        self.home = self.root / "home"
        self.home.mkdir()
        self.environment = dict(os.environ, GC_AUTO_HOME=str(self.home),
                                GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1")
        self.screen = self.root / "screen.md"

    # ---- helpers ------------------------------------------------------------

    def tool(self, *arguments):
        return subprocess.run(
            [sys.executable, "-I", "-B", str(Path(checker.__file__)), *arguments],
            cwd=self.root, env=self.environment, text=True, capture_output=True)

    def queue_tool(self, *arguments, queue=None, expect=0):
        if not HAS_PSUTIL:
            self.skipTest(PSUTIL_SKIP)
        if arguments and arguments[0] == "run" and not NON_MACOS and sandbox_problem():
            self.skipTest(sandbox_problem())  # Revision 16: `run` applies a sandbox-exec profile
        helper = AUTO_QUEUE_STANDIN if NON_MACOS else AUTO_QUEUE
        command = [sys.executable, "-I", "-B", str(helper), "--home", str(self.home)]
        if queue is not None:
            command += ["--queue", str(queue)]
        result = subprocess.run(command + list(arguments), cwd=self.root,
                                env=self.environment, text=True, capture_output=True)
        if expect is not None:
            self.assertEqual(result.returncode, expect,
                             f"auto_queue {arguments}: {result.stdout}{result.stderr}")
        return result

    def git(self, repo, *arguments):
        if shutil.which("git") is None:
            self.skipTest(GIT_SKIP)
        result = subprocess.run(
            ["git", "-C", str(repo), "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
             "-c", "commit.gpgsign=false", *arguments],
            env=self.environment, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.strip()

    def put(self, text):
        self.screen.write_text(text)
        return self.screen

    def put_bytes(self, data):
        self.screen.write_bytes(data)
        return self.screen

    def check(self, kind, screen, *options):
        return self.tool("check", kind, str(screen), *options)

    def assert_accepted(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def assert_refused(self, result):
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertTrue(result.stderr.startswith("ERROR: "), result.stderr)

    # ---- auto_start ---------------------------------------------------------

    def start_text(self, queue_id=None, extra_results=(), stop_line=None,
                   not_authorized="Nothing is merged to master and nothing is published."):
        stop_line = stop_line or f"{STOP_PYTHON} {CHECKER_AUTO_QUEUE} stop"
        lines = [f"AUTO START | {queue_id or self.QUEUE_ID}",
                 "JOBS",
                 "Job A: clarify one comment, needs the unit check.",
                 "Job B: fix one test, needs the unit check.",
                 "CHECKS",
                 "Each job runs its unit check in an isolated copy.",
                 "INSTALLATION",
                 "The shared test installation is used under one lock.",
                 "RESULTS",
                 "A report table is shown when the queue finishes.",
                 "A usage limit stops the queue; resume is manual with /gc auto resume.",
                 *extra_results,
                 "STOP",
                 "Reply stop, or run this command:",
                 stop_line,
                 "NOT AUTHORIZED",
                 not_authorized,
                 AUTO_START_ACTION]
        return "\n".join(lines) + "\n"

    def start_refused(self, text, *options):
        self.assert_refused(self.check("auto_start", self.put(text), *options))

    def test_auto_start_valid_screen_accepted(self):
        """Positive control for every auto_start refusal below."""
        self.assert_accepted(self.check("auto_start", self.put(self.start_text())))
        # NOT AUTHORIZED: case-insensitive, and `published` counts as publish.
        self.assert_accepted(self.check("auto_start", self.put(self.start_text(
            not_authorized="Nothing goes to MASTER; nothing is Published."))))
        self.assert_accepted(self.check("auto_start", self.put(self.start_text(
            not_authorized="No push to master and no publish step."))))

    def test_auto_start_queue_id_pattern(self):
        """Line 1 `AUTO START | <queue-id>`, id `[a-z0-9][a-z0-9-]{0,39}`."""
        for queue_id in ("q", "0", "a" * 40, "a-" + "b" * 38, "night-queue-2026-10-08"):
            with self.subTest(accepted=queue_id):
                self.assert_accepted(self.check("auto_start",
                                                self.put(self.start_text(queue_id))))
        for queue_id in ("a" * 41, "-night", "Night-1", "night_1", "night.1", "night 1",
                         "nïght"):
            with self.subTest(refused=queue_id):
                self.start_refused(self.start_text(queue_id))
        good = self.start_text()
        for label, first in (("report title", f"AUTO REPORT | {self.QUEUE_ID}"),
                             ("no separator", f"AUTO START {self.QUEUE_ID}"),
                             ("empty id", "AUTO START | "),
                             ("lower case title", f"auto start | {self.QUEUE_ID}"),
                             ("trailing text", f"AUTO START | {self.QUEUE_ID} | extra")):
            with self.subTest(line1=label):
                self.start_refused(first + "\n" + good.split("\n", 1)[1])

    def test_auto_start_headings_once_in_order_first_on_line_2(self):
        good = self.start_text()
        for heading in AUTO_START_HEADINGS:
            with self.subTest(missing=heading):
                lines = good.splitlines(keepends=True)
                self.start_refused("".join(l for l in lines if l != heading + "\n"))
            with self.subTest(repeated=heading):
                self.start_refused(good.replace(heading + "\n", heading + "\n" + "Note.\n"
                                                + heading + "\n", 1))
        with self.subTest(order="CHECKS before JOBS"):
            swapped = good.replace(
                "JOBS\nJob A: clarify one comment, needs the unit check.\n"
                "Job B: fix one test, needs the unit check.\n"
                "CHECKS\nEach job runs its unit check in an isolated copy.\n",
                "CHECKS\nEach job runs its unit check in an isolated copy.\n"
                "JOBS\nJob A: clarify one comment, needs the unit check.\n"
                "Job B: fix one test, needs the unit check.\n")
            self.assertNotEqual(swapped, good)
            self.start_refused(swapped)
        with self.subTest(first_heading="on line 3"):
            self.start_refused(good.replace("JOBS\n", "Prepared tonight.\nJOBS\n", 1))
        with self.subTest(heading="lower case"):
            self.start_refused(good.replace("CHECKS\n", "Checks\n", 1))

    def test_auto_start_blank_line_inside_section_refused(self):
        good = self.start_text()
        for label, altered in (
                ("empty", good.replace("Job B:", "\nJob B:", 1)),
                ("blank with spaces", good.replace("Job B:", "   \nJob B:", 1)),
                ("before a heading", good.replace("STOP\n", "\nSTOP\n", 1))):
            with self.subTest(form=label):
                self.assertNotEqual(altered, good)
                self.start_refused(altered)

    def test_auto_start_action_line(self):
        good = self.start_text()
        for label, altered in (
                ("other action", good.replace(AUTO_START_ACTION, "ACTION: APPROVE / STOP")),
                ("decision action",
                 good.replace(AUTO_START_ACTION, AUTO_REPORT_ACTION)),
                ("trailing text", good.replace(AUTO_START_ACTION, AUTO_START_ACTION + " now")),
                ("line after action", good + "postscript\n"),
                ("second ACTION line",
                 good.replace("A report table", "ACTION: STOP\nA report table", 1)),
                ("action missing", good.replace(AUTO_START_ACTION + "\n", ""))):
            with self.subTest(form=label):
                self.start_refused(altered)

    def test_auto_start_encoding_newlines_and_placeholders(self):
        good = self.start_text()
        for label, data in (
                ("CRLF", good.replace("\n", "\r\n").encode()),
                ("no final newline", good.rstrip("\n").encode()),
                ("not UTF-8", good.replace("Job A", "Job \xff").encode("latin-1")),
                ("angle placeholder", good.replace("one comment", "<comment>").encode()),
                ("brace placeholder", good.replace("one comment", "{{comment}}").encode())):
            with self.subTest(form=label):
                self.assert_refused(self.check("auto_start", self.put_bytes(data)))

    def test_auto_start_budget_28_lines_and_220_words(self):
        base = self.start_text()
        room = 28 - len(base.splitlines())
        self.assertGreater(room, 0)
        fill = [f"Filler line {n}." for n in range(room)]
        at_limit = self.start_text(extra_results=fill)
        self.assertEqual(len(at_limit.splitlines()), 28)
        self.assert_accepted(self.check("auto_start", self.put(at_limit)))
        over = self.start_text(extra_results=fill + ["One more."])
        self.assertEqual(len(over.splitlines()), 29)
        self.start_refused(over)
        # Words, counted over the whole screen split on whitespace.
        missing = 220 - word_count(base)
        self.assertGreater(missing, 0)
        words = self.start_text(extra_results=[" ".join(["word"] * missing)])
        self.assertEqual(word_count(words), 220)
        self.assertLessEqual(len(words.splitlines()), 28)
        self.assert_accepted(self.check("auto_start", self.put(words)))
        words = self.start_text(extra_results=[" ".join(["word"] * (missing + 1))])
        self.assertEqual(word_count(words), 221)
        self.start_refused(words)

    def test_auto_start_stop_section_needs_absolute_stop_command(self):
        """Revision 2 (S5): STOP contains exactly the stop command of <AQ>, the resolved
        auto_queue.py beside this screen_check.py; any other path, absolute or not, is
        refused. Revision 16: the command is `/ABS/libtbx.python -I -B <AQ> stop` (the helper
        needs psutil)."""
        self.assertTrue(CHECKER_AUTO_QUEUE.is_file())
        self.assert_accepted(self.check("auto_start", self.put(self.start_text(
            stop_line=f"{STOP_PYTHON} {CHECKER_AUTO_QUEUE} stop"))))
        unresolved = Path(checker.__file__).parent / ".." / ".." / "auto_queue.py"
        for label, line in (
                ("other absolute path", f"{STOP_PYTHON} /opt/guided_coding/auto_queue.py stop"),
                ("nonexistent absolute path", f"{STOP_PYTHON} /nonexistent/auto_queue.py stop"),
                ("same file, unresolved path", f"{STOP_PYTHON} {unresolved} stop"),
                ("trailing text", f"{STOP_PYTHON} {CHECKER_AUTO_QUEUE} stop now"),
                ("leading text", f"Run: {STOP_PYTHON} {CHECKER_AUTO_QUEUE} stop"),
                ("python instead of libtbx.python", f"python -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("bare libtbx.python (repair 1)", f"libtbx.python -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("relative path", f"{STOP_PYTHON} auto_queue.py stop"),
                ("dot relative path", f"{STOP_PYTHON} ./payload/auto_queue.py stop"),
                ("home tilde path", f"{STOP_PYTHON} ~/GuidedCoding/auto_queue.py stop"),
                ("no stop subcommand", f"{STOP_PYTHON} {CHECKER_AUTO_QUEUE} status"),
                ("other tool",
                 f"{STOP_PYTHON} {CHECKER_AUTO_QUEUE.with_name('screen_check.py')} stop"),
                ("no command", "Ask the Guide to stop the queue.")):
            with self.subTest(form=label):
                self.start_refused(self.start_text(stop_line=line))
        with self.subTest(form="command only in RESULTS, not in STOP"):
            self.start_refused(self.start_text(
                extra_results=[f"{STOP_PYTHON} {CHECKER_AUTO_QUEUE} stop"],
                stop_line="Ask the Guide to stop the queue."))

    def test_auto_start_not_authorized_mentions_master_and_publish(self):
        for label, text in (("no master", "Nothing is published and nothing is merged."),
                            ("no publish", "Nothing is merged to master."),
                            ("neither", "Nothing leaves this machine.")):
            with self.subTest(form=label):
                self.start_refused(self.start_text(not_authorized=text))
        with self.subTest(form="mentioned only outside NOT AUTHORIZED"):
            self.start_refused(self.start_text(
                extra_results=["Nothing is merged to master and nothing is published."],
                not_authorized="Nothing leaves this machine."))

    def test_auto_start_refuses_queue_evidence_reading_and_disposition(self):
        queue = self.root / "some-queue"
        queue.mkdir()
        other = self.root / "other.txt"
        other.write_text("x\n")
        screen = self.put(self.start_text())
        self.assert_accepted(self.check("auto_start", screen))
        for option, value in (("--queue", queue), ("--evidence", self.root),
                              ("--reading", other), ("--disposition", other)):
            with self.subTest(option=option):
                self.assert_refused(self.check("auto_start", screen, option, str(value)))

    def test_auto_start_present_shows_no_action_attention_line(self):
        screen = self.put(self.start_text())
        saved = screen.read_bytes()
        shown = self.tool("present", "auto_start", str(screen))
        self.assert_accepted(shown)
        display = shown.stdout
        self.assertTrue(display.startswith("**PLEASE READ — NO ACTION NEEDED**\n\n```\n"),
                        display)
        self.assertTrue(display.endswith("\n```\n"), display)
        self.assertIn(f"AUTO START | {self.QUEUE_ID}\n", display)
        for heading in AUTO_START_HEADINGS:
            self.assertIn("\n\n" + heading + "\n", display)
        self.assertIn(AUTO_START_ACTION, display)
        self.assertEqual(screen.read_bytes(), saved)
        refused = self.tool("present", "auto_start",
                            str(self.put(self.start_text(stop_line=f"{STOP_PYTHON} auto_queue.py stop"))))
        self.assert_refused(refused)
        self.assertNotIn("PLEASE READ", refused.stdout)

    # ---- auto_report --------------------------------------------------------

    def build_queue(self, queue_id=None, extra_ready=(), provisional=()):
        """A real queue: job a Ready for approval (one passing current run),
        job b Blocked; then each job in `extra_ready` Ready the same way, and a
        pending provisional choice recorded for each job in `provisional`.
        Returns the queue directory."""
        queue_id = queue_id or self.QUEUE_ID
        base = self.root / f"build-{queue_id}"
        repo = base / "repo"
        repo.mkdir(parents=True)
        self.git(repo, "init", "-q")
        (repo / "sample.py").write_text("# sample\n")
        self.git(repo, "add", "sample.py")
        self.git(repo, "commit", "-q", "-m", "sample")
        tree = self.git(repo, "rev-parse", "HEAD^{tree}")
        (base / "GRANT.md").write_text("Grant for the test queue.\n")
        jobs = [{"id": "a", "title": "Clarify comment", "requires": ["unit"], "depends_on": []},
                {"id": "b", "title": "Fix test", "requires": ["unit"], "depends_on": []}]
        jobs += [{"id": job, "title": f"Extra job {job}", "requires": ["unit"], "depends_on": []}
                 for job in extra_ready]
        (base / "JOBS.json").write_text(json.dumps(jobs) + "\n")
        (base / "criterion-1.txt").write_text("Criterion one.\n")
        created = self.queue_tool(
            "init", "--records", str(base / "records"), "--queue-id", queue_id,
            "--grant", str(base / "GRANT.md"), "--jobs", str(base / "JOBS.json"),
            "--lock-dir", str(base / "RUN.lock"), "--min-free-gib", "0",
            "--allowed-root", str(repo))
        queue = Path(created.stdout.strip().splitlines()[-1])
        self.assertTrue((queue / "QUEUE.json").is_file(), created.stdout)
        self.queue_tool("lock", "take", queue=queue)
        for job in ("a", *extra_ready):
            packet = self.make_ready(queue, base, repo, tree, job)
            if job == "a":
                self.packet_manifest = checker.digest((packet / "MANIFEST.sha256").read_bytes())
                self.queue_tool("set", "b", "Blocked", "--note", "needs a design decision",
                                queue=queue)
        for job in provisional:
            self.queue_tool("provisional", job, "--note", "used a smaller test map", queue=queue)
        self.queue_tool("lock", "release", queue=queue)
        self.queue_base = base
        return queue

    def make_ready(self, queue, base, repo, tree, job):
        """Criterion, candidate, one passing current run and a frozen packet; then
        Ready for approval. Returns the packet directory."""
        packet = base / f"packet-{job}" if job != "a" else base / "packet"
        packet.mkdir()
        self.queue_tool("criterion", job, "--file", str(base / "criterion-1.txt"), queue=queue)
        self.queue_tool("candidate", job, "--repo", str(repo), "--commit", "HEAD", queue=queue)
        # Revision 2: run needs the job Preparing/Testing, a registered worker
        # for it, and probe paths that exist; readiness needs the worker ended.
        self.queue_tool("set", job, "Preparing", queue=queue)
        self.queue_tool("worker", "start", job, "--id", f"worker-{job}", queue=queue)
        self.queue_tool("set", job, "Testing", queue=queue)
        self.assertTrue((repo / "sample.py").is_file())
        self.queue_tool("run", job, "unit", "--log", str(base / f"unit-{job}.log"),
                        "--probe", f"echo {repo / 'sample.py'}", "--", sys.executable, "-c",
                        "pass", queue=queue)
        self.queue_tool("worker", "end", job, "--id", f"worker-{job}", "--outcome", "finished",
                        queue=queue)
        # Revision 4: the packet is really frozen with this procedure's
        # screen_check.py; Approved would need --manifest <its SHA-256>.
        # AUTO_QUEUE Revision 11: the recorded candidate id and every candidate tree.
        status = json.loads(self.queue_tool("status", queue=queue).stdout)
        jobs = status["jobs"]
        entry = (jobs[job] if isinstance(jobs, dict)
                 else [item for item in jobs if item.get("id") == job][0])
        (packet / "CODE_IDENTITY.txt").write_text(
            f"candidate_id: {entry['candidate_id']}\n"
            + "".join(f"tested_tree: {item['tree']}\n" for item in entry["candidate"]))
        self.assertIn(f"tested_tree: {tree}\n", (packet / "CODE_IDENTITY.txt").read_text())
        # AUTO_QUEUE Revision 9 (Q1): a ready packet carries a valid screening record.
        (packet / "SCREENING.txt").write_text("test: unit | outside: none | included\n")
        frozen = self.tool("freeze", str(packet))
        self.assertEqual(frozen.returncode, 0, frozen.stdout + frozen.stderr)
        self.assertTrue((packet / "MANIFEST.sha256").is_file())
        self.queue_tool("set", job, "Ready for approval", "--packet", str(packet), queue=queue)
        return packet

    def report_lines(self, queue):
        report = self.queue_tool("report", queue=queue).stdout
        return [line for line in report.splitlines() if line.strip()]

    def ready_jobs(self, jobs_lines):
        """Revision 9: the ids of the rows whose step begins `approve for
        integration` or `decide P`, in report order (id = first word of cell 1)."""
        ready = []
        for line in jobs_lines:
            cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
            if line.startswith("|") and len(cells) == 3 and cells[2].startswith(READY_STEPS):
                ready.append(cells[0].split()[0])
        return ready

    def summary_block(self, job, title="", **changes):
        """A valid seven-line approval summary for `job`; `changes` maps a label
        (bug, fix, test, criterion, limits, approval) to its replacement line."""
        lines = [f"Job {job}:" + (f" {title}" if title else ""),
                 changes.get("bug", "Bug: the comment describes the wrong unit."),
                 changes.get("fix", "Fix: the comment now names the unit used."),
                 changes.get("test", "Test: fails before the fix, works after it."),
                 changes.get("criterion", "Criterion: the comment matches the code."),
                 changes.get("limits", "Limits: no behaviour change."),
                 changes.get("approval", APPROVAL_LINE)]
        return lines

    def summaries_for(self, jobs_lines):
        ready = self.ready_jobs(jobs_lines)
        if not ready:
            return ["None ready."]
        return [line for job in ready for line in self.summary_block(job, f"title of {job}")]

    def report_text(self, jobs_lines, queue_id=None, simple=("Job a passed its local checks.",
                                                             "Job b needs a design decision."),
                    action=AUTO_REPORT_ACTION, summaries=None):
        """Revision 9: APPROVAL SUMMARIES holds one block per ready job of
        `jobs_lines` (or `None ready.`) unless `summaries` gives its lines."""
        if summaries is None:
            summaries = self.summaries_for(jobs_lines)
        lines = [f"AUTO REPORT | {queue_id or self.QUEUE_ID}",
                 "JOBS", *jobs_lines,
                 "IN SIMPLE WORDS", *simple,
                 "APPROVAL SUMMARIES", *summaries,
                 "NEXT",
                 "Approve job a for integration, or inspect, revise or discard it.",
                 action]
        return "\n".join(lines) + "\n"

    def report_check(self, text, queue, *options):
        return self.check("auto_report", self.put(text), "--queue", str(queue), *options)

    def test_auto_report_exact_copy_accepted_and_altered_jobs_refused(self):
        queue = self.build_queue()
        rows = self.report_lines(queue)
        passed = [row for row in rows if "local checks passed" in row]
        blocked = [row for row in rows if "blocked" in row]
        self.assertEqual(len(passed), 1, rows)
        self.assertEqual(len(blocked), 1, rows)
        self.assert_accepted(self.report_check(self.report_text(rows), queue))
        with self.subTest(change="Result cell of the blocked job claims passed"):
            cells = blocked[0].split("|")
            cells[2] = " local checks passed; candidate and ticket saved "
            changed = [("|".join(cells) if row == blocked[0] else row) for row in rows]
            self.assertNotEqual(changed, rows)
            self.assert_refused(self.report_check(self.report_text(changed), queue))
        with self.subTest(change="Result cell reworded"):
            changed = [row.replace("needs a design decision", "needs a decision") for row in rows]
            self.assertNotEqual(changed, rows)
            self.assert_refused(self.report_check(self.report_text(changed), queue))
        with self.subTest(change="missing row"):
            self.assert_refused(self.report_check(
                self.report_text([row for row in rows if row != blocked[0]]), queue))
        with self.subTest(change="missing table header"):
            self.assert_refused(self.report_check(self.report_text(rows[1:]), queue))
        with self.subTest(change="extra row"):
            self.assert_refused(self.report_check(self.report_text(
                rows + ["| c Extra job | local checks passed | approve |"]), queue))
        with self.subTest(change="rows reordered"):
            reordered = rows[:2] + list(reversed(rows[2:]))
            self.assertNotEqual(reordered, rows)
            self.assert_refused(self.report_check(self.report_text(reordered), queue))
        with self.subTest(change="row with trailing space"):
            self.assert_refused(self.report_check(
                self.report_text([row + " " if row == passed[0] else row for row in rows]),
                queue))
        self.assert_accepted(self.report_check(self.report_text(rows), queue))

    def test_auto_report_requires_queue_and_refuses_other_options(self):
        queue = self.build_queue()
        text = self.report_text(self.report_lines(queue))
        self.assert_accepted(self.report_check(text, queue))
        with self.subTest(option="no --queue (while ACTIVE names this queue)"):
            # Required: the checker must not fall back to <home>/ACTIVE.
            self.assertTrue((self.home / "ACTIVE").is_file())
            self.assert_refused(self.check("auto_report", self.put(text)))
        with self.subTest(option="--queue names a missing directory"):
            self.assert_refused(self.report_check(text, self.root / "no-such-queue"))
        other = self.root / "other.txt"
        other.write_text("x\n")
        for option, value in (("--evidence", self.root), ("--reading", other),
                              ("--disposition", other)):
            with self.subTest(option=option):
                self.assert_refused(self.report_check(text, queue, option, str(value)))

    def test_auto_report_stale_run_old_passed_copy_refused_fresh_copy_accepted(self):
        queue = self.build_queue()
        old_rows = self.report_lines(queue)
        old_text = self.report_text(old_rows)
        self.assertTrue(any("local checks passed" in row for row in old_rows), old_rows)
        self.assert_accepted(self.report_check(old_text, queue))
        criterion = self.queue_base / "criterion-2.txt"
        criterion.write_text("Criterion two, changed.\n")
        self.queue_tool("criterion", "a", "--file", str(criterion), queue=queue)
        fresh_rows = self.report_lines(queue)
        self.assertNotEqual(fresh_rows, old_rows)
        self.assertFalse(any("passed" in row for row in fresh_rows), fresh_rows)
        self.assertTrue(any("stale" in row for row in fresh_rows), fresh_rows)
        self.assert_refused(self.report_check(old_text, queue))
        self.assert_accepted(self.report_check(self.report_text(
            fresh_rows, simple=("Job a must be checked again.", "Job b needs a decision.")),
            queue))

    def test_auto_report_jobs_are_the_non_blank_report_lines_after_a_stop(self):
        queue = self.build_queue()
        before = self.report_lines(queue)
        stop = self.queue_tool("stop", queue=queue, expect=None)
        self.assertIn(stop.returncode, (0, 10), stop.stdout + stop.stderr)
        report = self.queue_tool("report", queue=queue).stdout
        rows = [line for line in report.splitlines() if line.strip()]
        self.assertEqual(rows[0], "**Stopped**", report)
        self.assertNotEqual(rows, before)
        self.assert_accepted(self.report_check(self.report_text(rows), queue))
        with self.subTest(change="stop line omitted"):
            self.assert_refused(self.report_check(self.report_text(rows[1:]), queue))
        with self.subTest(change="pre-stop copy"):
            self.assert_refused(self.report_check(self.report_text(before), queue))

    def test_auto_report_stopped_word_needs_a_stopped_report(self):
        """Revision 2 (S6): `stopped` (whole word, any case) in IN SIMPLE WORDS
        or NEXT is refused unless the report starts with **Stopped**."""
        queue = self.build_queue()
        rows = self.report_lines(queue)
        self.assertFalse(rows[0].startswith("**"), rows)
        self.assert_accepted(self.report_check(self.report_text(rows), queue))
        with self.subTest(control="not a whole word"):
            self.assert_accepted(self.report_check(self.report_text(
                rows, simple=("Job a passed its local checks; work went on unstopped.",)), queue))
        for label, text in (
                ("IN SIMPLE WORDS", self.report_text(
                    rows, simple=("Job a passed its local checks.", "The queue stopped."))),
                ("upper case", self.report_text(
                    rows, simple=("Job a passed its local checks.", "The queue STOPPED."))),
                ("capitalized", self.report_text(
                    rows, simple=("Job a passed its local checks.", "Stopped: the queue."))),
                ("NEXT", self.report_text(rows).replace(
                    "Approve job a for integration,",
                    "The queue stopped; approve job a for integration,", 1))):
            with self.subTest(running_queue=label):
                self.assert_refused(self.report_check(text, queue))
        # Positive control: the same words once the report starts with **Stopped**.
        self.queue_tool("stop", queue=queue)
        rows = self.report_lines(queue)
        self.assertEqual(rows[0], "**Stopped**", rows)
        for label, text in (
                ("IN SIMPLE WORDS", self.report_text(
                    rows, simple=("Job a passed its local checks.", "The queue stopped."))),
                ("NEXT", self.report_text(rows).replace(
                    "Approve job a for integration,",
                    "The queue STOPPED; approve job a for integration,", 1))):
            with self.subTest(stopped_queue=label):
                self.assert_accepted(self.report_check(text, queue))

    def test_auto_report_stopped_word_refused_while_only_stopping(self):
        """Revision 2 (S6): a report starting **Stopping** does not allow
        `stopped`."""
        queue = self.build_queue()
        self.queue_tool("worker", "start", "b", "--id", "worker-b", queue=queue)
        self.queue_tool("stop", queue=queue, expect=10)
        rows = self.report_lines(queue)
        self.assertTrue(rows[0].startswith("**Stopping**"), rows)
        # A registered worker also suspends job a's readiness (AUTO_QUEUE
        # Revision 2, B4), so these summaries avoid `passed`.
        self.assert_accepted(self.report_check(self.report_text(
            rows, simple=("Job a waits for the worker to end.", "The queue is stopping.")),
            queue))
        self.assert_refused(self.report_check(self.report_text(
            rows, simple=("Job a waits for the worker to end.", "The queue stopped.")), queue))

    def test_auto_report_passed_word_needs_a_passed_row(self):
        """Revision 2 (S6): `passed` in IN SIMPLE WORDS or NEXT is refused
        unless some JOBS row contains `passed`."""
        queue = self.build_queue()
        rows = self.report_lines(queue)
        self.assertTrue(any("passed" in row for row in rows), rows)
        # Positive control: a row says passed, so the summary may say it.
        self.assert_accepted(self.report_check(self.report_text(rows), queue))
        self.assert_accepted(self.report_check(self.report_text(rows).replace(
            "Approve job a for integration,", "Job a passed; approve it for integration,", 1),
            queue))
        criterion = self.queue_base / "criterion-2.txt"
        criterion.write_text("Criterion two, changed.\n")
        self.queue_tool("criterion", "a", "--file", str(criterion), queue=queue)
        rows = self.report_lines(queue)
        self.assertFalse(any("passed" in row for row in rows), rows)
        neutral = ("Job a must be checked again.", "Job b needs a design decision.")
        self.assert_accepted(self.report_check(self.report_text(rows, simple=neutral), queue))
        for label, text in (
                ("IN SIMPLE WORDS", self.report_text(
                    rows, simple=("Job a passed its local checks.",
                                  "Job b needs a design decision."))),
                ("NEXT", self.report_text(rows, simple=neutral).replace(
                    "Approve job a for integration,",
                    "Job a passed earlier; approve it for integration,", 1))):
            with self.subTest(section=label):
                self.assert_refused(self.report_check(text, queue))

    def test_auto_report_jobs_include_test_installation_unavailable_line(self):
        """AUTO_QUEUE Revision 4 (installation availability): when the report
        starts with the `**Test installation unavailable:**` line, JOBS is
        still exactly the report's non-blank lines."""
        queue = self.build_queue()
        before = self.report_lines(queue)
        script = self.queue_base / "fingerprint.sh"
        script.write_text(stand_in(["FINGERPRINT x"]))
        script.chmod(0o755)
        self.queue_tool("install", "check", "--script", str(script), "--expect", "y",
                        queue=queue, expect=12)
        rows = self.report_lines(queue)
        self.assertTrue(rows[0].startswith("**Test installation unavailable:**"), rows)
        self.assertNotEqual(rows, before)
        neutral = ("The test installation does not match.", "Job b needs a design decision.")
        self.assert_accepted(self.report_check(self.report_text(rows, simple=neutral), queue))
        with self.subTest(change="installation line omitted"):
            self.assert_refused(self.report_check(self.report_text(rows[1:], simple=neutral),
                                                  queue))
        with self.subTest(change="pre-check copy"):
            self.assert_refused(self.report_check(self.report_text(before, simple=neutral),
                                                  queue))

    def test_auto_report_waiting_row_never_says_resumes_automatically(self):
        """AUTO_QUEUE Revision 4 (resume claims): a Waiting row does not say the
        queue resumes automatically; the copied report is accepted."""
        queue = self.build_queue()
        self.queue_tool("set", "b", "Preparing", queue=queue)
        self.queue_tool("wait", "b", "--until", "2026-10-09T06:00:00Z", queue=queue)
        rows = self.report_lines(queue)
        waiting = [row for row in rows if "| b " in row]
        self.assertEqual(len(waiting), 1, rows)
        self.assertIn("waiting", waiting[0].lower(), rows)
        report = self.queue_tool("report", queue=queue).stdout
        self.assertNotIn("resumes automatically", report.lower(), report)
        self.assert_accepted(self.report_check(self.report_text(
            rows, simple=("Job a passed its local checks.", "Job b waits for a reset.")),
            queue))

    def test_auto_report_form_rules(self):
        queue = self.build_queue()
        rows = self.report_lines(queue)
        good = self.report_text(rows)
        self.assert_accepted(self.report_check(good, queue))
        summaries = "APPROVAL SUMMARIES\n" + "".join(l + "\n" for l in self.summaries_for(rows))
        self.assertIn(summaries, good)
        cases = [
            ("start title", good.replace("AUTO REPORT |", "AUTO START |", 1)),
            ("bad queue id", self.report_text(rows, queue_id="Night_1")),
            ("first heading on line 3", good.replace("JOBS\n", "Overnight.\nJOBS\n", 1)),
            ("other action", self.report_text(rows, action=AUTO_START_ACTION)),
            ("action with trailing text", self.report_text(rows, action=AUTO_REPORT_ACTION + " x")),
            ("line after action", good + "postscript\n"),
            ("second ACTION line", self.report_text(
                rows, simple=("ACTION: APPROVE", "Job b needs a design decision."))),
            ("blank line in section", self.report_text(
                rows, simple=("Job a passed its local checks.", "",
                              "Job b needs a design decision."))),
            ("angle placeholder", self.report_text(rows, simple=("Job <id> passed.",))),
            ("brace placeholder", self.report_text(rows, simple=("Job {{id}} passed.",))),
            ("NEXT before IN SIMPLE WORDS", good.replace(
                "IN SIMPLE WORDS\nJob a passed its local checks.\nJob b needs a design decision.\n"
                + summaries +
                "NEXT\nApprove job a for integration, or inspect, revise or discard it.\n",
                "NEXT\nApprove job a for integration, or inspect, revise or discard it.\n"
                + summaries +
                "IN SIMPLE WORDS\nJob a passed its local checks.\nJob b needs a design decision.\n")),
            ("NEXT before APPROVAL SUMMARIES", good.replace(
                summaries + "NEXT\nApprove job a for integration, or inspect, revise or discard it.\n",
                "NEXT\nApprove job a for integration, or inspect, revise or discard it.\n"
                + summaries)),
            ("APPROVAL SUMMARIES before IN SIMPLE WORDS", good.replace(
                "IN SIMPLE WORDS\nJob a passed its local checks.\nJob b needs a design decision.\n"
                + summaries,
                summaries + "IN SIMPLE WORDS\nJob a passed its local checks.\n"
                "Job b needs a design decision.\n")),
        ]
        for heading in AUTO_REPORT_HEADINGS:
            lines = good.splitlines(keepends=True)
            cases.append((f"missing {heading}",
                          "".join(l for l in lines if l != heading + "\n")))
            cases.append((f"repeated {heading}", good.replace(
                AUTO_REPORT_ACTION, f"{heading}\nSee above.\n{AUTO_REPORT_ACTION}")))
        for label, altered in cases:
            with self.subTest(form=label):
                self.assertNotEqual(altered, good)
                self.assert_refused(self.report_check(altered, queue))
        for label, data in (("CRLF", good.replace("\n", "\r\n").encode()),
                            ("no final newline", good.rstrip("\n").encode()),
                            ("not UTF-8", good.replace("Job a passed", "Job \xff passed")
                             .encode("latin-1"))):
            with self.subTest(form=label):
                self.assert_refused(self.check("auto_report", self.put_bytes(data),
                                               "--queue", str(queue)))

    def test_auto_report_budget_90_lines_and_900_words(self):
        """Revision 9: at most 90 lines and 900 words (was 60 and 600)."""
        queue = self.build_queue()
        rows = self.report_lines(queue)
        base = self.report_text(rows)
        room = 90 - len(base.splitlines())
        self.assertGreater(room, 2)
        simple = ["Job a passed its local checks.", "Job b needs a design decision."]
        fill = simple + [f"Filler line {n}." for n in range(room)]
        at_limit = self.report_text(rows, simple=fill)
        self.assertEqual(len(at_limit.splitlines()), 90)
        self.assert_accepted(self.report_check(at_limit, queue))
        over = self.report_text(rows, simple=fill + ["One more."])
        self.assertEqual(len(over.splitlines()), 91)
        self.assert_refused(self.report_check(over, queue))
        missing = 900 - word_count(base)
        self.assertGreater(missing, 0)
        words = self.report_text(rows, simple=simple + [" ".join(["word"] * missing)])
        self.assertEqual(word_count(words), 900)
        self.assertLessEqual(len(words.splitlines()), 90)
        self.assert_accepted(self.report_check(words, queue))
        words = self.report_text(rows, simple=simple + [" ".join(["word"] * (missing + 1))])
        self.assertEqual(word_count(words), 901)
        self.assert_refused(self.report_check(words, queue))

    def test_auto_report_present_shows_action_needed_attention_line(self):
        queue = self.build_queue()
        rows = self.report_lines(queue)
        screen = self.put(self.report_text(rows))
        shown = self.tool("present", "auto_report", str(screen), "--queue", str(queue))
        self.assert_accepted(shown)
        display = shown.stdout
        self.assertTrue(display.startswith("**PLEASE READ — ACTION NEEDED**\n\n```\n"), display)
        self.assertTrue(display.endswith("\n```\n"), display)
        self.assertIn(f"AUTO REPORT | {self.QUEUE_ID}\n", display)
        for heading in AUTO_REPORT_HEADINGS:
            self.assertIn("\n\n" + heading + "\n", display)
        for row in rows:
            self.assertIn(row, display)
        self.assertIn(AUTO_REPORT_ACTION, display)
        missing_queue = self.tool("present", "auto_report", str(screen))
        self.assert_refused(missing_queue)
        self.assertNotIn("PLEASE READ", missing_queue.stdout)
        wrong = self.tool("present", "auto_report",
                          str(self.put(self.report_text(rows[:-1]))), "--queue", str(queue))
        self.assert_refused(wrong)
        self.assertNotIn("PLEASE READ", wrong.stdout)

    # ---- Revision 9: approval summaries, usage limit -------------------------

    def test_auto_start_must_mention_usage_limit_and_manual_resume(self):
        """Revision 9: auto_start mentions, anywhere in its body, `usage limit`
        (any case) and `/gc auto resume`; otherwise refused with the stated reason."""
        reason = "auto_start must say a usage limit stops the queue and resume is manual"
        good = self.start_text()
        line = "A usage limit stops the queue; resume is manual with /gc auto resume.\n"
        self.assertIn(line, good)
        for label, text in (
                ("neither", good.replace(line, "")),
                ("no usage limit", good.replace(line, "Resume is manual with /gc auto resume.\n")),
                ("no /gc auto resume", good.replace(line, "A usage limit stops the queue.\n")),
                ("resume command misspelt", good.replace("/gc auto resume", "/gc resume")),
                ("usage and limit apart", good.replace("usage limit", "usage-limit"))):
            with self.subTest(refused=label):
                self.assertNotEqual(text, good)
                result = self.check("auto_start", self.put(text))
                self.assert_refused(result)
                self.assertIn(reason, result.stderr)
        for label, text in (
                ("upper case", good.replace("usage limit", "USAGE LIMIT")),
                ("other sections", good.replace(line, "").replace(
                    "Job A: clarify one comment,", "Job A (stops at a Usage Limit): clarify one comment,")
                 .replace("Nothing is merged to master and nothing is published.",
                          "Nothing is merged to master and nothing is published; /gc auto resume "
                          "is manual."))):
            with self.subTest(accepted=label):
                self.assertNotEqual(text, good)
                self.assert_accepted(self.check("auto_start", self.put(text)))

    def summaries_refused(self, queue, rows, summaries, simple=None):
        """Refused, and (no verbatim reason in the spec) the reason names the
        APPROVAL SUMMARIES section or a job's block."""
        options = {} if simple is None else {"simple": simple}
        result = self.report_check(self.report_text(rows, summaries=summaries, **options), queue)
        self.assert_refused(result)
        self.assertTrue("APPROVAL SUMMARIES" in result.stderr
                        or re.search(r"\bjob \S+", result.stderr, re.IGNORECASE), result.stderr)

    def test_auto_report_one_summary_block_per_ready_job_in_report_order(self):
        """Revision 9: APPROVAL SUMMARIES has exactly one block per ready job in the
        order of the report rows. Refused: a missing block, a block for a job that
        is not ready (or not in the report), a duplicate, blocks out of order,
        extra lines, `None ready.` while a job is ready."""
        queue = self.build_queue(extra_ready=("c",))
        rows = self.report_lines(queue)
        self.assertEqual(self.ready_jobs(rows), ["a", "c"], rows)
        a, b, c = (self.summary_block(job, f"title of {job}") for job in ("a", "b", "c"))
        self.assert_accepted(self.report_check(self.report_text(rows, summaries=a + c), queue))
        with self.subTest(control="Job line without optional text"):
            self.assert_accepted(self.report_check(self.report_text(
                rows, summaries=self.summary_block("a") + self.summary_block("c")), queue))
        for label, summaries in (
                ("missing block for c", a),
                ("missing block for a", c),
                ("no blocks at all", []),
                ("block for blocked job b", a + b + c),
                ("block for job b instead of c", a + b),
                ("block for a job not in the report", a + c + self.summary_block("z")),
                ("duplicate block for a", a + a + c),
                ("duplicate block for c", a + c + c),
                ("out of order", c + a),
                ("extra line before blocks", ["Two jobs are ready."] + a + c),
                ("extra line between blocks", a + ["Next block:"] + c),
                ("extra line after blocks", a + c + ["That is all."]),
                ("None ready. while jobs are ready", ["None ready."]),
                ("None ready. with the blocks", ["None ready."] + a + c),
                ("block cut to six lines", a[:-1] + c),
                ("Job line names no job", ["Job:"] + a[1:] + c),
                ("Job line without the colon", ["Job a"] + a[1:] + c),
                ("literal `Job ID:` before the id", ["Job ID: a"] + a[1:] + c)):
            with self.subTest(refused=label):
                self.summaries_refused(queue, rows, summaries)

    def test_auto_report_decide_p_row_is_ready(self):
        """Revision 9: a row whose step begins `decide P` (pending provisional
        choice) is ready: its block is required; with it the screen is accepted."""
        queue = self.build_queue(extra_ready=("c",), provisional=("c",))
        rows = self.report_lines(queue)
        steps = {row.split("|")[1].split()[0]: row.split("|")[3].strip() for row in rows[2:]}
        self.assertTrue(steps["c"].startswith("decide P1"), rows)
        self.assertTrue(steps["a"].startswith("approve for integration"), rows)
        self.assertEqual(self.ready_jobs(rows), ["a", "c"], rows)
        a, c = self.summary_block("a"), self.summary_block("c")
        self.assert_accepted(self.report_check(self.report_text(rows, summaries=a + c), queue))
        self.summaries_refused(queue, rows, a)
        self.summaries_refused(queue, rows, c + a)

    def test_auto_report_summary_labels_present_ordered_and_non_empty(self):
        """Revision 9: after `Job ID:` a block has, in order, lines starting `Bug: `,
        `Fix: `, `Test: `, `Criterion: `, `Limits: `, `Approval: `, each with
        non-empty text; a missing, misordered or empty label is refused."""
        queue = self.build_queue()
        rows = self.report_lines(queue)
        self.assertEqual(self.ready_jobs(rows), ["a"], rows)
        good = self.summary_block("a", "Clarify comment")
        self.assert_accepted(self.report_check(self.report_text(rows, summaries=good), queue))
        for index, label in enumerate(SUMMARY_LABELS, 1):
            with self.subTest(missing=label):
                self.summaries_refused(queue, rows, good[:index] + good[index + 1:])
            with self.subTest(replaced_by_note=label):
                self.summaries_refused(queue, rows, good[:index] + ["Note: something else."]
                                       + good[index + 1:])
            for empty in (label, label.rstrip(), label + "   "):
                with self.subTest(empty=repr(empty)):
                    self.summaries_refused(queue, rows, good[:index] + [empty] + good[index + 1:])
            with self.subTest(no_space_after_colon=label):
                text = good[index][len(label):]
                self.summaries_refused(queue, rows, good[:index] + [label.rstrip() + text]
                                       + good[index + 1:])
        for first, second in ((1, 2), (3, 4), (5, 6), (1, 6)):
            with self.subTest(swapped=(SUMMARY_LABELS[first - 1], SUMMARY_LABELS[second - 1])):
                swapped = list(good)
                swapped[first], swapped[second] = swapped[second], swapped[first]
                self.summaries_refused(queue, rows, swapped)
        with self.subTest(job_line="after the labels"):
            self.summaries_refused(queue, rows, good[1:] + good[:1])
        with self.subTest(job_line="missing"):
            self.summaries_refused(queue, rows, good[1:])

    def test_auto_report_approval_line_needs_exact_packet_and_merges_nothing(self):
        """Revision 9: the Approval line contains `exact packet` and `merges nothing`
        (case-insensitive); lacking either is refused."""
        queue = self.build_queue()
        rows = self.report_lines(queue)
        for text in ("Approval: accepts this exact packet and merges nothing.",
                     "Approval: Accepts This EXACT PACKET; it MERGES NOTHING.",
                     "Approval: merges nothing and approves the exact packet only."):
            with self.subTest(accepted=text):
                self.assert_accepted(self.report_check(self.report_text(
                    rows, summaries=self.summary_block("a", approval=text)), queue))
        for text in ("Approval: accepts this packet and merges nothing.",
                     "Approval: accepts this exact packet.",
                     "Approval: accepts this exact packet and merges no code.",
                     "Approval: accepts the exact-packet and merges nothing.",
                     "Approval: accepts this exactly packed set and merges nothing.",
                     "Approval: yes."):
            with self.subTest(refused=text):
                self.summaries_refused(queue, rows, self.summary_block("a", approval=text))
        with self.subTest(refused="phrases only on another line"):
            self.summaries_refused(queue, rows, self.summary_block(
                "a", limits="Limits: exact packet; merges nothing.", approval="Approval: yes."))

    def test_auto_report_none_ready_when_no_job_is_ready(self):
        """Revision 9: with no ready job the section is the single line `None ready.`;
        a block, an empty section or other text is refused."""
        queue = self.build_queue()
        criterion = self.queue_base / "criterion-2.txt"
        criterion.write_text("Criterion two, changed.\n")
        self.queue_tool("criterion", "a", "--file", str(criterion), queue=queue)
        rows = self.report_lines(queue)
        self.assertEqual(self.ready_jobs(rows), [], rows)
        simple = ("Job a must be checked again.", "Job b needs a design decision.")
        self.assert_accepted(self.report_check(self.report_text(
            rows, simple=simple, summaries=["None ready."]), queue))
        for label, summaries in (
                ("block for the stale job", self.summary_block("a")),
                ("empty section", []),
                ("None ready. twice", ["None ready.", "None ready."]),
                ("None ready. with a note", ["None ready.", "Job a must be checked again."]),
                ("without the full stop", ["None ready"]),
                ("other words", ["Nothing is ready."]),
                ("trailing text", ["None ready. Check again later."])):
            with self.subTest(refused=label):
                self.summaries_refused(queue, rows, summaries, simple=simple)

    # ---- Revision 9 clarifications, prefixed stop command, re-read 8 -----------

    def test_auto_start_must_also_say_manual(self):
        """After re-read 8 (b) and after test round 2: auto_start also contains
        `manual` or `manually` as a whole word (any case), with `usage limit` and
        `/gc auto resume`; same refusal message."""
        reason = "auto_start must say a usage limit stops the queue and resume is manual"
        good = self.start_text()
        self.assertIn("resume is manual with", good)
        for label, text in (
                ("no manual", good.replace("resume is manual with", "resume with")),
                ("nonmanual alone", good.replace("resume is manual with",
                                                 "resume is nonmanual with")),
                ("manuals", good.replace("resume is manual with", "resume per manuals with")),
                ("manualy misspelt", good.replace("resume is manual with",
                                                  "resume manualy with"))):
            with self.subTest(refused=label):
                result = self.check("auto_start", self.put(text))
                self.assert_refused(result)
                self.assertIn(reason, result.stderr)
        for label, text in (("upper case", good.replace("manual", "MANUAL")),
                            ("manually", good.replace("resume is manual with",
                                                      "resume manually with")),
                            ("Manually capitalized", good.replace("resume is manual with",
                                                                  "resume Manually with")),
                            ("in parentheses", good.replace("resume is manual with",
                                                            "resume (manual) with")),
                            ("before a full stop", good.replace(
                                "resume is manual with /gc auto resume.",
                                "resume with /gc auto resume is manual.")),
                            ("elsewhere", good.replace("resume is manual with", "resume with")
                             .replace("Nothing is merged", "Manual steps only. Nothing is merged"))):
            with self.subTest(accepted=label):
                self.assert_accepted(self.check("auto_start", self.put(text)))

    PREFIX = 'TMPDIR="$(cd "${TMPDIR:-/tmp}" && pwd -P)"'

    def test_auto_start_isolated_stop_command_accepted_exactly(self):
        """Stop command form (option 1; supersedes the prefixed form). Revision 16: STOP
        accepts exactly `/ABS/libtbx.python -I -B AQ stop`, or the same with a build's
        /ABS/conda_base/bin/python3 in place of libtbx.python (repair 1: absolute paths only,
        so that the line works whatever PATH is; a bare libtbx.python is refused); the earlier
        `python3 AQ stop` and `python3 -I -B AQ stop` (a python3 that may lack psutil, so the
        stop would fail with exit 18), the TMPDIR-prefixed form and near-misses are refused
        with "auto_start needs the stop command of this procedure's auto_queue.py"."""
        reason = "auto_start needs the stop command of this procedure's auto_queue.py"
        isolated = f"{STOP_PYTHON} {CHECKER_AUTO_QUEUE} stop"
        conda = f"{STOP_BUILD}/conda_base/bin/python3"
        for label, line in (("absolute libtbx.python", isolated),
                            ("libtbx.python through a symbolic link",
                             f"{STOP_BUILD}/link/bin/libtbx.python -I -B {CHECKER_AUTO_QUEUE} stop"),
                            ("conda python3", f"{conda} -I -B {CHECKER_AUTO_QUEUE} stop"),
                            ("conda python3.11", f"{conda}.11 -I -B {CHECKER_AUTO_QUEUE} stop")):
            with self.subTest(accepted=label):
                self.assert_accepted(self.check("auto_start", self.put(
                    self.start_text(stop_line=line))))
        unresolved = Path(checker.__file__).parent / ".." / ".." / "auto_queue.py"
        for label, line in (
                ("python3 -I -B (Revision 16)", f"python3 -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("nonexistent libtbx.python (repair 2)",
                 f"/nonexistent/zz/libtbx.python -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("nonexistent conda python3 (repair 2)",
                 f"{STOP_BUILD}/conda_base/bin/python3.99 -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("non-executable libtbx.python (repair 2)",
                 f"{STOP_BUILD}/noexec/build/bin/libtbx.python -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("directory named libtbx.python (repair 2)",
                 f"{STOP_BUILD}/dir/bin/libtbx.python -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("python3 plain (Revision 16)", f"python3 {CHECKER_AUTO_QUEUE} stop"),
                ("libtbx.python plain", f"libtbx.python {CHECKER_AUTO_QUEUE} stop"),
                ("bare libtbx.python -I -B (repair 1)", f"libtbx.python -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("relative libtbx.python", f"build/bin/libtbx.python -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("relative conda python", f"conda_base/bin/python3 -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("other python by path", f"/usr/bin/python3 -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("TMPDIR prefix with -I -B", f"{self.PREFIX} {isolated}"),
                ("TMPDIR prefix, plain", f"{self.PREFIX} python3 {CHECKER_AUTO_QUEUE} stop"),
                ("-B -I order", isolated.replace("-I -B", "-B -I")),
                ("-I only", isolated.replace("-I -B", "-I")),
                ("-B only", isolated.replace("-I -B", "-B")),
                ("-IB combined", isolated.replace("-I -B", "-IB")),
                ("extra space after the interpreter",
                 isolated.replace("libtbx.python ", "libtbx.python  ", 1)),
                ("extra space between -I and -B", isolated.replace("-I -B", "-I  -B")),
                ("extra space before stop", isolated.replace(" stop", "  stop")),
                ("another path", f"{STOP_PYTHON} /opt/guided_coding/auto_queue.py stop"),
                ("unresolved path", f"{STOP_PYTHON} {unresolved} stop"),
                ("relative path", f"{STOP_PYTHON} auto_queue.py stop"),
                ("leading space", " " + isolated),
                ("trailing space", isolated + " "),
                ("python instead of libtbx.python", f"python -I -B {CHECKER_AUTO_QUEUE} stop"),
                ("env prefix", f"env {isolated}")):
            with self.subTest(refused=label):
                self.assertNotEqual(line, isolated)
                result = self.check("auto_start", self.put(self.start_text(stop_line=line)))
                self.assert_refused(result)
                self.assertIn(reason, result.stderr)

    def test_auto_report_names_psutil_when_the_helper_needs_it(self):
        """Revision 16: when the helper stops for lack of psutil (exit 18), auto_report is
        refused with the helper's message, not only "report failed". A copy of this
        screen_check.py beside a stand-in auto_queue.py that answers as the helper does
        without psutil; control: a stand-in that fails otherwise gives the old message."""
        for code, expected in ((18, "auto_queue.py report failed: REFUSED: requires psutil"),
                               (1, "auto_queue.py report failed")):
            with self.subTest(exit=code):
                kit = self.root / f"kit-{code}"
                (kit / "payload" / "tools").mkdir(parents=True)
                shutil.copy(checker.__file__, kit / "payload" / "tools" / "screen_check.py")
                (kit / "auto_queue.py").write_text(
                    "import sys\nprint('REFUSED: requires psutil: run with the cctbx conda "
                    "Python (libtbx.python) or pip install psutil', file=sys.stderr)\n"
                    "sys.exit(%d)\n" % code)
                queue = self.root / "queue"
                queue.mkdir(exist_ok=True)
                text = self.report_text(["| a | local checks passed | approve |"])
                result = subprocess.run(
                    [sys.executable, "-I", "-B", str(kit / "payload" / "tools" / "screen_check.py"),
                     "check", "auto_report", str(self.put(text)), "--queue", str(queue)],
                    cwd=self.root, env=self.environment, text=True, capture_output=True)
                self.assert_refused(result)
                self.assertEqual(result.stderr.strip(), "ERROR: " + expected
                                 + ("" if code == 1 else ": run with the cctbx conda Python "
                                    "(libtbx.python) or pip install psutil"))

    def test_auto_report_passed_and_stopped_allowed_in_approval_summaries(self):
        """Clarification: the revision 2 checks on `passed` and `stopped` apply to IN
        SIMPLE WORDS and NEXT only: in APPROVAL SUMMARIES both words are accepted even
        when no row says passed and the report is not Stopped."""
        queue = self.build_queue()
        rows = self.report_lines(queue)
        self.assertFalse(rows[0].startswith("**"), rows)
        block = self.summary_block("a", bug="Bug: the run stopped early and nothing passed.",
                                   test="Test: failed before, passed after; never stopped.")
        self.assert_accepted(self.report_check(self.report_text(rows, summaries=block), queue))
        # The same words in IN SIMPLE WORDS are still refused here (control for the
        # section boundary): `stopped` without a **Stopped** report.
        self.assert_refused(self.report_check(self.report_text(
            rows, summaries=block, simple=("Job a passed.", "The queue stopped.")), queue))

    def test_auto_report_labels_and_job_word_are_case_sensitive(self):
        """Clarification: labels are case-sensitive and the first line is the word
        `Job`, a space, the id and a colon: lowercase labels and `job a:` are refused."""
        queue = self.build_queue()
        rows = self.report_lines(queue)
        good = self.summary_block("a")
        self.assert_accepted(self.report_check(self.report_text(rows, summaries=good), queue))
        for index, label in enumerate(SUMMARY_LABELS, 1):
            for changed in (label.lower(), label.upper()):
                with self.subTest(label=changed):
                    self.summaries_refused(queue, rows, good[:index]
                                           + [changed + good[index][len(label):]]
                                           + good[index + 1:])
        for first in ("job a:", "JOB a:", "Job A:", "Job  a:", "Job a :", "Joba:"):
            with self.subTest(first=first):
                self.summaries_refused(queue, rows, [first] + good[1:])

    def test_auto_report_job_id_with_pattern_characters(self):
        """Revision 9: the id `a.b+` is matched literally: `Job a.b+:` is accepted;
        `Job aXb+:`, `Job a.bb:` and `Job a.b:` are refused."""
        queue = self.build_queue(extra_ready=("a.b+",))
        rows = self.report_lines(queue)
        self.assertEqual(self.ready_jobs(rows), ["a", "a.b+"], rows)
        a = self.summary_block("a")
        self.assert_accepted(self.report_check(self.report_text(
            rows, summaries=a + self.summary_block("a.b+", "pattern id")), queue))
        for other in ("aXb+", "a.bb", "a.b", "a.b++"):
            with self.subTest(block=other):
                self.summaries_refused(queue, rows, a + self.summary_block(other))

    def test_auto_report_only_exact_header_and_separator_rows_skipped(self):
        """After re-read 8 (a): only the exact header and separator rows are skipped;
        a ready job whose id is `Job` is found and needs its block. (A job id `---`
        is refused by init since AUTO_QUEUE "After test round 2", so it cannot occur.)"""
        queue = self.build_queue(extra_ready=("Job",))
        rows = self.report_lines(queue)
        self.assertEqual(self.ready_jobs(rows), ["a", "Job"], rows)
        a, job = self.summary_block("a"), self.summary_block("Job", "named Job")
        self.assert_accepted(self.report_check(self.report_text(rows, summaries=a + job), queue))
        self.summaries_refused(queue, rows, a)
        self.summaries_refused(queue, rows, job + a)

    def check_summaries_follow_rows(self, queue, rows, simple):
        """Accepted with the blocks the rows call for; refused with the opposite."""
        ready = self.ready_jobs(rows)
        self.assert_accepted(self.report_check(self.report_text(rows, simple=simple), queue))
        opposite = (["None ready."] if ready
                    else self.summary_block(rows[-1].split("|")[1].split()[0]))
        self.summaries_refused(queue, rows, opposite, simple=simple)
        return ready

    def test_auto_report_summaries_with_stop_and_installation_lines(self):
        """Revision 9 with the report's leading lines: under **Stopping**, **Stopped**
        and `**Test installation unavailable:**` the summaries still follow the rows
        (one block per ready row, else `None ready.`)."""
        queue = self.build_queue(queue_id="night-stopping")
        self.queue_tool("worker", "start", "b", "--id", "worker-b", queue=queue)
        self.queue_tool("stop", queue=queue, expect=10)
        rows = self.report_lines(queue)
        self.assertTrue(rows[0].startswith("**Stopping**"), rows)
        with self.subTest(lead="Stopping"):
            self.check_summaries_follow_rows(queue, rows, ("The queue is stopping.",))
        self.queue_tool("worker", "end", "b", "--id", "worker-b", "--outcome", "cancelled",
                        queue=queue)
        self.queue_tool("stop", queue=queue)
        rows = self.report_lines(queue)
        self.assertEqual(rows[0], "**Stopped**", rows)
        with self.subTest(lead="Stopped"):
            self.assertEqual(self.check_summaries_follow_rows(
                queue, rows, ("The queue stopped.",)), ["a"])
        queue = self.build_queue(queue_id="night-install")
        script = self.queue_base / "fingerprint.sh"
        script.write_text(stand_in(["FINGERPRINT x"]))
        script.chmod(0o755)
        self.queue_tool("install", "check", "--script", str(script), "--expect", "y",
                        queue=queue, expect=12)
        rows = self.report_lines(queue)
        self.assertTrue(rows[0].startswith("**Test installation unavailable:**"), rows)
        with self.subTest(lead="Test installation unavailable"):
            self.check_summaries_follow_rows(
                queue, rows, ("The test installation does not match.",))

    def test_existing_stop_kind_still_accepted_by_the_command_line(self):
        """Existing kinds behave as before (control for the CLI used here)."""
        stop = self.put("STOP | comment-1\nWHY\nThe intended meaning is uncertain.\n"
                        "PRESERVED\nThe isolated draft is saved; no integration.\n"
                        "NEEDED\nDecide intended meaning, or cancel this change.\n"
                        "ACTION: DECIDE / CANCEL\n")
        self.assert_accepted(self.check("stop", stop))


if __name__ == "__main__":
    unittest.main()
