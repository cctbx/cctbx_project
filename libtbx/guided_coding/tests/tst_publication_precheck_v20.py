"""Controls for the claims that publication_precheck.py makes (SPEC_A section 2).

The tool is only ever run through subprocess; nothing here imports it. Every
`check` call is bracketed by a snapshot of the work repository (HEAD, refs,
local config, working tree, FETCH_HEAD) and of the bare remote, which must be
identical before and after: the tool is read-only unless --fetch is named.
All Git state lives under a temporary directory, and the Developer's real Git
configuration is kept out with GIT_CONFIG_GLOBAL and GIT_CONFIG_NOSYSTEM.
"""

import itertools
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

TOOL = Path(__file__).resolve().parent.parent / "payload" / "tools" / "publication_precheck.py"
HEX = "0123456789abcdef0123456789abcdef01234567"
KEYS = ("repository", "remote", "remote_url", "base", "commit", "tree", "refspec")



def fake_git_script(real_git, calls, trigger, out, err, code):
    """Revision 15: the text of a Python test double for git (formerly /bin/sh). A call whose
    argument words include every word of `trigger` appends its words (space-joined) to the
    file `calls`, prints the lines `out` to stdout and `err` to stderr and exits `code`; every
    other call is replaced by the real git with the same arguments."""
    return ("#!%s\n" % sys.executable
            + "# test double: one intercepted git question; every other call is the real git\n"
            + "import os, sys\n"
            + "args = sys.argv[1:]\n"
            + "if all(word in args for word in %r):\n" % (tuple(trigger),)
            + "    with open(%r, 'a') as log:\n" % str(calls)
            + "        log.write(' '.join(args) + '\\n')\n"
            + "    for line in %r:\n" % (list(out),)
            + "        print(line)\n"
            + "    for line in %r:\n" % (list(err),)
            + "        print(line, file=sys.stderr)\n"
            + "    sys.exit(%d)\n" % code
            + "os.execv(%r, [%r] + args)\n" % (real_git, real_git))


def git_environment(home):
    """An environment in which git sees no personal or system configuration."""
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith("GIT_") and key != "XDG_CONFIG_HOME"}
    environment.update(GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1",
                       HOME=str(home), GIT_TERMINAL_PROMPT="0", LC_ALL="C")
    return environment


# Revision 16: `check` and `outgoing` read a git repository; without git on PATH their tests
# skip with this reason (`vet` needs no git), and NoGitRefusal checks the tool's refusal.
GIT_SKIP = ("git not on PATH: the check and outgoing tests need git repositories; the tool "
            "refuses those commands without git (NoGitRefusal)")


@unittest.skipIf(shutil.which("git") is None, GIT_SKIP)
class PrecheckFixture(unittest.TestCase):
    """A work repository one commit ahead of a bare remote it tracks."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        home = self.root / "home"
        home.mkdir()
        self.env = git_environment(home)
        self.bare = self.root / "remote.git"
        self.repo = self.root / "work"
        self.git("init", "-q", "--bare", str(self.bare), cwd=self.root)
        self.git("symbolic-ref", "HEAD", "refs/heads/master", cwd=self.bare)
        self.repo.mkdir()
        self.git("init", "-q", cwd=self.repo)
        self.git("symbolic-ref", "HEAD", "refs/heads/master", cwd=self.repo)
        self.configure_user(self.repo)
        self.git("remote", "add", "origin", str(self.bare), cwd=self.repo)
        self.base = self.commit_file(self.repo, "a.txt", "base\n", "base")
        self.git("push", "-q", "origin", "master", cwd=self.repo)
        self.git("fetch", "-q", "origin", cwd=self.repo)
        self.commit = self.commit_file(self.repo, "a.txt", "change\n", "change")
        self.tree = self.git("rev-parse", f"{self.commit}^{{tree}}", cwd=self.repo)
        self.assertEqual(self.git("rev-parse", "refs/remotes/origin/master", cwd=self.repo),
                         self.base)
        self.outgoing = self.root / "OUTGOING.txt"
        self.outgoing.write_text(self.block())

    def git(self, *words, cwd):
        result = subprocess.run(["git", *words], cwd=cwd, env=self.env,
                                text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.strip()

    def configure_user(self, repo):
        self.git("config", "user.name", "Test Writer", cwd=repo)
        self.git("config", "user.email", "writer@example.invalid", cwd=repo)

    def commit_file(self, repo, name, content, message):
        (repo / name).write_text(content)
        self.git("add", name, cwd=repo)
        self.git("commit", "-q", "-m", message, cwd=repo)
        return self.git("rev-parse", "HEAD", cwd=repo)

    def block(self, **fields):
        """One OUTGOING.txt block (SPEC_A 1.2) bound to the fixture repository."""
        values = dict(repository="phenix", remote="origin", remote_url=str(self.bare),
                      base=self.base, commit=self.commit, tree=self.tree)
        values.update(fields)
        values.setdefault("refspec", f"{values['commit']}:refs/heads/master")
        return "".join(f"{key}: {values[key]}\n" for key in KEYS)

    def snapshot(self):
        work = [self.git(*words, cwd=self.repo) for words in (
            ("rev-parse", "HEAD"), ("for-each-ref",), ("config", "--list", "--local"),
            ("status", "--porcelain", "--untracked-files=all"))]
        fetch_head = self.repo / ".git" / "FETCH_HEAD"
        work.append((fetch_head.read_bytes(), fetch_head.stat().st_mtime_ns)
                    if fetch_head.exists() else None)
        remote = [self.git(*words, cwd=self.bare) for words in (
            ("for-each-ref",), ("config", "--list", "--local"))]
        return work, remote

    def check(self, *options, fetch=False):
        """Run `check` and prove it changed nothing (with --fetch: nothing but
        remote-tracking refs and FETCH_HEAD)."""
        before = self.snapshot()
        result = subprocess.run(
            [sys.executable, "-I", "-B", str(TOOL), "check", str(self.repo),
             str(self.outgoing), *options],
            env=self.env, cwd=self.root, text=True, capture_output=True)
        after = self.snapshot()
        if not fetch:
            self.assertEqual(after, before, "check changed Git state")
        else:
            self.assertEqual(after[1], before[1], "fetch changed the remote")
            for index in (0, 2, 3):
                self.assertEqual(after[0][index], before[0][index])
            self.assertEqual([line for line in after[0][1].splitlines() if "refs/remotes/" not in line],
                             [line for line in before[0][1].splitlines() if "refs/remotes/" not in line])
        return result

    def advance_remote(self):
        """Another clone pushes a commit to master; the work repository is not told."""
        other = self.root / "other"
        self.git("clone", "-q", str(self.bare), str(other), cwd=self.root)
        self.configure_user(other)
        self.other = self.commit_file(other, "b.txt", "other\n", "other")
        self.git("push", "-q", "origin", "master", cwd=other)
        self.assertEqual(self.git("rev-parse", "refs/heads/master", cwd=self.bare), self.other)
        return self.other


class CheckCommand(PrecheckFixture):
    def test_one_block_passes_and_prints_the_command(self):
        """SPEC_A 2.1 step 6: PRECHECK OK line, then the exact push command on
        the next line; --repository may name the single block."""
        ok_line = f"PRECHECK OK phenix {self.commit} -> origin {self.commit}:refs/heads/master"
        command = (f"git -C {self.repo} -c push.followTags=false push --no-follow-tags "
                   f"origin {self.commit}:refs/heads/master")
        for options in ((), ("--repository", "phenix")):
            with self.subTest(options=options):
                result = self.check(*options)
                self.assertEqual(result.returncode, 0, result.stderr)
                lines = result.stdout.splitlines()
                self.assertIn(ok_line, lines)
                self.assertEqual(lines[lines.index(ok_line) + 1], command)
                self.assertNotIn("ERROR", result.stderr)

    def test_record_parse_refusals_match_screen_check(self):
        """SPEC_A 2.1 step 1: OUTGOING is parsed as in 1.2 with the same refusals."""
        for label, text, reason in (
                ("empty", "", "missing or empty publication evidence"),
                ("missing key", self.block().replace(f"tree: {self.tree}\n", ""),
                 "outgoing record needs repository, remote, remote_url, base, commit, tree and refspec"),
                ("not hex", self.block(base="A" * 40), "invalid Git identity in outgoing record"),
                ("bad refspec", self.block(refspec=f"{self.commit}:refs/heads/*"),
                 "outgoing refspec must be"),
                ("duplicate", self.block() + self.block(), "duplicate repository in outgoing record")):
            with self.subTest(record=label):
                self.outgoing.write_text(text)
                result = self.check()
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn("ERROR:", result.stderr)
                self.assertIn(reason, result.stderr)
        self.outgoing.write_bytes(b"repository: phenix\n\xff\xfe\n")
        result = self.check()
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("outgoing record is not UTF-8", result.stderr)

    def test_several_blocks_require_a_repository_name(self):
        """SPEC_A 2.1 step 1: two blocks without --repository are refused; the
        named block is checked; an unknown name cannot succeed."""
        self.outgoing.write_text(self.block() + "\n" + self.block(repository="cctbx"))
        result = self.check()
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("several repositories; name one with --repository", result.stderr)
        named = self.check("--repository", "phenix")
        self.assertEqual(named.returncode, 0, named.stderr)
        self.assertIn(f"PRECHECK OK phenix {self.commit}", named.stdout)
        unknown = self.check("--repository", "nowhere")
        self.assertEqual(unknown.returncode, 2, unknown.stdout)
        self.assertIn("ERROR:", unknown.stderr)

    def test_destination_must_equal_remote_url(self):
        """SPEC_A 2.1 step 2: `remote get-url --all --push` must print exactly
        one line equal to remote_url; the lines found are reported."""
        self.outgoing.write_text(self.block(remote_url="ssh://git@example.invalid/phenix.git"))
        result = self.check()
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("push destination differs from outgoing record", result.stderr)
        self.assertIn(str(self.bare), result.stderr)
        self.outgoing.write_text(self.block())
        # Two push URLs (both configured as pushurl) are two lines: step 2
        # refuses before the config step sees the pushurl setting.
        self.git("remote", "set-url", "--add", "--push", "origin", str(self.bare), cwd=self.repo)
        self.git("remote", "set-url", "--add", "--push", "origin", str(self.root / "x.git"),
                 cwd=self.repo)
        result = self.check()
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("push destination differs from outgoing record", result.stderr)
        self.git("config", "--unset-all", "remote.origin.pushurl", cwd=self.repo)
        self.assertEqual(self.check().returncode, 0)
        self.outgoing.write_text(self.block(remote="upstream"))
        result = self.check()
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("ERROR:", result.stderr)

    def test_forbidden_push_settings_are_refused(self):
        """SPEC_A 2.1 step 3: each listed setting in the repository's local
        config is refused, naming the setting and its origin."""
        for key, value in (
                ("url.https://example.invalid/.insteadOf", "ssh://nowhere.invalid/"),
                ("url.https://example.invalid/.pushInsteadOf", "ssh://nowhere.invalid/"),
                ("remote.origin.pushurl", str(self.bare)),
                ("remote.origin.mirror", "false"),
                ("remote.origin.mirror", "true"),
                ("remote.origin.push", "refs/heads/master:refs/heads/master"),
                ("remote.origin.tagOpt", "--no-tags"),
                ("push.followTags", "true")):
            with self.subTest(setting=key, value=value):
                self.git("config", key, value, cwd=self.repo)
                result = self.check()
                self.git("config", "--unset", key, cwd=self.repo)
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn("forbidden push setting: ", result.stderr)
                self.assertIn(key.lower(), result.stderr.lower())
                self.assertIn(".git/config", result.stderr)
        self.assertEqual(self.check().returncode, 0)

    def test_other_push_keys_are_notes_and_followtags_false_is_fine(self):
        """SPEC_A 2.1 step 3: push.followTags=false passes; another push.* key
        is noted, not refused."""
        self.git("config", "push.followTags", "false", cwd=self.repo)
        result = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.git("config", "push.default", "simple", cwd=self.repo)
        result = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("push.default", (result.stdout + result.stderr).lower())
        self.assertIn(f"PRECHECK OK phenix {self.commit}", result.stdout)

    def test_binding_mismatches_are_refused(self):
        """SPEC_A 2.1 step 4: wrong tree, wrong first parent, remote-tracking
        ref not at base, missing tracking ref, or unknown commit."""
        reason = "outgoing binding mismatch"
        base_tree = self.git("rev-parse", f"{self.base}^{{tree}}", cwd=self.repo)
        with self.subTest(binding="tree"):
            self.outgoing.write_text(self.block(tree=base_tree))
            result = self.check()
            self.assertEqual(result.returncode, 2, result.stdout)
            self.assertIn(reason, result.stderr)
        with self.subTest(binding="first parent"):
            third = self.commit_file(self.repo, "a.txt", "third\n", "third")
            third_tree = self.git("rev-parse", f"{third}^{{tree}}", cwd=self.repo)
            # base is the grandparent: origin/master still equals base, so
            # only the first-parent binding is wrong.
            self.outgoing.write_text(self.block(commit=third, tree=third_tree))
            result = self.check()
            self.assertEqual(result.returncode, 2, result.stdout)
            self.assertIn(reason, result.stderr)
        with self.subTest(binding="remote-tracking ref"):
            self.outgoing.write_text(self.block())
            self.git("update-ref", "refs/remotes/origin/master", self.commit, cwd=self.repo)
            result = self.check()
            self.git("update-ref", "refs/remotes/origin/master", self.base, cwd=self.repo)
            self.assertEqual(result.returncode, 2, result.stdout)
            self.assertIn(reason, result.stderr)
        with self.subTest(binding="no remote-tracking ref for the branch"):
            self.outgoing.write_text(self.block(refspec=f"{self.commit}:refs/heads/other"))
            result = self.check()
            self.assertEqual(result.returncode, 2, result.stdout)
            self.assertIn(reason, result.stderr)
        with self.subTest(binding="unknown commit"):
            self.outgoing.write_text(self.block(commit="f" * 40))
            result = self.check()
            self.assertEqual(result.returncode, 2, result.stdout)
            self.assertIn(reason, result.stderr)
        self.outgoing.write_text(self.block())
        self.assertEqual(self.check().returncode, 0)

    def test_dry_run_against_the_bare_remote_fast_forwards_one_ref(self):
        """SPEC_A 2.1 step 5: a fast-forward dry run of one branch passes and
        pushes nothing."""
        result = self.check("--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f"PRECHECK OK phenix {self.commit} -> origin "
                      f"{self.commit}:refs/heads/master", result.stdout)
        self.assertEqual(self.git("rev-parse", "refs/heads/master", cwd=self.bare), self.base)

    def test_dry_run_refused_when_the_remote_has_moved(self):
        """SPEC_A 2.1 steps 4-5: without --fetch the stale tracking ref still
        satisfies the bindings and no fetch happens; the dry run refuses the
        non-fast-forward; --fetch then exposes the binding mismatch."""
        moved = self.advance_remote()
        still = self.check()
        self.assertEqual(still.returncode, 0, still.stderr)
        self.assertEqual(self.git("rev-parse", "refs/remotes/origin/master", cwd=self.repo),
                         self.base)
        dry = self.check("--dry-run")
        self.assertEqual(dry.returncode, 2, dry.stdout)
        self.assertIn("ERROR:", dry.stderr)
        self.assertIn("dry run", dry.stderr)
        self.assertNotIn("PRECHECK OK", dry.stdout)
        self.assertEqual(self.git("rev-parse", "refs/heads/master", cwd=self.bare), moved)
        fetched = self.check("--fetch", fetch=True)
        self.assertEqual(fetched.returncode, 2, fetched.stdout)
        self.assertIn("outgoing binding mismatch", fetched.stderr)
        self.assertEqual(self.git("rev-parse", "refs/remotes/origin/master", cwd=self.repo),
                         moved)

    # ---- SPEC_A 6: refinements after the internal reviewer's reading --------

    def test_followtags_truthy_values_in_any_case_are_refused(self):
        """SPEC_A 6 (2.1 step 3): push.followTags set to yes, on or 1 (any
        case, like true) in the repository's local config is a forbidden push
        setting; false and no pass."""
        for value in ("yes", "on", "1", "On", "YES", "ON", "Yes", "true", "TRUE"):
            with self.subTest(value=value):
                self.git("config", "push.followTags", value, cwd=self.repo)
                result = self.check()
                self.git("config", "--unset", "push.followTags", cwd=self.repo)
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn("ERROR:", result.stderr)
                self.assertIn("forbidden push setting", result.stderr)
                self.assertIn("push.followtags", result.stderr.lower())
                self.assertNotIn("PRECHECK OK", result.stdout)
        for value in ("false", "no"):
            with self.subTest(value=value):
                self.git("config", "push.followTags", value, cwd=self.repo)
                result = self.check()
                self.git("config", "--unset", "push.followTags", cwd=self.repo)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(f"PRECHECK OK phenix {self.commit}", result.stdout)
                self.assertNotIn("forbidden push setting", result.stderr)
        self.assertEqual(self.check().returncode, 0)

    def test_dry_run_exit_status_is_checked_before_its_output(self):
        """SPEC_A 8 (2.1 step 5): Git's exit status is checked first. A fake
        git first on PATH answers the dry-run push with "To <url>", one
        plausible fast-forward line for refs/heads/master and "Done" but exits
        1: check --dry-run is refused with "dry run failed (exit 1)" and prints
        no PRECHECK OK. The same fake exiting 0 is the positive control. Every
        other git question is passed through to the real git."""
        real_git = shutil.which("git", path=self.env.get("PATH", os.defpath))
        self.assertTrue(real_git, "real git not found on PATH")
        fake_bin = self.root / "fakebin"
        fake_bin.mkdir()
        fake = fake_bin / "git"
        calls = self.root / "dry-run-calls.txt"

        def install(status):
            # Revision 15: a Python test double (formerly /bin/sh) with the same output.
            fake.write_text(fake_git_script(
                real_git, calls, ("push", "--dry-run"),
                [f"To {self.bare}",
                 f" \t{self.commit}:refs/heads/master\t{self.base[:7]}..{self.commit[:7]}",
                 "Done"], [f"fake git: dry run exits {status}"], status))
            fake.chmod(0o755)

        env = dict(self.env, PATH=str(fake_bin) + os.pathsep + self.env.get("PATH", os.defpath))

        def check_dry_run():
            before = self.snapshot()
            result = subprocess.run(
                [sys.executable, "-I", "-B", str(TOOL), "check", str(self.repo),
                 str(self.outgoing), "--dry-run"],
                env=env, cwd=self.root, text=True, capture_output=True)
            self.assertEqual(self.snapshot(), before, "check changed Git state")
            return result

        install(1)
        result = check_dry_run()
        self.assertTrue(calls.exists(), "the fake git never answered the dry-run push")
        recorded = calls.read_text()
        for word in ("push", "--no-follow-tags", "--dry-run", "--porcelain", "origin",
                     f"{self.commit}:refs/heads/master"):
            self.assertIn(word, recorded)
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("ERROR:", result.stderr)
        self.assertIn("dry run failed (exit 1)", result.stderr)
        self.assertNotIn("PRECHECK OK", result.stdout)
        calls.unlink()
        install(0)
        result = check_dry_run()
        self.assertTrue(calls.exists(), "the fake git never answered the dry-run push")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f"PRECHECK OK phenix {self.commit} -> origin "
                      f"{self.commit}:refs/heads/master", result.stdout)
        self.assertEqual(self.git("rev-parse", "refs/heads/master", cwd=self.bare), self.base)

    # ---- SPEC_A 10: every Git query's exit status is checked first ----------

    def fake_git(self, name, trigger, out, err, code):
        """Put a fake `git` first on PATH for the tool subprocess only. A call
        whose words include every word of `trigger` logs itself, prints the
        lines `out` and `err` and exits `code`; every other call is exec'd to the
        real git, so the test does not guess which other questions the tool
        asks. Returns (environment, log file). Revision 15: a Python test
        double (fake_git_script; formerly /bin/sh)."""
        real_git = shutil.which("git", path=self.env.get("PATH", os.defpath))
        self.assertTrue(real_git, "real git not found on PATH")
        fake_bin = self.root / f"fakebin-{name}"
        fake_bin.mkdir()
        calls = fake_bin / "intercepted.txt"
        script = fake_git_script(real_git, calls, trigger, out, err, code)
        (fake_bin / "git").write_text(script)
        (fake_bin / "git").chmod(0o755)
        env = dict(self.env, PATH=str(fake_bin) + os.pathsep + self.env.get("PATH", os.defpath))
        return env, calls

    def check_with(self, env, *options):
        """Run `check` under `env` inside the fixture's Git-state snapshot bracket."""
        before = self.snapshot()
        result = subprocess.run(
            [sys.executable, "-I", "-B", str(TOOL), "check", str(self.repo),
             str(self.outgoing), *options],
            env=env, cwd=self.root, text=True, capture_output=True)
        self.assertEqual(self.snapshot(), before, "check changed Git state")
        return result

    def test_configuration_query_failure_is_refused_before_its_output(self):
        """SPEC_A 10 (2.1 step 3): the config --show-origin --get-regexp query
        is accepted only with exit 0, or exit 1 with empty stdout and stderr.
        A fake git answering it with exit 128 (printing nothing, or one benign
        `push.default simple` line) makes check refuse with "git configuration
        query failed (exit 128)" and accept no NOTE; the real git with no
        matching settings (exit 1, empty output) is the positive control."""
        failure = "fatal: injected configuration query failure"
        trigger = ("config", "--show-origin", "--get-regexp")
        for label, out in (("nothing on stdout", []),
                           ("one benign line on stdout", ["file:fixture\tpush.default simple"])):
            with self.subTest(form=label):
                env, calls = self.fake_git(label.replace(" ", "-"), trigger, out, [failure], 128)
                result = self.check_with(env)
                self.assertTrue(calls.exists(), "the fake git never saw the configuration query")
                self.assertIn("--get-regexp", calls.read_text())
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn("ERROR:", result.stderr)
                self.assertIn("git configuration query failed (exit 128)", result.stderr)
                self.assertIn("injected configuration query failure", result.stderr)
                self.assertNotIn("PRECHECK OK", result.stdout)
                self.assertNotIn("NOTE", result.stdout)
                self.assertNotIn("push.default", result.stdout)
        # positive control: the real git has no push.* setting in the fixture
        probe = subprocess.run(["git", "-C", str(self.repo), "config", "--show-origin",
                                "--get-regexp", r"^push\."],
                               env=self.env, text=True, capture_output=True)
        self.assertEqual((probe.returncode, probe.stdout, probe.stderr), (1, "", ""))
        control = self.check()
        self.assertEqual(control.returncode, 0, control.stderr)
        self.assertIn(f"PRECHECK OK phenix {self.commit}", control.stdout)

    def test_binding_query_failure_is_refused_despite_correct_output(self):
        """SPEC_A 10 (2.1 step 4): a fake git that answers rev-parse --verify
        -q <commit>^{commit} with the correct commit on stdout but exits 128
        makes check refuse with "outgoing binding mismatch" and "(exit 128)";
        the real git is the positive control."""
        failure = "fatal: injected object lookup failure"
        env, calls = self.fake_git(
            "binding", ("rev-parse", "--verify", "-q", f"{self.commit}^{{commit}}"),
            [self.commit], [failure], 128)
        result = self.check_with(env)
        self.assertTrue(calls.exists(), "the fake git never saw the commit query")
        self.assertIn(f"{self.commit}^{{commit}}", calls.read_text())
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("ERROR:", result.stderr)
        self.assertIn("outgoing binding mismatch", result.stderr)
        self.assertIn("(exit 128)", result.stderr)
        self.assertNotIn("PRECHECK OK", result.stdout)
        control = self.check()
        self.assertEqual(control.returncode, 0, control.stderr)
        self.assertIn(f"PRECHECK OK phenix {self.commit}", control.stdout)


class OutgoingCommand(PrecheckFixture):
    """Checkpoint 3 (item 17): `outgoing` writes the strict one-commit OUTGOING.txt
    block from the live repository, checks it with the same parser and the `check`
    comparisons, and refuses a chain, a wrong parent, a missing remote, several push
    URLs and a dirty tree. Local bare repositories only."""

    def outgoing_tool(self, *options, commit=None):
        before = self.snapshot()
        result = subprocess.run(
            [sys.executable, "-I", "-B", str(TOOL), "outgoing", str(self.repo),
             "--commit", commit or self.commit, "--repository", "phenix", *options],
            env=self.env, cwd=self.root, text=True, capture_output=True)
        if "--fetch" not in options:
            self.assertEqual(self.snapshot(), before, "outgoing changed Git state")
        return result

    def refused(self, result, reason):
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("ERROR:", result.stderr)
        self.assertIn(reason, result.stderr)
        self.assertEqual(result.stdout, "")

    def test_generated_block_is_exact_parses_and_passes_check_dry_run(self):
        result = self.outgoing_tool()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, self.block())
        self.assertIn(f"OUTGOING OK phenix {self.commit}", result.stderr)
        self.outgoing.write_text(result.stdout)
        checked = self.check("--dry-run")
        self.assertEqual(checked.returncode, 0, checked.stderr)
        self.assertIn(f"PRECHECK OK phenix {self.commit}", checked.stdout)
        written = self.root / "written" / "OUTGOING.txt"
        written.parent.mkdir()
        result = self.outgoing_tool("--output", str(written))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(written.read_text(), self.block())
        again = self.outgoing_tool("--output", str(written))
        self.assertEqual(again.returncode, 2, again.stdout)
        self.assertEqual(written.read_text(), self.block())
        with self.subTest(fetch=True):
            fetched = self.outgoing_tool("--fetch")
            self.assertEqual(fetched.returncode, 0, fetched.stderr)
            self.assertEqual(fetched.stdout, self.block())

    def test_push_setting_note_goes_to_stderr_and_stdout_passes_check(self):
        """Review CP3 finding 1: with a harmless push setting (push.default=simple,
        common in ~/.gitconfig) the standard output is still exactly the block, the
        NOTE goes to standard error, and the redirected file passes check --dry-run."""
        self.git("config", "push.default", "simple", cwd=self.repo)
        self.git("config", "push.autoSetupRemote", "true", cwd=self.repo)
        result = self.outgoing_tool()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, self.block())
        self.assertNotIn("NOTE", result.stdout)
        self.assertIn("NOTE push setting present: push.default=simple", result.stderr)
        self.outgoing.write_text(result.stdout)
        checked = self.check("--dry-run")
        self.assertEqual(checked.returncode, 0, checked.stderr)
        self.assertIn(f"PRECHECK OK phenix {self.commit}", checked.stdout)
        self.assertIn("NOTE push setting present: push.default=simple", checked.stdout)

    def test_two_repositories_appended_to_one_file_pass_check(self):
        """Review CP45 finding 4: the documented route (WORKER.md, PUBLICATION_PACKET.md)
        appends each repository's standard output to one OUTGOING.txt; both blocks then
        pass `check --dry-run` for their repository. `--output` cannot add the second."""
        bare2, repo2 = self.root / "second.git", self.root / "second"
        self.git("init", "-q", "--bare", str(bare2), cwd=self.root)
        self.git("symbolic-ref", "HEAD", "refs/heads/master", cwd=bare2)
        repo2.mkdir()
        self.git("init", "-q", cwd=repo2)
        self.git("symbolic-ref", "HEAD", "refs/heads/master", cwd=repo2)
        self.configure_user(repo2)
        self.git("remote", "add", "origin", str(bare2), cwd=repo2)
        base2 = self.commit_file(repo2, "b.txt", "base\n", "base")
        self.git("push", "-q", "origin", "master", cwd=repo2)
        commit2 = self.commit_file(repo2, "b.txt", "change\n", "change")
        packet = self.root / "packet"
        packet.mkdir()
        outgoing = packet / "OUTGOING.txt"
        for repo, commit, name in ((self.repo, self.commit, "phenix"), (repo2, commit2, "second")):
            result = subprocess.run(
                [sys.executable, "-I", "-B", str(TOOL), "outgoing", str(repo), "--commit", commit,
                 "--repository", name, "--remote", "origin", "--branch", "master", "--fetch"],
                env=self.env, cwd=self.root, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with open(outgoing, "a", encoding="utf-8") as stream:  # the shell's >>
                stream.write(result.stdout)
        text = outgoing.read_text()
        self.assertTrue(text.startswith(self.block()), text)
        self.assertIn("repository: second\n", text)
        self.assertIn(f"base: {base2}\ncommit: {commit2}\n", text)
        for repo, commit, name in ((self.repo, self.commit, "phenix"), (repo2, commit2, "second")):
            checked = subprocess.run(
                [sys.executable, "-I", "-B", str(TOOL), "check", str(repo), str(outgoing),
                 "--repository", name, "--fetch", "--dry-run"],
                env=self.env, cwd=self.root, text=True, capture_output=True)
            self.assertEqual(checked.returncode, 0, checked.stderr)
            self.assertIn(f"PRECHECK OK {name} {commit}", checked.stdout)
        refused = subprocess.run(
            [sys.executable, "-I", "-B", str(TOOL), "outgoing", str(repo2), "--commit", commit2,
             "--repository", "second", "--output", str(outgoing)],
            env=self.env, cwd=self.root, text=True, capture_output=True)
        self.assertEqual(refused.returncode, 2, refused.stdout)
        self.assertEqual(outgoing.read_text(), text)

    def test_two_commit_chain_is_refused(self):
        second = self.commit_file(self.repo, "a.txt", "change 2\n", "second")
        self.refused(self.outgoing_tool(commit=second), "is a chain of 2 commits")

    def test_wrong_parent_is_refused(self):
        self.advance_remote()
        self.refused(self.outgoing_tool("--fetch"), "is not the remote tip")
        self.refused(self.outgoing_tool(commit=self.base), "has 0 parents")

    def test_missing_remote_is_refused(self):
        self.refused(self.outgoing_tool("--remote", "nowhere"), "remote not found: nowhere")
        self.refused(self.outgoing_tool("--branch", "other"),
                     "remote-tracking branch refs/remotes/origin/other not found")

    def test_several_push_urls_are_refused(self):
        self.git("remote", "set-url", "--add", "--push", "origin", str(self.bare), cwd=self.repo)
        self.git("remote", "set-url", "--add", "--push", "origin", str(self.root / "elsewhere.git"),
                 cwd=self.repo)
        self.refused(self.outgoing_tool(), "has 2 push URLs")

    def test_dirty_tree_is_refused(self):
        (self.repo / "a.txt").write_text("uncommitted\n")
        self.refused(self.outgoing_tool(), "working tree is not clean")
        (self.repo / "a.txt").write_text("change\n")
        (self.repo / "new.txt").write_text("untracked\n")
        self.refused(self.outgoing_tool(), "working tree is not clean")


class NoPushVariant(PrecheckFixture):
    """Revision 17: the no-push publication variant, for a project whose rules forbid the
    agent any `git push` (for example a project AGENTS.md: "Do not run `git push` in any form", read as including
    dry runs). The agent
    runs `check` without --dry-run and `vet`s the exact commands; the Developer runs them;
    `confirm` then reads the remote with `git ls-remote` (read-only: no push and no fetch)
    and confirms the branch is at the outgoing commit and the tag state is unchanged."""

    def tool(self, *words, env=None):
        """Run the tool and prove it changed no Git state, here or on the remote."""
        before = self.snapshot()
        result = subprocess.run([sys.executable, "-I", "-B", str(TOOL), *words],
                                env=env or self.env, cwd=self.root, text=True,
                                capture_output=True)
        self.assertEqual(self.snapshot(), before, f"{words[0]} changed Git state")
        return result

    def confirm(self, *options, env=None):
        return self.tool("confirm", str(self.repo), str(self.outgoing), *options, env=env)

    def push_as_developer(self, line):
        """The Developer runs the printed command himself (here: the test, outside the tool)."""
        words = line.split()
        self.assertEqual(words[:3], ["git", "-C", str(self.repo)])
        self.git(*words[3:], cwd=self.repo)

    def logging_git(self, trigger):
        """A git that logs every call whose words include `trigger` and runs the real git."""
        real_git = shutil.which("git", path=self.env.get("PATH", os.defpath))
        fake_bin = self.root / f"logbin-{'-'.join(trigger)}"
        fake_bin.mkdir()
        calls = fake_bin / "calls.txt"
        script = ("#!%s\nimport os, sys\nargs = sys.argv[1:]\n" % sys.executable
                  + "if any(word in args for word in %r):\n" % (tuple(trigger),)
                  + "    open(%r, 'a').write(' '.join(args) + '\\n')\n" % str(calls)
                  + "os.execv(%r, [%r] + args)\n" % (real_git, real_git))
        (fake_bin / "git").write_text(script)
        (fake_bin / "git").chmod(0o755)
        return dict(self.env, PATH=str(fake_bin) + os.pathsep + self.env.get("PATH", "")), calls

    def tag_state(self, result):
        line = next(x for x in result.stdout.splitlines() if x.startswith("REMOTE TAGS "))
        return line.split()[-1]

    def test_whole_variant_check_vet_developer_push_confirm(self):
        """The variant end to end: confirm --before passes and prints the tag state; check
        without --dry-run prints the push command, which vet accepts, as it does the dry-run
        form the Developer may run first; the agent's tools never call `git push` or `git
        fetch` (a logging git); after the Developer's push, confirm passes with the tag state
        from before, and confirm --before now refuses. Repair 1: after the push --tags is
        required."""
        env, calls = self.logging_git(("push", "fetch"))
        before = self.confirm("--before", env=env)
        self.assertEqual(before.returncode, 0, before.stderr)
        state = self.tag_state(before)
        self.assertTrue(state.startswith("0/"), state)
        self.assertIn("CONFIRMED before the push", before.stdout)
        not_yet = self.confirm("--tags", state, env=env)
        self.assertEqual(not_yet.returncode, 2, not_yet.stdout)
        self.assertIn("the push has not happened", not_yet.stderr)
        checked = self.tool("check", str(self.repo), str(self.outgoing), env=env)
        self.assertEqual(checked.returncode, 0, checked.stderr)
        self.assertNotIn("DRY RUN", checked.stdout)
        command = checked.stdout.splitlines()[-1]
        self.assertFalse(calls.exists(), "the agent's tools called git push or fetch: "
                         + (calls.read_text() if calls.exists() else ""))
        dry = command.replace("--no-follow-tags", "--no-follow-tags --dry-run --porcelain")
        for line in (command, dry):
            with self.subTest(vet=line):
                vetted = self.tool("vet", "--", *line.split())
                self.assertEqual(vetted.returncode, 0, vetted.stderr)
                self.assertIn("VET OK", vetted.stdout)
        self.push_as_developer(command)
        self.assertEqual(self.git("rev-parse", "refs/heads/master", cwd=self.bare), self.commit)
        no_tags = self.confirm(env=env)
        self.assertEqual(no_tags.returncode, 2, no_tags.stdout)
        self.assertIn("give --tags", no_tags.stderr)
        after = self.confirm("--tags", state, env=env)
        self.assertEqual(after.returncode, 0, after.stderr)
        self.assertIn(f"CONFIRMED phenix {self.commit} is origin refs/heads/master",
                      after.stdout)
        self.assertEqual(self.tag_state(after), state)
        late = self.confirm("--before", env=env)
        self.assertEqual(late.returncode, 2, late.stdout)
        self.assertIn("not the recorded base", late.stderr)
        self.assertFalse(calls.exists(), "confirm called git push or fetch")

    def test_confirm_refuses_a_changed_tag_set(self):
        """Negative controls: a tag pushed with the commit changes the state, and so does
        one tag removed and another added (the count unchanged, repair 1); the actual state
        passes."""
        self.git("tag", "v-old", self.base, cwd=self.repo)
        self.git("push", "-q", "origin", "refs/tags/v-old", cwd=self.repo)
        state = self.tag_state(self.confirm("--before"))
        self.assertTrue(state.startswith("1/"), state)
        self.git("tag", "v-test", cwd=self.repo)
        self.git("push", "-q", "origin", f"{self.commit}:refs/heads/master",
                 "refs/tags/v-test", cwd=self.repo)
        tagged = self.confirm("--tags", state)
        self.assertEqual(tagged.returncode, 2, tagged.stdout)
        self.assertIn("a tag was added, removed or moved", tagged.stderr)
        self.git("push", "-q", "origin", ":refs/tags/v-old", cwd=self.repo)
        swapped = self.confirm("--tags", state)
        self.assertEqual(swapped.returncode, 2, swapped.stdout)
        self.assertIn("a tag was added, removed or moved", swapped.stderr)
        current = subprocess.run([sys.executable, "-I", "-B", str(TOOL), "confirm", str(self.repo),
                                  str(self.outgoing), "--tags", "x"], env=self.env,
                                 text=True, capture_output=True)
        actual = current.stderr.split("remote tag state is ")[1].split(",")[0]
        self.assertTrue(actual.startswith("1/") and actual != state, actual)
        self.assertEqual(self.confirm("--tags", actual).returncode, 0)

    def test_confirm_refuses_a_remote_moved_elsewhere(self):
        """Negative control: a remote moved by someone else (not to the outgoing commit) is
        refused after and before the push."""
        other = self.advance_remote()
        for options, reason in ((("--tags", "0/unused"), "not the outgoing commit"),
                                (("--before",), "not the recorded base")):
            with self.subTest(options=options):
                result = self.confirm(*options)
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn(reason, result.stderr)
                self.assertIn(other, result.stderr)

    def test_confirm_after_the_push_refuses_a_deleted_branch(self):
        """Repair 1 of e9421261 (reading 1, finding 2): after the push, a remote branch that
        is absent (deleted) is refused, never confirmed; the same remote with the branch
        restored to the outgoing commit is the positive control."""
        state = self.tag_state(self.confirm("--before"))
        self.git("push", "-q", "origin", f"{self.commit}:refs/heads/master", cwd=self.repo)
        self.git("update-ref", "-d", "refs/heads/master", cwd=self.bare)
        result = self.confirm("--tags", state)
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("refs/heads/master is absent, not the outgoing commit", result.stderr)
        self.assertNotIn("CONFIRMED", result.stdout)
        self.git("update-ref", "refs/heads/master", self.commit, cwd=self.bare)
        self.assertEqual(self.confirm("--tags", state).returncode, 0)

    def env_logging_git(self, name):
        """A git that records GIT_SSH_COMMAND and GIT_TERMINAL_PROMPT of every ls-remote call
        and runs the real git."""
        real_git = shutil.which("git", path=self.env.get("PATH", os.defpath))
        fake_bin = self.root / f"envbin-{name}"
        fake_bin.mkdir()
        calls = fake_bin / "env.txt"
        script = ("#!%s\nimport os, sys\nargs = sys.argv[1:]\n" % sys.executable
                  + "if 'ls-remote' in args:\n"
                  + "    open(%r, 'a').write('%%s|%%s\\n' %% (os.environ.get('GIT_SSH_COMMAND'), "
                    "os.environ.get('GIT_TERMINAL_PROMPT')))\n" % str(calls)
                  + "os.execv(%r, [%r] + args)\n" % (real_git, real_git))
        (fake_bin / "git").write_text(script)
        (fake_bin / "git").chmod(0o755)
        return fake_bin, calls

    def test_confirm_runs_git_without_prompts_and_ssh_in_batch_mode(self):
        """Repair 1 of e9421261 (reading 1, finding 3): every `git ls-remote` that confirm
        runs has GIT_TERMINAL_PROMPT=0 and a GIT_SSH_COMMAND ending in `-o BatchMode=yes`,
        appended to the ssh command git would otherwise use: plain ssh; the repository's
        core.sshCommand; a GIT_SSH_COMMAND already set; a GIT_SSH program (quoted)."""
        cases = (("plain", {}, None, "ssh -o BatchMode=yes"),
                 ("core.sshCommand", {}, "ssh -i /keys/id_test", "ssh -i /keys/id_test -o BatchMode=yes"),
                 ("GIT_SSH_COMMAND", {"GIT_SSH_COMMAND": "ssh -p 2222"}, "ssh -i /ignored",
                  "ssh -p 2222 -o BatchMode=yes"),
                 ("GIT_SSH", {"GIT_SSH": "/opt/my ssh"}, None, "'/opt/my ssh' -o BatchMode=yes"))
        for label, extra, config, expected in cases:
            with self.subTest(source=label):
                if config:
                    self.git("config", "core.sshCommand", config, cwd=self.repo)
                fake_bin, calls = self.env_logging_git(label.replace(".", "-"))
                env = dict(self.env, PATH=str(fake_bin) + os.pathsep + self.env.get("PATH", ""),
                           **extra)
                result = self.confirm("--before", env=env)
                if config:
                    self.git("config", "--unset", "core.sshCommand", cwd=self.repo)
                self.assertEqual(result.returncode, 0, result.stderr)
                lines = calls.read_text().splitlines()
                self.assertEqual(len(lines), 2, lines)  # the branch and the tags
                self.assertEqual(set(lines), {f"{expected}|0"}, lines)

    def test_confirm_refuses_another_destination_and_a_missing_branch(self):
        """Negative controls: a push URL that differs from the record is refused before the
        remote is read; a branch the remote does not have is reported absent."""
        self.outgoing.write_text(self.block(remote_url=str(self.root / "elsewhere.git")))
        result = self.confirm("--before")
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("push destination differs from outgoing record", result.stderr)
        self.outgoing.write_text(self.block(refspec=f"{self.commit}:refs/heads/absent"))
        result = self.confirm("--before")
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("refs/heads/absent is absent", result.stderr)


class NoGitRefusal(unittest.TestCase):
    """Revision 16: with no git on PATH, `check` and `outgoing` (revision 17: and `confirm`) refuse (exit 2) with the
    cause before doing anything; `vet` still works. PATH holds only an empty directory."""

    def test_check_and_outgoing_refuse_without_git_vet_works(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(os.path.realpath(temp))
            empty = root / "empty-bin"
            empty.mkdir()
            env = dict(git_environment(root), PATH=str(empty))
            (root / "OUTGOING.txt").write_text("repository: x\n")
            for words in (["check", str(root), str(root / "OUTGOING.txt")],
                          ["outgoing", str(root), "--commit", HEX],
                          ["confirm", str(root), str(root / "OUTGOING.txt")]):  # revision 17
                with self.subTest(command=words[0]):
                    result = subprocess.run([sys.executable, "-I", "-B", str(TOOL), *words],
                                            env=env, cwd=root, text=True, capture_output=True)
                    self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                    self.assertEqual(result.stderr.strip(),
                                     f"ERROR: requires git on PATH: `{words[0]}` reads the "
                                     "repository with git (`vet` needs no git)")
                    self.assertEqual(result.stdout, "")
            result = subprocess.run([sys.executable, "-I", "-B", str(TOOL), "vet", "--", "git",
                                     "-C", "X", "-c", "push.followTags=false", "push",
                                     "--no-follow-tags", "origin", f"{HEX}:refs/heads/master"],
                                    env=env, cwd=root, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("VET OK", result.stdout)


class VetCommand(unittest.TestCase):
    REFSPEC = f"{HEX}:refs/heads/master"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.env = git_environment(Path(self.temp.name))

    def vet(self, *words):
        return subprocess.run([sys.executable, "-I", "-B", str(TOOL), "vet", "--", *words],
                              env=self.env, cwd=self.temp.name, text=True, capture_output=True)

    def canonical(self, repo="X", refspec=None):
        return ["git", "-C", repo, "-c", "push.followTags=false", "push",
                "--no-follow-tags", "origin", refspec or self.REFSPEC]

    def test_canonical_command_is_accepted(self):
        """SPEC_A 2.2: the exact command printed by `check` passes with VET OK."""
        for repo in ("X", self.temp.name, "/Users/dev/unix/PHENIX/modules/phenix"):
            with self.subTest(repo=repo):
                result = self.vet(*self.canonical(repo))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("VET OK", result.stdout)
        with self.subTest(refspec="nested branch"):
            result = self.vet(*self.canonical(refspec=f"{HEX}:refs/heads/release/1.0"))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("VET OK", result.stdout)
        with self.subTest(form="without the -c pair"):
            result = self.vet("git", "-C", "X", "push", "--no-follow-tags", "origin", self.REFSPEC)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("VET OK", result.stdout)

    def test_forbidden_words_are_refused(self):
        """SPEC_A 2.2: each listed word, and the -c push.followTags=true pair in
        any case, is refused as a forbidden push form."""
        for word in ("--tags", "--follow-tags", "--mirror", "--all", "--force", "-f",
                     "--force-with-lease", f"--force-with-lease=master:{HEX}", "--prune"):
            with self.subTest(word=word):
                command = self.canonical()
                command.insert(command.index("push") + 1, word)
                result = self.vet(*command)
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn(f"forbidden push form: {word}", result.stderr)
                self.assertNotIn("VET OK", result.stdout)
        for pair in (("-c", "push.followTags=true"), ("-c", "PUSH.FOLLOWTAGS=TRUE"),
                     ("-c", "push.followtags=True")):
            with self.subTest(pair=pair):
                result = self.vet("git", "-C", "X", *pair, "push", "--no-follow-tags",
                                  "origin", self.REFSPEC)
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn("forbidden push form", result.stderr)
                self.assertNotIn("VET OK", result.stdout)

    def test_refspec_forms_are_refused(self):
        """SPEC_A 2.2: a refspec with *, an empty side, any bare ref word that is
        not <40hex>:refs/heads/<branch>, or a second refspec is refused."""
        for refspec in (f"{HEX}:refs/heads/*", ":refs/heads/master", f"{HEX}:", "master",
                        "HEAD:refs/heads/master", "refs/heads/master:refs/heads/master",
                        f"{HEX}:refs/tags/v1", f"{HEX}:master", HEX, "refs/heads/master",
                        f"{HEX[:39]}:refs/heads/master"):
            with self.subTest(refspec=refspec):
                result = self.vet(*self.canonical(refspec=refspec))
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertTrue(result.stderr.startswith("ERROR:"), result.stderr)
                self.assertNotIn("VET OK", result.stdout)
        with self.subTest(refspec="two refspecs"):
            result = self.vet(*self.canonical(), f"{HEX}:refs/heads/other")
            self.assertEqual(result.returncode, 2, result.stdout)
            self.assertTrue(result.stderr.startswith("ERROR:"), result.stderr)
            self.assertNotIn("VET OK", result.stdout)

    def test_missing_no_follow_tags_is_refused(self):
        """SPEC_A 2.2: --no-follow-tags is required."""
        result = self.vet("git", "-C", "X", "-c", "push.followTags=false", "push",
                          "origin", self.REFSPEC)
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertIn("push command must disable tag following", result.stderr)
        self.assertNotIn("VET OK", result.stdout)

    # ---- SPEC_A 6: refinements after the internal reviewer's reading --------

    def test_bare_branch_beside_the_refspec_is_refused(self):
        """SPEC_A 6 (2.2), refusal forms per SPEC_A 8: words after push are
        read positionally (first non-option word = remote, every later one =
        refspec), so `origin feature <40hex>:refs/heads/master` is refused:
        exactly one refspec, forbidden push form naming the bare branch word,
        or the documented-shape message."""
        for label, words in (
                ("plain", ("git", "push", "--no-follow-tags", "origin", "feature", self.REFSPEC)),
                ("canonical prefix", (*self.canonical()[:-1], "feature", self.REFSPEC)),
                ("branch after the refspec", (*self.canonical(), "feature"))):
            with self.subTest(form=label):
                self.assert_refused(self.vet(*words), "feature", "exactly one refspec")

    def test_remote_must_be_a_name_not_a_url_or_path(self):
        """SPEC_A 6 (2.2), refusal forms per SPEC_A 8: the remote operand must
        be a name; a word containing `:` or `/` (an scp-style or https URL, an
        absolute or relative path) is refused."""
        for remote in ("git@github.com:cctbx/cctbx_project.git",
                       "https://github.com/cctbx/cctbx_project.git",
                       "ssh://git@example.invalid/phenix.git",
                       "/tmp/remote.git", "../remote.git"):
            with self.subTest(remote=remote):
                command = self.canonical()
                command[command.index("origin")] = remote
                self.assert_refused(self.vet(*command), remote)

    def test_forbidden_c_settings_are_refused_in_any_case(self):
        """SPEC_A 6 (2.2), refusal forms per SPEC_A 8: a -c setting whose key
        matches url.*.insteadof, url.*.pushinsteadof, remote.*.pushurl,
        remote.*.mirror, remote.*.push or remote.*.tagopt (any case) is refused,
        even beside the canonical -c push.followTags=false pair (under the
        SPEC_A 8 whitelist no other -c is admitted at all)."""
        for setting in (
                "url.https://example.invalid/.insteadOf=ssh://nowhere.invalid/",
                "url.https://example.invalid/.pushInsteadOf=ssh://nowhere.invalid/",
                "URL.https://example.invalid/.INSTEADOF=ssh://nowhere.invalid/",
                "Url.x.PushInsteadOf=y",
                "remote.origin.pushurl=/tmp/elsewhere.git",
                "REMOTE.ORIGIN.PUSHURL=/tmp/elsewhere.git",
                "remote.upstream.pushurl=/tmp/elsewhere.git",
                "remote.origin.mirror=true", "remote.origin.mirror=false",
                "Remote.Origin.Mirror=no",
                "remote.origin.push=refs/heads/master:refs/heads/master",
                "remote.origin.PUSH=refs/heads/master",
                "remote.origin.tagOpt=--tags", "remote.origin.tagopt=--no-tags",
                "REMOTE.ORIGIN.TAGOPT=--tags"):
            with self.subTest(setting=setting):
                self.assert_refused(
                    self.vet("git", "-C", "X", "-c", "push.followTags=false", "-c", setting,
                             "push", "--no-follow-tags", "origin", self.REFSPEC),
                    setting, "-c")

    def test_followtags_truthy_values_are_refused_and_false_is_accepted(self):
        """SPEC_A 6 (2.2), refusal forms per SPEC_A 8: -c push.followTags with
        true, yes, on or 1 in any case (key or value) is refused;
        -c push.followTags=false is still accepted."""
        pairs = [("push.followTags", value)
                 for value in ("true", "yes", "on", "1", "TRUE", "Yes", "ON", "True")]
        pairs += [("PUSH.FOLLOWTAGS", "yes"), ("push.followtags", "On"), ("PUSH.FOLLOWTAGS", "1")]
        for key, value in pairs:
            with self.subTest(key=key, value=value):
                self.assert_refused(
                    self.vet("git", "-C", "X", "-c", f"{key}={value}", "push",
                             "--no-follow-tags", "origin", self.REFSPEC),
                    f"{key}={value}", "-c")
        result = self.vet("git", "-C", "X", "-c", "push.followTags=false", "push",
                          "--no-follow-tags", "origin", self.REFSPEC)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("VET OK", result.stdout)
        self.assertNotIn("ERROR", result.stderr)

    def test_delete_and_transport_words_are_refused(self):
        """SPEC_A 6 (2.2), refusal forms per SPEC_A 8: --delete, -d,
        --receive-pack, --exec, and words starting with --receive-pack= or
        --exec= are refused wherever they stand after push."""
        words = ("--delete", "-d", "--receive-pack", "--exec",
                 "--receive-pack=/usr/bin/git-receive-pack", "--exec=/usr/bin/git-receive-pack",
                 "--receive-pack=", "--exec=x")
        for word in words:
            with self.subTest(word=word, position="after push"):
                command = self.canonical()
                command.insert(command.index("push") + 1, word)
                self.assert_refused(self.vet(*command), word)
        for word in ("--delete", "-d", "--exec=x"):
            with self.subTest(word=word, position="after the refspec"):
                self.assert_refused(self.vet(*self.canonical(), word), word)
        with self.subTest(word="--receive-pack", position="with its own argument"):
            command = self.canonical()
            at = command.index("push") + 1
            command[at:at] = ["--receive-pack", "/usr/bin/git-receive-pack"]
            self.assert_refused(self.vet(*command), "--receive-pack", "/usr/bin/git-receive-pack")

    def test_canonical_command_still_accepted_after_the_refinements(self):
        """SPEC_A 6 (2.2) positive control: `git -C <dir> -c push.followTags=false
        push --no-follow-tags origin <40hex>:refs/heads/master` is VET OK; the
        -C directory (which contains `/`) is an option argument, not the
        remote operand."""
        for repo in ("X", self.temp.name, "/Users/dev/unix/PHENIX/modules/phenix"):
            with self.subTest(repo=repo):
                result = self.vet(*self.canonical(repo))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("VET OK", result.stdout)
                self.assertNotIn("ERROR", result.stderr)

    def test_valueless_c_key_means_true(self):
        """SPEC_A 7 (2.2), narrowed by SPEC_A 8: a `-c <key>` word without "="
        means true to Git, so a valueless -c push.followTags (any case) and a
        valueless forbidden key (remote.*.mirror, remote.*.pushurl,
        url.*.insteadOf, ...) are refused; -c push.followTags=false stays
        accepted. SPEC_A 7 still accepted an unrelated valueless key such as
        -c color.ui; under the SPEC_A 8 whitelist it is refused as well."""
        for key in ("push.followTags", "PUSH.FOLLOWTAGS", "push.followtags",
                    "remote.origin.mirror", "remote.origin.pushurl", "url.x.insteadOf",
                    "url.x.pushInsteadOf", "remote.origin.push", "remote.origin.tagOpt",
                    "REMOTE.ORIGIN.MIRROR", "Url.X.InsteadOf"):
            with self.subTest(key=key):
                self.assert_refused(self.vet("git", "-C", "X", "-c", key, "push",
                                             "--no-follow-tags", "origin", self.REFSPEC),
                                    key, "-c")
        with self.subTest(key="push.followTags", position="beside the canonical pair"):
            self.assert_refused(self.vet("git", "-C", "X", "-c", "push.followTags=false", "-c",
                                         "push.followTags", "push", "--no-follow-tags", "origin",
                                         self.REFSPEC), "push.followTags", "-c")
        result = self.vet(*self.canonical())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("VET OK", result.stdout)
        self.assertNotIn("ERROR", result.stderr)
        # SPEC_A 7 accepted an unrelated valueless key; the SPEC_A 8 whitelist
        # admits only -c push.followTags=false, so these are refused now.
        for label, words in (
                ("valueless color.ui alone",
                 ["git", "-C", "X", "-c", "color.ui", "push", "--no-follow-tags", "origin",
                  self.REFSPEC]),
                ("valueless color.ui beside the canonical pair",
                 ["git", "-C", "X", "-c", "push.followTags=false", "-c", "color.ui", "push",
                  "--no-follow-tags", "origin", self.REFSPEC])):
            with self.subTest(form=label):
                self.assert_refused(self.vet(*words), "color.ui", "-c")

    # ---- SPEC_A 8: the vet is a whitelist of the documented shape -----------

    def assert_refused(self, result, *words):
        """Exit 2, no VET OK, and an ERROR: line in one of the SPEC_A 8 refusal
        forms: "forbidden push form" naming one of `words`, or "push command
        must have the documented shape". A word containing a space is a message
        phrase accepted on its own (SPEC_A 6: "exactly one refspec"); a word
        without any letter or digit (`--`, `.`, `..`) must appear as a whole
        token, so it cannot match inside --no-follow-tags or a sentence."""
        self.assertEqual(result.returncode, 2, result.stdout)
        self.assertNotIn("VET OK", result.stdout)
        errors = [line for line in result.stderr.splitlines() if line.startswith("ERROR:")]
        self.assertTrue(errors, result.stderr)
        phrases = ["documented shape"] + [word for word in words if " " in word]
        names = [word for word in words if " " not in word]

        def named(line, word):
            if not any(ch.isalnum() for ch in word):
                tokens = line.replace("'", " ").replace('"', " ").replace("`", " ").split()
                return word in tokens
            return word in line

        def refusal(line):
            return (any(phrase in line for phrase in phrases)
                    or ("forbidden push form" in line and any(named(line, w) for w in names)))

        self.assertTrue(any(refusal(line) for line in errors), (words, result.stderr))

    def test_whitelist_refuses_every_word_outside_the_documented_shape(self):
        """SPEC_A 8 (2.2): any word outside the documented shape, wherever it
        stands, is refused (exit 2, ERROR: on stderr, no VET OK); the ERROR
        line names the refused word or says "documented shape"."""
        base = self.canonical()
        at_push = base.index("push")

        def after_push(*words):
            return [*base[:at_push + 1], *words, *base[at_push + 1:]]

        def before_push(*words):
            return [*base[:at_push], *words, *base[at_push:]]

        second = f"{HEX}:refs/heads/other"
        cases = (
            ("--ta", after_push("--ta"), ("--ta",)),
            ("--tags", after_push("--tags"), ("--tags",)),
            ("-vf", after_push("-vf"), ("-vf",)),
            ("-f", after_push("-f"), ("-f",)),
            ("-v", after_push("-v"), ("-v",)),
            ("--config-env before push", before_push("--config-env=remote.origin.pushurl=x"),
             ("--config-env=remote.origin.pushurl=x",)),
            ("--receive-pack=/x", after_push("--receive-pack=/x"), ("--receive-pack=/x",)),
            ("-c color.ui after push", after_push("-c", "color.ui"), ("color.ui", "-c")),
            ("-o ci.skip", after_push("-o", "ci.skip"), ("ci.skip", "-o")),
            ("--repo x", after_push("--repo", "x"), ("--repo",)),
            ("second -C /y", [*base[:3], "-C", "/y", *base[3:]], ("/y", "-C")),
            ("-c push.followTags=true before push", before_push("-c", "push.followTags=true"),
             ("push.followTags=true", "-c")),
            ("-c color.ui=auto before push", before_push("-c", "color.ui=auto"),
             ("color.ui=auto", "-c")),
            ("-c push.followTags=FALSE (value not exactly false)",
             [*base[:4], "push.followTags=FALSE", *base[5:]], ("push.followTags=FALSE", "-c")),
            ("-- separator before the operands", [*base[:-2], "--", *base[-2:]], ("--",)),
            ("three operands", [*base, second], (second, "exactly one refspec")),
        )
        for label, command, words in cases:
            with self.subTest(form=label):
                self.assert_refused(self.vet(*command), *words)

    def test_whitelist_accepts_exactly_the_documented_shape(self):
        """SPEC_A 8 (2.2): VET OK for the documented shape only: optional git,
        optional -C <path>, optional -c push.followTags=false (key case-
        insensitive, value exactly false), push, option words from
        {--no-follow-tags, --dry-run, --porcelain} in any order with
        --no-follow-tags present, then a remote name and one
        <40hex>:refs/heads/<branch> refspec."""
        base = self.canonical(self.temp.name)
        at_push = base.index("push")
        forms = {
            "canonical": base,
            "without git": base[1:],
            "without -C": [base[0], *base[3:]],
            "without -c": [*base[:3], *base[at_push:]],
            "without git, -C and -c": base[at_push:],
            "key case PUSH.FOLLOWTAGS=false": [*base[:4], "PUSH.FOLLOWTAGS=false", *base[5:]],
        }
        options = ("--no-follow-tags", "--dry-run", "--porcelain")
        for order in itertools.permutations(options):
            forms["options " + " ".join(order)] = [*base[:at_push + 1], *order, *base[at_push + 2:]]
        forms["options --no-follow-tags --dry-run"] = [
            *base[:at_push + 2], "--dry-run", *base[at_push + 2:]]
        forms["options --porcelain --no-follow-tags"] = [
            *base[:at_push + 1], "--porcelain", *base[at_push + 1:]]
        for label, command in forms.items():
            with self.subTest(form=label):
                result = self.vet(*command)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("VET OK", result.stdout)
                self.assertNotIn("ERROR", result.stderr)

    def test_remote_grammar_and_exact_word_order(self):
        """SPEC_A 9 (2.2): the remote operand must match
        [A-Za-z0-9_][A-Za-z0-9._-]* and be neither . nor ..; the documented
        word order is enforced exactly (-C before -c, the -C argument not
        starting with -, options before the operands, each option at most
        once). Legal remote names pass the vet with the canonical shape (check
        binds them to the configured remote later)."""
        base = self.canonical()
        at_remote = base.index("origin")
        for remote in (".", "..", "origin=x", "origin;id", "-origin"):
            with self.subTest(remote=remote):
                command = list(base)
                command[at_remote] = remote
                self.assert_refused(self.vet(*command), remote)
        with self.subTest(remote="ori gin (two words: three operands)"):
            command = [*base[:at_remote], "ori", "gin", *base[at_remote + 1:]]
            self.assert_refused(self.vet(*command), "ori", "gin", "exactly one refspec")
        for label, command, words in (
                ("-c before -C",
                 ["git", "-c", "push.followTags=false", "-C", "/x", "push", "--no-follow-tags",
                  "origin", self.REFSPEC], ("-c", "-C", "push.followTags=false", "/x")),
                ("option after an operand",
                 ["git", "push", "--no-follow-tags", "origin", self.REFSPEC, "--dry-run"],
                 ("--dry-run",)),
                ("repeated option",
                 ["git", "push", "--no-follow-tags", "--no-follow-tags", "origin", self.REFSPEC],
                 ("--no-follow-tags",)),
                ("-C argument starting with -",
                 ["git", "-C", "-x", "push", "--no-follow-tags", "origin", self.REFSPEC],
                 ("-x", "-C"))):
            with self.subTest(form=label):
                self.assert_refused(self.vet(*command), *words)
        for remote in ("origin", "upstream-2", "my.remote", "origin."):
            with self.subTest(remote=remote):
                command = list(base)
                command[at_remote] = remote
                result = self.vet(*command)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("VET OK", result.stdout)
                self.assertNotIn("ERROR", result.stderr)
        with self.subTest(form="options in another order, all before the operands"):
            result = self.vet("git", "-C", "/x", "-c", "push.followTags=false", "push",
                              "--dry-run", "--no-follow-tags", "--porcelain", "origin",
                              self.REFSPEC)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("VET OK", result.stdout)
            self.assertNotIn("ERROR", result.stderr)


if __name__ == "__main__":
    unittest.main()
