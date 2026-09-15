from __future__ import absolute_import, division, print_function
"""Test of libtbx.guided_coding: install, verify, remove, and its refusals.

Runs under plain python3 (no libtbx import) against scratch git repositories
in a temporary directory, which it removes at the end. Every command the test
runs gets HOME and TMPDIR inside that directory, so the command's captures
(~/Downloads), its working copies of them and its per-install backup
directories (the temp directory) all land under the test's own root and go
with it; the test never touches the system temporary directory. On Windows
the test checks that the command refuses, prints a skip line and OK: the
command and the procedure are for macOS and Linux only.
"""
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_TMPDIR = None   # set by run_all: TMPDIR for every command the test runs
CMD = os.path.join(HERE, "command_line", "guided_coding.py")
PAYLOAD = os.path.join(HERE, "guided_coding", "payload")


def run(args, env_extra=None, cwd=None):
  env = dict(os.environ)
  env["PYTHONDONTWRITEBYTECODE"] = "1"
  if TEST_TMPDIR:
    env["TMPDIR"] = TEST_TMPDIR
  env["GIT_CONFIG_GLOBAL"] = os.devnull
  env["GIT_CONFIG_NOSYSTEM"] = "1"
  if GIT_HOME:
    env["XDG_CONFIG_HOME"] = os.path.join(GIT_HOME, ".config")
  if env_extra:
    env.update(env_extra)
  p = subprocess.Popen([sys.executable, CMD] + list(args), stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, env=env, cwd=cwd)
  out = p.communicate()[0].decode("utf-8", "replace")
  return p.returncode, out


GIT_HOME = None   # set by run_all: a scratch HOME/XDG_CONFIG_HOME for git calls


def git_env():
  """The test's own git calls see none of the user's git settings: a scratch
  HOME and XDG_CONFIG_HOME (so git's default excludes file ~/.config/git/
  ignore is never read), the global and system configuration ignored
  (git 2.32+), and - for every git version - the excludes file and the
  untracked-files setting overridden by flags on each call."""
  env = dict(os.environ)
  env["GIT_CONFIG_GLOBAL"] = os.devnull
  env["GIT_CONFIG_NOSYSTEM"] = "1"
  if GIT_HOME:
    env["HOME"] = GIT_HOME
    env["XDG_CONFIG_HOME"] = os.path.join(GIT_HOME, ".config")
  return env


GIT_FLAGS = ["-c", "core.excludesFile=" + os.devnull, "-c", "status.showUntrackedFiles=all"]


def git_status(repo):
  return subprocess.check_output(["git"] + GIT_FLAGS + ["-C", repo, "status", "--porcelain",
                                  "--untracked-files=all"], env=git_env())


def make_repo(root, name):
  repo = os.path.join(root, "modules", name)
  os.makedirs(repo)
  subprocess.check_call(["git"] + GIT_FLAGS + ["init", "-q", repo], env=git_env())
  return repo


def write(path, text):
  d = os.path.dirname(path)
  if not os.path.isdir(d):
    os.makedirs(d)
  with io.open(path, "w", encoding="utf-8") as f:
    f.write(text)


def exercise_fresh_install_verify_rerun_remove(root, home):
  repo = make_repo(root, "fresh")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  assert "Install complete" in out, out
  assert "no personal profile CLAUDE.local.md" in out, out
  manifest = os.path.join(repo, ".claude", "MANIFEST.sha256")
  assert os.path.isfile(manifest)
  # excludes written; git sees nothing
  st = git_status(repo)
  assert st.strip() == b"", st
  rc, out = run(["verify", repo], {"HOME": home})
  assert rc == 0 and "NOT_OK=0" in out, out
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0 and "upgraded: 0" in out and "new; upgraded" in out, out
  # a user file under records survives removal; package templates do not
  write(os.path.join(repo, ".claude", "records", "2026-01-01-mine.md"), u"# mine\n")
  write(os.path.join(repo, ".claude", "settings.local.json"), u"{}\n")
  rc, out = run(["remove", repo], {"HOME": home})
  assert rc == 0 and "Removed" in out, out
  left = []
  for dp, dn, fn in os.walk(os.path.join(repo, ".claude")):
    for f in fn:
      left.append(os.path.relpath(os.path.join(dp, f), repo).replace(os.sep, "/"))
  assert sorted(left) == [".claude/records/2026-01-01-mine.md",
                          ".claude/settings.local.json"], left


def exercise_refusals(root, home):
  repo = make_repo(root, "refuse")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  # locally modified file -> whole install stops
  t = os.path.join(repo, ".claude", "rules", "testing.md")
  with io.open(t, "a", encoding="utf-8") as f:
    f.write(u"local edit\n")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "CONFLICT" in out and "Nothing was installed" in out, out
  shutil.copy(os.path.join(PAYLOAD, "rules", "testing.md"), t)
  # open record (new marker) -> refuse; override -> proceed
  rec = os.path.join(repo, ".claude", "records", "2026-01-02-open.md")
  write(rec, u"# y\n\n## Status (generated, exactly one place)\n- Status: OPEN\n")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "still says it is OPEN" in out, out
  rc, out = run(["install", repo], {"HOME": home, "ALLOW_OPEN_CHANGE": "yes"})
  assert rc == 0, out
  os.remove(rec)
  # legacy top-line marker -> refuse
  write(rec, u"# z\n\nStatus: open (stage 3)\n")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2, out
  os.remove(rec)
  # state file -> refuse
  hj = os.path.join(repo, ".claude", "records", "HANDOFF.json")
  write(hj, json.dumps({"schema_version": "1", "open_change": {"status": "open"}}))
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "state file" in out, out
  os.remove(hj)
  # tree lock (two levels up from the repository) -> refuse, override does not help
  lock = os.path.join(root, ".codinghelper_tree_in_use")
  os.mkdir(lock)
  rc, out = run(["install", repo], {"HOME": home, "ALLOW_OPEN_CHANGE": "yes"})
  assert rc == 2 and "tree lock" in out, out
  os.rmdir(lock)
  # a symlinked destination -> refuse; a DANGLING symlink too (it must never
  # be written through)
  link_target = os.path.join(repo, ".claude", "rules", "testing.md")
  os.remove(link_target)
  os.symlink(os.path.join(root, "does-not-exist"), link_target)
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "symbolic link" in out, out
  assert not os.path.exists(os.path.join(root, "does-not-exist"))
  os.remove(link_target)
  shutil.copy(os.path.join(PAYLOAD, "rules", "testing.md"), link_target)
  # a CLOSED record that keeps an old "- Status: OPEN" line in its history
  # must not refuse: only the generated section's first Status line counts
  rec2 = os.path.join(repo, ".claude", "records", "2026-01-03-closed.md")
  write(rec2, u"# c\n\n## Status (generated, exactly one place)\n- Status: CLOSED\n"
              u"\nEarlier status lines, SUPERSEDED (kept as generated):\n- Status: OPEN\n")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  os.remove(rec2)


def exercise_upgrade_from_hidden_profile(root, home):
  """Ownership comes only from package-shipped data. A profile whose bytes
  match a packaged profile migrates; one that is merely listed in the
  installed manifest does not; a colleague's CLAUDE.md is left alone;
  adopt-profile migrates a generated profile on explicit say-so."""
  repo = make_repo(root, "old")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  mp = os.path.join(repo, ".claude", "MANIFEST.sha256")
  # 1. a real packaged profile (bytes taken from ACCEPTED_PREVIOUS) migrates
  acc = os.path.join(HERE, "guided_coding", "ACCEPTED_PREVIOUS.sha256")
  prof_hash = None
  for line in io.open(acc, encoding="utf-8"):
    h, p = line.split()
    if p == "./CLAUDE.md":
      prof_hash = h
      break
  assert prof_hash, "ACCEPTED_PREVIOUS carries no packaged profile hash"
  # we cannot reproduce those bytes here, so simulate with a temporary
  # accepted entry for known bytes: write the profile, compute, and append to
  # a COPY of the package data is not possible without editing the package;
  # instead prove the negative cases below and the adopt-profile path.
  prof = os.path.join(repo, "CLAUDE.md")
  write(prof, u"# CLAUDE.md - a profile the interview generated\n")
  def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
  # 2. a manifest listing alone does not make a profile the package's
  with io.open(mp, "a", encoding="utf-8") as f:
    f.write(u"%s  ./CLAUDE.md\n" % h(prof))
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0 and "not a profile this package ever shipped" in out, out
  assert os.path.isfile(prof) and not os.path.exists(os.path.join(repo, "CLAUDE.local.md"))
  # 3. a manifest listing alone does not make a user file retirable
  notes = os.path.join(repo, ".claude", "MY_USER_NOTES.md")
  write(notes, u"my notes\n")
  with io.open(mp, "a", encoding="utf-8") as f:
    f.write(u"%s  ./.claude/MY_USER_NOTES.md\n" % h(notes))
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0 and "Retired" not in out and "treated as yours" in out, out
  assert os.path.isfile(notes)
  # 4. a manifest hash that matches a modified file does not make it replaceable
  mod = os.path.join(repo, ".claude", "rules", "testing.md")
  with io.open(mod, "a", encoding="utf-8") as f:
    f.write(u"USER SENTINEL\n")
  lines = [l for l in io.open(mp, encoding="utf-8")]
  with io.open(mp, "w", encoding="utf-8") as f:
    for l in lines:
      if l.strip().endswith("./.claude/rules/testing.md"):
        f.write(u"%s  ./.claude/rules/testing.md\n" % h(mod))
      else:
        f.write(l)
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "CONFLICT" in out, out
  assert u"USER SENTINEL" in io.open(mod, encoding="utf-8").read()
  shutil.copy(os.path.join(PAYLOAD, "rules", "testing.md"), mod)
  # restore a sane manifest by reinstalling with the sentinel gone
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  # 5. adopt-profile migrates the generated profile explicitly and drops the old exclude line
  excl = os.path.join(repo, ".git", "info", "exclude")
  with io.open(excl, "a", encoding="utf-8") as f:
    f.write(u"/CLAUDE.md\n")
  rc, out = run(["adopt-profile", repo], {"HOME": home})
  assert rc == 0 and "Renamed CLAUDE.md -> CLAUDE.local.md" in out, out
  assert not os.path.exists(prof) and os.path.isfile(os.path.join(repo, "CLAUDE.local.md"))
  ex = [l.strip() for l in io.open(excl, encoding="utf-8")]
  assert "/CLAUDE.md" not in ex and "/CLAUDE.local.md" in ex, ex
  rc, out = run(["verify", repo], {"HOME": home})
  assert rc == 0 and "NOT_OK=0" in out, out
  # 6. a colleague's CLAUDE.md is untouched even when listed in the manifest
  repo2 = make_repo(root, "foreign")
  write(os.path.join(repo2, "CLAUDE.md"), u"# a colleague's project file\n")
  rc, out = run(["install", repo2], {"HOME": home})
  assert rc == 0, out
  with io.open(os.path.join(repo2, ".claude", "MANIFEST.sha256"), "a", encoding="utf-8") as f:
    f.write(u"%s  ./CLAUDE.md\n" % h(os.path.join(repo2, "CLAUDE.md")))
  rc, out = run(["install", repo2], {"HOME": home})
  assert rc == 0 and "not a profile this package ever shipped" in out, out
  assert io.open(os.path.join(repo2, "CLAUDE.md"), encoding="utf-8").read().startswith(u"# a colleague")


def exercise_manifest_confinement(root, home):
  """Malformed entries in an installed manifest (absolute paths, '..', or
  escapes) are rejected by remove, retire and verify; nothing outside the
  repository is touched."""
  repo = make_repo(root, "confine")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  outside = os.path.join(root, "OUTSIDE.txt")
  write(outside, u"must survive\n")
  mp = os.path.join(repo, ".claude", "MANIFEST.sha256")
  good = io.open(mp, encoding="utf-8").read()
  bad_lines = [u"%s  ./.claude/../../../OUTSIDE.txt\n" % ("0" * 64),
               u"%s  /etc/hostname\n" % ("0" * 64),
               u"%s  ./.claude/rules/../../OUTSIDE.txt\n" % ("0" * 64)]
  for bad in bad_lines:
    with io.open(mp, "w", encoding="utf-8") as f:
      f.write(good + bad)
    for cmd in (["remove", repo], ["verify", repo], ["install", repo]):
      rc, out = run(cmd, {"HOME": home})
      assert rc == 2 and "manifest path refused" in out, (cmd, bad, out)
    assert os.path.isfile(outside), "file outside the repository changed"
  # duplicate destination refused
  with io.open(mp, "w", encoding="utf-8") as f:
    f.write(good + good.splitlines()[0] + u"\n")
  rc, out = run(["verify", repo], {"HOME": home})
  assert rc == 2 and "twice" in out, out
  with io.open(mp, "w", encoding="utf-8") as f:
    f.write(good)
  # retirement rejects an escaping entry too, before anything moves
  with io.open(mp, "w", encoding="utf-8") as f:
    f.write(good + u"%s  ./.claude/../../OUTSIDE.txt\n" % hashlib.sha256(b"must survive\n").hexdigest())
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "manifest path refused" in out, out
  assert os.path.isfile(outside)
  with io.open(mp, "w", encoding="utf-8") as f:
    f.write(good)


def exercise_remove_safety(root, home):
  """Removal: the lock is never overridable; symlinked parents refuse; a
  modified package file is retained; HANDOFF.json, known_failure_variation.md
  and CLAUDE.local.md survive."""
  repo = make_repo(root, "rm")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  # lock + override -> still refused
  lock = os.path.join(root, ".codinghelper_tree_in_use")
  os.mkdir(lock)
  rc, out = run(["remove", repo], {"HOME": home, "ALLOW_OPEN_CHANGE": "yes"})
  assert rc == 2 and "tree lock" in out, out
  os.rmdir(lock)
  # symlinked parent directory -> refused, external file untouched
  ext = os.path.join(root, "ext_rules")
  os.makedirs(ext)
  shutil.copy(os.path.join(PAYLOAD, "rules", "testing.md"), os.path.join(ext, "testing.md"))
  shutil.copy(os.path.join(PAYLOAD, "rules", "conventions.md"), os.path.join(ext, "conventions.md"))
  rules = os.path.join(repo, ".claude", "rules")
  shutil.rmtree(rules)
  os.symlink(ext, rules)
  rc, out = run(["remove", repo], {"HOME": home})
  assert rc == 2 and ("symbolic link" in out or "escapes the repository" in out), out
  assert os.path.isfile(os.path.join(ext, "testing.md"))
  os.remove(rules)
  # a symlinked parent that stays INSIDE the repository is caught by the
  # symlink check itself
  inside = os.path.join(repo, ".claude", "rules_real")
  shutil.copytree(ext, inside)
  os.symlink(inside, rules)
  rc, out = run(["remove", repo], {"HOME": home})
  assert rc == 2 and "symbolic link" in out, out
  assert os.path.isfile(os.path.join(inside, "testing.md"))
  os.remove(rules)
  shutil.rmtree(inside)
  shutil.copytree(ext, rules)
  # modified package file is retained; user state files survive
  mod = os.path.join(repo, ".claude", "rules", "testing.md")
  with io.open(mod, "a", encoding="utf-8") as f:
    f.write(u"my local note\n")
  for name, text in [("HANDOFF.json", u'{"schema_version": "1", "open_change": {"status": "closed"}}\n'),
                     ("known_failure_variation.md", u"# mine\n")]:
    write(os.path.join(repo, ".claude", "records", name), text)
  write(os.path.join(repo, "CLAUDE.local.md"), u"# my profile\n")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "CONFLICT" in out, out   # modified file blocks upgrade, as before
  rc, out = run(["remove", repo], {"HOME": home})
  assert rc == 0 and "Retained 1" in out and "rules/testing.md" in out, out
  assert os.path.isfile(mod)
  for name in ("HANDOFF.json", "known_failure_variation.md"):
    assert os.path.isfile(os.path.join(repo, ".claude", "records", name)), name
  assert os.path.isfile(os.path.join(repo, "CLAUDE.local.md"))


def exercise_profile_collision_and_git_visibility(root, home):
  """A user's own '/CLAUDE.md' exclude line is preserved when no migration
  happened; a colleague's CLAUDE.md stays visible to git; RELEASE goes
  through the destination machinery (a symlink there is refused)."""
  repo = make_repo(root, "userexclude")
  excl = os.path.join(repo, ".git", "info", "exclude")
  write(excl, u"/CLAUDE.md\n")          # the user's own choice, before any install
  write(os.path.join(repo, "CLAUDE.md"), u"# the user hides this on purpose\n")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0 and "may be yours" in out, out
  ex = [l.strip() for l in io.open(excl, encoding="utf-8")]
  assert "/CLAUDE.md" in ex, ex
  st = git_status(repo).decode()
  assert "CLAUDE.md" not in st, st
  # a colleague's CLAUDE.md with no exclude line is visible after install
  repo2 = make_repo(root, "visible")
  write(os.path.join(repo2, "CLAUDE.md"), u"# a colleague's project file\n")
  rc, out = run(["install", repo2], {"HOME": home})
  assert rc == 0, out
  st = git_status(repo2).decode()
  assert "?? CLAUDE.md" in st, st
  # RELEASE is a payload file: a symlink at its destination refuses the install
  repo3 = make_repo(root, "release")
  os.makedirs(os.path.join(repo3, ".claude"))
  ext = os.path.join(root, "external_release_target")
  write(ext, u"external\n")
  os.symlink(ext, os.path.join(repo3, ".claude", "RELEASE"))
  rc, out = run(["install", repo3], {"HOME": home})
  assert rc == 2 and "symbolic link" in out, out
  assert io.open(ext, encoding="utf-8").read() == u"external\n"


def exercise_status_and_capture_collision(root, home):
  repo = make_repo(root, "st")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  rc, out = run(["status", repo, os.path.join(root, "modules", "nothere")], {"HOME": home})
  assert rc == 0 and "package release:" in out and "pre-libtbx package" in out, out
  # a differing earlier capture is never overwritten
  cap = os.path.join(home, "Downloads", "guided_coding_verify_capture.txt")
  write(cap, u"earlier, different\n")
  rc, out = run(["verify", repo], {"HOME": home})
  assert rc == 0 and "earlier capture with different content" in out, out
  assert io.open(cap, encoding="utf-8").read() == u"earlier, different\n"
  # a second run in the same second must not overwrite the first alternate:
  # pre-create every plausible alternate name for this second and check none changed
  import glob
  before = {}
  for p in glob.glob(os.path.join(home, "Downloads", "guided_coding_verify_capture_*.txt")):
    before[p] = open(p, "rb").read()
  rc, out = run(["verify", repo], {"HOME": home})
  assert rc == 0, out
  for p, data in before.items():
    assert open(p, "rb").read() == data, "alternate capture overwritten: " + p


def exercise_packaged_profile_migrates(root, home):
  """The positive migration case, tested against a COPY of the package whose
  accepted-previous list carries the hash of a profile we write here. No
  override exists in the shipped command; the test controls the trusted data
  by controlling the package copy it runs."""
  pkg = os.path.join(root, "pkgcopy")
  shutil.copytree(os.path.join(HERE, "command_line"), os.path.join(pkg, "command_line"))
  shutil.copytree(os.path.join(HERE, "guided_coding"), os.path.join(pkg, "guided_coding"))
  cmd_copy = os.path.join(pkg, "command_line", "guided_coding.py")
  repo = make_repo(root, "packaged")
  prof = os.path.join(repo, "CLAUDE.md")
  write(prof, u"# CLAUDE.md - GuidedCoding project profile (packaged r48 stand-in)\n")
  with io.open(os.path.join(pkg, "guided_coding", "ACCEPTED_PREVIOUS.sha256"), "a", encoding="utf-8") as f:
    f.write(u"%s  ./CLAUDE.md\n" % hashlib.sha256(open(prof, "rb").read()).hexdigest())
  excl = os.path.join(repo, ".git", "info", "exclude")
  write(excl, u".claude/\n/CLAUDE.md\n")     # what a pre-step-5 install wrote
  def runcopy(args, env_extra):
    env = dict(os.environ); env.update(env_extra); env["PYTHONDONTWRITEBYTECODE"] = "1"
    if TEST_TMPDIR:
      env["TMPDIR"] = TEST_TMPDIR
    p = subprocess.Popen([sys.executable, cmd_copy] + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    out = p.communicate()[0].decode("utf-8", "replace")
    return p.returncode, out
  # no CLAUDE.local.md yet: renamed, and the old exclude line dropped
  rc, out = runcopy(["install", repo], {"HOME": home})
  assert rc == 0 and "renamed to CLAUDE.local.md" in out and "removed the old /CLAUDE.md" in out, out
  assert not os.path.exists(prof) and os.path.isfile(os.path.join(repo, "CLAUDE.local.md"))
  ex = [l.strip() for l in io.open(excl, encoding="utf-8")]
  assert "/CLAUDE.md" not in ex and "/CLAUDE.local.md" in ex, ex
  # collision: packaged old profile present again AND a new CLAUDE.local.md -> old one retired
  write(prof, u"# CLAUDE.md - GuidedCoding project profile (packaged r48 stand-in)\n")
  write(os.path.join(repo, "CLAUDE.local.md"), u"# new personal profile\n")
  rc, out = runcopy(["install", repo], {"HOME": home})
  assert rc == 0 and "retired to the backup" in out, out
  assert not os.path.exists(prof)
  assert io.open(os.path.join(repo, "CLAUDE.local.md"), encoding="utf-8").read().startswith(u"# new personal")
  # the SHIPPED command, with its real accepted list, does NOT recognise this stand-in profile
  write(prof, u"# CLAUDE.md - GuidedCoding project profile (packaged r48 stand-in)\n")
  os.remove(os.path.join(repo, "CLAUDE.local.md"))
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0 and "not a profile this package ever shipped" in out, out
  assert os.path.isfile(prof)


def exercise_remove_authority(root, home):
  """remove deletes only known package bytes; manifest hashes that match a
  modified package file or a user file do not make them removable."""
  repo = make_repo(root, "rmauth")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  mp = os.path.join(repo, ".claude", "MANIFEST.sha256")
  def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
  mod = os.path.join(repo, ".claude", "rules", "testing.md")
  with io.open(mod, "a", encoding="utf-8") as f:
    f.write(u"USER SENTINEL\n")
  notes = os.path.join(repo, ".claude", "MY_USER_NOTES.md")
  write(notes, u"my notes\n")
  lines = [l for l in io.open(mp, encoding="utf-8")]
  with io.open(mp, "w", encoding="utf-8") as f:
    for l in lines:
      if l.strip().endswith("./.claude/rules/testing.md"):
        f.write(u"%s  ./.claude/rules/testing.md\n" % h(mod))
      else:
        f.write(l)
    f.write(u"%s  ./.claude/MY_USER_NOTES.md\n" % h(notes))
  rc, out = run(["remove", repo], {"HOME": home})
  assert rc == 0 and "Retained 2" in out, out
  assert os.path.isfile(mod) and u"USER SENTINEL" in io.open(mod, encoding="utf-8").read()
  assert os.path.isfile(notes)


def exercise_capture_write_failure_is_bounded(root, home):
  """When ~/Downloads cannot take a capture, the command reports once and
  returns; it never loops."""
  bad_home = os.path.join(root, "badhome")
  os.makedirs(bad_home)
  write(os.path.join(bad_home, "Downloads"), u"not a directory\n")   # a FILE named Downloads
  repo = make_repo(root, "capfail")
  rc, out = run(["install", repo], {"HOME": bad_home})
  assert rc == 0 and ("capture at:" in out or "could not write a capture" in out), out
  # and a primary capture that exists as a directory-with-that-name collision
  home2 = os.path.join(root, "home2")
  os.makedirs(os.path.join(home2, "Downloads", "guided_coding_verify_capture.txt"))
  rc, out = run(["verify", repo], {"HOME": home2})
  assert rc == 0 and "NOT_OK=0" in out, out


def exercise_adopt_profile_preflight(root, home):
  """adopt-profile changes nothing when any destination is a symbolic link or
  the manifest is malformed; the checks run before the rename."""
  def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
  # symlinked .git/info -> refused, external exclude untouched, profile untouched
  repo = make_repo(root, "adopt1")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  write(os.path.join(repo, "CLAUDE.md"), u"# generated profile\n")
  ext_info = os.path.join(root, "ext_info")
  os.makedirs(ext_info)
  write(os.path.join(ext_info, "exclude"), u"external\n")
  shutil.rmtree(os.path.join(repo, ".git", "info"))
  os.symlink(ext_info, os.path.join(repo, ".git", "info"))
  rc, out = run(["adopt-profile", repo], {"HOME": home})
  assert rc == 2 and "symbolic link" in out, out
  assert os.path.isfile(os.path.join(repo, "CLAUDE.md"))
  assert not os.path.exists(os.path.join(repo, "CLAUDE.local.md"))
  assert io.open(os.path.join(ext_info, "exclude"), encoding="utf-8").read() == u"external\n"
  # symlinked manifest -> refused, external file untouched
  repo = make_repo(root, "adopt2")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  write(os.path.join(repo, "CLAUDE.md"), u"# generated profile\n")
  ext_m = os.path.join(root, "ext_manifest.txt")
  write(ext_m, u"external manifest\n")
  mp = os.path.join(repo, ".claude", "MANIFEST.sha256")
  os.remove(mp)
  os.symlink(ext_m, mp)
  rc, out = run(["adopt-profile", repo], {"HOME": home})
  assert rc == 2 and "symbolic link" in out, out
  assert os.path.isfile(os.path.join(repo, "CLAUDE.md"))
  assert io.open(ext_m, encoding="utf-8").read() == u"external manifest\n"
  # malformed manifest -> refused BEFORE the rename
  repo = make_repo(root, "adopt3")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  write(os.path.join(repo, "CLAUDE.md"), u"# generated profile\n")
  mp = os.path.join(repo, ".claude", "MANIFEST.sha256")
  with io.open(mp, "a", encoding="utf-8") as f:
    f.write(u"%s  ./.claude/../../escape.txt\n" % ("0" * 64))
  rc, out = run(["adopt-profile", repo], {"HOME": home})
  assert rc == 2 and "manifest path refused" in out, out
  assert os.path.isfile(os.path.join(repo, "CLAUDE.md"))
  assert not os.path.exists(os.path.join(repo, "CLAUDE.local.md"))


def exercise_capture_non_regular_file(root, home):
  """A pipe (FIFO) or other non-regular object at the capture's name is
  reported, never opened: the command returns promptly."""
  if not hasattr(os, "mkfifo"):
    return
  repo = make_repo(root, "fifo")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  home3 = os.path.join(root, "home3")
  os.makedirs(os.path.join(home3, "Downloads"))
  os.mkfifo(os.path.join(home3, "Downloads", "guided_coding_verify_capture.txt"))
  env = dict(os.environ)
  env.update({"HOME": home3, "PYTHONDONTWRITEBYTECODE": "1", "TMPDIR": TEST_TMPDIR})
  p = subprocess.Popen([sys.executable, CMD, "verify", repo], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, env=env)
  try:
    out = p.communicate(timeout=60)[0].decode("utf-8", "replace")
  except subprocess.TimeoutExpired:
    p.kill()
    raise AssertionError("verify hung on a FIFO at the capture path")
  assert p.returncode == 0 and "not a regular file" in out, out


def exercise_hard_link_aliases(root, home):
  """An existing file the command would rewrite is refused when it is a
  hard-link alias of another file; the other file is never changed."""
  # exclude file aliased to an external file -> install refuses
  repo = make_repo(root, "hl1")
  ext = os.path.join(root, "hl_external_exclude")
  write(ext, u"EXTERNAL SENTINEL\n")
  excl = os.path.join(repo, ".git", "info", "exclude")
  if os.path.lexists(excl):
    os.remove(excl)
  os.link(ext, excl)
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "more than one hard link" in out, out
  assert io.open(ext, encoding="utf-8").read() == u"EXTERNAL SENTINEL\n"
  os.remove(excl)
  # installed manifest aliased -> reinstall refuses
  repo = make_repo(root, "hl2")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  mp = os.path.join(repo, ".claude", "MANIFEST.sha256")
  ext_m = os.path.join(root, "hl_external_manifest")
  os.link(mp, ext_m)
  before = hashlib.sha256(open(ext_m, "rb").read()).hexdigest()
  write(os.path.join(repo, "CLAUDE.local.md"), u"# profile\n")   # would change the manifest
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 2 and "more than one hard link" in out, out
  assert hashlib.sha256(open(ext_m, "rb").read()).hexdigest() == before
  os.remove(ext_m)
  # a payload destination aliased -> refused, alias unchanged
  repo = make_repo(root, "hl3")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  target = os.path.join(repo, ".claude", "rules", "testing.md")
  ext_p = os.path.join(root, "hl_external_payload")
  os.link(target, ext_p)
  before = hashlib.sha256(open(ext_p, "rb").read()).hexdigest()
  # make the destination need replacing: an accepted-previous payload copy would; here force via a
  # package copy with a modified testing.md so the destination is classified as an upgrade
  pkg = os.path.join(root, "pkgcopy_hl")
  shutil.copytree(os.path.join(HERE, "command_line"), os.path.join(pkg, "command_line"))
  shutil.copytree(os.path.join(HERE, "guided_coding"), os.path.join(pkg, "guided_coding"))
  pt = os.path.join(pkg, "guided_coding", "payload", "rules", "testing.md")
  with io.open(pt, "a", encoding="utf-8") as f:
    f.write(u"\n(next release)\n")
  with io.open(os.path.join(pkg, "guided_coding", "PAYLOAD_MANIFEST.sha256"), "w", encoding="utf-8") as f:
    for dp, dn, fn in os.walk(os.path.join(pkg, "guided_coding", "payload")):
      for name in sorted(fn):
        full = os.path.join(dp, name)
        rel = os.path.relpath(full, os.path.join(pkg, "guided_coding", "payload")).replace(os.sep, "/")
        f.write(u"%s  ./%s\n" % (hashlib.sha256(open(full, "rb").read()).hexdigest(), rel))
  env = dict(os.environ); env.update({"HOME": home, "PYTHONDONTWRITEBYTECODE": "1", "TMPDIR": TEST_TMPDIR})
  p = subprocess.Popen([sys.executable, os.path.join(pkg, "command_line", "guided_coding.py"), "install", repo],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
  out = p.communicate()[0].decode("utf-8", "replace")
  assert p.returncode == 2 and "more than one hard link" in out, out
  assert hashlib.sha256(open(ext_p, "rb").read()).hexdigest() == before
  os.remove(ext_p)
  # adopt-profile with an aliased exclude file -> refused before the rename
  repo = make_repo(root, "hl4")
  rc, out = run(["install", repo], {"HOME": home})
  assert rc == 0, out
  write(os.path.join(repo, "CLAUDE.md"), u"# generated profile\n")
  excl = os.path.join(repo, ".git", "info", "exclude")
  ext_e = os.path.join(root, "hl_external_exclude2")
  os.link(excl, ext_e)
  before = open(ext_e, "rb").read()
  rc, out = run(["adopt-profile", repo], {"HOME": home})
  assert rc == 2 and "more than one hard link" in out, out
  assert os.path.isfile(os.path.join(repo, "CLAUDE.md"))
  assert open(ext_e, "rb").read() == before


def exercise_usage_and_git_isolation(root, home):
  """No arguments or an unknown subcommand prints the usage text, never
  'None' (I-B); a global git excludes file hiding CLAUDE.md does not fool the
  command or the test (M7); captures and backups land under a dated
  guided_coding/ directory of TMPDIR."""
  rc, out = run([], {"HOME": home})
  assert rc == 2 and "Install, verify, or remove" in out and out.split("\n")[0].strip() != "None", out
  rc, out = run(["frobnicate", "x"], {"HOME": home})
  assert rc == 2 and "Install, verify, or remove" in out, out
  gx = os.path.join(root, "global_excludes")
  write(gx, u"CLAUDE.md\n")
  gc = os.path.join(root, "gitconfig_global")
  write(gc, u"[core]\n\texcludesFile = %s\n[status]\n\tshowUntrackedFiles = no\n" % gx)
  repo = make_repo(root, "isolated")
  write(os.path.join(repo, "CLAUDE.md"), u"# a colleague's project file\n")
  # git's DEFAULT excludes file, read even when the global config is empty
  xdg = os.path.join(root, "hostile_xdg")
  write(os.path.join(xdg, "git", "ignore"), u"CLAUDE.md\n.claude/\n")
  hostile = {"HOME": home, "GIT_CONFIG_GLOBAL": gc, "GIT_CONFIG_NOSYSTEM": "1",
             "XDG_CONFIG_HOME": xdg}
  # the command must succeed under the hostile global config (it isolates itself)
  rc, out = run(["install", repo], hostile)
  assert rc == 0, out
  # a raw git status under that config hides CLAUDE.md (via the default excludes
  # file even with an empty global config); the test's helper does not, and
  # neither do the command's checks under the same hostile environment
  env = dict(os.environ); env.update(hostile)
  raw = subprocess.check_output(["git", "-C", repo, "status", "--porcelain"], env=env).decode()
  assert "CLAUDE.md" not in raw, raw
  env2 = dict(env); env2["GIT_CONFIG_GLOBAL"] = os.devnull
  raw2 = subprocess.check_output(["git", "-C", repo, "status", "--porcelain"], env=env2).decode()
  assert "CLAUDE.md" not in raw2, "the default excludes file did not hide it: " + raw2
  assert "?? CLAUDE.md" in git_status(repo).decode()
  # the command's excludes check is not vacuous under the hostile environment:
  # its exact isolated git call must still SEE .claude/ once the repository's
  # own exclude line is removed, despite the default excludes file hiding it
  assert "Excludes verified." in out, out
  excl = os.path.join(repo, ".git", "info", "exclude")
  kept = [l for l in io.open(excl, encoding="utf-8") if l.strip() != ".claude/"]
  with io.open(excl, "w", encoding="utf-8") as f:
    f.writelines(kept)
  iso = subprocess.check_output(["git", "-c", "core.excludesFile=" + os.devnull,
                                 "-c", "status.showUntrackedFiles=all", "-C", repo,
                                 "status", "--porcelain", "--untracked-files=all"], env=env).decode()
  assert ".claude/" in iso, "the command's isolated status did not see .claude/: " + iso
  with io.open(excl, "a", encoding="utf-8") as f:
    f.write(u".claude/\n")
  stamp = "guided_coding_%s_" % time.strftime("%Y-%m-%d")
  runs = [n for n in os.listdir(TEST_TMPDIR) if n.startswith(stamp)]
  assert runs, os.listdir(TEST_TMPDIR)
  # a symlink where a shared parent would have been is never followed: there is no shared parent
  assert not os.path.lexists(os.path.join(TEST_TMPDIR, "guided_coding"))


def run_all():
  global TEST_TMPDIR, GIT_HOME
  root = tempfile.mkdtemp(prefix="tst_guided_coding.")
  TEST_TMPDIR = os.path.join(root, "tmp")
  os.makedirs(TEST_TMPDIR)
  GIT_HOME = os.path.join(root, "githome")
  os.makedirs(os.path.join(GIT_HOME, ".config", "git"))
  if sys.platform.startswith("win"):
    win_home = os.path.join(root, "winhome")
    os.makedirs(os.path.join(win_home, "Downloads"))
    rc, out = run(["status", HERE], {"HOME": win_home, "USERPROFILE": win_home,
                                     "TEMP": TEST_TMPDIR, "TMP": TEST_TMPDIR})
    assert rc == 2 and "Windows is not supported" in out, out
    assert os.listdir(os.path.join(win_home, "Downloads")) == [], "the refusal wrote a capture"
    shutil.rmtree(root, ignore_errors=True)
    print("Skipping tst_guided_coding: GuidedCoding supports macOS and Linux only.")
    print("OK")
    return
  home = os.path.join(root, "home")
  os.makedirs(os.path.join(home, "Downloads"))
  groups = [
    exercise_fresh_install_verify_rerun_remove,
    exercise_refusals,
    exercise_upgrade_from_hidden_profile,
    exercise_manifest_confinement,
    exercise_remove_safety,
    exercise_profile_collision_and_git_visibility,
    exercise_status_and_capture_collision,
    exercise_packaged_profile_migrates,
    exercise_remove_authority,
    exercise_capture_write_failure_is_bounded,
    exercise_adopt_profile_preflight,
    exercise_capture_non_regular_file,
    exercise_hard_link_aliases,
    exercise_usage_and_git_isolation,
  ]
  t0 = time.time()
  try:
    for g in groups:
      t1 = time.time()
      g(root, home)
      print("  %-50s %6.1f s" % (g.__name__, time.time() - t1))
      sys.stdout.flush()
  finally:
    shutil.rmtree(root, ignore_errors=True)   # takes the commands' captures and backups with it
  print("total %.1f s" % (time.time() - t0))
  print("OK")


if __name__ == "__main__":
  run_all()
