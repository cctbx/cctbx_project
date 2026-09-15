"""Install, verify, or remove the GuidedCoding procedure in a repository.
(libtbx.guided_coding, package release 1.1)

Usage:
  libtbx.guided_coding install <repository-path>   place or upgrade the
      procedure files under <repository>/.claude/ and record them in
      <repository>/.claude/MANIFEST.sha256. The personal profile
      CLAUDE.local.md is NEVER written by this command: INSTALL.md's
      interview writes it, and this command only records its hash.
  libtbx.guided_coding verify <repository-path>    check every installed
      file against the manifest.
  libtbx.guided_coding remove <repository-path>    package-owned removal:
      delete only the files the manifest lists as the package's own, then
      the manifest. Records, settings.local.json, HANDOFF.json,
      known_failure_variation.md and CLAUDE.local.md stay.
  libtbx.guided_coding status <repository-path>... report the installed
      release in each repository against this package's release.
  libtbx.guided_coding adopt-profile <repository-path>  on your explicit
      say-so, rename a CLAUDE.md that IS your GuidedCoding profile (one the
      interview generated or you revised, so its bytes are not a packaged
      profile the command can recognise) to CLAUDE.local.md, drop the old
      /CLAUDE.md exclude line, and record it in the manifest.

Rules, each earned in use:
  - The whole install stops before placing a single file if any destination
    is a symbolic link, if an installed file matches neither this release
    nor any previously delivered one (locally modified), or if a change is
    open in the repository (tree lock, a record whose generated Status
    section's FIRST "- Status:" line says OPEN, a pre-r40 top-line
    "Status: open", or HANDOFF.json saying open). The lock is looked for
    two directories above the repository ($PHENIX/.codinghelper_tree_in_use
    for $PHENIX/modules/<repo>); a repository elsewhere is not yet supported.
  - Installation is not transactional: if it fails halfway through, the
    files already placed stay placed. Every replaced file was backed up
    first (the backup path is in the capture), so nothing is lost, but the
    repository can be left between two releases until install is re-run.
  - Platforms: macOS and Linux. On Windows the command refuses to run before
    creating anything; the procedure it installs is not supported there.
  - The command's own git call ignores every user-level git setting that
    could hide a file: the global and system configuration (GIT_CONFIG_GLOBAL,
    GIT_CONFIG_NOSYSTEM, git 2.32+) AND, for every git version, the excludes
    file and the untracked-files setting, overridden on the command line
    (-c core.excludesFile=/dev/null, --untracked-files=all). git's default
    excludes file ~/.config/git/ignore is covered by the -c override.
  - Each run's capture and backups live in one private directory created
    atomically under the temp directory, named guided_coding_<date>_<random>. ALLOW_OPEN_CHANGE=yes in the environment overrides a stale
    record, never the lock.
  - Every replaced file is backed up; every placed file is re-hashed
    afterwards; the capture goes to ~/Downloads without overwriting a
    differing file.
  - Any existing regular file the command would rewrite - a payload
    destination, the manifest, the exclude file, a profile - is refused when
    it has more than one hard link, since rewriting it would change every
    other name for the same bytes. Deleting (remove) is not affected: unlink
    removes one name only.
  - Ownership is decided ONLY by package-shipped data (this payload's
    manifest and ACCEPTED_PREVIOUS.sha256, the hashes of everything ever
    delivered, profiles included). The installed manifest is an ordinary
    editable file: it says which paths to look at, never whether bytes are
    the package's.
  - The old hidden profile CLAUDE.md is renamed to CLAUDE.local.md only when
    its bytes are a profile the package shipped; any other CLAUDE.md is left
    alone and visible to git, and adopt-profile exists for the developer to
    migrate a generated or revised profile explicitly.

This module has no libtbx imports on purpose: it runs the same way under
libtbx.python and under plain python3, and its own test exercises it under
python3 against scratch repositories.
"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME libtbx.guided_coding
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

PACKAGE_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "guided_coding"))
PAYLOAD_DIR = os.path.join(PACKAGE_DIR, "payload")
PAYLOAD_MANIFEST = os.path.join(PACKAGE_DIR, "PAYLOAD_MANIFEST.sha256")
ACCEPTED_PREVIOUS = os.path.join(PACKAGE_DIR, "ACCEPTED_PREVIOUS.sha256")
RELEASE_FILE = os.path.join(PAYLOAD_DIR, "RELEASE")

LOCK_NAME = ".codinghelper_tree_in_use"   # kept: old and new sessions must
                                          # exclude each other by one name
PROFILE_LOCAL = "CLAUDE.local.md"
PROFILE_OLD = "CLAUDE.md"
EXCLUDE_LINES = [".claude/", "/" + PROFILE_LOCAL]
RETIRED_EXCLUDE_LINES = ["/" + PROFILE_OLD]   # written by pre-step-5 installs; removed on install
NEVER_REMOVE = set(["records", "settings.local.json", "worktrees"])


def sha256_file(path):
  h = hashlib.sha256()
  with open(path, "rb") as f:
    for block in iter(lambda: f.read(1 << 16), b""):
      h.update(block)
  return h.hexdigest()


def read_manifest(path):
  """Return a list of (hash, relative_path) from a manifest file."""
  out = []
  if not os.path.isfile(path):
    return out
  for line in io.open(path, encoding="utf-8"):
    line = line.rstrip("\n")
    if not line.strip():
      continue
    h, p = line.split(None, 1)
    out.append((h, p.strip()))
  return out


def safe_rel(repo, p):
  """Return the validated relative path for a manifest entry, or raise Stop.

  A manifest is a file on disk and could be edited; nothing read from it may
  point outside the repository or through '..'. Accepted: a relative path
  with no absolute prefix, no '..' component, no empty component, no
  backslash, whose real parent directory lies inside the repository."""
  rel = rel_of(p)
  if not rel or rel.startswith("/") or "\\" in rel or rel.startswith("~"):
    raise Stop("manifest path refused (absolute or malformed): %r" % p)
  parts = rel.split("/")
  if any(c in ("", ".", "..") for c in parts):
    raise Stop("manifest path refused (empty, '.' or '..' component): %r" % p)
  target = os.path.join(repo, rel)
  root = os.path.realpath(repo)
  parent_real = os.path.realpath(os.path.dirname(target))
  if parent_real != root and not parent_real.startswith(root + os.sep):
    raise Stop("manifest path refused (escapes the repository): %r" % p)
  return rel


def validated_manifest(repo, path):
  """(hash, rel) pairs from a manifest, every path confined, no duplicates."""
  out, seen = [], set()
  for h, p in read_manifest(path):
    rel = safe_rel(repo, p)
    if rel in seen:
      raise Stop("manifest lists %s twice" % rel)
    seen.add(rel)
    if not re.match(r"^[0-9a-f]{64}$", h):
      raise Stop("manifest hash malformed for %s" % rel)
    out.append((h, rel))
  return out


def release_string():
  if os.path.isfile(RELEASE_FILE):
    return io.open(RELEASE_FILE, encoding="utf-8").read().strip()
  return "(no RELEASE file)"


def rel_of(p):
  """Strip one leading './' from a manifest path; never strip dots from names."""
  return p[2:] if p.startswith("./") else p


_WORK_DIR = None


def work_dir():
  """One private directory per process for this command's capture and
  backups, created atomically by mkdtemp directly under the temp directory
  (TMPDIR), with the date in its name so it is easy to find and to remove.
  No shared or predictable parent is ever created or followed."""
  global _WORK_DIR
  if _WORK_DIR is None:
    _WORK_DIR = tempfile.mkdtemp(prefix="guided_coding_%s_" % time.strftime("%Y-%m-%d"))
  return _WORK_DIR


class Log(object):
  def __init__(self, name):
    self.lines = []
    fd, self.path = tempfile.mkstemp(prefix=name + ".", suffix=".txt",
                                     dir=work_dir())
    os.close(fd)
    self.name = name

  def __call__(self, msg):
    print(msg)
    self.lines.append(msg)
    with io.open(self.path, "a", encoding="utf-8") as f:
      f.write(msg + u"\n")

  def finish(self):
    """Copy the capture to ~/Downloads without ever overwriting a differing
    file, and without ever raising: the capture is evidence, and a failure
    to file it is reported, not thrown."""
    import errno
    try:
      downloads = os.path.expanduser("~/Downloads")
      if not os.path.isdir(downloads):
        print("capture at: %s" % self.path)
        return
      dest = os.path.join(downloads, self.name + "_capture.txt")
      import stat as _stat
      if os.path.lexists(dest):
        st = os.lstat(dest)
        if not _stat.S_ISREG(st.st_mode):
          print("NOTE: %s exists and is not a regular file (link, directory, "
                "pipe or device); capture left at %s" % (dest, self.path))
          return
      data = open(self.path, "rb").read()
      if os.path.lexists(dest):
        if open(dest, "rb").read() == data:
          print("Capture already in ~/Downloads and identical: %s" % dest)
          return
        stamp = time.strftime("%H%M%S")
        for n in range(1000):
          suffix = stamp if n == 0 else "%s_%d" % (stamp, n)
          cand = os.path.join(downloads, "%s_capture_%s.txt" % (self.name, suffix))
          try:
            fd = os.open(cand, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
          except OSError as e:
            if e.errno == errno.EEXIST:
              continue
            print("NOTE: could not write a capture into ~/Downloads (%s); "
                  "the capture is at %s" % (e, self.path))
            return
          with os.fdopen(fd, "wb") as out:
            out.write(data)
          print("NOTE: an earlier capture with different content exists; "
                "this run saved as %s" % cand)
          return
        print("NOTE: too many captures with this name in ~/Downloads; the "
              "capture is at %s" % self.path)
        return
      try:
        fd = os.open(dest, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
      except OSError as e:
        print("NOTE: could not write a capture into ~/Downloads (%s); "
              "the capture is at %s" % (e, self.path))
        return
      with os.fdopen(fd, "wb") as out:
        out.write(data)
      print("Capture copied to: %s" % dest)
    except Exception as e:   # never let the capture step raise
      print("NOTE: capture could not be filed (%s); it is at %s" % (e, self.path))


class Stop(Exception):
  pass


def nolink_path(repo, target):
  """Every component of target below repo must be a non-symlink."""
  rel = os.path.relpath(target, repo)
  acc = repo
  for c in rel.split(os.sep):
    acc = os.path.join(acc, c)
    if os.path.islink(acc):
      return False
  return True


def multiply_linked(path):
  """True if path exists as a regular file with more than one hard link -
  rewriting it would also change every other name for the same bytes."""
  try:
    st = os.lstat(path)
  except OSError:
    return False
  import stat as _stat
  return _stat.S_ISREG(st.st_mode) and st.st_nlink > 1


GIT_ISOLATION_FLAGS = ["-c", "core.excludesFile=" + os.devnull,
                       "-c", "status.showUntrackedFiles=all"]


def git_env():
  """Environment for the command's own git calls: the user's global and the
  system git configuration are ignored (git 2.32+ honours these variables;
  older gits ignore them, which is why the FLAGS below carry the same
  isolation for every git version)."""
  env = dict(os.environ)
  env["GIT_CONFIG_GLOBAL"] = os.devnull
  env["GIT_CONFIG_NOSYSTEM"] = "1"
  return env


def is_git_repo(repo):
  return os.path.isdir(os.path.join(repo, ".git"))


def open_change_signals(repo):
  """Return a list of reasons a change appears to be open in repo."""
  reasons = []
  lock = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(repo))),
                      LOCK_NAME)
  if os.path.isdir(lock):
    reasons.append(("lock", lock))
  records = os.path.join(repo, ".claude", "records")
  if os.path.isdir(records):
    for name in sorted(os.listdir(records)):
      if not (name.startswith("20") and name.endswith(".md")):
        continue
      path = os.path.join(records, name)
      try:
        text = io.open(path, encoding="utf-8", errors="replace").read()
      except (IOError, OSError):
        continue
      head = "\n".join(text.split("\n")[:20])
      is_open = bool(re.search(r"^Status: *open", head, re.I | re.M))
      m = re.search(r"^## Status \(generated.*?$", text, re.I | re.M)
      if m:
        section = text[m.end():]
        first = re.search(r"^- Status: *(\w+)", section, re.I | re.M)
        if first and first.group(1).upper() == "OPEN":
          is_open = True
      if is_open:
        reasons.append(("record", path))
    hj = os.path.join(records, "HANDOFF.json")
    if os.path.isfile(hj):
      try:
        d = json.load(io.open(hj, encoding="utf-8"))
        if str(d.get("open_change", {}).get("status", "")).lower() == "open":
          reasons.append(("state file", hj))
      except (ValueError, IOError, OSError):
        pass
  return reasons


def accepted_hashes():
  """Hashes of every previously delivered payload file, by relative path."""
  acc = {}
  for h, p in read_manifest(ACCEPTED_PREVIOUS):
    p = rel_of(p)
    if p.startswith(".claude/"):
      p = p[len(".claude/"):]
    elif p != PROFILE_OLD:
      continue
    acc.setdefault(p, set()).add(h)
  return acc


_BACKUP = {}


def backup_dir_for(repo):
  return _BACKUP[repo]


def claude_dir_for(repo):
  return os.path.join(repo, ".claude")


def cmd_install(repo, log):
  repo = os.path.abspath(repo)
  log("=== GuidedCoding install: %s ===" % release_string())
  log("Date: %s" % time.strftime("%Y-%m-%d %H:%M:%S %Z"))
  log("Target repository: %s" % repo)
  if not is_git_repo(repo):
    raise Stop("%s is not a git repository (no .git directory)." % repo)
  payload = [(h, rel_of(q)) for h, q in read_manifest(PAYLOAD_MANIFEST)]
  for h, q in payload:
    safe_rel(claude_dir_for(repo), q)   # payload paths confined too
  if not payload:
    raise Stop("this package has no payload manifest at %s" % PAYLOAD_MANIFEST)
  # Package self-check.
  for h, p in payload:
    src = os.path.join(PAYLOAD_DIR, p)
    if not os.path.isfile(src) or sha256_file(src) != h:
      raise Stop("package payload does not match its own manifest: %s" % p)
  log("Package payload matches its manifest (%d files)." % len(payload))

  # Refuse under an open change.
  signals = open_change_signals(repo)
  lock_signals = [s for s in signals if s[0] == "lock"]
  other = [s for s in signals if s[0] != "lock"]
  if lock_signals:
    raise Stop("the tree lock %s exists - a change is running or was "
               "interrupted. Close it (or /wrap it) first. Nothing was "
               "installed." % lock_signals[0][1])
  if other and os.environ.get("ALLOW_OPEN_CHANGE", "") != "yes":
    log("A change in this repository still says it is OPEN:")
    for kind, path in other:
      log("  %s: %s" % (kind, path))
    raise Stop("installing while a change is open is not allowed. Close it "
               "with /wrap first. If the record is stale and you are certain "
               "nothing is running, re-run with ALLOW_OPEN_CHANGE=yes. "
               "Nothing was installed.")

  claude_dir = os.path.join(repo, ".claude")
  installed_manifest_path = os.path.join(claude_dir, "MANIFEST.sha256")
  installed = dict((rel, h) for h, rel in
                   validated_manifest(repo, installed_manifest_path))

  # Symlink refusal on every destination.
  links = []
  for h, p in payload:
    t = os.path.join(claude_dir, p)
    if os.path.lexists(t) and (os.path.islink(t) or not nolink_path(repo, t)):
      links.append(t)
  for x in [claude_dir, installed_manifest_path,
            os.path.join(repo, ".git"), os.path.join(repo, ".git", "info"),
            os.path.join(repo, ".git", "info", "exclude")]:
    if os.path.islink(x):
      links.append(x)
  if links:
    for x in links:
      log("CONFLICT: %s (or a parent) is a symbolic link." % x)
    raise Stop("%d symlinked destination(s). Nothing was installed." % len(links))
  aliased = []
  for h, p in payload:
    t = os.path.join(claude_dir, p)
    if os.path.isfile(t) and multiply_linked(t):
      aliased.append(t)
  for x in [installed_manifest_path, os.path.join(repo, ".git", "info", "exclude"),
            os.path.join(repo, PROFILE_OLD), os.path.join(repo, PROFILE_LOCAL)]:
    if os.path.isfile(x) and multiply_linked(x):
      aliased.append(x)
  if aliased:
    for x in aliased:
      log("CONFLICT: %s has more than one hard link; rewriting it would change "
          "another file too." % x)
    raise Stop("%d hard-linked destination(s). Nothing was installed." % len(aliased))

  # Classify.
  acc = accepted_hashes()
  new, upgrade, same, conflicts = [], [], [], []
  for h, p in payload:
    t = os.path.join(claude_dir, p)
    if not os.path.lexists(t):
      new.append(p)
    else:
      th = sha256_file(t)
      if th == h:
        same.append(p)
      elif th in acc.get(p, set()):
        upgrade.append(p)   # a known package byte-string from a prior release
      else:
        conflicts.append((p, th))
  if conflicts:
    for p, th in conflicts:
      log("CONFLICT: .claude/%s is neither this release nor a previously "
          "delivered one (locally modified; sha256 %s)." % (p, th[:12]))
    raise Stop("%d locally modified file(s). Nothing was installed or "
               "replaced. Send this capture to the Guide." % len(conflicts))

  backup = tempfile.mkdtemp(prefix="guided_coding_backup.", dir=work_dir())
  _BACKUP[repo] = backup

  # Old hidden profile -> CLAUDE.local.md, only when it is ours.
  old_profile = os.path.join(repo, PROFILE_OLD)
  local_profile = os.path.join(repo, PROFILE_LOCAL)
  migrated_profile = False
  if os.path.isfile(old_profile) and not os.path.islink(old_profile):
    ours = sha256_file(old_profile) in acc.get(PROFILE_OLD, set())
    if ours and not os.path.exists(local_profile):
      migrated_profile = True
      os.rename(old_profile, local_profile)
      log("Profile: %s was this package's hidden profile; renamed to %s "
          "(the package no longer uses a file named CLAUDE.md)."
          % (PROFILE_OLD, PROFILE_LOCAL))
    elif ours:
      migrated_profile = True
      shutil.move(old_profile, os.path.join(backup_dir_for(repo), PROFILE_OLD))
      log("Profile: %s was this package's OLD hidden profile and %s already "
          "exists; the old one is retired to the backup so Claude Code does "
          "not load both." % (PROFILE_OLD, PROFILE_LOCAL))
    else:
      log("Profile: a CLAUDE.md exists whose bytes are not a profile this "
          "package ever shipped; left untouched. If it IS your GuidedCoding "
          "profile (for example one the INSTALL.md interview generated or you "
          "revised), run: libtbx.guided_coding adopt-profile <repository> - "
          "that renames it to %s on your explicit say-so. Otherwise it is a "
          "project file and stays as it is." % PROFILE_LOCAL)

  # Place.
  for p in new + upgrade:
    src = os.path.join(PAYLOAD_DIR, p)
    t = os.path.join(claude_dir, p)
    d = os.path.dirname(t)
    if not os.path.isdir(d):
      os.makedirs(d)
    if os.path.exists(t):
      bd = os.path.dirname(os.path.join(backup, p))
      if not os.path.isdir(bd):
        os.makedirs(bd)
      shutil.copy2(t, os.path.join(backup, p))
    shutil.copy2(src, t)
  # Retire files the old manifest listed that this payload no longer ships,
  # but only when they are known package bytes (never a user's file).
  payload_paths = set(p for h, p in payload)
  retired = []
  for old_rel, old_h in installed.items():
    if not old_rel.startswith(".claude/"):
      continue
    inner = old_rel[len(".claude/"):]
    if inner in payload_paths or inner == "MANIFEST.sha256":
      continue
    t = os.path.join(repo, old_rel)
    if os.path.islink(t) or not nolink_path(repo, t):
      log("Retire skipped (symbolic link in path): %s" % old_rel)
      continue
    if os.path.isfile(t) and sha256_file(t) in acc.get(inner, set()):
      bd = os.path.dirname(os.path.join(backup, inner))
      if not os.path.isdir(bd):
        os.makedirs(bd)
      shutil.move(t, os.path.join(backup, inner))
      retired.append(inner)
  if retired:
    log("Retired %d file(s) this release no longer ships (copies kept in the "
        "backup): %s" % (len(retired), ", ".join(sorted(retired))))
  unknown_listed = [r for r, hh in installed.items()
                    if r.startswith(".claude/") and r[len(".claude/"):] not in payload_paths
                    and r[len(".claude/"):] != "MANIFEST.sha256" and r[len(".claude/"):] not in retired
                    and os.path.isfile(os.path.join(repo, r))]
  if unknown_listed:
    log("Left in place (listed in the old manifest but not known package bytes, "
        "so treated as yours): %s" % ", ".join(sorted(unknown_listed)))
  # Manifest: payload + the personal profile if present.
  lines = [(h, "./.claude/" + p) for h, p in payload]
  if os.path.isfile(local_profile):
    lines.append((sha256_file(local_profile), "./" + PROFILE_LOCAL))
  if os.path.exists(installed_manifest_path):
    shutil.copy2(installed_manifest_path,
                 os.path.join(backup, "MANIFEST.sha256.prev"))
  with io.open(installed_manifest_path, "w", encoding="utf-8") as f:
    for h, p in lines:
      f.write(u"%s  %s\n" % (h, p))
  log("Installed: %d new; upgraded: %d (old copies kept at %s); unchanged: %d."
      % (len(new), len(upgrade), backup, len(same)))

  # Excludes.
  excl = os.path.join(repo, ".git", "info", "exclude")
  if not os.path.isdir(os.path.dirname(excl)):
    os.makedirs(os.path.dirname(excl))
  have = []
  if os.path.isfile(excl):
    have = [l.rstrip("\n") for l in io.open(excl, encoding="utf-8")]
  if migrated_profile:
    dropped = [l for l in have if l in RETIRED_EXCLUDE_LINES]
    kept_lines = [l for l in have if l not in RETIRED_EXCLUDE_LINES]
  else:
    dropped = []
    kept_lines = list(have)
    if any(l in RETIRED_EXCLUDE_LINES for l in have):
      log("Excludes: your .git/info/exclude has a /CLAUDE.md line. This "
          "release did not write it and does not remove it, because it may "
          "be yours; if an earlier GuidedCoding install wrote it, "
          "adopt-profile removes it, or delete the line yourself.")
  for line in EXCLUDE_LINES:
    if line not in kept_lines:
      kept_lines.append(line)
  with io.open(excl, "w", encoding="utf-8") as f:
    for l in kept_lines:
      f.write(l + u"\n")
  if dropped:
    log("Excludes: removed the old %s line(s) written by earlier installs - "
        "a project CLAUDE.md is no longer hidden from git." % ", ".join(dropped))
  try:
    out = subprocess.check_output(["git"] + GIT_ISOLATION_FLAGS +
                                  ["-C", repo, "status", "--porcelain",
                                   "--untracked-files=all"], env=git_env())
    if re.search(rb"^\?\? (\.claude|CLAUDE\.local\.md)", out, re.M):
      raise Stop("'.claude/' or CLAUDE.local.md appears in git status despite "
                 "the excludes - send this capture to the Guide.")
  except subprocess.CalledProcessError:
    log("NOTE: git status could not be run here; excludes were written but "
        "not confirmed.")
  log("Excludes verified.")

  # Post-install verification of the written manifest.
  bad = []
  for h, p in validated_manifest(repo, installed_manifest_path):
    t = os.path.join(repo, p)
    if not os.path.isfile(t) or sha256_file(t) != h:
      bad.append(p)
  if bad:
    for p in bad:
      log("POST-INSTALL MISMATCH: %s" % p)
    raise Stop("the installed manifest does not verify.")
  log("Post-install verification: every file matches the manifest "
      "(%d lines)." % len(lines))
  log("")
  log("=== Install complete ===")
  log("Your files are here:")
  for h, p in lines:
    log("  %s" % os.path.join(repo, rel_of(p)))
  log("  %s" % installed_manifest_path)
  if not os.path.isfile(local_profile):
    log("NOTE: no personal profile %s in this repository yet. A Claude Code "
        "session told \"install guided_coding from libtbx\" runs the INSTALL.md "
        "interview and writes it; re-run this install afterwards so the "
        "manifest records it." % PROFILE_LOCAL)
  log("Old copies of upgraded files (safe to ignore): %s" % backup)


def cmd_verify(repo, log):
  repo = os.path.abspath(repo)
  mp = os.path.join(repo, ".claude", "MANIFEST.sha256")
  entries = validated_manifest(repo, mp)
  if not entries:
    raise Stop("no installed manifest at %s" % mp)
  ok, bad = 0, []
  for h, p in entries:
    t = os.path.join(repo, p)
    if os.path.isfile(t) and sha256_file(t) == h:
      ok += 1
    else:
      bad.append(p)
  for p in bad:
    log("NOT OK: %s" % p)
  log("verify %s: OK=%d NOT_OK=%d (installed release: %s)"
      % (repo, ok, len(bad), installed_release(repo)))
  if bad:
    raise Stop("%d file(s) do not match the manifest." % len(bad))


def installed_release(repo):
  rm = os.path.join(repo, ".claude", "RELEASE")
  if os.path.isfile(rm):
    return io.open(rm, encoding="utf-8").read().strip()
  return "(pre-libtbx package; see CLAUDE.local.md's 'this release' line)"


def cmd_remove(repo, log):
  repo = os.path.abspath(repo)
  claude_dir = os.path.join(repo, ".claude")
  mp = os.path.join(claude_dir, "MANIFEST.sha256")
  entries = validated_manifest(repo, mp)
  if not entries:
    raise Stop("no installed manifest at %s - nothing to remove" % mp)
  signals = open_change_signals(repo)
  lock_signals = [x for x in signals if x[0] == "lock"]
  other = [x for x in signals if x[0] != "lock"]
  if lock_signals:
    raise Stop("the tree lock %s exists; removal refused, and the override "
               "does not apply to the lock." % lock_signals[0][1])
  if other and os.environ.get("ALLOW_OPEN_CHANGE", "") != "yes":
    for kind, path in other:
      log("  %s: %s" % (kind, path))
    raise Stop("a change appears open; removal refused. Close it first, or "
               "re-run with ALLOW_OPEN_CHANGE=yes if the record is stale.")
  # Preflight: every path a regular file reached without any symlink component.
  for h, rel in entries:
    t = os.path.join(repo, rel)
    if os.path.lexists(t) and (os.path.islink(t) or not nolink_path(repo, t)):
      raise Stop("%s is, or lies under, a symbolic link; removal refused." % rel)
  # Deletion authority is package-shipped data ONLY: a file is removed when
  # its on-disk bytes are this payload's or a previously delivered release's.
  # The installed manifest (editable) only says which paths to look at; its
  # hashes are never the reason to delete anything.
  acc = accepted_hashes()
  payload_hashes = {}
  for h, q in read_manifest(PAYLOAD_MANIFEST):
    payload_hashes.setdefault(rel_of(q), set()).add(h)
  removed, retained, kept = 0, [], []
  for h, rel in entries:
    if rel == PROFILE_LOCAL or not rel.startswith(".claude/"):
      kept.append(rel)
      continue
    t = os.path.join(repo, rel)
    if not os.path.isfile(t):
      continue
    inner = rel[len(".claude/"):]
    known = acc.get(inner, set()) | payload_hashes.get(inner, set())
    if sha256_file(t) not in known:
      retained.append(rel)   # not known package bytes: yours, not ours to delete
      continue
    os.remove(t)
    removed += 1
  t = os.path.join(claude_dir, "MANIFEST.sha256")
  if os.path.isfile(t) and not os.path.islink(t):
    os.remove(t)
  for dp, dn, fn in os.walk(claude_dir, topdown=False):
    if dp == claude_dir or os.path.islink(dp):
      continue
    if os.path.relpath(dp, claude_dir).split(os.sep)[0] in NEVER_REMOVE:
      continue
    if not os.listdir(dp):
      os.rmdir(dp)
  log("Removed %d package files from %s (including the package's record "
      "templates); kept: %s, every non-package file under .claude/records/, "
      "and .claude/settings.local.json."
      % (removed, claude_dir, ", ".join(kept) or "nothing else"))
  if retained:
    log("Retained %d manifest-listed file(s) whose bytes are not a release "
        "this package ever shipped - they are yours, not the package's: %s"
        % (len(retained), ", ".join(retained)))
  log("Not removed on purpose: %s, your change records, HANDOFF.json, "
      "known_failure_variation.md, settings.local.json, and the "
      ".git/info/exclude lines (harmless)." % PROFILE_LOCAL)


def cmd_adopt_profile(repo, log):
  """On the developer's explicit say-so, rename CLAUDE.md to CLAUDE.local.md,
  drop the old /CLAUDE.md exclude line, and record the profile in the
  manifest. Every destination is checked BEFORE anything is changed: no
  symbolic link at or above any file this touches, the installed manifest
  valid, no change open, no CLAUDE.local.md already present."""
  repo = os.path.abspath(repo)
  old_profile = os.path.join(repo, PROFILE_OLD)
  local_profile = os.path.join(repo, PROFILE_LOCAL)
  excl = os.path.join(repo, ".git", "info", "exclude")
  mp = os.path.join(repo, ".claude", "MANIFEST.sha256")
  # --- preflight: nothing below mutates ---
  if not is_git_repo(repo):
    raise Stop("%s is not a git repository." % repo)
  if not os.path.lexists(old_profile) or os.path.islink(old_profile) \
     or not os.path.isfile(old_profile):
    raise Stop("no ordinary file %s in %s" % (PROFILE_OLD, repo))
  if os.path.lexists(local_profile):
    raise Stop("%s already exists; nothing renamed. Merge them yourself." % PROFILE_LOCAL)
  for x in [os.path.join(repo, ".git"), os.path.join(repo, ".git", "info"),
            os.path.join(repo, ".claude"), mp, excl, old_profile]:
    if os.path.islink(x):
      raise Stop("%s is a symbolic link; adopt-profile refuses to write "
                 "through links. Nothing was changed." % x)
  if os.path.lexists(excl) and not os.path.isfile(excl):
    raise Stop("%s exists and is not a regular file. Nothing was changed." % excl)
  for x in [excl, mp, old_profile]:
    if os.path.isfile(x) and multiply_linked(x):
      raise Stop("%s has more than one hard link; adopt-profile refuses to "
                 "rewrite an aliased file. Nothing was changed." % x)
  entries = []
  if os.path.lexists(mp):
    if not os.path.isfile(mp):
      raise Stop("%s exists and is not a regular file. Nothing was changed." % mp)
    entries = validated_manifest(repo, mp)   # raises Stop on any bad entry
  signals = open_change_signals(repo)
  if signals:
    raise Stop("a change appears open in this repository; adopt the profile "
               "after it has closed. Nothing was changed.")
  # --- all preconditions hold: mutate, in the order least harmful to interrupt ---
  os.rename(old_profile, local_profile)
  log("Renamed %s -> %s on your instruction." % (PROFILE_OLD, PROFILE_LOCAL))
  have = []
  if os.path.isfile(excl):
    have = [l.rstrip("\n") for l in io.open(excl, encoding="utf-8")]
  kept = [l for l in have if l not in RETIRED_EXCLUDE_LINES]
  for line in EXCLUDE_LINES:
    if line not in kept:
      kept.append(line)
  d = os.path.dirname(excl)
  if not os.path.isdir(d):
    os.makedirs(d)
  with io.open(excl, "w", encoding="utf-8") as f:
    for l in kept:
      f.write(l + u"\n")
  log("Excludes updated: /CLAUDE.md removed, /CLAUDE.local.md and .claude/ present.")
  if os.path.isfile(mp):
    entries = [(h, r) for h, r in entries if r not in (PROFILE_OLD, PROFILE_LOCAL)]
    entries.append((sha256_file(local_profile), PROFILE_LOCAL))
    with io.open(mp, "w", encoding="utf-8") as f:
      for h, r in entries:
        f.write(u"%s  ./%s\n" % (h, r))
    log("Manifest updated to record %s." % PROFILE_LOCAL)


def cmd_status(repos, log):
  log("package release: %s" % release_string())
  for repo in repos:
    repo = os.path.abspath(repo)
    log("%s: %s" % (repo, installed_release(repo)))


def main(argv):
  if len(argv) < 2 or argv[0] not in ("install", "verify", "remove", "status", "adopt-profile"):
    print(__doc__)
    return 2
  cmd = argv[0]
  if sys.platform.startswith("win"):
    print("STOP: GuidedCoding supports macOS and Linux; Windows is not supported.")
    return 2
  log = Log("guided_coding_" + cmd)
  try:
    if cmd == "install":
      cmd_install(argv[1], log)
    elif cmd == "verify":
      cmd_verify(argv[1], log)
    elif cmd == "remove":
      cmd_remove(argv[1], log)
    elif cmd == "adopt-profile":
      cmd_adopt_profile(argv[1], log)
    else:
      cmd_status(argv[1:], log)
    rc = 0
  except Stop as e:
    log("STOP: %s" % e)
    log("Nothing further was done.")
    rc = 2
  except Exception as e:   # an unexpected error is still captured, then re-raised
    import traceback
    log("ERROR: %s" % e)
    log(traceback.format_exc())
    log.finish()
    raise
  log.finish()
  return rc


if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
