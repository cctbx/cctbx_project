"""Tests for the chat logging helpers: session-log path generation and the
API-key redaction helper used before any payload is written to disk."""

import os
import shutil
import tempfile
from pathlib import Path

from libtbx.utils import format_cpu_times

from qttbx.widgets.chat.logging_setup import (
  LOG_KEEP, _prune_old_logs, open_raw_log, open_session_log, redact_secrets)


def exercise_open_session_log_creates_dir_and_writes():
  tmp = Path(tempfile.mkdtemp())
  try:
    log, path = open_session_log(chat_root=tmp)
    try:
      print("hello", file=log)
      log.flush()
      content = Path(path).read_text()
      assert "hello" in content
    finally:
      log.close()
    assert (tmp / "logs").is_dir()
  finally:
    shutil.rmtree(tmp)


def exercise_redact_secrets_handles_anthropic_keys():
  s = 'x-api-key: sk-ant-abc123XYZ'
  assert "sk-ant-" not in redact_secrets(s)
  assert "<REDACTED>" in redact_secrets(s)


def exercise_redact_secrets_handles_authorization_header():
  s = 'Authorization: Bearer eyJabcdefghijklmnop'
  out = redact_secrets(s)
  assert "Bearer" in out
  assert "eyJabcdefghijklmnop" not in out
  assert "<REDACTED>" in out


def exercise_redact_secrets_passes_clean_text():
  assert redact_secrets("nothing to hide") == "nothing to hide"


def exercise_redact_secrets_is_idempotent():
  s = ('x-api-key: sk-ant-abc123 '
       'Authorization: Bearer eyJxyz '
       '"api_key": "sk-ant-def456"')
  once = redact_secrets(s)
  twice = redact_secrets(once)
  assert once == twice, (once, twice)


def exercise_prune_keeps_only_log_keep_most_recent():
  tmp = Path(tempfile.mkdtemp())
  try:
    log_dir = tmp / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    # Seed 10 pre-existing logs with distinct mtimes (oldest first).
    seeded = []
    base_mtime = 1_000_000_000.0
    for i in range(10):
      p = log_dir / ("chat-seed-%02d.log" % i)
      p.write_text("seed %d" % i)
      mtime = base_mtime + i  # later index => more recent
      os.utime(p, (mtime, mtime))
      seeded.append((p, mtime))
    log, path = open_session_log(chat_root=tmp)
    try:
      surviving = sorted(log_dir.glob("chat-*.log"))
      # The newly created log should exist.
      assert path.exists()
      # Pre-existing survivors: LOG_KEEP most-recent by mtime.
      expected_survivors = {p for p, _ in seeded[-LOG_KEEP:]}
      actual_seeded = {p for p in surviving if p.name.startswith("chat-seed-")}
      assert actual_seeded == expected_survivors, (
        actual_seeded, expected_survivors)
      # Total files on disk: LOG_KEEP pre-existing + 1 newly created.
      assert len(surviving) == LOG_KEEP + 1, surviving
    finally:
      log.close()
  finally:
    shutil.rmtree(tmp)


def exercise_prune_survives_unstattable_entry():
  """A globbed log whose stat() fails mid-prune (e.g. a TOCTOU delete between
  glob and stat) must not crash _prune_old_logs: the unstattable entry sorts as
  oldest and the guarded unlink handles it. Regression for the sort key that
  used an unguarded p.stat()."""
  import pathlib
  tmp = Path(tempfile.mkdtemp())
  try:
    log_dir = tmp / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    for i in range(LOG_KEEP + 2):
      (log_dir / ("coot-%02d.log" % i)).write_text("x")
    (log_dir / "coot-vanish.log").write_text("x")   # stat faked to fail below
    real_stat = pathlib.Path.stat
    def fake_stat(self, *args, **kwargs):
      if self.name == "coot-vanish.log":
        raise FileNotFoundError(2, "No such file", str(self))
      return real_stat(self, *args, **kwargs)
    pathlib.Path.stat = fake_stat
    try:
      _prune_old_logs(log_dir, "coot-*.log")          # must not raise
    finally:
      pathlib.Path.stat = real_stat
    remaining = sorted(p.name for p in log_dir.glob("coot-*.log"))
    assert len(remaining) == LOG_KEEP, remaining
    assert "coot-vanish.log" not in remaining, remaining
  finally:
    shutil.rmtree(tmp)


def exercise_open_raw_log_rotates_and_is_raw():
  """open_raw_log returns a REAL file handle (has fileno, unlike the redacting
  session/debug logs) under chat_root/logs, writes through un-redacted, and
  rotates <prefix>-*.log to LOG_KEEP -- the unified machinery for the Coot
  subprocess-capture log, so callers don't reinvent path/prune/open."""
  tmp = Path(tempfile.mkdtemp())
  try:
    log_dir = tmp / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    for i in range(LOG_KEEP + 3):
      (log_dir / ("coot-seed-%02d.log" % i)).write_text("x")
    fh, path = open_raw_log(tmp, "coot")
    try:
      assert path.parent == log_dir, path
      assert path.name.startswith("coot-") and path.name.endswith(".log"), path
      fh.write("sk-ant-not-redacted\n")    # raw log: no redaction wrapper
      fh.flush()
      fh.fileno()                          # a real fd (redacting wrapper raises)
    finally:
      fh.close()
    assert path.read_text() == "sk-ant-not-redacted\n", path.read_text()
    # Pre-existing coot-*.log pruned to LOG_KEEP (the new file is created after
    # pruning, so it isn't counted against the seeds).
    seeds = [p for p in log_dir.glob("coot-*.log") if "seed" in p.name]
    assert len(seeds) == LOG_KEEP, [p.name for p in seeds]
  finally:
    shutil.rmtree(tmp)


def exercise_log_dir_and_files_are_owner_only():
  """On a shared project dir the chat logs (only best-effort-redacted, and the
  raw coot log un-redacted) must not be world-readable: the logs dir is created
  0700 and each log file 0600 via an explicit chmod, not the ambient umask.
  POSIX-only (Windows perms don't map)."""
  if os.name != "posix":
    return
  import stat
  tmp = Path(tempfile.mkdtemp())
  try:
    log, path = open_session_log(chat_root=tmp)
    try:
      log_dir = tmp / "logs"
      assert stat.S_IMODE(log_dir.stat().st_mode) == 0o700, \
        oct(stat.S_IMODE(log_dir.stat().st_mode))
      assert stat.S_IMODE(Path(path).stat().st_mode) == 0o600, \
        oct(stat.S_IMODE(Path(path).stat().st_mode))
    finally:
      log.close()
    # The raw (un-redacted) coot capture log is tightened too.
    fh, rawpath = open_raw_log(tmp, "coot")
    try:
      assert stat.S_IMODE(Path(rawpath).stat().st_mode) == 0o600, \
        oct(stat.S_IMODE(Path(rawpath).stat().st_mode))
    finally:
      fh.close()
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_open_session_log_creates_dir_and_writes()
  exercise_log_dir_and_files_are_owner_only()
  exercise_redact_secrets_handles_anthropic_keys()
  exercise_redact_secrets_handles_authorization_header()
  exercise_redact_secrets_passes_clean_text()
  exercise_redact_secrets_is_idempotent()
  exercise_prune_keeps_only_log_keep_most_recent()
  exercise_prune_survives_unstattable_entry()
  exercise_open_raw_log_rotates_and_is_raw()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
