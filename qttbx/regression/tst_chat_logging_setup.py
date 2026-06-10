"""Tests for the chat logging helpers: session-log path generation and the
API-key redaction helper used before any payload is written to disk."""

import os
import shutil
import tempfile
from pathlib import Path

from libtbx.utils import format_cpu_times

from qttbx.widgets.chat.logging_setup import (
  LOG_KEEP, open_session_log, redact_secrets, session_log_path)


def exercise_session_log_path_format():
  tmp = Path(tempfile.mkdtemp())
  try:
    p = session_log_path(chat_root=tmp)
    assert p.parent == tmp / "logs"
    assert p.name.startswith("chat-")
    assert p.name.endswith(".log")
  finally:
    shutil.rmtree(tmp)


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


def exercise():
  exercise_session_log_path_format()
  exercise_open_session_log_creates_dir_and_writes()
  exercise_redact_secrets_handles_anthropic_keys()
  exercise_redact_secrets_handles_authorization_header()
  exercise_redact_secrets_passes_clean_text()
  exercise_redact_secrets_is_idempotent()
  exercise_prune_keeps_only_log_keep_most_recent()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
