"""Diagnostics generator test (pure-Python; no Qt needed)."""

import shutil
import tempfile
from pathlib import Path

from libtbx.utils import format_cpu_times

from qttbx.widgets.chat.diagnostics import build_diagnostics


class _MiniProfile:
  def __init__(self, name="phenix_expert", model="claude-opus-4-7"):
    self.name = name
    self.model = model


class _MiniStorage:
  def __init__(self, project_dir, chat_root):
    self.project_dir = project_dir
    self.chat_root = chat_root


def exercise_dump_contains_basic_fields():
  tmp = Path(tempfile.mkdtemp())
  try:
    log = tmp / "logs" / "chat-x.log"
    log.parent.mkdir(parents=True)
    log.write_text("line1\nline2\nline3\n")
    text = build_diagnostics(
      profile=_MiniProfile(),
      storage=_MiniStorage(project_dir=tmp, chat_root=tmp),
      log_path=log,
      log_tail_lines=10)
    assert "phenix_expert" in text
    assert "claude-opus-4-7" in text
    assert "line3" in text
    assert str(tmp) in text
  finally:
    shutil.rmtree(tmp)


def exercise_no_log_file_is_ok():
  tmp = Path(tempfile.mkdtemp())
  try:
    text = build_diagnostics(
      profile=_MiniProfile(),
      storage=_MiniStorage(project_dir=tmp, chat_root=tmp),
      log_path=tmp / "missing.log")
    assert "phenix_expert" in text
  finally:
    shutil.rmtree(tmp)


def exercise_non_utf8_log_bytes_are_readable():
  """Regression: the session-log tail is read with a strict utf-8 decode, but
  a log written under a legacy non-utf-8 locale can hold bytes that aren't
  valid utf-8. A strict read raises UnicodeDecodeError -- a ValueError, NOT the
  OSError the tail's except clause catches -- so the whole dump would crash.
  The read must tolerate undecodable bytes (errors='replace')."""
  tmp = Path(tempfile.mkdtemp())
  try:
    log = tmp / "logs" / "chat-x.log"
    log.parent.mkdir(parents=True)
    log.write_bytes(b"line1\nlegacy \xff\xfe bytes\nline3\n")
    text = build_diagnostics(
      profile=_MiniProfile(),
      storage=_MiniStorage(project_dir=tmp, chat_root=tmp),
      log_path=log,
      log_tail_lines=10)
    assert "phenix_expert" in text
    assert "line3" in text
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_dump_contains_basic_fields()
  exercise_no_log_file_is_ok()
  exercise_non_utf8_log_bytes_are_readable()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
