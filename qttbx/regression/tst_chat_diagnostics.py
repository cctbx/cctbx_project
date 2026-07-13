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


def exercise_log_tail_does_not_read_whole_file():
  """The log tail must be produced with collections.deque(fh, maxlen=N), NOT
  fh.readlines()[-N:] which slurps a possibly-huge log fully into memory on the
  GUI thread. Inject a file object that is line-iterable but whose readlines()
  raises: with the deque tail the dump still builds and holds exactly the last
  N lines; the old readlines() path would raise here instead."""
  from qttbx.widgets.chat import diagnostics as diag

  class _NoReadlinesFile:
    """Line-iterable stand-in for an open log whose readlines() is sabotaged,
    so any reliance on readlines() surfaces instead of silently passing."""
    def __init__(self, lines):
      self._lines = list(lines)
    def __enter__(self):
      return self
    def __exit__(self, *exc):
      return False
    def __iter__(self):
      return iter(self._lines)
    def readlines(self):
      raise AssertionError("must tail via deque iteration, not readlines()")

  all_lines = ["line%d\n" % i for i in range(1000)]

  def fake_open(path, *args, **kwargs):
    return _NoReadlinesFile(all_lines)

  diag.open = fake_open
  try:
    text = diag.build_diagnostics(
      profile=_MiniProfile(),
      storage=_MiniStorage(project_dir="/x", chat_root="/x"),
      log_path="/does/not/matter.log",
      log_tail_lines=5)
  finally:
    del diag.open

  # Exactly the last 5 lines (995..999) are present ...
  for i in range(995, 1000):
    assert "line%d" % i in text, (i, text)
  # ... and lines before the tail window are not.
  assert "line994\n" not in text, text
  assert "line0\n" not in text, text


def exercise():
  exercise_dump_contains_basic_fields()
  exercise_no_log_file_is_ok()
  exercise_non_utf8_log_bytes_are_readable()
  exercise_log_tail_does_not_read_whole_file()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
