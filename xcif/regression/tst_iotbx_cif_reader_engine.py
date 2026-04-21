"""Tests for the engine= kwarg on iotbx.cif.reader.

Covers:
- Unknown engine values are rejected with ValueError before any file I/O.
- Known engines ("ucif", "xcif") are accepted.
- engine="xcif" with file_path dispatches to xcif_ext.parse_file (zero-copy
  mmap) rather than reading the whole file into a Python string.
"""
from __future__ import absolute_import, division, print_function

import os
import tempfile

import iotbx.cif
import xcif_ext


def exercise_engine_validated_before_io():
  # A bad engine value must be rejected before the reader touches the
  # filesystem. Use a path that does not exist: if engine validation
  # happens first we get ValueError from the engine check; if file I/O
  # happens first we get an OSError/IOError/Sorry from smart_open.
  bogus_path = os.path.join(
    tempfile.gettempdir(),
    "definitely_not_a_real_cif_file_%d.cif" % os.getpid())
  # Make sure it really does not exist
  if os.path.exists(bogus_path):
    os.remove(bogus_path)
  try:
    iotbx.cif.reader(file_path=bogus_path, engine="not_a_real_engine")
  except ValueError as e:
    msg = str(e)
    assert "engine" in msg, \
      "expected engine-related ValueError, got: %s" % msg
    return
  except Exception as e:
    raise AssertionError(
      "engine validation should happen before file I/O; "
      "got %s instead of ValueError: %s" % (type(e).__name__, e))
  raise AssertionError("expected ValueError for bad engine, no exception raised")


def exercise_known_engines_accepted():
  input_string = "data_t\n_k v\n"
  for engine in ("ucif", "xcif"):
    reader = iotbx.cif.reader(input_string=input_string, engine=engine)
    model = reader.model()
    assert "t" in model, "block not parsed with engine=%s" % engine
    assert model["t"]["_k"] == "v", \
      "pair not parsed with engine=%s: got %r" % (engine, model["t"]["_k"])


_SAMPLE_CIF = """\
data_t
_k v
loop_
_a.id
_a.val
 1 x
 2 y
"""


def _write_sample(tmpdir):
  p = os.path.join(tmpdir, "sample.cif")
  with open(p, "w") as f:
    f.write(_SAMPLE_CIF)
  return p


def exercise_file_path_uses_parse_file():
  # Record which xcif_ext entry point the reader invokes.
  real_parse = xcif_ext.parse
  real_parse_file = xcif_ext.parse_file
  calls = {"parse": 0, "parse_file": 0}

  def _tracking_parse(text, **kw):
    calls["parse"] += 1
    return real_parse(text, **kw)

  def _tracking_parse_file(path, **kw):
    calls["parse_file"] += 1
    return real_parse_file(path, **kw)

  xcif_ext.parse = _tracking_parse
  xcif_ext.parse_file = _tracking_parse_file
  try:
    tmpdir = tempfile.mkdtemp(prefix="xcif_reader_item5_")
    try:
      path = _write_sample(tmpdir)

      # input_string + xcif: parse() is used, parse_file() is NOT.
      iotbx.cif.reader(input_string=_SAMPLE_CIF, engine="xcif").model()
      assert calls["parse"] == 1, \
        "input_string path should call parse: %r" % calls
      assert calls["parse_file"] == 0, \
        "input_string path must NOT call parse_file: %r" % calls

      # file_path + xcif: parse_file() is used, parse() is NOT.
      calls["parse"] = 0
      calls["parse_file"] = 0
      iotbx.cif.reader(file_path=path, engine="xcif").model()
      assert calls["parse_file"] == 1, \
        "file_path path should call parse_file: %r" % calls
      assert calls["parse"] == 0, \
        "file_path path must NOT call parse (avoid Python-string copy): %r" % (
          calls,)
    finally:
      import shutil
      shutil.rmtree(tmpdir, ignore_errors=True)
  finally:
    xcif_ext.parse = real_parse
    xcif_ext.parse_file = real_parse_file


def exercise_file_path_matches_input_string():
  tmpdir = tempfile.mkdtemp(prefix="xcif_reader_item5_")
  try:
    path = _write_sample(tmpdir)
    m_str = iotbx.cif.reader(
      input_string=_SAMPLE_CIF, engine="xcif").model()
    m_file = iotbx.cif.reader(file_path=path, engine="xcif").model()
    # Compare block set, pair items, loop tag+content.
    assert list(m_str.keys()) == list(m_file.keys())
    for name in m_str:
      b_s = m_str[name]
      b_f = m_file[name]
      assert list(b_s._set) == list(b_f._set), \
        "source-order drift between input_string and file_path"
      assert dict(b_s._items) == dict(b_f._items)
      assert set(b_s.loops.keys()) == set(b_f.loops.keys())
      for ln in b_s.loops:
        lp_s = b_s.loops[ln]
        lp_f = b_f.loops[ln]
        assert set(lp_s.keys()) == set(lp_f.keys())
        for tag in lp_s.keys():
          assert list(lp_s[tag]) == list(lp_f[tag])
  finally:
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)


def run():
  exercise_engine_validated_before_io()
  exercise_known_engines_accepted()
  exercise_file_path_uses_parse_file()
  exercise_file_path_matches_input_string()
  print("OK")


if __name__ == "__main__":
  run()
