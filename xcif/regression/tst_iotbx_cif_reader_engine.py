"""Tests for the engine= kwarg on iotbx.cif.reader.

Covers:
- Unknown engine values are rejected with ValueError before any file I/O.
- Known engines ("ucif", "xcif") are accepted.
"""
from __future__ import absolute_import, division, print_function

import os
import tempfile

import iotbx.cif


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


def run():
  exercise_engine_validated_before_io()
  exercise_known_engines_accepted()
  print("OK")


if __name__ == "__main__":
  run()
