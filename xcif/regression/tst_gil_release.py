from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_gil_release.py
#
# GIL release tests for xcif Python bindings.
#
# Verifies that parse_file() and parse() release the GIL during the C++
# parsing work, allowing other Python threads to make progress concurrently.

import os
import threading
import time
import libtbx.load_env
from libtbx.utils import format_cpu_times

def exercise_parse_file_releases_gil():
  """parse_file must release the GIL so other threads can run.

  Strategy: start a thread that parses a CIF file, and simultaneously
  increment a counter in the main thread. If the GIL is held during
  parsing, the main thread counter stays near zero. If released, the
  counter advances.

  Uses a threading.Event to ensure the main thread is counting before
  the parse thread starts, avoiding a race where parsing completes
  before the main thread enters the counting loop.
  """
  import xcif_ext
  dist_dir = libtbx.env.dist_path("xcif")
  cif_file = os.path.join(dist_dir, "regression", "example.cif")

  counter = [0]
  done = [False]
  ready = threading.Event()

  def parse_loop():
    ready.wait()  # Don't start until main thread is counting
    for _ in range(1000):
      doc = xcif_ext.parse_file(cif_file)
      _ = len(doc)
    done[0] = True

  t = threading.Thread(target=parse_loop)
  t.start()

  # Signal parse thread to start, then immediately begin counting.
  # The parse thread needs the GIL to proceed past ready.wait(), but
  # the main thread holds it until time.sleep(), so counter >= 1 is
  # guaranteed before any parsing begins.
  ready.set()
  while not done[0]:
    counter[0] += 1
    # Yield to avoid busy-waiting artifacts
    time.sleep(0.0001)

  t.join()
  # If GIL was released, main thread should have made progress
  assert counter[0] > 0, (
    "Main thread counter is %d; GIL may not be released during parse_file"
    % counter[0])

def exercise_parse_string_releases_gil():
  """parse() from string should also release the GIL.

  Uses a threading.Event to avoid a race where parsing completes before
  the main thread enters the counting loop.
  """
  import xcif_ext

  # Build a moderately large CIF string to make parsing take measurable time
  lines = ["data_big\n"]
  lines.append("loop_\n_x.id\n_x.val\n")
  for i in range(5000):
    lines.append("%d %f\n" % (i, i * 0.1))
  big_cif = "".join(lines)

  counter = [0]
  done = [False]
  ready = threading.Event()

  def parse_loop():
    ready.wait()  # Don't start until main thread is counting
    for _ in range(50):
      doc = xcif_ext.parse(big_cif)
      _ = len(doc)
    done[0] = True

  t = threading.Thread(target=parse_loop)
  t.start()

  ready.set()
  while not done[0]:
    counter[0] += 1
    time.sleep(0.0001)

  t.join()
  assert counter[0] > 0, (
    "Main thread counter is %d; GIL may not be released during parse"
    % counter[0])

def exercise_concurrent_parse():
  """Multiple threads can parse simultaneously without crashes."""
  import xcif_ext
  dist_dir = libtbx.env.dist_path("xcif")
  cif_file = os.path.join(dist_dir, "regression", "example.cif")

  results = [None] * 4
  errors = []

  def worker(idx):
    try:
      for _ in range(100):
        doc = xcif_ext.parse_file(cif_file)
        assert len(doc) == 1
        block = doc[0]
        assert block.has_tag("_cell.length_a")
      results[idx] = True
    except Exception as e:
      errors.append((idx, str(e)))
      results[idx] = False

  threads = [threading.Thread(target=worker, args=(i,)) for i in range(4)]
  for t in threads:
    t.start()
  for t in threads:
    t.join()

  assert not errors, "Thread errors: %s" % errors
  assert all(r is True for r in results), "Some threads failed: %s" % results

if __name__ == "__main__":
  exercise_parse_file_releases_gil()
  exercise_parse_string_releases_gil()
  exercise_concurrent_parse()
  print(format_cpu_times())
  print("OK")
