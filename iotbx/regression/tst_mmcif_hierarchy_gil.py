"""
build_mmcif_hierarchy releases the GIL for the C++ build
(iotbx/pdb/mmcif_hierarchy_bpl.cpp): another Python thread makes progress
while a hierarchy is built.
"""
import threading
import time
from cctbx.array_family import flex
from iotbx.pdb import atom_site_columns, build_mmcif_hierarchy

N_ROWS = 300000

def synthetic_columns(n):
  """n ATOM rows in one chain, ten atoms per residue."""
  cols = atom_site_columns()
  cols.type_symbol = flex.std_string(["C"] * n)
  cols.atom_name = flex.std_string(["C%d" % (i % 10) for i in range(n)])
  cols.alt_id = flex.std_string(["."] * n)
  cols.label_asym_id = flex.std_string(["A"] * n)
  cols.auth_asym_id = flex.std_string(["A"] * n)
  cols.comp_id = flex.std_string(["ALA"] * n)
  cols.seq_id = flex.std_string([str(i // 10 + 1) for i in range(n)])
  cols.cartn_x = flex.std_string(["1.0"] * n)
  cols.cartn_y = flex.std_string(["2.0"] * n)
  cols.cartn_z = flex.std_string(["3.0"] * n)
  cols.b_iso = flex.std_string(["10.0"] * n)
  cols.occupancy = flex.std_string(["1.0"] * n)
  return cols

def exercise_build_releases_gil():
  """
  One build in a worker thread while the main thread wakes every 0.1 ms
  and records the time. With the GIL held for the build the main thread's
  longest pause is the whole build; with it released, the sleep
  granularity (16 ms on Windows before Python 3.11).
  """
  cols = synthetic_columns(N_ROWS)
  done = [False]
  result = {}
  ready = threading.Event()
  def build():
    ready.wait()
    try:
      t0 = time.perf_counter()
      result["n_atoms"] = build_mmcif_hierarchy(cols).atoms_size()
      result["build_time"] = time.perf_counter() - t0
    except Exception as e:
      result["error"] = e
    finally:
      done[0] = True
  t = threading.Thread(target=build)
  t.start()
  ready.set()
  ticks = [time.perf_counter()]
  while not done[0]:
    time.sleep(0.0001)
    ticks.append(time.perf_counter())
  t.join()
  assert "error" not in result, result["error"]
  assert result["n_atoms"] == N_ROWS, result["n_atoms"]
  longest_pause = max(b - a for a, b in zip(ticks, ticks[1:]))
  assert longest_pause < 0.5 * result["build_time"], (
    "main thread paused %.3f s during a %.3f s build; GIL not released?"
    % (longest_pause, result["build_time"]))

def exercise_concurrent_builds():
  """Four threads build at once; every hierarchy is complete."""
  cols = synthetic_columns(N_ROWS // 10)
  results = [None] * 4
  errors = []
  def build(i):
    try:
      results[i] = build_mmcif_hierarchy(cols).atoms_size()
    except Exception as e:
      errors.append((i, str(e)))
  threads = [threading.Thread(target=build, args=(i,)) for i in range(4)]
  for t in threads:
    t.start()
  for t in threads:
    t.join()
  assert not errors, errors
  assert results == [N_ROWS // 10] * 4, results

def run():
  exercise_build_releases_gil()
  exercise_concurrent_builds()

if __name__ == "__main__":
  run()
  print("OK")
