from __future__ import absolute_import, division, print_function
import os
import subprocess
import sys
import iotbx.pdb
from mmtbx.pair_interaction import pair_interaction

def fragments():
  """
  Fragments (cluster + buffer) computed one after another in the same process.
  Element order in get_interactions follows set() iteration order, which
  depends on string hashing, so the result must not depend on what was
  computed before.
  """
  file_name = os.path.join(pair_interaction.qr_unit_tests_data, "2lvr.pdb")
  ph = iotbx.pdb.input(file_name=file_name).construct_hierarchy()
  clusters = [[1, 2],
              [3, 4, 5, 13, 14, 15, 16, 17, 18, 19, 20, 21],
              [6, 7, 8, 9, 10, 11, 12, 22, 23, 26, 31],
              [24, 25, 27, 28, 29, 30]]
  expected = [
    [1, 2, 3, 4],
    [1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
     23, 24, 25],
    [3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
     31],
    [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]]
  for cluster, e in zip(clusters, expected):
    fragment = pair_interaction.run(ph.deep_copy(), core=cluster)[2]
    assert sorted(fragment) == e, [cluster, sorted(fragment)]

def run():
  if pair_interaction.qr_unit_tests_data is None: return # skip test
  # Hash seeds that gave wrong fragments when wave functions were kept from
  # previous calls.
  for seed in ["11", "15"]:
    env = dict(os.environ, PYTHONHASHSEED=seed)
    rc = subprocess.call([sys.executable, __file__, "fragments"], env=env)
    assert rc == 0, "PYTHONHASHSEED=%s" % seed

if __name__=="__main__":
  if sys.argv[1:] == ["fragments"]:
    fragments()
  else:
    run()
