from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from mmtbx.secondary_structure.build.tst_2 import tst_01_start_lines
import libtbx.load_env
import os.path
import time

def exercise_02(prefix="tst_mi_map_test_02"):
  """
  Simple run to a completion with reference map. With SS annotations.
  """
  h_records = """\
HELIX    1   1 PRO A    3  ALA A   21  1                                  19
HELIX    2   2 ARG A   23  GLN A   44  1                                  22
"""
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(h_records)
  pdb_file.write(tst_01_start_lines)
  pdb_file.close()
  cmd = " ".join([
      "phenix.model_idealization",
      "use_map_for_reference=True",
      "loop_idealization.number_of_ccd_trials=1",
      "number_of_refinement_cycles=1",
      "n_macro=1",
      "debug=True",
      "%s_start.pdb" % prefix,
      ">%s.log" % prefix])
  print(cmd)
  assert not easy_run.call(cmd)
  res_log = open("%s.log" % prefix, "r")
  log_lines = res_log.readlines()
  for l in [
      # "Replacing ss-elements with ideal ones:\n",
      "  Minimizing...\n",
      "Using map as reference\n",
      # "Ramachandran outliers:      0.00      0.00      0.00      0.00      0.00\n",
      "All done.\n"]:
    assert l in log_lines, "'%s' not in log file." % l
  res_log.close()
  assert os.path.isfile("%s_start.pdb_all_idealized.pdb" % prefix)
  res_pdb = open("%s_start.pdb_all_idealized.pdb" % prefix, "r")
  res_pdb_lines = res_pdb.readlines()
  res_pdb.close()
  for l in [
      "HELIX    1   1 PRO A    3  ALA A   21  1                                  19\n",
      "HELIX    2   2 ARG A   23  GLN A   44  1                                  22\n",
      ]:
    assert l in res_pdb_lines, "'%s' not in pdb file." % l

if (__name__ == "__main__"):
  t0 = time.time()
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping: probe not configured")
  else:
    exercise_02()
  print("Time: %.2f" % (time.time() - t0))
  print("OK")
