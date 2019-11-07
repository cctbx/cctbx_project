from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from mmtbx.secondary_structure.build.tst_2 import tst_01_start_lines
import libtbx.load_env
import os.path
import time

def exercise_01(prefix="tst_mi_map_test_01"):
  """
  Simple run to a completion with reference map. no SS annotations.
  """
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(tst_01_start_lines)
  pdb_file.close()
  cmd = " ".join([
      "phenix.model_idealization",
      "%s_start.pdb" % prefix,
      "use_map_for_reference=True",
      "number_of_refinement_cycles=1",
      "run_minimization_first=False",
      "loop_idealization.number_of_ccd_trials=1",
      "n_macro=1",
      "debug=True",
      ">%s.log" % prefix])
  print(cmd)
  assert not easy_run.call(cmd)
  res_log = open("%s.log" % prefix, "r")
  log_lines = res_log.readlines()
  for l in [
      # "Secondary structure substitution step will be skipped\n",
      "  Minimizing...\n",
      "Using map as reference\n",
      # "Ramachandran outliers:      0.00      0.00      0.00      0.00      0.00\n",
      "All done.\n"]:
    assert l in log_lines, "'%s' not in log file." % l
  res_log.close()
  # assert os.path.isfile("%s_start.pdb_idealized.pdb" % prefix)

if (__name__ == "__main__"):
  t0 = time.time()
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping: probe not configured")
  else:
    exercise_01()
  print("Time: %.2f" % (time.time() - t0))
  print("OK")
