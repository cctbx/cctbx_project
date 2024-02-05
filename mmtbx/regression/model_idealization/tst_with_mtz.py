from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os.path
import time
import iotbx
from mmtbx.utils import fmodel_from_xray_structure
from mmtbx.secondary_structure.build.tst_2 import tst_01_start_lines
import mmtbx.programs.fmodel

def exercise_01(prefix="tst_mi_mtz_01"):
  """
  Simple run to a completion with reference map. no SS annotations.
  """
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(tst_01_start_lines)
  pdb_file.close()

  pdb_inp = iotbx.pdb.input(file_name="%s_start.pdb" % prefix)
  pdb_h = pdb_inp.construct_hierarchy()
  xrs = pdb_h.extract_xray_structure(crystal_symmetry=pdb_inp.crystal_symmetry())

  params = mmtbx.programs.fmodel.master_phil.extract()
  params.high_resolution = 3
  params.low_resolution = 20
  params.output.label="FOBS"
  params.output.type="real"
  params.r_free_flags_fraction=0.1
  mmtbx.utils.fmodel_from_xray_structure(
    xray_structure = xrs,
    f_obs          = None,
    add_sigmas     = params.add_sigmas,
    params         = params,
    twin_law       = params.twin_law,
    twin_fraction  = params.twin_fraction,
    out            = None).write_to_file(file_name = "%s_start.mtz" % prefix,
      obs_type=params.output.obs_type)

  cmd = " ".join([
      "phenix.model_idealization",
      "%s_start.pdb" % prefix,
      "%s_start.mtz" % prefix,
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
      "Processing input hkl file...\n",
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
