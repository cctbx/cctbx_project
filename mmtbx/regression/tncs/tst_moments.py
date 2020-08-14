from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex # import dependency
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_ncs_ext")
import iotbx.pdb
from mmtbx.ncs import tncs
import libtbx.load_env
import os, time


def exercise_01(reflections_per_bin=150):
  """
  tncs_epsfac calculation with radius refinement
  """
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/tncs/model_3.pdb",
    test=os.path.isfile)
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  xray_structure = pdb_inp.xray_structure_simple()
  xray_structure = xray_structure.set_b_iso(value=10)
  f_obs = abs(xray_structure.structure_factors(d_min=2).f_calc())
  f_obs = f_obs.set_observation_type_xray_amplitude()
  result = tncs.compute_eps_factor(
    f_obs               = f_obs,
    pdb_hierarchy       = pdb_inp.construct_hierarchy(),
    reflections_per_bin = reflections_per_bin)
  result.show_summary()
  M2 = f_obs.second_moments_centric_acentric(
    reflections_per_bin = reflections_per_bin)
  M2_corr = f_obs.second_moments_centric_acentric(
    reflections_per_bin = reflections_per_bin,
    eps_fac = result.epsfac)
  print("%s %4.2f"%M2[0], "%s %4.2f"%M2_corr[0])

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_01()
  print("Time: %6.3f"%(time.time()-t0))
  print("OK")
