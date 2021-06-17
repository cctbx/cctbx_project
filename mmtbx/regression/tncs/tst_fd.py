from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex # import dependency
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_ncs_ext")
import iotbx.pdb
from scitbx.array_family import flex
from mmtbx.ncs import tncs
import libtbx.load_env
import os, time

def f_obs_and_tncs_pairs_from_pdb(file_name, reflections_per_bin):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/tncs/%s"%file_name,
    test=os.path.isfile)
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  xray_structure = pdb_inp.xray_structure_simple()
  xray_structure = xray_structure.set_b_iso(value=10)
  f_obs = abs(xray_structure.structure_factors(d_min=2.0).f_calc())
  f_obs.set_sigmas(sigmas = flex.double(f_obs.data().size(), 0.0))
  reflections_per_bin = min(f_obs.data().size(), reflections_per_bin)
  f_obs.setup_binner(reflections_per_bin = reflections_per_bin)
  print(f_obs.binner().n_bins_used())
  ncs_pairs = tncs.groups(
    pdb_hierarchy    = pdb_inp.construct_hierarchy(),
    crystal_symmetry = f_obs.crystal_symmetry()).ncs_pairs
  tncs.initialize_rho_mn(
    ncs_pairs       = ncs_pairs,
    d_spacings_data = f_obs.d_spacings().data(),
    binner          = f_obs.binner())
  return f_obs, ncs_pairs

def exercise_00(file_name, reflections_per_bin=150):
  """
  Finite differences test for radii.
  """
  f_obs, ncs_pairs = f_obs_and_tncs_pairs_from_pdb(file_name = file_name,
    reflections_per_bin = reflections_per_bin)
  tncs.finite_differences_grad_radius(ncs_pairs=ncs_pairs,
    f_obs=f_obs, reflections_per_bin=reflections_per_bin, tolerance=1.)

def exercise_01(file_name, reflections_per_bin=5000):
  """
  Finite differences test for rho_mn.
  """
  f_obs, ncs_pairs = f_obs_and_tncs_pairs_from_pdb(file_name = file_name,
    reflections_per_bin = reflections_per_bin)
  tncs.finite_differences_rho_mn(ncs_pairs=ncs_pairs,
    f_obs=f_obs, reflections_per_bin=reflections_per_bin, tolerance=1.)

def run():
  for file_name in ["model_2.pdb", "model_4.pdb"]:
    print(file_name)
    exercise_00(file_name = file_name)
    exercise_01(file_name = file_name)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.3f"%(time.time()-t0))
  print("OK")
