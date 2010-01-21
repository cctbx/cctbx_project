import libtbx.load_env
import sys, os, math, time
from mmtbx import utils
from libtbx.test_utils import approx_equal, run_command
from libtbx.utils import show_times_at_exit
from iotbx import reflection_file_reader
import iotbx.pdb
from scitbx.array_family import flex
import mmtbx.f_model


def exercise_f_model_option_default(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "phe_e.pdb")
  high_resolution = 3.0
  cmd = " ".join(['phenix.fmodel',
                  '%s'%file_name,
                  'high_resolution=%s'%str(high_resolution)])
  run_command(command=cmd, verbose=verbose)
  reflection_file = reflection_file_reader.any_reflection_file(
    file_name          = os.path.basename(file_name)+".mtz",
    ensure_read_access = False)
  miller_arrays = reflection_file.as_miller_arrays()
  assert len(miller_arrays) == 1
  ma = miller_arrays[0]
  assert "%s"%ma.info() == \
    "%s:FMODEL,PHIFMODEL"%(os.path.basename(file_name)+".mtz")
  xray_structure = iotbx.pdb.input(file_name=file_name).xray_structure_simple()
  xray_structure.scattering_type_registry(table="n_gaussian",d_min =ma.d_min())
  assert ma.is_complex_array()
  assert approx_equal(ma.d_min(), high_resolution, 0.1)
  fmodel = mmtbx.f_model.manager(
    xray_structure = xray_structure,
    r_free_flags   = ma.generate_r_free_flags(fraction = 0.1),
    target_name    = "ls_wunit_k1",
    f_obs          = abs(ma))
  r1 = abs(ma).data().min_max_mean().as_tuple()
  r2 = abs(fmodel.f_model()).data().min_max_mean().as_tuple()
  r3 = abs(fmodel.f_calc()).data().min_max_mean().as_tuple()
  r4 = abs(fmodel.f_obs).data().min_max_mean().as_tuple()
  assert approx_equal(r1, r2, 1.e-5)
  assert approx_equal(r3, r4, 1.e-5)
  assert approx_equal(r1, r4, 1.e-5)

def exercise_f_model_option_custom(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "enk_gbr.pdb")
  high_resolution = 2.0
  format = "cns"
  type = "real"
  label = "Fobs"
  low_resolution = 6.0
  algorithm = "direct"
  par = (0.35,60,3,[1,2,-3,0,0,0])
  par_str = 'k_sol=%s b_sol=%s scale=%s b_cart="%s"'%(par[0],
    par[1], par[2], " ".join([str(i) for i in par[3]]).strip())
  for type in ["real", "complex"]:
    for table in ["wk1995", "neutron"]:
      cmd = " ".join(
        ['phenix.fmodel',
         '%s'%file_name,
         'high_resolution=%s'%str(high_resolution),
         'format=%s'%format,
         'type=%s'%type,
         'label=%s'%label,
         'low_resolution=%s'%str(low_resolution),
         'algorithm=%s'%algorithm,
         'scattering_table=%s'%table,
         '%s'%par_str])
      run_command(command=cmd, verbose=verbose)
      xray_structure=iotbx.pdb.input(file_name=file_name).xray_structure_simple()
      reflection_file = reflection_file_reader.any_reflection_file(
        file_name          = os.path.basename(file_name)+".hkl",
        ensure_read_access = False)
      miller_arrays = reflection_file.as_miller_arrays(crystal_symmetry =
        xray_structure.crystal_symmetry())
      assert len(miller_arrays) == 1
      ma = miller_arrays[0]
      if(table == "neutron"):
        xray_structure.switch_to_neutron_scattering_dictionary()
      else:
        xray_structure.scattering_type_registry(table = table,d_min = ma.d_min())
      if(type == "real"): assert ma.is_real_array()
      if(type == "complex"): assert ma.is_complex_array()
      d_max, d_min = ma.d_max_min()
      assert approx_equal(d_min, high_resolution, 0.1)
      assert approx_equal(d_max, low_resolution, 0.1)
      assert "%s"%ma.info() == "%s:FOBS"%(os.path.basename(file_name)+".hkl")
      sf_calc_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
      sf_calc_params.algorithm = algorithm
      fmodel = mmtbx.f_model.manager(
        xray_structure = xray_structure,
        sf_and_grads_accuracy_params = sf_calc_params,
        r_free_flags   = ma.generate_r_free_flags(fraction = 0.1),
        target_name    = "ml",
        f_obs          = abs(ma))
      fmodel.update_solvent_and_scale(verbose = -1)
      # tolerances MUST be small, otherwise ring a bell
      assert approx_equal(fmodel.r_work(),        0, 1.e-4)
      assert approx_equal(fmodel.r_free(),        0, 1.e-4)
      assert approx_equal(fmodel.k_sol(),    par[0], 1.e-2)
      assert approx_equal(fmodel.b_sol(),    par[1], 1.e-2)
      assert approx_equal(fmodel.scale_k1(), par[2], 1.e-4)
      assert approx_equal(fmodel.b_cart(),   par[3], 1.e-3)

def exercise_01(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "t.pdb")
  high_resolution = 3.0
  cmd = " ".join(['phenix.fmodel',
                  '%s'%file_name,
                  'high_resolution=%s'%str(high_resolution)])
  result = run_command(command=cmd, verbose=verbose, sorry_expected = True,
    join_stdout_stderr = True)
  sorry_found = False
  for line in result.stdout_lines:
    if(line == "Sorry: CRYST1 record in input PDB file is incomplete or missing."):
      sorry_found = True
  assert sorry_found == True


def exercise(args):
  if ("--show-everything" in args):
    verbose = 2
  elif ("--verbose" in args):
    verbose = 1
  else:
    verbose = 0
  pdb_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb", test=os.path.isdir)
  if (pdb_dir is None):
    print "Skipping exercise(): input files not available"
    return
  eargs = {"pdb_dir": pdb_dir, "verbose": verbose}
  exercise_f_model_option_default(**eargs)
  exercise_f_model_option_custom(**eargs)
  exercise_01(**eargs)
  print "OK"

if (__name__ == "__main__"):
  show_times_at_exit()
  exercise(sys.argv[1:])
