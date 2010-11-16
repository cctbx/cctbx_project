from iotbx import pdb
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx import easy_run
import libtbx.load_env
import random
import sys, os

def calculate_fobs(resolution   = 1.0,
                   algorithm = "direct"):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/enk_gor.pdb", test=os.path.isfile)
  xray_structure = pdb.input(file_name=pdb_file).xray_structure_simple()
  xray_structure.scattering_type_registry(table = "wk1995")
  f_calc = xray_structure.structure_factors(
                                        d_min          = resolution,
                                        anomalous_flag = False,
                                        cos_sin_table  = False,
                                        algorithm      = algorithm).f_calc()
  f_calc = abs(f_calc.structure_factors_from_scatterers(
                                     xray_structure = xray_structure).f_calc())
  r_free_flags = f_calc.generate_r_free_flags(fraction = 0.01,
                                              max_free = 200000)

  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "f_obs")

  mtz_dataset.add_miller_array(miller_array      = r_free_flags,
                               column_root_label = "TEST")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "enk_gor.mtz")

def exercise_1(hkl = "enk_gor.mtz"):
  par_str = """
refinement {
  refine {
    strategy = individual_sites rigid_body individual_adp group_adp tls \
               *occupancies group_anomalous
    occupancies {
      constrained_group {
        selection = resseq 1
      }
      constrained_group {
        selection = resseq 2
      }
      constrained_group {
        selection = resseq 3
      }
      constrained_group {
        selection = resseq 4
      }
      constrained_group {
        selection = resseq 5
      }
    }
  }
  main {
    bulk_solvent_and_scale = False
    target = *ml mlhl ml_sad ls
    scattering_table = *wk1995 it1992 n_gaussian neutron
    occupancy_max = 100
    occupancy_min = -100
    fake_f_obs = true
  }
  fake_f_obs {
    scattering_table = *wk1995 it1992 n_gaussian neutron
    structure_factors_accuracy {
      algorithm = fft *direct
      cos_sin_table = false
    }
  }
  modify_start_model {
    occupancies {
      randomize = true
    }
  }
  group_occupancy {
    run_finite_differences_test = True
  }
  ls_target_names {
    target_name = ls_wunit_k1 ls_wunit_k2 *ls_wunit_kunit ls_wunit_k1_fixed \
                  ls_wunit_k1ask3_fixed ls_wexp_k1 ls_wexp_k2 ls_wexp_kunit \
                  ls_wff_k1 ls_wff_k2 ls_wff_kunit ls_wff_k1_fixed \
                  ls_wff_k1ask3_fixed lsm_kunit lsm_k1 lsm_k2 lsm_k1_fixed \
                  lsm_k1ask3_fixed
  }
  structure_factors_and_gradients_accuracy {
    algorithm = fft *direct
    cos_sin_table = false
  }
}
"""
  pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/enk_gor.pdb", test=os.path.isfile)
  par_file = open("tst_group_occ_exercise_1.params", "w").write(par_str)
  cmd = " ".join(["phenix.refine", pdb, hkl, "tst_group_occ_exercise_1.params",
                  "--overwrite", "output.prefix=tst_group_occ_exercise_1",
                  "wu=0"])
  easy_run.call(cmd)
  r_work_final, r_free_final, r_work_start, r_free_start = [None,]*4
  for line in open("tst_group_occ_exercise_1_001.pdb","r").readlines():
    if(line.startswith("REMARK Start: r_work =")):
      line = line.split()
      r_work_start = float(line[4])
      r_free_start = float(line[7])
    elif(line.startswith("REMARK Final: r_work =")):
      line = line.split()
      r_work_final = float(line[4])
      r_free_final = float(line[7])
  assert approx_equal(r_work_final, 0.0)
  assert approx_equal(r_free_final, 0.0)
  assert r_work_start > 0.5
  assert r_free_start > 0.5

def run(args):
  random_seed = None
  forever = False
  for arg in args:
    if (arg.startswith("--random_seed=")):
      random_seed = int(arg.split("=", 1)[1])
    if (arg.startswith("--forever")):
      forever = True
  if forever: assert random_seed is None
  if random_seed is not None: assert forever == False
  while 1:
     if (random_seed is None):
       random_seed = flex.get_random_seed()
     print "random_seed:", random_seed
     sys.stdout.flush()
     random.seed(random_seed)
     flex.set_random_seed(value=random_seed)
     calculate_fobs()
     exercise_1()
     sys.stdout.flush()
     if(not forever):
       print "random_seed last used:", random_seed
       break

if (__name__ == "__main__"):
  run(sys.argv[1:])
