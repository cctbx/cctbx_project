from libtbx.test_utils import approx_equal
import libtbx.load_env
import libtbx.path
import time
import math
import sys, os
from iotbx import pdb
import iotbx.pdb.interpretation
from cctbx.array_family import flex

def remove_refinement_files():
  sys.stdout.flush()
  os.system("rm -rf *_refine_data.mtz *.map  *.eff *.def *.geo *_coeffs.mtz")

def evaluate(pdb_file,
             cycle,
             rw_tol,
             rf_tol,
             ksol,
             ksol_tol,
             bsol,
             bsol_tol,
             u,
             u_tol,
             pher_tol,
             fom_tol,
             alpha_tol,
             beta_tol,
             n_water,
             n_water_tol):
  to_evaluate = str(cycle)+"_adp"
  st = []
  for line in open(pdb_file).read().splitlines():
    if(line.count(to_evaluate)):
      st.append(line.split())

  # R-factors
  print "R-factors= ", float(st[0][2]), float(st[0][3])
  assert approx_equal(float(st[0][2]), 0.0, rw_tol)
  assert approx_equal(float(st[0][3]), 0.0, rf_tol)
  # ksol, bsol, uaniso
  print "ksol, bsol, uaniso= ", float(st[2][2]),float(st[2][3])
  assert approx_equal(float(st[2][2]),ksol,ksol_tol)
  assert approx_equal(float(st[2][3]),bsol,bsol_tol)
  assert approx_equal((float(st[2][4]),float(st[2][5]),float(st[2][6]),\
                float(st[2][7]),float(st[2][8]),float(st[2][9])), u,u_tol)
  # phase error, fom, alpha, beta
  print "phase error, fom, alpha, beta= ", float(st[3][2]),float(st[3][3]),float(st[3][4]),float(st[3][5])
  assert approx_equal(float(st[3][2]),0.1,pher_tol)
  assert approx_equal(float(st[3][3]),1.0,fom_tol)
  assert approx_equal(float(st[3][4]),1.0,alpha_tol)
  assert approx_equal(float(st[3][5]),0.0,beta_tol)
  # Number_of_waters =
  print "Number_of_waters= ", int(st[9][2])
  assert approx_equal(int(st[9][2]),n_water,n_water_tol)


def exercise_1(pdb_file):
  pdb = libtbx.env.find_in_repositories(
        relative_path="regression/pdb/"+pdb_file, test=os.path.isfile)
  hkl = libtbx.env.find_in_repositories(
        relative_path="regression/reflection_files/lysozyme.pdb_fcalc_fft_wk1995_resolution_1.8.cv", \
        test=os.path.isfile)
  print pdb,hkl
  output_file_prefix = pdb_file[:-4]
  opt1= "output.prefix="+output_file_prefix+" main.max_number_of_iterations=25 scattering_table=wk1995"
  opt2= "main.number_of_macro_cycles=3 main.ordered_solvent=True apply_back_trace_of_b_cart=True "
  opt3 = "--overwrite target_weights.wxc_scale=3.0 target_weights.wxu_scale=3.0 target_weights.shake_sites=false "
  cmd = " ".join(["phenix.refine", pdb, hkl, opt1, opt2, opt3])
  print cmd
  print
  sys.stdout.flush()
  os.system(cmd)
  evaluate(pdb_file = output_file_prefix+"_001.pdb",
           cycle    = 3,
           rw_tol   = 0.007,
           rf_tol   = 0.0095,
           ksol     = 0.0,
           ksol_tol = 0.000,
           bsol     = 0.0,
           bsol_tol = 0.000,
           u        = [0.0,0.0,0.0,0.0,0.0,0.0],
           u_tol    = 0.01,
           pher_tol = 1.3,
           fom_tol  = 0.05,
           alpha_tol= 0.05,
           beta_tol = 120.0,
           n_water  = 186,
           n_water_tol = 0)
  remove_refinement_files()

def run():
  t1 = time.time()
  pdb_files = ["lysozyme_nohoh.pdb", "lysozyme_nohoh_plus6H.pdb"]
  f_calcs = []
  for pdb_file in pdb_files:
      exercise_1(pdb_file = pdb_file)
      xray_structure = pdb.interpretation.stage_1(file_name =
                             pdb_file[:-4]+"_001.pdb").extract_xray_structure()
      f_calc = xray_structure.structure_factors(algorithm = "fft",
                                                d_min     = 1.5).f_calc()
      f_calcs.append(f_calc)
  assert len(f_calcs) == 2
  assert approx_equal(f_calcs[0].data(), f_calcs[1].data())
  assert approx_equal(flex.mean(flex.abs(f_calcs[0].data())),
                      flex.mean(flex.abs(f_calcs[1].data())))
  t2 = time.time()
  print "OK. Total time: %.2f" % abs(t2-t1)

if (__name__ == "__main__"):
  run()
  remove_refinement_files()
