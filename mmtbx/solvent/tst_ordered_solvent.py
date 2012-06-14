
from libtbx.test_utils import approx_equal
from libtbx import easy_run
import libtbx.load_env
import libtbx.path
import time
import sys, os
from iotbx import pdb
from cctbx.array_family import flex

def evaluate(pdb_file,
             rw_tol,
             rf_tol,
             n_water,
             n_water_tol):
  evaluated = 0
  start_looking = False
  print pdb_file
  cntr = 0
  for line in open(pdb_file).read().splitlines():
    if(line.startswith("REMARK Final: r_work =")):
      rwork = float(line.split()[4])
      rfree = float(line.split()[7])
      assert approx_equal(rwork, 0.0, rw_tol)
      assert approx_equal(rfree, 0.0, rf_tol)
      evaluated += 1
    if(line.startswith("REMARK  stage  number of ordered solvent")):
      start_looking = True
    if(start_looking and line.startswith("REMARK      3_bss:") and cntr > 2):
      cntr += 1
      n_waters_located = int(line.split()[2])
      assert approx_equal(n_waters_located, n_water, n_water_tol)
      evaluated += 1
  assert evaluated >= 0

def exercise_1 (args) :
  pdb_file, kick = args
  rcs = []
  pdb = libtbx.env.find_in_repositories(
        relative_path="phenix_regression/pdb/"+pdb_file, test=os.path.isfile)
  hkl = libtbx.env.find_in_repositories(
        relative_path="phenix_regression/reflection_files/lysozyme.pdb_fcalc_fft_wk1995_resolution_1.8.cv", \
        test=os.path.isfile)
  print pdb,hkl
  base = os.path.basename(os.path.splitext(pdb_file)[0])
  suffix = "_wat_pik_no_h"
  if (kick) :
    suffix += "_kick"
  new_pdb = base + suffix + ".pdb"
  if(pdb_file == "lysozyme_nohoh.pdb"):
    rcs.append(easy_run.call("phenix.pdbtools %s output.file_name=%s"%(
      pdb, new_pdb)))
  else:
    rcs.append(easy_run.call(
      "phenix.pdbtools %s remove='element H' output.file_name=%s"%(
      pdb, new_pdb)))
  output_file_prefix = pdb_file[:-4]
  opt1= "output.prefix="+output_file_prefix+" main.max_number_of_iterations=25 main.scattering_table=wk1995"
  opt2= "main.number_of_macro_cycles=3 main.ordered_solvent=True "
  opt3 = "--overwrite target_weights.wxc_scale=3.0 target_weights.wxu_scale=3.0 target_weights.shake_sites=false "
  opt4 = " ordered_solvent.mode=every_macro_cycle main.nqh_flips=False optimize_mask=False"
  if (not kick) :
    cmd = " ".join(["phenix.refine", new_pdb, hkl, opt1, opt2, opt3, opt4,
      "> /dev/null"])
    print cmd
    print
    sys.stdout.flush()
    rcs.append(easy_run.call(cmd))
    evaluate(pdb_file = output_file_prefix+"_001.pdb",
             rw_tol   = 0.0155,
             rf_tol   = 0.0155,
             n_water  = 186,
             n_water_tol = 0)
  else :
    opt5 = "ordered_solvent.use_kick_maps=True output.serial=2"
    cmd = " ".join(["phenix.refine", new_pdb, hkl, opt1, opt2, opt3, opt4, opt5,
      "> /dev/null"])
    print cmd
    print
    rcs.append(easy_run.call(cmd))
    rcs.append(easy_run.call("rm -rf %s"%new_pdb))
  return (rcs == [ 0 for x in range(len(rcs)) ])

def run(args=()):
  mp = False
  if ("--parallel" in args) or ("--mp" in args) :
    mp = True
  t1 = time.time()
  pdb_files = ["lysozyme_nohoh.pdb", "lysozyme_nohoh_plus6H.pdb"]
  args = []
  for fn in pdb_files :
    args.append((fn, False))
    args.append((fn, True))
  if (mp) :
    try :
      import multiprocessing
    except ImportError :
      mp = False
    else :
      pool = multiprocessing.Pool(processes=len(args))
      results = pool.map(exercise_1, args)
      pool.close()
      assert (not False in results)
  if (not mp) :
    for arg_list in args :
      assert (exercise_1(arg_list))
  f_calcs = []
  for pdb_file in pdb_files:
    xray_structure = pdb.input(
      file_name=pdb_file[:-4]+"_001.pdb").xray_structure_simple()
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
  run(sys.argv[1:])
