from iotbx import pdb
from cctbx.array_family import flex
import time, math,os,sys
from libtbx.test_utils import approx_equal
from libtbx import introspection
import libtbx.load_env


def calculate_fobs(resolution   = 1.0,
                   sf_algorithm = "direct"):
  pdb_file = libtbx.env.find_in_repositories(
                              relative_path="regression/pdb/enk_gbr.pdb", test=os.path.isfile)
  xray_structure = pdb.as_xray_structure(pdb_file)
  xray_structure.scattering_type_registry(table = "wk1995")
  f_calc = xray_structure.structure_factors(
                                        d_min          = resolution,
                                        anomalous_flag = False,
                                        algorithm      = sf_algorithm).f_calc()
  f_calc = abs(f_calc.structure_factors_from_scatterers(
                                     xray_structure = xray_structure).f_calc())
  r_free_flags = f_calc.generate_r_free_flags(fraction = 0.01,
                                              max_free = 200000)

  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "FOBS")

  mtz_dataset.add_miller_array(miller_array      = r_free_flags,
                               column_root_label = "TEST")

  sigmas = r_free_flags.array(data = flex.double(r_free_flags.data().size(),1))
  mtz_dataset.add_miller_array(miller_array      = sigmas,
                               column_root_label = "SIGMA")

  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "enk_gbr.mtz")

def exercise_1(hkl = "enk_gbr.mtz"):
  pdb = libtbx.env.find_in_repositories(
                              relative_path="regression/pdb/enk_gbr_e.pdb", test=os.path.isfile)
  opt0= "main.number_of_macro_cycles=3 main.refine_site=false"
  opt1= "main.target=ls main.refine_adp_individual=false"
  opt2= "output.write_map_file=false output.write_geo_file=true output.write_def_file=false"
  opt3= "output.write_eff_file=false output.write_refined_mtz_file=false"
  opt4= "main.group_bfactor_refinement=true group_b_iso.one_residue_one_group=true"
  opt5= "main.sf_algorithm=direct scattering_table=wk1995 --overwrite"
  opt6= "refinement.input.xray_data.labels=FOBS main.bulk_solvent_and_scale=false  output.prefix=ref1"
  cmd = " ".join(["phenix.refine", pdb, hkl, opt0, opt1, opt2, opt3, opt4, opt5, opt6])
  os.system(cmd)

def exercise_2(pdb = "enk_gbr_e.pdb", hkl = "enk_gbr.mtz"):
  pdb = libtbx.env.find_in_repositories(
                              relative_path="regression/pdb/enk_gbr_e.pdb", test=os.path.isfile)
  opt0= "main.number_of_macro_cycles=3 main.refine_site=false"
  opt1= "main.target=ls main.refine_adp_individual=false"
  opt2= "output.write_map_file=false output.write_geo_file=true output.write_def_file=false"
  opt3= "output.write_eff_file=false output.write_refined_mtz_file=false"
  opt4= "main.group_bfactor_refinement=true group_b_iso.one_residue_one_group=false"
  opt5= "main.sf_algorithm=direct scattering_table=wk1995 --overwrite"
  opt6= "refinement.input.xray_data.labels=FOBS main.bulk_solvent_and_scale=false"
  opt7= "group_b_iso.selection="+""""chain A" """ +" group_b_iso.selection="+""""chain B" """
  opt8= "group_b_iso.selection="+""""chain C" """+" group_b_iso.selection="+""""chain D" """+" output.prefix=ref2"
  cmd = " ".join(["phenix.refine", pdb, hkl, opt0, opt1, opt2, opt3, opt4, opt5, opt6, opt7, opt8])
  os.system(cmd)

def check_result():
  for st in open("ref1_001.pdb","r").read().splitlines():
    if(st.count("REMARK Final: r_work =")==1):
       st = st.split()
       r1 = float(st[4])
  for st in open("ref2_001.pdb","r").read().splitlines():
    if(st.count("REMARK Final: r_work =")==1):
       st = st.split()
       r2 = float(st[4])
  assert approx_equal(r1,r2, 0.001)
  assert r1 < 0.01 and r2 < 0.01


def run():
  calculate_fobs()
  exercise_1()
  exercise_2()
  check_result()
  os.system("rm -rf *.geo *.mtz *.log ref1* ref2*")

if (__name__ == "__main__"):
  run()
