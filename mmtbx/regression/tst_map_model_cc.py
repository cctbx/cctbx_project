from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
from libtbx import easy_run
import time
import mmtbx.maps.map_model_cc
from libtbx import easy_pickle

pdb_str = """
CRYST1   14.755   15.940   19.523  90.00  90.00  90.00 P 1
ATOM      1  N   ALA L 139       5.232   6.337   5.079  1.00 25.00           N
ATOM      2  CA  ALA L 139       5.342   5.476   6.250  1.00 25.00           C
ATOM      3  C   ALA L 139       6.027   6.204   7.402  1.00 25.00           C
ATOM      4  O   ALA L 139       7.033   6.886   7.206  1.00 25.00           O
ATOM      5  N   TRP L 140       5.477   6.054   8.602  1.00 25.00           N
ATOM      6  CA  TRP L 140       6.034   6.696   9.786  1.00 25.00           C
ATOM      7  C   TRP L 140       6.303   5.675  10.888  1.00 25.00           C
ATOM      8  O   TRP L 140       5.382   5.017  11.372  1.00 25.00           O
ATOM      9  CB  TRP L 140       5.092   7.787  10.298  1.00 30.00           C
ATOM     10  CG  TRP L 140       4.833   8.873   9.299  1.00 30.00           C
ATOM     11  CD1 TRP L 140       3.836   8.910   8.369  1.00 30.00           C
ATOM     12  CD2 TRP L 140       5.584  10.082   9.131  1.00 30.00           C
ATOM     13  NE1 TRP L 140       3.919  10.066   7.631  1.00 30.00           N
ATOM     14  CE2 TRP L 140       4.984  10.803   8.079  1.00 30.00           C
ATOM     15  CE3 TRP L 140       6.705  10.624   9.766  1.00 30.00           C
ATOM     16  CZ2 TRP L 140       5.467  12.037   7.649  1.00 30.00           C
ATOM     17  CZ3 TRP L 140       7.183  11.849   9.338  1.00 30.00           C
ATOM     18  CH2 TRP L 140       6.565  12.542   8.290  1.00 30.00           C
ATOM     19  N   ALA L 141       7.567   5.544  11.283  1.00 25.00           N
ATOM     20  CA  ALA L 141       8.649   6.331  10.700  1.00 25.00           C
ATOM     21  C   ALA L 141       9.725   5.427  10.110  1.00 25.00           C
ATOM     22  O   ALA L 141       9.623   4.991   8.963  1.00 25.00           O
TER
END
"""

def write_ccp4_map(map_data, cs, file_name):
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
    file_name=file_name,
    unit_cell=cs.unit_cell(),
    space_group=cs.space_group(),
    map_data=map_data,
    labels=flex.std_string([""]))

def run(prefix="tst_map_model_cc"):
  """
  Make sure it works with map having origin != (0,0,0)
  """
  # original (zero-origin) map and model
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  pdb_inp.write_pdb_file(file_name="%s_%s"%(prefix,"orig.pdb"))
  ph = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  fc = xrs.structure_factors(d_min=1.5).f_calc()
  fft_map = fc.fft_map(resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  assert map_data.all() == (40, 45, 54)
  assert map_data.origin() == (0, 0, 0)
  assert map_data.focus() == (40, 45, 54)
  write_ccp4_map(map_data=map_data, cs=xrs.crystal_symmetry(),
    file_name="%s_%s"%(prefix,"orig.ccp4"))
  # shift origin of the map
  g = flex.grid((-20,-25,-27), (20,20,27))
  map_data.reshape(g)
  assert map_data.all() == (40, 45, 54)
  assert map_data.origin() == (-20, -25, -27)
  assert map_data.focus() == (20, 20, 27)
  write_ccp4_map(map_data=map_data, cs=xrs.crystal_symmetry(),
    file_name="%s_%s"%(prefix,"shifted.ccp4"))
  # apply same shift to the model
  a,b,c = xrs.crystal_symmetry().unit_cell().parameters()[:3]
  N = map_data.all()
  O=map_data.origin()
  sites_cart = ph.atoms().extract_xyz()
  sx,sy,sz = a/N[0]*O[0], b/N[1]*O[1], c/N[2]*O[2]
  sites_cart_shifted = sites_cart+\
    flex.vec3_double(sites_cart.size(), [sx,sy,sz])
  ph.atoms().set_xyz(sites_cart_shifted)
  ph.write_pdb_file(file_name="%s_%s"%(prefix,"shifted.pdb"))
  # run phenix.real_space_refine
  checked = 0
  cmd = " ".join([
    "phenix.map_model_cc",
    "force",
    "%s_shifted.pdb"%prefix,
    "%s_shifted.ccp4"%prefix,
    "resolution=1.5",
    "> %s.zlog"%prefix
  ])
  print(cmd)
  assert not easy_run.call(cmd)
  # check results
  fo = open("%s.zlog"%prefix,"r")
  for l in fo.readlines():
    if(l.startswith("  CC_mask  :")):
      cc = float(l.split()[2])
      assert cc>0.98
      checked+=1
  fo.close()
  assert checked==1
  # Exercise corresponding library function
  params = mmtbx.maps.map_model_cc.master_params().extract()
  params.map_model_cc.resolution=1.5
  task_obj = mmtbx.maps.map_model_cc.map_model_cc(
    map_data         = map_data,
    pdb_hierarchy    = ph,
    crystal_symmetry = xrs.crystal_symmetry(),
    params           = params.map_model_cc)
  task_obj.validate()
  task_obj.run()
  result = task_obj.get_results()
  assert result.cc_mask  >0.98
  assert result.cc_peaks >0.98
  assert result.cc_volume>0.98
  assert result.cc_side_chain.cc>0.98
  assert result.cc_main_chain.cc>0.98
  easy_pickle.dump(
    file_name = "map_model_cc_test.pkl",
    obj       = result)
  easy_pickle.load("map_model_cc_test.pkl")

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.3f"%(time.time()-t0))
  print("OK")
