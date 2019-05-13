from __future__ import division
from __future__ import print_function
from mmtbx.wwpdb import standard_geometry_cif
import libtbx.load_env
import sys, os

def exercise(args):
  assert len(args) == 0
  std_geo_cif = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wwpdb/standard_geometry.cif",
    test=os.path.isfile,
    optional=False)
  std_comps = standard_geometry_cif.process_chem_comps(file_name=std_geo_cif)
  chem_comp = std_comps["LEU"]
  assert chem_comp.comp_id == "LEU"
  assert len(chem_comp.atoms) == 24
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
