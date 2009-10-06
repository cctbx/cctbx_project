from mmtbx.rotamer.sidechain_angles import SidechainAngles
from mmtbx.rotamer import rotamer_eval
from phenix.validation.rotalyze import rotalyze
from iotbx import pdb
import libtbx.load_env
import sys, os

def run():
  #load PDB file to be tested
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb13gs.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_regression(): input pdb (pdb1jxt.ent) not available"
    return
  
  #build PDB hierarchy to iterate through
  pdb_io = pdb.input(file_name=regression_pdb)
  hierarchy = pdb_io.construct_hierarchy()
  
  #initialize classes needed for rotamer check
  sa = SidechainAngles(False)
  r = rotamer_eval.RotamerEval()
  rot = rotalyze()
  
  #interate through hierarchy and evaluate each sidechain 
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        all_dict = rot.construct_complete_sidechain(rg)
        is_outlier, value = rot.evaluate_residue(rg, sa, r, all_dict)

if __name__ == "__main__":
  run()