from __future__ import division

import os,time 
from mmtbx.command_line.place_model_in_map import run

import libtbx.load_env

pdb_file = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/place_model_in_map/std_full.pdb",
    test=os.path.isfile)
map_file= libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/place_model_in_map/poor_map_2mFo-DFc.ccp4",
    test=os.path.isfile)

def tst_01():

  from mmtbx.command_line.place_model_in_map import run
  args=[pdb_file,map_file,"resolution=3"]
  file_name,cc_to_map=run(args)
  print "File: %s, CC to map: %5.2f " %(file_name,cc_to_map)
  assert file_name=="placed_model.pdb"
  assert cc_to_map>=0.75 and cc_to_map <=0.78


if (__name__ == "__main__"):
  t0=time.time()
  tst_01()
  print "Time: %6.4f"%(time.time()-t0)
  print "test 1 OK"
