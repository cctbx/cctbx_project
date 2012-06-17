
import libtbx.load_env
from libtbx.utils import null_out
import os

def exercise () :
  if (not libtbx.env.has_module("phenix_regression")) :
    print "phenix_regression not configured, skipping test"
    return False
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1x9t.pdb.gz",
    test=os.path.isfile)
  from iotbx.command_line import pdb_expand_ncs
  from iotbx.file_reader import any_file
  pdb_out = pdb_expand_ncs.run(args=[pdb_file], out=null_out())
  assert (os.path.isfile(pdb_out))
  pdb_in = any_file(pdb_out, force_type="pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  chains = hierarchy.models()[0].chains()
  assert (len(chains) == 180)

if (__name__ == "__main__") :
  if (exercise()) : print "OK"
