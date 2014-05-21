#(jEdit options) :folding=explicit:collapseFolds=1:

from __future__ import division
from libtbx.utils import null_out
import libtbx.load_env
import sys, os

#{{{ exercise_rna_validate
def exercise_rna_validate():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb2goz_refmac_tls.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_regression(): input pdb (pdb2goz_refmac_tls.ent) not available"
    return
  from mmtbx.command_line import rna_validate
  rv = rna_validate.run(args=[regression_pdb], out=null_out())
  assert len(rv.puckers.results) == 2
  assert len(rv.bonds.results) == 2
  assert len(rv.angles.results) == 14
  assert len(rv.suites.results) == 4
  from iotbx.file_reader import any_file
  pdb_in = any_file(regression_pdb)
  import mmtbx.validation.rna_validate
  result = mmtbx.validation.rna_validate.rna_validation(
    pdb_hierarchy=pdb_in.file_object.construct_hierarchy(),
    geometry_restraints_manager=None,
    params=None)

def run():
  verbose = "--verbose" in sys.argv[1:]
  if (not libtbx.env.has_module(name="suitename")):
    print \
      "Skipping exercise_rna_validate():" \
      " phenix not available"
  else:
    exercise_rna_validate()
    print "OK"

if (__name__ == "__main__"):
  run()
