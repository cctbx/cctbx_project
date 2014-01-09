from __future__ import division
#(jEdit options) :folding=explicit:collapseFolds=1:
from mmtbx.validation.rna_validate import rna_validate
from iotbx import pdb
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
  pdb_io = pdb.input(file_name=regression_pdb)
  rv=rna_validate()
  rv.analyze_pdb(pdb_io=pdb_io)
  assert len(rv.pucker_outliers) == 2
  assert len(rv.bond_outliers) == 2
  assert len(rv.angle_outliers) == 0
  assert len(rv.suite_validation) == 4

def run():
  verbose = "--verbose" in sys.argv[1:]
  if (not libtbx.env.has_module(name="phenix")):
    print \
      "Skipping exercise_rna_validate():" \
      " phenix not available"
  else:
    exercise_rna_validate()
    print "OK"

if (__name__ == "__main__"):
  run()
