"""Interpret TLS in terms of rotations and translations"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.tls_analysis

import sys
from libtbx.utils import Sorry
import time
import mmtbx.tls.analysis

legend = """phenix.tls_analysis:
  Given PDB file with TLS records analyze each (T,L,S) triplet and interpret it
  in terms of parameters of elemental motions (rotations and translations).

Citation:
  From deep TLS validation to ensembles of atomic models built from elemental
  motions
  Urzhumtsev,A., Afonine,P.V., Van Benschoten,A.H., Fraser,J.S. & Adams,P.D.
  Acta Cryst. (2015). D71

How to run:
  phenix.tls_analysis model.pdb"""

def run(args):
  t0 = time.time()
  print("-"*79)
  print(legend)
  print("-"*79)
  if(len(args) != 1): raise Sorry("PDB file must be provided.")
  mmtbx.tls.analysis.cmd_driver(pdb_file_name = args[0])
  print("Time: %-10.3f"%(time.time()-t0))

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

