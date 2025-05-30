"""Summarize B-factor statistics from a model"""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb.b_factor_stats

import iotbx.pdb
from cctbx import adptbx
from cctbx.array_family import flex
import libtbx.utils
import sys

def run(args):
  for file_name in args:
    print("File name:", file_name)
    try:
      pdb_inp = iotbx.pdb.input(file_name=file_name)
    except KeyboardInterrupt: raise
    except Exception:
      libtbx.utils.format_exception()
    isotropic_b_factors = flex.double()
    all_eigenvalues = flex.double()
    for atom in pdb_inp.atoms():
      if (atom.uij == (-1,-1,-1,-1,-1,-1)):
        isotropic_b_factors.append(atom.b)
      else:
        all_eigenvalues.extend(flex.double(adptbx.eigenvalues(atom.uij)))
    all_eigenvalues *= adptbx.u_as_b(1)
    print("Number of isotropic atoms:  ", isotropic_b_factors.size())
    print("Number of anisotropic atoms:", all_eigenvalues.size() // 3)
    if (isotropic_b_factors.size() != 0):
      print("Histogram of isotropic B-factors:")
      flex.histogram(data=isotropic_b_factors, n_slots=10).show(
        prefix="  ", format_cutoffs="%7.2f")
    if (all_eigenvalues.size() != 0):
      print("Histogram of eigenvalues of anisotropic B-factors:")
      flex.histogram(data=all_eigenvalues, n_slots=10).show(
        prefix="  ", format_cutoffs="%7.2f")
    print()

if (__name__ == "__main__"):
  run(sys.argv[1:])
