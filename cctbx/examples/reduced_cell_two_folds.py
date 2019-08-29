"""
Enumeration of the 81 2-fold symmetry operations possible for reduced cells.
Show the matrix elements and the axis directions in direct space and
reciprocal space.
"""
from __future__ import absolute_import, division, print_function

from cctbx import sgtbx
from cctbx.array_family import flex

def run():
  for elements in flex.nested_loop([-1]*9,[1+1]*9):
    r = sgtbx.rot_mx(elements)
    if (r.determinant() != 1): continue
    if (not r.multiply(r).is_unit_mx()): continue
    if (r.is_unit_mx()): continue
    print(elements, r.info().ev(), r.transpose().info().ev())
  print("OK")

if (__name__ == "__main__"):
  run()
