"""Tools for interfacing with Pymol
"""
from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args

class pml_stick(object):

  def __init__(self, begin, end, colors=None, width=None):
    if (colors is None):
      colors = [[1,0,0]]*2
    if (width is None):
      width = 0.1
    adopt_init_args(self, locals())

def pml_write(f, label, sticks):
  print("from pymol.cgo import *", file=f)
  print("from pymol import cmd", file=f)
  print("obj = [", file=f)
  for stick in sticks:
    print("CYLINDER,", "%.6g, %.6g, %.6g," % tuple(stick.begin), end=' ', file=f)
    print("%.6g, %.6g, %.6g," % tuple(stick.end), end=' ', file=f)
    print("%.6g," % stick.width, end=' ', file=f)
    print("%.6g, %.6g, %.6g," % tuple(stick.colors[0]), end=' ', file=f)
    print("%.6g, %.6g, %.6g," % tuple(stick.colors[1]), file=f)
  print("]", file=f)
  print('cmd.load_cgo(obj, "%s")' % label, file=f)
