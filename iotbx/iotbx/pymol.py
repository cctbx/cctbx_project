from scitbx.python_utils.misc import adopt_init_args

class pml_stick:

  def __init__(self, begin, end, colors=None, width=None):
    if (colors is None):
      colors = [[1,0,0]]*2
    if (width is None):
      width = 0.1
    adopt_init_args(self, locals())

def pml_write(f, label, sticks):
  print >> f, "from pymol.cgo import *"
  print >> f, "from pymol import cmd"
  print >> f, "obj = ["
  for stick in sticks:
    print >> f, "CYLINDER,", "%.6g, %.6g, %.6g," % tuple(stick.begin),
    print >> f, "%.6g, %.6g, %.6g," % tuple(stick.end),
    print >> f, "%.6g," % stick.width,
    print >> f, "%.6g, %.6g, %.6g," % tuple(stick.colors[0]),
    print >> f, "%.6g, %.6g, %.6g," % tuple(stick.colors[1])
  print >> f, "]"
  print >> f, 'cmd.load_cgo(obj, "%s")' % label
