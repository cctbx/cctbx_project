from __future__ import absolute_import, division, print_function
try:
  import OpenGL  # implicit import
  from OpenGL.GLU import *
except ImportError:
  import boost_adaptbx.boost.python as bp
  ext = bp.import_ext("gltbx_glu_ext")
  from gltbx_glu_ext import *
