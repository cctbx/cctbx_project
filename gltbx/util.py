import scitbx.array_family.flex

import boost.python
ext = boost.python.import_ext("gltbx_util_ext")
from gltbx_util_ext import *

def show_versions():
  from gltbx import gl
  from gltbx import glu
  if (hasattr(gl, "GL_VENDOR")):
    print "GL_VENDOR:", gl.glGetString(gl.GL_VENDOR)
  if (hasattr(gl, "GL_RENDERER")):
    print "GL_RENDERER:", gl.glGetString(gl.GL_RENDERER)
  if (hasattr(gl, "GL_VERSION")):
    print "GL_VERSION:", gl.glGetString(gl.GL_VERSION)
  if (hasattr(gl, "GL_EXTENSIONS")):
    print "GL_EXTENSIONS:", gl.glGetString(gl.GL_EXTENSIONS)
  if (hasattr(glu, "GLU_VERSION")):
    print "GLU_VERSION:", glu.gluGetString(glu.GLU_VERSION)
  if (hasattr(glu, "GLU_EXTENSIONS")):
    print "GLU_EXTENSIONS:", glu.gluGetString(glu.GLU_EXTENSIONS)
