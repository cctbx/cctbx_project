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


class version(object):

  _shared_state = {}

  def __init__(self):
    from gltbx import gl
    import re
    self.__dict__ = self._shared_state
    if not self._shared_state:
      vers_pat = re.compile("^((\d+)\.(\d+))(?:\.(\d+))?(?: (.*))?$")
      m = vers_pat.search(gl.glGetString(gl.GL_VERSION))
      self.__dict__.update(dict(zip(
        ["principal", "major_number", "minor_number", "release_number",
         "vendor_info"], m.groups())))
      self.principal = float(self.principal)


class rescale_normals(object):

  _shared_state = {}

  GL_RESCALE_NORMAL = 0x803A

  def __init__(self, fallback_to_normalize=False):
    self.__dict__ = self._shared_state
    if not self._shared_state:
      self.has_rescale_normal = version().principal >= 1.2
    self.fallback_to_normalize = fallback_to_normalize

  def enable(self, flag=True):
    from gltbx import gl
    if self.has_rescale_normal:
      mode = self.GL_RESCALE_NORMAL
    else:
      assert self.fallback_to_normalize,\
             "Rescale normals only available from OpenGL 1.2 onward"
      mode = gl.GL_NORMALIZE
    if flag:
      gl.glEnable(mode)
    else:
      gl.glDisable(mode)

  def is_enabled(self):
    from gltbx import gl
    return self.has_rescale_normal and gl.glIsEnabled(self.GL_RESCALE_NORMAL)
