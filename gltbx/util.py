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

  _shared_state = None

  def __init__(self):
    from gltbx import gl
    import re
    if self.__class__._shared_state is None:
      vers_pat = re.compile("^((\d+)\.(\d+))(?:\.(\d+))? (.*)$")
      m = vers_pat.search(gl.glGetString(gl.GL_VERSION))
      state = dict(zip(
        ["principal", "major_number", "minor_number", "release_number",
         "vendor_info"], m.groups()))
      state['principal'] = float(state['principal'])
      self.__class__._shared_state = state
    self.__dict__ = self.__class__._shared_state


class normalizing_normals(object):

  _shared_state = {}

  GL_RESCALE_NORMAL = 0x803A

  def __init__(self):
    self.__dict__ = self._shared_state

  def enable_rescale(self, shall_rescale):
    from gltbx.gl import *
    if shall_rescale:
      func = glEnable
    else:
      func = glDisable
    if version().principal >= 1.2:
      func(self.GL_RESCALE_NORMAL)
    else:
      func(GL_NORMALIZE)

  def enable_normalize(self, shall_normalize):
    from gltbx.gl import *
    if shall_normalize:
      glEnable(GL_NORMALIZE)
    else:
      glDisable(GL_NORMALIZE)

  def is_rescale_enabled(self):
    from gltbx.gl import *
    return glIsEnabled(self.GL_RESCALE_NORMAL)

  def is_normalize_enabled(self):
    from gltbx.gl import *
    return glIsEnabled(GL_NORMALIZE)
