from __future__ import absolute_import, division, print_function
import scitbx.array_family.flex # import dependency
import scitbx.matrix
from libtbx import easy_run
import boost_adaptbx.boost.python as bp
from six.moves import zip
ext = bp.import_ext("gltbx_util_ext")
from gltbx_util_ext import *
import re
import sys

def handle_error():
  '''
  Windows will sometimes throw extra errors that can be ignored
  This function will ignore errors with "invalid" in the message (e.g.
  GL_INVALID_ENUM, GL_INVALID_OPERATION) on Windows.
  '''
  try:
    ext.handle_error()
  except RuntimeError as e:
    if (sys.platform == 'win32') and ('invalid' in repr(e)):
      pass
    else:
      raise

def show_versions():
  from gltbx import gl
  from gltbx import glu
  if (hasattr(gl, "GL_VENDOR")):
    print("GL_VENDOR:", gl.glGetString(gl.GL_VENDOR))
  if (hasattr(gl, "GL_RENDERER")):
    print("GL_RENDERER:", gl.glGetString(gl.GL_RENDERER))
  if (hasattr(gl, "GL_VERSION")):
    print("GL_VERSION:", gl.glGetString(gl.GL_VERSION))
  if (hasattr(gl, "GL_EXTENSIONS")):
    print("GL_EXTENSIONS:", gl.glGetString(gl.GL_EXTENSIONS))
  if (hasattr(glu, "GLU_VERSION")):
    print("GLU_VERSION:", glu.gluGetString(glu.GLU_VERSION))
  if (hasattr(glu, "GLU_EXTENSIONS")):
    print("GLU_EXTENSIONS:", glu.gluGetString(glu.GLU_EXTENSIONS))

# this is essential for Linux - if the X server does not support GLX,
# attempting to use OpenGL will crash the entire program.  this usually
# only happens with remote display on Windows. . .
def check_glx_availability():
  glxerr = easy_run.fully_buffered("glxinfo -b").stderr_lines
  for line in glxerr :
    if re.search('extension "GLX" missing', line):
      return False
  return True

class version(object):

  _shared_state = {}

  def __init__(self):
    from gltbx import gl
    from libtbx.utils import to_str
    import re
    self.__dict__ = self._shared_state
    if not self._shared_state:
      vers_pat = re.compile(r"^((\d+)\.(\d+))(?:\.(\d+))?(?: (.*))?$")
      m = vers_pat.search(to_str(gl.glGetString(gl.GL_VERSION)))
      self.__dict__.update(dict(zip(
        ["principal", "major_number", "minor_number", "release_number",
         "vendor_info"], m.groups())))
      self.principal = float(self.principal)

class extensions(set):

  def __init__(self):
    from gltbx import gl
    ext = gl.glGetString(gl.GL_EXTENSIONS)
    super(extensions, self).__init__(ext.split())

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

def modelview_matrix_as_rt():
  return scitbx.matrix.rt((
    extract_rotation_from_gl_modelview_matrix(),
    extract_translation_from_gl_modelview_matrix()))

def augment_3x3(m):
  a,b,c,d,e,f,g,h,i = m
  return (a,b,c,0,
          d,e,f,0,
          g,h,i,0,
          0,0,0,1)
