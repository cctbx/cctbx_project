from __future__ import absolute_import, division, print_function
import scitbx.array_family.flex # import dependency
import time

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("gltbx_viewer_utils_ext")
from gltbx_viewer_utils_ext import *

def read_pixels_to_str(x, y, width, height):
  ext = bp.import_ext("gltbx_gl_ext")
  from gltbx_gl_ext import glPixelStorei, glReadPixels, \
    GL_PACK_ALIGNMENT, GL_RGB, GL_UNSIGNED_BYTE
  # from gltbx.gl import glPixelStorei, glReadPixels, \
  #   GL_PACK_ALIGNMENT, GL_RGB, GL_UNSIGNED_BYTE
  glPixelStorei(GL_PACK_ALIGNMENT, 1)
  pixels = []
  glReadPixels(
    x=x, y=y, width=width, height=height,
    format=GL_RGB, type=GL_UNSIGNED_BYTE,
    pixels=pixels)
  return pixels[0]

def read_pixels_to_pil_image(x, y, width, height):
  try:
    import PIL.Image
  except ImportError:
    return None
  mode = "RGB"
  size = (width, height)
  data = read_pixels_to_str(x=x, y=y, width=width, height=height)
  decoder_name = "raw"
  raw_mode = "RGB"
  stride = 0
  orientation = -1
  return PIL.Image.frombytes(
    mode, size, data, decoder_name, raw_mode, stride, orientation)

class fps_monitor(object):
  def __init__(self):
    self._t_start = time.time()
    self._n = 0

  def update(self):
    self._n += 1
    if (self._n % 10 == 0):
      t_curr = time.time()
      t_elapsed = t_curr - self._t_start
      self._t_start = t_curr
      print("%.2f fps" % (10 / t_elapsed))
      self._n = 0
