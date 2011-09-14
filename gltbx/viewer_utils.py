import scitbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("gltbx_viewer_utils_ext")
from gltbx_viewer_utils_ext import *

def read_pixels_to_str(x, y, width, height):
  from gltbx.gl import glPixelStorei, glReadPixels, \
    GL_PACK_ALIGNMENT, GL_RGB, GL_UNSIGNED_BYTE
  glPixelStorei(GL_PACK_ALIGNMENT, 1)
  pixels = []
  glReadPixels(
    x=0, y=0, width=width, height=height,
    format=GL_RGB, type=GL_UNSIGNED_BYTE,
    pixels=pixels)
  return pixels[0]

def read_pixels_to_pil_image(x, y, width, height):
  try:
    import Image
  except ImportError:
    return None
  mode = "RGB"
  size = (width, height)
  data = read_pixels_to_str(x=x, y=y, width=width, height=height)
  decoder_name = "raw"
  raw_mode = "RGB"
  stride = 0
  orientation = -1
  return Image.fromstring(
    mode, size, data, decoder_name, raw_mode, stride, orientation)
