import boost.python
ext = boost.python.import_ext("gltbx_fonts_ext")
from gltbx_fonts_ext import *

ucs_bitmap_8x13 = ucs_bitmap(short_name="8x13")
ucs_bitmap_9x15 = ucs_bitmap(short_name="9x15")
ucs_bitmap_10x20 = ucs_bitmap(short_name="10x20")

class _ucs_bitmap(boost.python.injector, ucs_bitmap):

  def render_text(self, position, text, relative_line_spacing=1.0):
    from gltbx.gl import glRasterPos2f, glBitmap
    line_spacing = round(self.height() * relative_line_spacing)
    for i,string in enumerate(text.splitlines()):
      glRasterPos2f(*position)
      glBitmap(0, 0, 0.0, 0.0, 0.0, -i*line_spacing, "")
      self.render_string(string=string.expandtabs())
