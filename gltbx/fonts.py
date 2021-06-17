from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("gltbx_fonts_ext")
from gltbx_fonts_ext import *

ucs_bitmap_8x13 = ucs_bitmap(short_name="8x13")
ucs_bitmap_9x15 = ucs_bitmap(short_name="9x15")
ucs_bitmap_10x20 = ucs_bitmap(short_name="10x20")

@bp.inject_into(ucs_bitmap)
class _():

  def render_text(self, position, text, relative_line_spacing=1.0,
                  use_3d_position=False):
    from gltbx.gl import glRasterPos2f, glRasterPos3f, glBitmap
    if use_3d_position: glRasterPos = glRasterPos3f
    else: glRasterPos = glRasterPos2f
    line_spacing = round(self.height() * relative_line_spacing)
    for i,string in enumerate(text.splitlines()):
      glRasterPos(*position)
      glBitmap(0, 0, 0.0, 0.0, 0.0, -i*line_spacing, b"")
      self.render_string(string=string.expandtabs())
