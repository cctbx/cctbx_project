from __future__ import absolute_import, division, print_function
from gltbx import fonts

def exercise():
  for short_name,full_name,width,height,xorig,yorig in [
      ("8x13", "-Misc-Fixed-Medium-R-Normal--13-120-75-75-C-80-ISO10646-1",
        8,13,0,-2),
      ("9x15", "-Misc-Fixed-Medium-R-Normal--15-140-75-75-C-90-ISO10646-1",
        9,15,0,-3),
      ("10x20", "-Misc-Fixed-Medium-R-Normal--20-200-75-75-C-100-ISO10646-1",
        10,20,0,-4)]:
    bitmap = fonts.ucs_bitmap(short_name=short_name)
    assert bitmap.short_name() == short_name
    assert bitmap.full_name() == full_name
    assert bitmap.width() == width
    assert bitmap.height() == height
    assert bitmap.xorig() == xorig
    assert bitmap.yorig() == yorig
  print("OK")

if (__name__ == "__main__"):
  exercise()
