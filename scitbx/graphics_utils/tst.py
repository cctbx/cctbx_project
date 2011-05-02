
from scitbx import graphics_utils
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

def exercise () :
  g = graphics_utils.make_rainbow_gradient(nbins=21)
  assert approx_equal(g[0], (0.0, 0.0, 1.0))
  assert approx_equal(g[10], (0.0, 1.0, 0.0))
  assert approx_equal(g[-1], (1.0, 0.0, 0.0))
  sel = flex.bool([True] * 10)
  sel[5] = False
  r = graphics_utils.color_rainbow(
    selection=sel,
    color_all=False)
  assert approx_equal(r[0], (0.0,0.0,1.0))
  assert approx_equal(r[-1], (1.0,0.0,0.0))
  assert approx_equal(r[4], (0.0,1.0,0.0))
  r2 = graphics_utils.color_rainbow(
    selection=sel,
    color_all=True)
  assert approx_equal(r2[0], (0.0,0.0,1.0))
  assert approx_equal(r2[-1], (1.0,0.0,0.0))
  assert approx_equal(r2[6], (2/3.,1.0,0.0))
  b = flex.double([4.0,5.2,1.7,6.9,9.5,24.3])
  c = graphics_utils.color_by_property(
    properties=b,
    selection=flex.bool(b.size(), True))
  assert approx_equal(c[2], (0.0,0.0,1.0))
  c2 = graphics_utils.scale_selected_colors(
    input_colors=c,
    selection=flex.bool(c.size(), True),
    scale=0.9)
  assert approx_equal(c2[2], (0.0,0.0,0.9))

if (__name__ == "__main__") :
  exercise()
  print "OK"
