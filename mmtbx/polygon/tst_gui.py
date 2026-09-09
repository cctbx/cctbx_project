from __future__ import absolute_import, division, print_function

import mmtbx.polygon.output
from six.moves import cPickle as pickle
import six

class test_renderer(mmtbx.polygon.output.renderer):
  def draw_box(self, out, points, color):
    pass

  def draw_solid_line(self, out, start, end, color):
    pass

  def draw_dashed_line(self, out, start, end, color):
    pass

  def draw_labels(self, out, label, min, max, value, pos, angle):
    pass

def exercise():
  stats = {
    "r_work" : 0.25,
    "r_free" : 0.28,
    "adp_mean_all" : 20.0,
    "bond_rmsd" : 0.02,
    "angle_rmsd" : 1.8,
    "clashscore" : 20.0
  }
  data = mmtbx.polygon.output.get_basic_histogram_data(d_min=2.5)
  s = pickle.dumps(data)
  histograms = mmtbx.polygon.output.convert_histogram_data(data)
  assert (len(histograms) == 6)
  for stat_key, histogram in six.iteritems(histograms):
    bins = [ n for n in histogram.slots() ]
    #print "%-16s : %s" % (stat_key, " ".join([ "%5d" % n for n in bins ]))
  renderer = test_renderer(
    histogram_data=data,
    structure_stats=stats)
  renderer.draw(out=None)
  exercise_color_models(renderer)
  print("OK")

def _assert_rgb_ints(color, where):
  # wxPython 4 (Phoenix) wx.Colour/wx.Pen/wx.Brush reject float components,
  # so every colour handed to the GUI must be ints
  assert len(color) == 3, (where, color)
  for c in color:
    assert isinstance(c, int), (where, color)
    assert 0 <= c <= 255, (where, color)

def exercise_color_models(renderer):
  from mmtbx.polygon.output import hsv2rgb
  for h in range(0, 361, 30):
    for s in (0.0, 0.25, 0.5, 1.0):
      for v in (0.0, 0.5, 1.0):
        _assert_rgb_ints(hsv2rgb(h, s, v), ("hsv2rgb", h, s, v))
  for model_name in ("original", "rainbow", "rmb", "blue", "red", "gray"):
    for relative in (True, False):
      renderer.set_color_model(model_name, relative_scaling=relative)
      colors, cutoffs = renderer.get_color_key()
      assert len(colors) > 0
      for color in colors:
        _assert_rgb_ints(color, (model_name, relative, "key"))
      for value in (0.0, 0.3, 0.5, 1.0):
        _assert_rgb_ints(renderer.colors.get_bin_color(value),
          (model_name, "bin", value))
      renderer.draw(out=None)

if __name__ == "__main__" :
  exercise()
