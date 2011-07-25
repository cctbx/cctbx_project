
from mmtbx.polygon import gui as polygon_gui
import cPickle

class test_renderer (polygon_gui.renderer) :
  def draw_box (self, out, points, color) :
    pass

  def draw_solid_line (self, out, start, end, color) :
    pass

  def draw_dashed_line (self, out, start, end, color) :
    pass

  def draw_labels (self, out, label, min, max, value, pos, angle) :
    pass

def exercise () :
  stats = {
    "r_work" : 0.25,
    "r_free" : 0.28,
    "adp_mean_all" : 20.0,
    "bond_rmsd" : 0.02,
    "angle_rmsd" : 1.8,
    "clashscore" : 20.0
  }
  data = polygon_gui.get_basic_histogram_data(d_min=2.5)
  s = cPickle.dumps(data)
  histograms = polygon_gui.convert_histogram_data(data)
  assert (len(histograms) == 6)
  for stat_key, histogram in histograms.iteritems() :
    bins = [ n for n in histogram.slots() ]
    #print "%-16s : %s" % (stat_key, " ".join([ "%5d" % n for n in bins ]))
  renderer = test_renderer(
    histogram_data=data,
    structure_stats=stats)
  renderer.draw(out=None)
  print "OK"

if __name__ == "__main__" :
  exercise()
