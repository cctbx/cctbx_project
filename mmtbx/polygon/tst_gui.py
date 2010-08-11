
from mmtbx.polygon import gui as polygon_gui
import cPickle

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
  print "OK"

if __name__ == "__main__" :
  exercise()
