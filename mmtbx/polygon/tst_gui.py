
from mmtbx.polygon import gui as polygon_gui
import cPickle

def exercise () :
  stats={"r_work_pdb" : 0.25,
         "r_free_pdb" : 0.28,
         "adp_mean" : 20.0,
         "bonds_rmsd" : 0.02,
         "angles_rmsd" : 1.8  }
  data = polygon_gui.get_histogram_data(d_min=2.5)
  s = cPickle.dumps(data)
  histograms = polygon_gui.convert_histogram_data(data)
  assert len(histograms) == 5
  for stat_key, histogram in histograms.iteritems() :
    bins = [ n for n in histogram.slots() ]
    #print "%-16s : %s" % (stat_key, " ".join([ "%5d" % n for n in bins ]))
  print "OK"

if __name__ == "__main__" :
  pass
#  exercise()
