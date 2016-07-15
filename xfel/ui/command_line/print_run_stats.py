from __future__ import division

'''
Author      : Young, I.D.
Created     : 07/14/2016
Last Changed: 07/14/2016
Description : XFEL UI plot real-time run stats
'''

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.stats import HitrateStats
import sys
from xfel.ui.command_line.plot_run_stats import phil_scope

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception, e:
      raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  print "Printing results for trial", params.trial, "using a hit cutoff of", params.hit_cutoff, "reflections"
  print
  print "Run    N Hits   (%) N Indexed   (%) N Frames"

  hit_total = 0
  indexed_total = 0
  overall_total = 0

  app = xfel_db_application(params)
  for run_no in params.run:
    timestamps, n_strong, average_intensity, average_sigma, average_i_sigi = HitrateStats(app, run_no, params.trial, params.rungroup)()
    n_hit = (n_strong >= params.hit_cutoff).count(True)
    n_indexed = (average_i_sigi > 0).count(True)
    n_total = len(timestamps)
    print "% 4d  % 7d % 5.1f   % 7d % 5.1f  % 7d" % (run_no, n_hit, 100*n_hit/n_total, n_indexed, 100*n_indexed/n_total, n_total)

    hit_total += n_hit
    indexed_total += n_indexed
    overall_total += n_total

  if len(params.run) > 1:
    print "-" * 80
    print "Total % 7d % 5.1f   % 7d % 5.1f  % 7d" % (hit_total, 100*hit_total/overall_total, indexed_total, 100*indexed_total/overall_total, overall_total)

if __name__ == "__main__":
  run(sys.argv[1:])
