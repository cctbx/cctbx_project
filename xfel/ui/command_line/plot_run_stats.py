from __future__ import absolute_import, division, print_function

'''
Author      : Young, I.D.
Created     : 07/14/2016
Last Changed: 07/14/2016
Description : XFEL UI plot real-time run stats
'''

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui import master_phil_str, db_phil_str
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.stats import HitrateStats
from xfel.ui.components.run_stats_plotter import plot_multirun_stats
import sys

phil_str = """
  run = None
    .type = str
    .multiple = True
  trial = None
    .type = int
  rungroup = None
    .type = int
  d_min = 2.0
    .type = float
    .help = High resolution bin for the I/sig(I) plot and per-run statistics.
  n_strong_cutoff = 16
    .type = int
    .help = Number of strong spots to consider an image a hit.
  i_sigi_cutoff = 1
    .type = float
    .help = Avg. I/sig(I) in a bin to reach the cutoff for producing a spot (low or high res) in the third plot.
  run_tags = None
    .type = str
    .multiple = True
    .help = Tags to be applied as labels to the runs.
  minimalist = False
    .type = bool
    .help = Generate final plot without run tags, per-run text summaries or vertical lines between runs.
  compress_runs = True
    .type = bool
    .help = When plotting multiple runs, adjust timestamps so there is no blank space between them.
    .help = This mode is not compatible with fetching events from timestamps.
  title = None
    .type = str
    .help = Plot title.
"""
phil_scope = parse(phil_str + master_phil_str + db_phil_str)

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception as e:
      raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()

  app = xfel_db_application(params)
  runs = []
  all_results = []
  if params.rungroup is None:
    assert len(params.run) == 0
    trial = app.get_trial(trial_number = params.trial)
    for rungroup in trial.rungroups:
      for run in rungroup.runs:
        stats = HitrateStats(app, run.run, trial.trial, rungroup.id, params.d_min)()
        if len(stats[0]) > 0:
          runs.append(run.run)
          all_results.append(stats)
  else:
    for run_no in params.run:
      runs.append(run_no)
      all_results.append(HitrateStats(app, run_no, params.trial, params.rungroup, params.d_min)())
  plot_multirun_stats(all_results, runs, params.d_min, n_strong_cutoff=params.n_strong_cutoff, \
    i_sigi_cutoff=params.i_sigi_cutoff, run_tags=params.run_tags, title=params.title, \
    minimalist=params.minimalist, interactive=True, compress_runs=params.compress_runs)

if __name__ == "__main__":
  run(sys.argv[1:])
