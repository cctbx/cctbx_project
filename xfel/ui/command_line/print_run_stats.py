'''
Author      : Young, I.D.
Created     : 07/14/2016
Last Changed: 07/14/2016
Description : XFEL UI plot real-time run stats
'''
from __future__ import absolute_import, division, print_function

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.stats import HitrateStats
import sys
from xfel.ui.command_line.plot_run_stats import phil_scope
from scitbx.array_family import flex
from six.moves import zip

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception as e:
      raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  print("Printing results for trial", params.trial, "using a hit cutoff of", params.n_strong_cutoff, "reflections")
  print()
  print("                 Run  N Drop Hits   (%)   N Hits   (%) N Indexed   (%) N Lattices N High qual   (%)  %HQR   N Frames")

  drop_hits_total = 0
  hit_total = 0
  indexed_total = 0
  lattices_total = 0
  high_quality_total = 0
  overall_total = 0

  ratio_cutoff = 1

  app = xfel_db_application(params)

  if params.run is None or len(params.run) == 0:
    trial = app.get_trial(trial_number=params.trial)
    runs = []
    run_ids = []
    rungroups = []
    for rg in trial.rungroups:
      for run in rg.runs:
        if run.id in run_ids: continue
        if params.run_tags:
          s1 = set(params.run_tags)
          s2 = set(tag.name for tag in run.tags)
          if not s1.intersection(s2): continue

        runs.append(run.run)
        run_ids.append(run.id)
        rungroups.append(rg.id)
  else:
    runs = params.run
    assert params.rungroup is not None
    rungroups = [params.rungroup] * len(runs)

  for run_no, rungroup_id in sorted(zip(runs, rungroups), key=lambda x: x[0]):
    if params.rungroup and params.rungroup != rungroup_id: continue
    try:
      timestamps, two_theta_low, two_theta_high, n_strong, average_i_sigi, n_lattices = HitrateStats(app, run_no, params.trial, rungroup_id, params.d_min, params.i_sigi_cutoff)()
    except Exception as e:
      print("Couldn't get run", run_no)
      continue

    n_hit = (n_strong >= params.n_strong_cutoff).count(True)
    n_indexed = (n_lattices > 0).count(True)
    n_lattices = flex.sum(n_lattices)
    n_total = len(timestamps)
    n_high_quality = ((average_i_sigi > 0) & (n_strong >= params.n_strong_cutoff)).count(True)

    invalid = (two_theta_low <= 0) or (two_theta_high < 0) # <= to prevent /0
    numerator = two_theta_high.set_selected(invalid, 0)
    denominator = two_theta_low.set_selected(two_theta_low == 0, 1) # prevent /0
    drop_ratios = numerator/denominator
    drop_hits = drop_ratios >= ratio_cutoff
    n_drop_hits = drop_hits.count(True)

    try:
      print("% 20s      % 7d % 5.1f  % 7d % 5.1f   % 7d % 5.1f    % 7d     % 7d % 5.1f % 5.1f    % 7d " % (run_no, n_drop_hits, 100*n_drop_hits/n_total, n_hit, 100*n_hit/n_total, n_indexed, 100*n_indexed/n_total, n_lattices, n_high_quality, 100*n_high_quality/n_total, 100*n_high_quality/n_indexed, n_total))
    except ZeroDivisionError:
      print("% 20s      % 7d % 5.1f  % 7d % 5.1f   % 7d % 5.1f    % 7d     % 7d % 5.1f % 5.1f    % 7d " % (run_no, n_drop_hits, 0, n_hit, 0, n_indexed, 0, n_lattices, n_high_quality, 0, 0, n_total))

    drop_hits_total += n_drop_hits
    hit_total += n_hit
    indexed_total += n_indexed
    lattices_total += n_lattices
    high_quality_total += n_high_quality
    overall_total += n_total

  if len(runs) > 1:
    print("-" * 80)
    try:
        print("Total                     % 7d % 5.1f  % 7d % 5.1f   % 7d % 5.1f    % 7d     % 7d % 5.1f % 5.1f    % 7d " % (drop_hits_total, 100*drop_hits_total/overall_total, hit_total, 100*hit_total/overall_total, indexed_total, 100*indexed_total/overall_total, lattices_total, high_quality_total, 100*high_quality_total/overall_total, 100*high_quality_total/indexed_total, overall_total))
    except ZeroDivisionError:
      if overall_total == 0:
        print("Total                     % 7d % 5.1f  % 7d % 5.1f   % 7d % 5.1f    % 7d     % 7d % 5.1f % 5.1f    % 7d " % (drop_hits_total, 0, hit_total, 0, indexed_total, 0, lattices_total, 0, 0, 0, overall_total))
      else:
        print("Total                     % 7d % 5.1f  % 7d % 5.1f   % 7d % 5.1f    % 7d     % 7d % 5.1f % 5.1f    % 7d " % (drop_hits_total, 100*drop_hits_total/overall_total, hit_total, 100*hit_total/overall_total, indexed_total, 100*indexed_total/overall_total, lattices_total, high_quality_total, 100*high_quality_total/overall_total, 0, overall_total))

if __name__ == "__main__":
  run(sys.argv[1:])
