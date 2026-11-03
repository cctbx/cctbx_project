'''
Author      : Young, I.D.
Created     : 07/14/2016
Last Changed: 07/14/2016
Description : XFEL UI print per-run hit, indexing and high quality counts for a trial
'''
from __future__ import absolute_import, division, print_function

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui import load_cached_settings, master_phil_str, db_phil_str
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.stats import HitrateStats
from xfel.ui.components.run_stats_plotter import drop_hit_selection, high_quality_selection
import sys
from scitbx.array_family import flex
from six.moves import zip

phil_str = """
  trial = None
    .type = int
  run = None
    .type = str
    .multiple = True
    .help = "Runs to print; requires rungroup. Default is every run in the trial's active run groups."
  rungroup = None
    .type = int
    .help = "With run, the run group the runs were processed in. Otherwise, print only this run group."
  d_min = 2.0
    .type = float
    .help = "Resolution an indexed image must reach to count as high quality."
  n_strong_cutoff = 16
    .type = int
    .help = "Number of strong spots to consider an image a hit."
  ratio_cutoff = 1
    .type = float
    .help = "Ratio of the radial average at the high two theta to the low one above which an image shows solvent."
  i_sigi_cutoff = 1
    .type = float
    .help = "Avg. I/sig(I) a resolution bin must reach to count towards an image's resolution."
  run_tags = None
    .type = str
    .multiple = True
    .help = "Print only runs with these tags."
  tag_selection_mode = *union intersection
    .type = choice
    .help = "Keep runs with any (union) or all (intersection) of run_tags."
"""
phil_scope = parse(phil_str + master_phil_str + db_phil_str)

def percent(n, d):
  ''' 100 * n / d, or 0 when d is 0 '''
  return 100 * n / d if d else 0

def format_row(label, n_drop_hits, n_hit, n_indexed, n_lattices, n_high_quality, n_total, status=""):
  return "% 20s      % 7d % 5.1f  % 7d % 5.1f   % 7d % 5.1f    % 7d     % 7d % 5.1f % 5.1f    % 7d %s" % (
    label, n_drop_hits, percent(n_drop_hits, n_total), n_hit, percent(n_hit, n_total),
    n_indexed, percent(n_indexed, n_total), n_lattices,
    n_high_quality, percent(n_high_quality, n_total), percent(n_high_quality, n_indexed), n_total, status)

def run_sort_key(run_no):
  ''' Numeric order for integer run names, string order for the rest '''
  try:
    return (0, int(run_no), "")
  except ValueError:
    return (1, 0, run_no)

def run(args):
  user_phil = []
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception as e:
      raise Sorry("Unrecognized argument %s"%arg)

  scope = load_cached_settings(scope=phil_scope, extract=False)
  params = scope.fetch(sources=user_phil).extract()
  if params.trial is None:
    raise Sorry("Specify a trial, e.g. trial=0")
  print("Printing results for trial", params.trial, "using a hit cutoff of", params.n_strong_cutoff, "strong spots")
  print()
  print("                 Run  N Drop Hits   (%)   N Hits   (%) N Indexed   (%) N Lattices N High qual   (%)  %HQR   N Frames Status")

  drop_hits_total = 0
  hit_total = 0
  indexed_total = 0
  lattices_total = 0
  high_quality_total = 0
  overall_total = 0

  app = xfel_db_application(params)
  # All jobs, not only those of active trials and run groups, so a run's
  # status is found whatever the state of its trial.
  jobs = app.get_all_jobs()

  trial = app.get_trial(trial_number=params.trial)
  # The trial's run group list holds only its active run groups.
  if params.rungroup is not None and params.rungroup not in [rg.id for rg in trial.rungroups]:
    raise Sorry("Run group %d is not an active run group of trial %d" % (params.rungroup, params.trial))

  if params.run is None or len(params.run) == 0:
    runs = []
    run_ids = []
    rungroups = []
    for rg in trial.rungroups:
      for run in rg.runs:
        if run.id in run_ids: continue
        if params.run_tags:
          requested = set(params.run_tags)
          present = set(tag.name for tag in run.tags)
          if params.tag_selection_mode == "union":
            if not requested & present: continue
          elif not requested <= present: continue

        runs.append(run.run)
        run_ids.append(run.id)
        rungroups.append(rg.id)
  else:
    runs = params.run
    if params.rungroup is None:
      raise Sorry("Specify the run group of the given runs, e.g. rungroup=1")
    rungroups = [params.rungroup] * len(runs)

  for run_no, rungroup_id in sorted(zip(runs, rungroups), key=lambda x: run_sort_key(x[0])):
    if params.rungroup and params.rungroup != rungroup_id: continue
    try:
      timestamps, two_theta_low, two_theta_high, n_strong, resolutions, n_lattices, wavelengths = HitrateStats(app, run_no, params.trial, rungroup_id, params.d_min, params.i_sigi_cutoff)()
    except Exception as e:
      print("Couldn't get run %s: %s" % (run_no, e))
      continue

    n_hit = (n_strong >= params.n_strong_cutoff).count(True)
    n_indexed = (n_lattices > 0).count(True)
    n_lattices = flex.sum(n_lattices)
    n_total = len(timestamps)
    n_high_quality = high_quality_selection(resolutions, params.d_min).count(True)
    n_drop_hits = drop_hit_selection(two_theta_low, two_theta_high, params.ratio_cutoff)[1].count(True)

    status = "Not found"
    for job in jobs:
      if job.run and job.run.run == run_no and job.rungroup and job.rungroup.id == rungroup_id and \
         job.trial and job.trial.trial == params.trial:
        status = job.status
        break

    print(format_row(run_no, n_drop_hits, n_hit, n_indexed, n_lattices, n_high_quality, n_total, status))

    drop_hits_total += n_drop_hits
    hit_total += n_hit
    indexed_total += n_indexed
    lattices_total += n_lattices
    high_quality_total += n_high_quality
    overall_total += n_total

  if len(runs) > 1:
    print("-" * 80)
    print(format_row("Total", drop_hits_total, hit_total, indexed_total, lattices_total, high_quality_total, overall_total))

if __name__ == "__main__":
  run(sys.argv[1:])
