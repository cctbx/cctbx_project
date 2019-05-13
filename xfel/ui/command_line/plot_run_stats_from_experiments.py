from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.plot_run_stats_from_experiments
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui.components.run_stats_plotter import plot_multirun_stats
import sys, os
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.integration.stills_significance_filter import SignificanceFilter, phil_scope as sf_scope

"""
Script to analyze the results of dials.stills_process and plot statisitics over time using the xfel gui.

Example command:
cctbx.xfel.plot_run_stats_from_experiments . d_min=2.0
"""

phil_str = """
  d_min = 2.0
    .type = float
    .help = High resolution bin for the I/sig(I) plot and per-run statistics.
  n_strong_cutoff = 40
    .type = int
    .help = Number of strong spots to consider an image a hit.
  i_sigi_cutoff = 1
    .type = float
    .help = Avg. I/sig(I) in a bin to reach the cutoff for producing a spot (low or high res) in the third plot.
  run_tags = None
    .type = str
    .multiple = True
    .help = Tags to be applied as labels to the runs.
  run_tags_from_filenames = True
    .type = bool
    .help = Attempt to find the run number in the filenames supplied.
  minimalist = False
    .type = bool
    .help = Generate final plot without run tags, per-run text summaries or vertical lines between runs.
  title = None
    .type = str
    .help = Plot title.
"""
phil_scope = parse(phil_str)

def run(args):
  user_phil = []
  input_dirs = []
  input_paths = []
  for arg in args:
    if os.path.isdir(arg):
      input_dirs.append(arg)
      continue
    elif os.path.exists(arg):
      input_paths.append(arg)
      continue
    try:
      user_phil.append(parse(arg))
    except Exception as e:
      raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  sf_params = sf_scope.extract()
  sf_params.significance_filter.isigi_cutoff = params.i_sigi_cutoff
  sf_params.significance_filter.d_min = params.d_min

  def get_paths(dirname):
    absolute = lambda name: os.path.join(dirname, name)
    names = os.listdir(dirname)
    return map(absolute, names)

  files_dict = {dirname:get_paths(dirname) for dirname in input_dirs}
  if params.run_tags_from_filenames:
    for path in input_paths:
      filename = os.path.basename(path)
      try:
        run = int(filename.split("idx-")[1].split("-")[0].split("run")[1])
      except Exception:
        run = int(filename.split("idx-")[1].split("-")[0].split("r")[1])
      except Exception:
        run = None
      try:
        files_dict[run].append(path)
      except KeyError:
        files_dict[run] = [path]
  else:
    files_dict[None] = input_paths
  all_results = []
  runs = []

  # iterate through grouped file paths and look for processing results
  for run in sorted(files_dict):
    files = files_dict[run]
    if len(files) == 0: continue
    runs.append(run)
    timestamps = flex.double()
    two_theta_low = flex.double()
    two_theta_high = flex.double()
    n_strong = flex.int()
    resolutions = flex.double()
    n_lattices = flex.int()

    for i, path in enumerate(sorted(files)):
      root = os.path.dirname(path)
      filename = os.path.basename(path)
      extension = os.path.splitext(filename)[1]
      split_fn = filename.split('_')
      if extension not in ['.pickle', '.mpack'] or len(split_fn) <= 0 or not split_fn[-1].startswith('strong'):
        continue
      base = os.path.join(root, "_".join(split_fn[:-1]))
      print(filename)
      strong_name = base + "_strong%s"%extension
      if not os.path.exists(strong_name):
        print("Couldn't log %s, strong%s not found"%(filename, exension))
        continue

      # Read the spotfinding results
      strong = flex.reflection_table.from_file(strong_name)
      print("N strong reflections: %d"%len(strong))

      timestamps.append(i)
      n_strong.append(len(strong))
      two_theta_low.append(0)
      two_theta_high.append(0)

      # Read indexing results if possible
      experiments_name = base + "_integrated_experiments.json"
      indexed_name = base + "_integrated%s"%extension
      if not os.path.exists(experiments_name) or not os.path.exists(indexed_name):
        print("Frame didn't index")
        resolutions.append(0)
        n_lattices.append(0)
        continue

      experiments = ExperimentListFactory.from_json_file(experiments_name, check_format=False)
      n_lattices.append(len(experiments))
      reflections = flex.reflection_table.from_file(indexed_name)
      reflections = reflections.select(reflections['intensity.sum.value'] > 0) # positive reflections only
      best_d_min = None
      for expt_id, experiment in enumerate(experiments):
        refls = reflections.select(reflections['id'] == expt_id)
        refls['id'] = flex.int(len(refls), 0)
        sig_filter = SignificanceFilter(sf_params)
        sig_filter(experiments[expt_id:expt_id+1], refls)
        if best_d_min is None or sig_filter.best_d_min < best_d_min:
          best_d_min = sig_filter.best_d_min
      resolutions.append(best_d_min or 0)

    all_results.append((timestamps, two_theta_low, two_theta_high, n_strong, resolutions, n_lattices))

  plot_multirun_stats(all_results, runs, params.d_min, n_strong_cutoff=params.n_strong_cutoff, \
    i_sigi_cutoff=params.i_sigi_cutoff, run_tags=params.run_tags, \
    minimalist=params.minimalist, interactive=True, compress_runs=True, \
    title=params.title)


if __name__ == "__main__":
  run(sys.argv[1:])
