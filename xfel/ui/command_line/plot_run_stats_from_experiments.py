from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.plot_run_stats_from_experiments
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui.components.run_stats_plotter import plot_multirun_stats
import sys, os
from scitbx.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx import easy_pickle

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
phil_scope = parse(phil_str)

def run(args):
  user_phil = []
  input_dirs = []
  for arg in args:
    if os.path.isdir(arg):
      input_dirs.append(arg)
      continue
    try:
      user_phil.append(parse(arg))
    except Exception, e:
      raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()

  all_results = []
  runs = []
  # iterate through supplied directories and look for processing results
  for dir_number, root in enumerate(input_dirs):
    runs.append(dir_number)
    timestamps = flex.double()
    two_theta_low = flex.double()
    two_theta_high = flex.double()
    n_strong = flex.int()
    average_i_sigi = flex.double()
    n_lattices = flex.int()

    for i, filename in enumerate(sorted(os.listdir(root))):
      split_fn = filename.split('_')
      if len(split_fn) <= 0 or split_fn[-1] != "datablock.json":
        continue
      base = os.path.join(root, "_".join(split_fn[:-1]))
      print filename
      strong_name = base + "_strong.pickle"
      if not os.path.exists(strong_name):
        print "Couldn't log %s, strong pickle not found"%filename
        continue

      # Read the spotfinding results
      strong = easy_pickle.load(strong_name)
      print "N strong reflections: %d"%len(strong)

      timestamps.append(i)
      n_strong.append(len(strong))
      two_theta_low.append(0)
      two_theta_high.append(0)

      # Read indexing results if possible
      experiments_name = base + "_integrated_experiments.json"
      indexed_name = base + "_integrated.pickle"
      if not os.path.exists(experiments_name) or not os.path.exists(indexed_name):
        print "Frame didn't index"
        average_i_sigi.append(0)
        n_lattices.append(0)
        continue

      experiments = ExperimentListFactory.from_json_file(experiments_name, check_format=False)
      n_lattices.append(len(experiments))
      reflections = easy_pickle.load(indexed_name)
      reflections = reflections.select(reflections['intensity.sum.value'] > 0) # positive reflections only
      best_avg_i_sigi = 0
      for expt_id, experiment in enumerate(experiments):
        crystal = experiment.crystal
        uc = crystal.get_unit_cell()
        refls = reflections.select(reflections['id'] == expt_id)
        d = uc.d(reflections['miller_index'])

        # Bin the reflections possibly present in this space group/cell so that we can report average I/sigma
        # in the highest requested bin
        from cctbx.crystal import symmetry
        cs = symmetry(unit_cell = uc, space_group_info=crystal.get_space_group().info())
        mset = cs.build_miller_set(anomalous_flag=False, d_min=params.d_min)
        binner = mset.setup_binner(n_bins=10)
        d_max, d_min = binner.bin_d_range(binner.range_used()[-1]) # highest res
        refls = refls.select((d <= d_max) & (d > d_min))
        n_refls = len(refls)

        avg_i_sigi = flex.mean(refls['intensity.sum.value'] /
                               flex.sqrt(refls['intensity.sum.variance'])) if n_refls > 0 else 0
        if avg_i_sigi > best_avg_i_sigi:
          best_avg_i_sigi = avg_i_sigi
        print best_avg_i_sigi, d_max, d_min, len(refls)
      average_i_sigi.append(best_avg_i_sigi)

    all_results.append((timestamps, two_theta_low, two_theta_high, n_strong, average_i_sigi, n_lattices))

  plot_multirun_stats(all_results, runs, params.d_min, n_strong_cutoff=params.n_strong_cutoff, \
    i_sigi_cutoff=params.i_sigi_cutoff, run_tags=params.run_tags, \
    minimalist=params.minimalist, interactive=True, compress_runs=params.compress_runs, \
    title=params.title)

if __name__ == "__main__":
  run(sys.argv[1:])
