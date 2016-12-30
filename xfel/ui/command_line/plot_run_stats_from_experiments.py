from __future__ import division

from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.ui.components.run_stats_plotter import plot_multirun_stats
import sys, os
from scitbx.array_family import flex
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from libtbx import easy_pickle

"""
Script to analyze the results of dials.stills_process and plot statisitics over time using the xfel gui.

Example command:
libtbx.python plot_run_stats_from_experiments.py . d_min=2.0
"""

phil_str = """
  hit_cutoff = 30
    .type = int
    .help = Number of reflections to consider an image a hit. Estimate by looking at plot of strong reflections/image.
  d_min = None
    .type = float
    .help = Highest resolution to consider for I/sigI plot
  compress_runs = True
    .type = bool
    .help = When plotting multiple runs, adjust timestamps so there is no blank space between them.
    .help = Thise mode is not compatible with fetching events from timestamps.
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
    average_i_sigi_low = flex.double()
    average_i_sigi_high = flex.double()

    for i, filename in enumerate(sorted(os.listdir(root))):
      split_fn = filename.split('_')
      if len(split_fn) <= 0 or split_fn[-1] != "datablock.json":
        continue
      base = "_".join(split_fn[:-1])
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
      experiments_name = base + "_refined_experiments.json"
      indexed_name = base + "_indexed.pickle"
      if not os.path.exists(experiments_name) or not os.path.exists(indexed_name):
        print "Frame didn't index"
        average_i_sigi_low.append(0)
        average_i_sigi_high.append(0)
        continue

      experiments = ExperimentListFactory.from_json_file(experiments_name, check_format=False)
      reflections = easy_pickle.load(indexed_name)
      reflections = reflections.select(reflections['intensity.sum.value'] > 0) # positive reflections only
      assert len(experiments) == 1 and len(experiments.crystals()) == 1
      crystal = experiments.crystals()[0]
      uc = crystal.get_unit_cell()
      d = uc.d(reflections['miller_index'])

      # Bin the reflections possibly present in this space group/cell so that we can report average I/sigma
      # in the highest and lowest resolution cell
      from cctbx.crystal import symmetry
      cs = symmetry(unit_cell = uc, space_group_info=crystal.get_space_group().info())
      mset = cs.build_miller_set(anomalous_flag=False, d_min=params.d_min)
      binner = mset.setup_binner(n_bins=10)
      for j, i_sigi in zip([0,-1], [average_i_sigi_low, average_i_sigi_high]):
        d_max, d_min = binner.bin_d_range(binner.range_used()[j])
        refls = reflections.select((d <= d_max) & (d > d_min))
        n_refls = len(refls)

        avg_i_sigi = flex.mean(refls['intensity.sum.value'] /
                               flex.sqrt(refls['intensity.sum.variance'])) if n_refls > 0 else 0
        i_sigi.append(avg_i_sigi)

    all_results.append((timestamps, two_theta_low, two_theta_high, n_strong, average_i_sigi_low, average_i_sigi_high))

  plot_multirun_stats(all_results, runs, params.d_min, n_strong_cutoff=params.hit_cutoff, \
    interactive=True, compress_runs=params.compress_runs)

if __name__ == "__main__":
  run(sys.argv[1:])
