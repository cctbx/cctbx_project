from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.plot_uc_cloud_from_experiments
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from libtbx.phil import parse

help_message = """
Plot a cloud of unit cell dimensions from stills. Provide either a combined_experiments.json file or a specify individual .json files on the command line.
"""

phil_str = """
  iqr_ratio = 1.5
    .type = float
    .help = Interquartile range multiplier for outlier rejection. Use None to disable outlier rejection.
  ranges = None
    .type = floats(6)
    .help = Lower and upper bounds for the ranges to display for each of the a, b and c axes
"""
phil_scope = parse(phil_str)

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name

    self.tag = None
    self.reference_detector = None

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_experiments=True,
      check_format=False,
      epilog=help_message
      )

  def run(self):
    '''Execute the script.'''
    from dials.util.options import flatten_experiments
    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    if all([len(e.data) == 1 for e in params.input.experiments]):
      experiments_list = [flatten_experiments(params.input.experiments)]
    else:
      experiments_list = [e.data for e in params.input.experiments]

    info_list = []
    for experiments in experiments_list:
      info = []
      for experiment in experiments:
        a, b, c, alpha, beta, gamma = experiment.crystal.get_unit_cell().parameters()

        info.append({'a':a,
                     'b':b,
                     'c':c,
                     'alpha':alpha,
                     'beta':beta,
                     'gamma':gamma,
                     'n_img':0})
      info_list.append(info)
    import xfel.ui.components.xfel_gui_plotter as pltr
    plotter = pltr.PopUpCharts()
    plotter.plot_uc_histogram(info_list=info_list, legend_list=[""]*len(experiments_list), iqr_ratio = params.iqr_ratio, ranges = params.ranges)
    if len(experiments_list) == 1:
      plotter.plot_uc_3Dplot(info, iqr_ratio = params.iqr_ratio)
    plotter.plt.show()

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
