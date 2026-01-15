from __future__ import absolute_import, division, print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.plot_uc_cloud_from_experiments
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from dials.util import show_mail_on_error
from iotbx.phil import parse
from cctbx import crystal

help_message = """
Plot a cloud of unit cell dimensions from stills. Provide either a combined.expt
file or a specify individual .expt files on the command line. To generate an overlay of
multiple plots (similar to grouping by run tag in the XFEL GUI), provide multiple
combined.expt files named as ${tag}_combined_*.expt and set extract_tags to True in
the phil scope.
"""

phil_str = """
  iqr_ratio = None
    .type = float
    .help = Interquartile range multiplier for outlier rejection. Use None to disable outlier rejection.
  ranges = None
    .type = floats(6)
    .help = Lower and upper bounds for the ranges to display for each of the a, b and c axes
  angle_ranges = None
    .type = floats(6)
    .help = Lower and upper bounds for the ranges to display for each of the cell angles
  extract_tags = False
    .type = bool
    .help = Extract tags from the names of multiple combined.expt filenames and use
    .help = these tags to label multiple groups of experiments.
  combine_all_input = False
    .type = bool
    .help = Combine all experiment lists to generate a single unit cell histogram. Useful for plotting
    .help = data from multiple runs outside the cctbx.xfel GUI.
  title = None
    .type = str
    .help = Title for the plot
  write_gnuplot_cloud = False
    .type = bool
    .help = If True, write a list of cells such that gnuplot can read them. Start gnuplot then use \
            splot "uccloud.dat".
  output.image_file = None
    .type = str
  scale = *linear log
    .type = choice
    .help = Z-axis scale for 2d histograms
  move_near {
    unit_cell = None
      .type = unit_cell
    space_group = P1
      .type = space_group
    relative_length_tolerance = 0.05
      .type = float
    absolute_angle_tolerance_deg = 3
      .type = float
    enable = True
      .type = bool
  }

"""
phil_scope = parse(phil_str)

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import ArgumentParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name

    self.tag = None
    self.reference_detector = None

    # Create the parser
    self.parser = ArgumentParser(
      usage=usage,
      phil=phil_scope,
      read_experiments=True,
      check_format=False,
      epilog=help_message
      )

  def run(self):
    '''Execute the script.'''
    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    if params.write_gnuplot_cloud:
      gnuplot = open("uccloud.dat", 'w')

    if params and params.move_near.enable and params.move_near.unit_cell is not None:
      # Create reference symmetry object
      self.ref_sym = crystal.symmetry(
        unit_cell=params.move_near.unit_cell,
        space_group_symbol=str(params.move_near.space_group)
      )


    def get_info(experiment, params=None):
      # Get the crystal symmetry
      crystal_sym = experiment.crystal

      # Apply move_near transformation if enabled
      if params and params.move_near.enable and params.move_near.unit_cell is not None:
        # Create reference symmetry object
        # Create experiment's symmetry object
        exp_sym = crystal.symmetry(
          unit_cell=crystal_sym.get_unit_cell(),
          space_group_symbol=str(crystal_sym.get_space_group().type().lookup_symbol())
        )
        # Find nearest setting
        moved_sym = self.ref_sym.nearest_setting(
          exp_sym,
          length_tolerance=params.move_near.relative_length_tolerance,
          angle_tolerance=params.move_near.absolute_angle_tolerance_deg
        )
        a, b, c, alpha, beta, gamma = moved_sym.unit_cell().parameters()
      else:
        a, b, c, alpha, beta, gamma = crystal_sym.get_unit_cell().parameters()

      return {'a':a,
              'b':b,
              'c':c,
              'alpha':alpha,
              'beta':beta,
              'gamma':gamma,
              'n_img':0}
    experiments_list = [e.data for e in params.input.experiments]
    if params.extract_tags:
      import os
      experiments_tags = [os.path.splitext(os.path.basename(f.filename))[0] for f in params.input.experiments]
      info_list = []
      for experiments in experiments_list:
        infos = []
        for experiment in experiments:
          info = get_info(experiment, params=params)
          infos.append(info)
          if params.write_gnuplot_cloud:
            gnuplot.write("% 3.10f % 3.10f % 3.10f\n"%(info['a'], info['b'], info['c']))
        info_list.append(infos)
    else:
      experiments_tags = [str(i) for i in range(len(experiments_list))]
      info_list = []
      for experiments in experiments_list:
        infos = []
        for experiment in experiments:
          info = get_info(experiment, params=params)
          infos.append(info)
          if params.write_gnuplot_cloud:
            gnuplot.write("% 3.10f % 3.10f % 3.10f\n"%(info['a'], info['b'], info['c']))
        info_list.append(infos)
    if params.combine_all_input:
      info_list = [[info for infos in info_list for info in infos]]
      experiments_tags = ["combined"]

    if params.write_gnuplot_cloud:
      gnuplot.close()

    import xfel.ui.components.xfel_gui_plotter as pltr
    interactive = params.output.image_file is None
    plotter = pltr.PopUpCharts(interactive=interactive)
    plotter.plot_uc_histogram(
      info_list=info_list,
      legend_list=experiments_tags,
      iqr_ratio = params.iqr_ratio,
      ranges = params.ranges,
      angle_ranges = params.angle_ranges,
      title = params.title,
      image_fname = params.output.image_file,
      hist_scale = params.scale)
    plotter.plt.show()

if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
