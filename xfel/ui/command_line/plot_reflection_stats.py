from __future__ import absolute_import, division, print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.plot_reflection_stats
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
from matplotlib import pyplot as plt
from xfel.command_line.detector_residuals import setup_stats
from dials.array_family import flex
from dials.util import show_mail_on_error
import math
from libtbx.phil import parse
from scitbx.math import five_number_summary

help_message = """
Two plots:
deltaXY: plot of the difference between observed and predicted reflections in microns vs. two theta. Only reflections with I/sigI >= 5 are included.
I/sigI: plot of I/sigI for all reflections vs. two theta.

Two theta bins with fewer than 10 reflections are ommited.

For each plot, the median value per two theta bin is plotted.  The shaded area comprised the upper and lower quartiles.
N reflections per two theta bin is also plotted.

Example:
cctbx.xfel.plot_reflection_stats integrated1.expt integrated1.refl integrated2.expt integrated2.refl
"""

phil_str = """
save_pdf = False
  .type = bool
  .help = If true, save results as pdf
tag = None
  .type = str
  .help = Title for the plots
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

    # Create the parser
    self.parser = ArgumentParser(
      usage=usage,
      phil=phil_scope,
      check_format=False,
      read_experiments=True,
      read_reflections=True,
      epilog=help_message
      )

  def run(self):
    '''Execute the script.'''
    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    n_bins = 10
    arbitrary_padding = 1
    legend = []

    # local container for axes
    figures = {}

    for plot in ['deltaXY','isigi']:
      fig, ax1 = plt.subplots()
      fig.suptitle(plot)
      ax2 = ax1.twinx()
      figures[plot] = {'fig':fig, 'ax1': ax1, 'ax2': ax2}

    def plotit(reflections, experiments):
      """
      Make the plots for a set of reflections and experiments.
      """
      beam = experiments.beams()[0] # only used to compute resolution of 2theta
      reflections = reflections.select(reflections['intensity.sum.variance'] > 0)

      # Setup up deltaXY and two theta bins
      reflections['difference_vector_norms'] = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()
      reflections = setup_stats(experiments, reflections, two_theta_only=True) # add two theta to reflection table
      sorted_two_theta = flex.sorted(reflections['two_theta_obs'])
      bin_low = [sorted_two_theta[int((len(sorted_two_theta)/n_bins) * i)] for i in range(n_bins)]
      bin_high = [bin_low[i+1] for i in range(n_bins-1)]
      bin_high.append(sorted_two_theta[-1]+arbitrary_padding)

      x_centers = flex.double()
      n_refls = flex.int()
      rmsds = flex.double()
      p25r = flex.double()
      p50r = flex.double()
      p75r = flex.double()
      p25i = flex.double()
      p50i = flex.double()
      p75i = flex.double()
      print("# 2theta Res N dXY IsigI")

      # Compute stats for each bin
      for i in range(n_bins):
        refls = reflections.select((reflections['two_theta_obs'] >= bin_low[i]) & (reflections['two_theta_obs'] < bin_high[i]))
        # Only compute deltaXY stats on reflections with I/sigI at least 5
        i_sigi = refls['intensity.sum.value']/flex.sqrt(refls['intensity.sum.variance'])
        refls = refls.select(i_sigi >= 5)
        n = len(refls)
        if n < 10: continue
        min_r, q1_r, med_r, q3_r, max_r = five_number_summary(1000*refls['difference_vector_norms'])

        n_refls.append(n)

        rmsds_ = 1000*math.sqrt(flex.sum_sq(refls['difference_vector_norms'])/n)

        min_i, q1_i, med_i, q3_i, max_i = five_number_summary(i_sigi)
        p25i.append(q1_i)
        p50i.append(med_i)
        p75i.append(q3_i)
        # x_center
        c = ((bin_high[i]-bin_low[i])/2) + bin_low[i]
        # resolution
        d = beam.get_wavelength()/(2*math.sin(math.pi*c/(2*180)))
        x_centers.append(c)
        rmsds.append(rmsds_)
        print("%d % 5.1f % 5.1f % 8d %.1f %.1f"%(i, c, d, n, med_r, med_i))
        p25r.append(q1_r)
        p50r.append(med_r)
        p75r.append(q3_r)

      # After binning, plot the results
      for plot in figures:
        ax1 = figures[plot]['ax1']
        ax2 = figures[plot]['ax2']
        if plot == 'isigi':
          line, = ax1.plot(x_centers.as_numpy_array(), p50i.as_numpy_array(), '-')
          line.set_label('Median')
          ax1.fill_between(x_centers.as_numpy_array(), p25i.as_numpy_array(), p75i.as_numpy_array(),
            interpolate = True, alpha = 0.50, color = line.get_color())
          line, = ax2.plot(x_centers.as_numpy_array(), n_refls.as_numpy_array(), '-', color = line.get_color())
          line.set_label('Median')
        elif plot == 'deltaXY':
          line, = ax1.plot(x_centers.as_numpy_array(), p50r.as_numpy_array(), '-')
          line.set_label('Median')
          ax1.fill_between(x_centers.as_numpy_array(), p25r.as_numpy_array(), p75r.as_numpy_array(),
            interpolate = True, alpha = 0.50, color = line.get_color())
          line, = ax2.plot(x_centers.as_numpy_array(), n_refls.as_numpy_array(), '-', color = line.get_color())
          line.set_label('Median')
        ax1.legend()
        ax2.legend()

    assert len(params.input.experiments) == len(params.input.reflections)

    # Plotit!
    for i in range(len(params.input.experiments)):
      plotit(params.input.reflections[i].data, params.input.experiments[i].data)

    # Set up labels
    for plot in figures:
      fig = figures[plot]['fig']
      ax1 = figures[plot]['ax1']
      ax2 = figures[plot]['ax2']
      if plot == 'isigi':
        ax1.set_ylabel("I/sigI")
        ax2.set_ylabel("N reflections")
        ax1.set_xlabel("Two theta (degrees)")
      elif plot == 'deltaXY':
        ax1.set_ylabel(r"$\Delta$XY")
        ax2.set_ylabel("N reflections")
        ax1.set_xlabel("Two theta (degrees)")
      if params.tag is not None:
        pass

    if params.save_pdf:
      from matplotlib.backends.backend_pdf import PdfPages
      pp = PdfPages('reflection_stats.pdf')
      for i in plt.get_fignums():
        pp.savefig(plt.figure(i), dpi=300)
      pp.close()
    else:
      plt.show()

if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
