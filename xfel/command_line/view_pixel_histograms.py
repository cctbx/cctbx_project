# LIBTBX_SET_DISPATCHER_NAME cxi.pixel_histograms
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import sys

from libtbx import easy_pickle
from libtbx.option_parser import option_parser
from libtbx.utils import frange
from scitbx.array_family import flex
from scitbx.math import curve_fitting

from matplotlib import pyplot

from xfel.cxi.cspad_ana import cspad_tbx

def run(args):
  assert len(args) > 0
  command_line = (option_parser()
                  .option("--roi",
                          type="string",
                          help="Region of interest for summing up histograms"
                          "from neighbouring pixels.")
                  .option("--log_scale",
                          action="store_true",
                          default=False,
                          help="Draw y-axis on a log scale.")
                  .option("--normalise",
                          action="store_true",
                          default=False,
                          help="Normalise by number of member images.")
                  ).process(args=args)
  log_scale = command_line.options.log_scale
  roi = cspad_tbx.getOptROI(command_line.options.roi)
  normalise = command_line.options.normalise
  args = command_line.args

  path = args[0]
  window_title = path
  d = easy_pickle.load(path)
  args = args[1:]

  if roi is not None:

    summed_hist = None
    for i in range(roi[2], roi[3]):
      for j in range(roi[0], roi[1]):
        if summed_hist is None:
          summed_hist = d[(i,j)]
        else:
          summed_hist.update_from_histogram(d[(i,j)])

    title = str(roi)
    plot(summed_hist, window_title=window_title, title=title,
         log_scale=log_scale, normalise=normalise)
    return

  if len(args) > 0:
    for arg in args:
      pixel = eval(arg)
      assert isinstance(pixel, tuple)
      assert len(pixel) == 2
      hist = d[pixel]
      print pixel
      title = str(pixel)
      plot(hist, window_title=window_title, title=title,
           log_scale=log_scale, normalise=normalise)
  else:
    for pixel, hist in sorted(d.items()):
      print pixel
      title = str(pixel)
      plot(hist, window_title=window_title, title=title,
           log_scale=log_scale, normalise=normalise)

def plot(hist, window_title=None, title=None, log_scale=False, normalise=True):
  if log_scale:
    pyplot.yscale("log")

  #fig = pyplot.figure()
  #fig.canvas.set_window_title(window_title)
  slots = hist.slots().as_double()
  if normalise:
    normalisation = (flex.sum(slots) + hist.n_out_of_slot_range()) / 1e5
    print "normalising by factor: ", normalisation
    slots /= normalisation
  pyplot.bar(hist.slot_centers(), slots, width=hist.slot_width())
  pyplot.xlim(hist.data_min(), hist.data_max())
  pyplot.suptitle(title)
  pyplot.show()
  return

  starting_gaussians = [
    curve_fitting.gaussian(scale=10, mu=0, sigma=1),
    #curve_fitting.gaussian(scale=1, mu=4, sigma=1),
    #curve_fitting.gaussian(scale=50, mu=200, sigma=2.1),
    #curve_fitting.gaussian(scale=20, mu=16, sigma=2.1),
    #curve_fitting.gaussian(scale=10, mu=24, sigma=2.1),
  ]
  x = hist.slot_centers()
  y = hist.slots().as_double()

  import scitbx.lbfgs
  termination_params = None
  #termination_params = scitbx.lbfgs.termination_parameters(max_iterations=0)
  #fit = curve_fitting.single_gaussian_fit(x, y)
  #print fit.scale, fit.mu, fit.sigma
  #scale = fit.scale
  #mu = fit.mu
  #sigma = fit.sigma
  ##print scale, mu, sigma
  ##x = flex.double(frange(hist_min, hist_max, step=0.1))
  #x = hist.slot_centers()
  #y = scale * flex.exp(-flex.pow2(x - mu) / (2 * sigma**2))
  #y = hist.slots().as_double() - y
  ##pyplot.plot(x,y)


  weights = y
  #weights = None
  fit = curve_fitting.gaussian_fit(x, y, starting_gaussians)#, weights=weights)
  #fit = curve_fitting.single_gaussian_fit(x, y)
  #fit.gaussians = [curve_fitting.gaussian(fit.scale, fit.mu, fit.sigma)]
  fit_gaussians = [
    curve_fitting.gaussian(scale=8100, mu=0.1, sigma=0.4),
    #curve_fitting.gaussian(scale=220, mu=8.5, sigma=1.7),
    #curve_fitting.gaussian(scale=80, mu=17, sigma=1.7),
    #curve_fitting.gaussian(scale=40, mu=25.5, sigma=1.7),
  ]
  x = flex.double(frange(hist.data_min(), hist.data_max(), step=0.1))
  y_calc = flex.double(x.size())
  #for g in fit.gaussians:
  for g in fit_gaussians:
    print "%.3f %.3f %.3f %.3f" %(g.scale, g.mu, g.sigma, g.area())
    scale = g.scale
    mu = g.mu
    sigma = g.sigma
    print scale, mu, sigma
    y = scale * flex.exp(-flex.pow2(x - mu) / (2 * sigma**2))
    y_calc += y
    pyplot.plot(x, y)
  pyplot.plot(x, y_calc)
  pyplot.ylim(ymin=0)
  pyplot.show()


if __name__ == '__main__':
  run(sys.argv[1:])
