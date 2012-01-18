# LIBTBX_SET_DISPATCHER_NAME cxi.pixel_histograms
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

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
                  .option("--save",
                          action="store_true",
                          default=False,
                          help="Save each plot as a png.")
                  .option("--start",
                          type="string",
                          help="Starting pixel coordinates")
                  .option("--fit_gaussians",
                          action="store_true",
                          default=False,
                          help="Fit gaussians to the peaks.")
                  ).process(args=args)
  log_scale = command_line.options.log_scale
  fit_gaussians = command_line.options.fit_gaussians
  roi = cspad_tbx.getOptROI(command_line.options.roi)
  normalise = command_line.options.normalise
  save_image = command_line.options.save
  starting_pixel = command_line.options.start
  if starting_pixel is not None:
    starting_pixel = eval(starting_pixel)
    assert isinstance(starting_pixel, tuple)
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
          summed_hist.update(d[(i,j)])

    title = str(roi)
    plot(hist, window_title=window_title, title=title,log_scale=log_scale,
         normalise=normalise, save_image=save_image, fit_gaussians=fit_gaussians)
    return

  if len(args) > 0:
    for arg in args:
      pixel = eval(arg)
      assert isinstance(pixel, tuple)
      assert len(pixel) == 2
      hist = d[pixel]
      print pixel
      title = str(pixel)
      plot(hist, window_title=window_title, title=title,log_scale=log_scale,
           normalise=normalise, save_image=save_image, fit_gaussians=fit_gaussians)
  else:
    pixels = sorted(d.keys())
    if starting_pixel is not None:
      pixels = pixels[pixels.index(starting_pixel):]
    for pixel in pixels:
      hist = d[pixel]
      print pixel
      title = str(pixel)
      plot(hist, window_title=window_title, title=title,log_scale=log_scale,
           normalise=normalise, save_image=save_image, fit_gaussians=fit_gaussians)

def plot(hist, window_title=None, title=None,
         log_scale=False, normalise=False, save_image=False,
         fit_gaussians=False):
  pyplot.clf() # clear current figure
  if log_scale:
    pyplot.yscale("log")

  slots = hist.slots().as_double()
  if normalise:
    normalisation = (flex.sum(slots) + hist.n_out_of_slot_range()) / 1e5
    print "normalising by factor: ", normalisation
    slots /= normalisation
  bins, data = hist_outline(hist)
  pyplot.plot(bins, data, '-k')
  #pyplot.bar(hist.slot_centers()-0.5*hist.slot_width(), slots, width=hist.slot_width())
  pyplot.xlim(hist.data_min(), hist.data_max())
  pyplot.suptitle(title)
  data_min = min([slot.low_cutoff for slot in hist.slot_infos() if slot.n > 0])
  data_max = max([slot.low_cutoff for slot in hist.slot_infos() if slot.n > 0])
  pyplot.xlim(data_min, data_max)

  if fit_gaussians:

    import scitbx.lbfgs

    n_gaussians_to_fit = 2

    fitted_gaussians = []

    lower_fit_thresholds = [-1000, 20, 50]
    upper_fit_thresholds = [20, 40, 70]
    estimated_peak_positions = [0, 30, 60]

    x = hist.slot_centers()
    y = hist.slots().as_double()

    for i in range(n_gaussians_to_fit):
      lower_threshold = lower_fit_thresholds[i]
      upper_threshold = upper_fit_thresholds[i]
      mean = estimated_peak_positions[i]

      fit = peak_fit(hist, lower_threshold, upper_threshold, mean)

      #diff_peak_position = estimated_peak_positions[i] - fit.functions[0].params[1]
      #if abs(diff_peak_position) > 1:
        #lower_threshold += diff_peak_position
        #upper_threshold += diff_peak_position
        #mean += diff_peak_position
        #print "Refitting peak: delta peak position: %.2f" %diff_peak_position
        #print fit.functions[0].params
        #fit = peak_fit(hist, lower_threshold, upper_threshold, mean)
        #print fit.functions[0].params

      fitted_gaussians += fit.functions

    y_calc = flex.double(x.size(), 0)
    for g in fitted_gaussians:
      #print "%.3f %.3f %.3f %.3f" %(g.a, g.b, g.c, g.area())
      #scale = g.a
      #mu = g.b
      #sigma = g.c
      #print scale, mu, sigma
      #y = scale * flex.exp(-flex.pow2(x - mu) / (2 * sigma**2))
      print g.params
      y = g(x)
      y_calc += y
      pyplot.plot(x, y)
    pyplot.plot(x, y_calc)
    if log_scale:
      pyplot.ylim(ymin=0.1)
  pyplot.show()

  if save_image:
    pyplot.savefig("%s.png" %title)
  else:
    pyplot.show()
  return

def peak_fit(hist, lower_threshold, upper_threshold, mean):
  starting_gaussians = [curve_fitting.gaussian(a=100, b=mean, c=5)]
  lower_slot = 0
  for slot in hist.slot_centers():
    lower_slot += 1
    if slot > lower_threshold: break
  upper_slot = 0
  for slot in hist.slot_centers():
    upper_slot += 1
    if slot > upper_threshold: break

  x = hist.slot_centers()
  y = hist.slots().as_double()
  if 1:
    fit = curve_fitting.lbfgs_minimiser(
      starting_gaussians, x[lower_slot:upper_slot], y[lower_slot:upper_slot])
  else:
    fit = curve_fitting.cma_es_minimiser(
      starting_gaussians, x[lower_slot:upper_slot], y[lower_slot:upper_slot])
  return fit



def hist_outline(hist):

  step_size = hist.slot_width()
  half_step_size = 0.5 * step_size
  n_slots = len(hist.slots())

  bins = flex.double(n_slots * 2 + 2, 0)
  data = flex.double(n_slots * 2 + 2, 0)
  for i in range(n_slots):
    bins[2 * i + 1] = hist.slot_centers()[i] - half_step_size
    bins[2 * i + 2] = hist.slot_centers()[i] + half_step_size
    data[2 * i + 1] = hist.slots()[i]
    data[2 * i + 2] = hist.slots()[i]

  bins[0] = bins[1] - step_size
  bins[-1] = bins[-2] + step_size
  data[0] = 0
  data[-1] = 0

  return (bins, data)



if __name__ == '__main__':
  run(sys.argv[1:])
