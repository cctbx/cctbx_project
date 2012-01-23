# LIBTBX_SET_DISPATCHER_NAME cxi.pixel_histograms
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import sys

from libtbx import easy_pickle
from libtbx.option_parser import option_parser
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

  pixels = None
  if len(args) > 0:
    pixels = [eval(arg) for arg in args]
    for pixel in pixels:
      assert isinstance(pixel, tuple)
      assert len(pixel) == 2

  #if roi is not None:

    #summed_hist = None
    #for i in range(roi[2], roi[3]):
      #for j in range(roi[0], roi[1]):
        #if summed_hist is None:
          #summed_hist = d[(i,j)]
        #else:
          #summed_hist.update(d[(i,j)])

    #title = str(roi)
    #plot(hist, window_title=window_title, title=title,log_scale=log_scale,
         #normalise=normalise, save_image=save_image, fit_gaussians=fit_gaussians)
    #return

  histograms = pixel_histograms(d)
  histograms.plot(
    pixels=pixels, starting_pixel=starting_pixel, fit_gaussians=fit_gaussians,
    window_title=window_title, log_scale=log_scale, save_image=save_image)


class pixel_histograms(object):

  def __init__(self, histograms):
    self.histograms = histograms

  def plot(self, pixels=None, starting_pixel=None, fit_gaussians=True,
           window_title=None, log_scale=False, save_image=False):
    normalise=False # XXX
    assert [pixels, starting_pixel].count(None) > 0
    if pixels is None:
      pixels = sorted(self.histograms.keys())
      if starting_pixel is not None:
        pixels = pixels[pixels.index(starting_pixel):]
    for pixel in pixels:
      hist = self.histograms[pixel]
      print pixel
      title = str(pixel)
      if fit_gaussians:
        pyplot.subplot(211)
      self.plot_one_histogram(
        hist, window_title=window_title, title=title,log_scale=log_scale,
        normalise=normalise, save_image=save_image)
      if fit_gaussians:
        self.plot_gaussians(pixel, log_scale=log_scale)

      if save_image:
        pyplot.savefig("%s.png" %title)
      else:
        pyplot.show()

  def plot_gaussians(self, pixel, log_scale=False):
    if log_scale:
      pyplot.ylim(ymin=0.1)
    hist = self.histograms[pixel]
    gaussians = self.fit_one_histogram(pixel)
    x = hist.slot_centers()
    y_calc = flex.double(x.size(), 0)
    for g in gaussians:
      print g.params
      y = g(x)
      y_calc += y
      pyplot.plot(x, y)
    #print "Peak height ratio: %.2f" %(gaussians[0].params[0]/gaussians[1].params[0])
    params = gaussians[0].params
    g = curve_fitting.gaussian(params[0], 0, params[2])
    pyplot.plot(x, g(x))
    pyplot.plot(x, y_calc)
    xlim = pyplot.xlim() # store for reuse below
    pyplot.subplot(212)
    pyplot.plot(x, hist.slots().as_double() - y_calc)
    pyplot.xlim(xlim)

  def plot_one_histogram(self, histogram,
                         window_title=None, title=None,
                         log_scale=False, normalise=False, save_image=False):
    if log_scale:
      pyplot.yscale("log")

    slots = histogram.slots().as_double()
    if normalise:
      normalisation = (flex.sum(slots) + histogram.n_out_of_slot_range()) / 1e5
      print "normalising by factor: ", normalisation
      slots /= normalisation
    bins, data = hist_outline(histogram)
    pyplot.plot(bins, data, '-k')
    #pyplot.bar(hist.slot_centers()-0.5*hist.slot_width(), slots, width=hist.slot_width())
    pyplot.xlim(histogram.data_min(), histogram.data_max())
    pyplot.suptitle(title)
    data_min = min([slot.low_cutoff for slot in histogram.slot_infos() if slot.n > 0])
    data_max = max([slot.low_cutoff for slot in histogram.slot_infos() if slot.n > 0])
    pyplot.xlim(data_min, data_max)

  def fit_one_histogram(self, pixel):
    histogram = self.histograms[pixel]
    n_gaussians_to_fit = 2
    fitted_gaussians = []

    # XXX These parameter should be configurable
    lower_fit_thresholds = [-1000, 20, 50]
    upper_fit_thresholds = [20, 40, 70]
    estimated_peak_positions = [0, 30, 60]

    slot_centers = histogram.slot_centers()
    slots = histogram.slots().as_double()

    for i in range(n_gaussians_to_fit):
      lower_threshold = lower_fit_thresholds[i]
      upper_threshold = upper_fit_thresholds[i]
      mean = estimated_peak_positions[i]
      #max_freq = 0
      #max_slot_i = 0
      #if i > 0: # XXX
        #for j, slot_center in enumerate(slot_centers):
          #if slot_center < lower_threshold: continue
          #if slot_center > upper_threshold: break
          #max_freq = max(max_freq, slots[j])
          #if max_freq == slots[j]:
            #max_slot_i = j
        #mean = slot_centers[max_slot_i]
        #lower_threshold = mean - 8
        #upper_threshold = mean + 8
      fit = self.single_peak_fit(histogram, lower_threshold, upper_threshold, mean)
      fitted_gaussians += fit.functions

    return fitted_gaussians

  def single_peak_fit(self, hist, lower_threshold, upper_threshold, mean):
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
      sigma = abs(fit.functions[0].params[2])
      if sigma < 1 or sigma > 10:
        print "using cma_es:", sigma
        fit = curve_fitting.cma_es_minimiser(
          starting_gaussians, x[lower_slot:upper_slot], y[lower_slot:upper_slot])
    else:
      fit = curve_fitting.cma_es_minimiser(
        starting_gaussians, x[lower_slot:upper_slot], y[lower_slot:upper_slot])
    return fit

  def pixels(self):
    for pixel in self.histograms.keys():
      yield pixel


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
