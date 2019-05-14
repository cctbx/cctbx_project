from __future__ import absolute_import, division, print_function
import sys

from matplotlib import pyplot

import iotbx.phil
from scitbx.array_family import flex

from scitbx import smoothing
from xfel.command_line import smooth_spectrum

master_phil_str = """\
x_offsets = None
  .type = floats
legend = None
  .type = strings
  .help = Labels for each plot.
bg_range = 100, 200
  .type = floats(size=2)
  .help = Range in pixels to align background of spectra.
plot_range = None
  .type = floats(size=2)
  .help = Range in pixels to plot
smoothing {
  method = savitzky_golay fourier_filter
    .type = choice
  savitzky_golay {
    half_window = 16
      .type = int
      .help = The number of values either side of a data point to use.
    degree = 4
      .type = int
      .help = The degree of polynomial to fit to data points.
  }
  fourier_filter_cutoff = None
    .type = int
    .help = Cutoff frequency for Fourier filtering.
}
"""

def run(args):
  master_phil = iotbx.phil.parse(master_phil_str)
  processed = iotbx.phil.process_command_line(
    args=args, master_string=master_phil_str)
  args = processed.remaining_args
  work_params = processed.work.extract()

  x_offsets = work_params.x_offsets
  bg_range_min, bg_range_max = work_params.bg_range
  if work_params.plot_range is not None:
    x_min, x_max = work_params.plot_range
  else:
    x_min, x_max = (0, 385)

  print(bg_range_min, bg_range_max)
  if x_offsets is None:
    x_offsets = [0]*len(args)
  legend = work_params.legend
  linewidth = 2
  fontsize = 26
  xy_pairs = []
  colours = ["cornflowerblue", "darkmagenta", "darkgreen", "black", "red", "blue", "pink"]
  colours[2] = "orangered"
  colours[1] = "olivedrab"
  min_background = 1e16
  #x_min, x_max = (0, 391)
  #x_min, x_max = (0, 360)
  #x_min, x_max = (200, 360)

  for i, filename in enumerate(args):
    print(filename)
    f = open(filename, 'rb')
    x, y = zip(*[line.split() for line in f.readlines() if not line.startswith("#")])
    x = flex.double(flex.std_string(x))
    y = flex.double(flex.std_string(y))

    if work_params.smoothing.method is not None:
      savitzky_golay_half_window = work_params.smoothing.savitzky_golay.half_window
      savitzky_golay_degree = work_params.smoothing.savitzky_golay.degree
      fourier_cutoff = work_params.smoothing.fourier_filter_cutoff

      method = work_params.smoothing.method
      if method == "fourier_filter":
        assert work_params.smoothing.fourier_filter_cutoff is not None

      if method == "savitzky_golay":
        x, y = smoothing.savitzky_golay_filter(
          x, y, savitzky_golay_half_window, savitzky_golay_degree)

      elif method == "fourier_filter":
        x, y = smooth_spectrum.fourier_filter(x, y, cutoff_frequency=fourier_cutoff)


    x += x_offsets[i]
    y = y.select((x <= x_max) & (x > 0))
    x = x.select((x <= x_max) & (x > 0))
    bg_sel = (x > bg_range_min) & (x < bg_range_max)
    xy_pairs.append((x,y))
    min_background = min(min_background, flex.mean(y.select(bg_sel))/flex.max(y))
    y -= min_background
    print("Peak maximum at: %i" %int(x[flex.max_index(y)]))
  for i, filename in enumerate(args):
    if legend is None:
      label = filename
    else:
      print(legend)
      assert len(legend) == len(args)
      label = legend[i]
    x, y = xy_pairs[i]
    if i == -1:
      x, y = interpolate(x, y)
      x, y = savitzky_golay_filter(x, y)
    #if i == 0:
      #y -= 10
    bg_sel = (x > bg_range_min) & (x < bg_range_max)
    y -= (flex.mean(y.select(bg_sel)) - min_background*flex.max(y))
    #y -= flex.min(y)
    y_min = flex.min(y.select(bg_sel))
    if i == -2:
      y += 0.2 * flex.max(y)
    print("minimum at: %i" %int(x[flex.min_index(y)]), flex.min(y))
    #print "fwhm: %.2f" %full_width_half_max(x, y)
    y /= flex.max(y)
    if len(colours) > i:
      pyplot.plot(x, y, label=label, linewidth=linewidth, color=colours[i])
    else:
      pyplot.plot(x, y, label=label, linewidth=linewidth)
  pyplot.ylabel("Intensity", fontsize=fontsize)
  pyplot.xlabel("Pixel column", fontsize=fontsize)
  if i > 0:
    # For some reason the line below causes a floating point error if we only
    # have one plot (i.e. i==0)
    legend = pyplot.legend(loc=2)
    for t in legend.get_texts():
      t.set_fontsize(fontsize)
  axes = pyplot.axes()
  for tick in axes.xaxis.get_ticklabels():
    tick.set_fontsize(20)
  for tick in axes.yaxis.get_ticklabels():
    tick.set_fontsize(20)
  pyplot.ylim(0,1)
  pyplot.xlim(x_min, x_max)
  ax = pyplot.axes()
  #ax.xaxis.set_minor_locator(pyplot.MultipleLocator(5))
  #ax.yaxis.set_major_locator(pyplot.MultipleLocator(0.1))
  #ax.yaxis.set_minor_locator(pyplot.MultipleLocator(0.05))
  pyplot.show()

def full_width_half_max(x, y):
  y = y/flex.max(y)
  perm = flex.sort_permutation(x)
  y = y.select(perm)
  x = x.select(perm)
  x_lower = None
  x_upper = None
  for x_i, y_i in zip(x, y):
    if x_lower is None:
      if y_i >= 0.5:
        x_lower = x_i
    elif x_upper is None:
      if y_i <= 0.5:
        x_upper = x_i
    else:
      break
  return (x_upper - x_lower)


def estimate_signal_to_noise(x, y):
  raise
  if 1:
    x, y = interpolate(x, y)
    #x, y_tr = fourier_filter(x, y)
    x, y_tr = savitzky_golay_filter(x, y)
    noise = y - y_tr
  else:

    from scitbx.math import chebyshev_polynome
    from scitbx.math import chebyshev_lsq_fit

    x_obs, y_obs = x, y
    w_obs = flex.double(y_obs.size(), 1)
    w_obs[0] = 1e16
    w_obs[-1] = 1e16
    ## determining the number of terms takes much, much longer than the fit
    n_terms = chebyshev_lsq_fit.cross_validate_to_determine_number_of_terms(
      x_obs, y_obs, w_obs,
      min_terms=2, max_terms=30,
      n_goes=20, n_free=20)
    #n_terms = 7
    print("n_terms:", n_terms)
    fit = chebyshev_lsq_fit.chebyshev_lsq_fit(n_terms, x_obs, y_obs, w_obs)
    fit_funct = chebyshev_polynome(
      n_terms, fit.low_limit, fit.high_limit, fit.coefs)
    y_fitted = fit_funct.f(x)
    y_tr = y_fitted
    n = y_tr.size()
    noise = y - y_tr


  noise_sq = flex.pow2(noise)
  from xfel.command_line.view_pixel_histograms import sliding_average
  #sigma_sq = sliding_average(noise_sq, n=31)
  sigma_sq = sliding_average(noise_sq, n=15)
  #sigma_sq = sliding_average(sigma_sq)
  #signal_to_noise = y/flex.sqrt(sigma_sq)
  import math
  signal_to_noise = y/math.sqrt(flex.mean(noise_sq[50:200]))
  #pyplot.plot(noise)
  #pyplot.plot(x,y)
  #pyplot.show()
  offset = 0.2 * flex.max(y)
  offset = 0
  pyplot.plot(x, y, linewidth=2)
  pyplot.plot(x, offset+y_tr, linewidth=2)
  pyplot.show()
  pyplot.plot(x, noise, linewidth=2)
  #pyplot.plot(x, flex.sqrt(sigma_sq), linewidth=2)
  #ax2 = pyplot.twinx()
  #ax2.plot(x, y)
  pyplot.show()
  pyplot.plot(x[:375], signal_to_noise[:375])
  #pyplot.xlim(
  #ax2 = pyplot.twinx()
  #ax2.plot(x, y)
  pyplot.show()

if __name__ == '__main__':
  run(sys.argv[1:])
