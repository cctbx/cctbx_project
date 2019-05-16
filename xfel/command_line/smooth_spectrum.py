from __future__ import absolute_import, division, print_function
import math
import os
import sys

#from libtbx.option_parser import option_parser
import iotbx.phil
from scitbx.array_family import flex
from scitbx import smoothing
from six.moves import range
from six.moves import zip

def pyplot_label_axes(xlabel="Pixel column", ylabel="Intensity", fontsize=20):
  from matplotlib import pyplot
  pyplot.ylabel("Intensity", fontsize=fontsize)
  pyplot.xlabel("Pixel column", fontsize=fontsize)
  axes = pyplot.axes()
  for tick in axes.xaxis.get_ticklabels():
    tick.set_fontsize(20)
  for tick in axes.yaxis.get_ticklabels():
    tick.set_fontsize(20)


master_phil_str = """\
smoothing {
  method = *savitzky_golay fourier_filter
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

master_phil = iotbx.phil.parse(master_phil_str)


def run(args):
  processed = iotbx.phil.process_command_line(
    args=args, master_string=master_phil_str)
  args = processed.remaining_args
  work_params = processed.work.extract().smoothing
  savitzky_golay_half_window = work_params.savitzky_golay.half_window
  savitzky_golay_degree = work_params.savitzky_golay.degree
  fourier_cutoff = work_params.fourier_filter_cutoff

  #assert (fourier_cutoff is not None or
          #[savitzky_golay_degree, savitzky_golay_half_window].count(None) == 0)

  method = work_params.method
  if method == "fourier_filter":
    assert work_params.fourier_filter_cutoff is not None

  for i, filename in enumerate(args):
    print(filename)
    f = open(filename, 'rb')
    x, y = zip(*[line.split() for line in f.readlines() if not line.startswith("#")])
    x = flex.double(flex.std_string(x))
    y = flex.double(flex.std_string(y))
    x = x[:-2]
    y = y[:-2]
    x_orig = x.deep_copy()
    y_orig = y.deep_copy()

    x, y = interpolate(x, y)

    if method == "savitzky_golay":
      x, y_smoothed = smoothing.savitzky_golay_filter(
        x, y, savitzky_golay_half_window, savitzky_golay_degree)

    elif method == "fourier_filter":
      x, y_smoothed = fourier_filter(x, y, cutoff_frequency=fourier_cutoff)

    from matplotlib import pyplot
    fontsize = 20
    pyplot.plot(x_orig, y_orig, color='black', linestyle='dotted', linewidth=2)
    #pyplot.plot(x_orig, y_orig, color='black', linewidth=2)
    pyplot.plot(x, y_smoothed, linewidth=2, color='red')
    pyplot_label_axes()
    pyplot.show()
    #pyplot.plot(x[:385], (y - y_smoothed)[:385], linewidth=2)
    #pyplot_label_axes()
    #pyplot.show()

    filename = os.path.join(os.path.dirname(filename), "smoothed_spectrum.txt")
    f = open(filename, "wb")
    print("\n".join(["%i %f" %(xi, yi)
                           for xi, yi in zip(x, y_smoothed)]), file=f)
    f.close()
    print("Smoothed spectrum written to %s" %filename)

    x_interp_size = x.size()
    for i, x_i in enumerate(reversed(x)):
      if x_i not in x_orig:
        assert x[x_interp_size - i - 1] == x_i
        del x[x_interp_size - i - 1]
        del y[x_interp_size - i - 1]
        del y_smoothed[x_interp_size - i - 1]

    x = x[10:-10]
    y = y[10:-10]
    y_smoothed = y_smoothed[10:-10]

    signal_to_noise = estimate_signal_to_noise(x, y, y_smoothed, plot=False)

    #pyplot.plot(x[:375], signal_to_noise[:375])
    #pyplot.show()


def interpolate(x, y, half_window=10):
  perm = flex.sort_permutation(x)
  x = x.select(perm)
  y = y.select(perm)
  x_all = flex.double()
  y_all = flex.double()
  for i in range(x.size()):
    x_all.append(x[i])
    y_all.append(y[i])
    if i < x.size()-1 and (x[i+1] - x[i]) > 1:
      window_left = min(half_window, i)
      window_right = min(half_window, x.size() - i)
      x_ = x[i-window_left:i+window_right]
      y_ = y[i-window_left:i+window_right]
      from scitbx.math import curve_fitting
      # fit a 2nd order polynomial through the missing points
      polynomial = curve_fitting.univariate_polynomial(1, 1)
      fit = curve_fitting.lbfgs_minimiser([polynomial], x_,y_).functions[0]
      missing_x = flex.double(range(int(x[i]+1), int(x[i+1])))
      x_all.extend(missing_x)
      y_all.extend(fit(missing_x))

  perm = flex.sort_permutation(x_all)
  x_all = x_all.select(perm)
  y_all = y_all.select(perm)
  return x_all, y_all


def fourier_filter(x, y, cutoff_frequency):
  assert cutoff_frequency < len(y)
  from scitbx import fftpack
  fft = fftpack.real_to_complex(len(y))
  n = fft.n_real()
  m = fft.m_real()
  y_tr = y.deep_copy()
  y_tr.extend(flex.double(m-n, 0))
  fft.forward(y_tr)
  for i in range(cutoff_frequency, m):
    y_tr[i] = 0
  fft.backward(y_tr)
  y_tr = y_tr[:n]
  scale = 1 / fft.n_real()
  y_tr *= scale
  return x, y_tr[:n]


def estimate_signal_to_noise(x, y_noisy, y_smoothed, plot=False):
  """Estimate noise in spectra by subtracting a smoothed spectrum from the
     original noisy unsmoothed spectrum.

     See:
       The extraction of signal to noise values in x-ray absorption spectroscopy
       A. J. Dent, P. C. Stephenson, and G. N. Greaves
       Rev. Sci. Instrum. 63, 856 (1992); https://doi.org/10.1063/1.1142627
  """
  noise = y_noisy - y_smoothed
  noise_sq = flex.pow2(noise)
  from xfel.command_line.view_pixel_histograms import sliding_average
  sigma_sq = sliding_average(noise_sq, n=31)
  sigma_sq = smoothing.savitzky_golay_filter(
    x.as_double(), flex.pow2(noise), half_window=20, degree=1)[1]
  sigma_sq.set_selected(sigma_sq <= 0, flex.mean(sigma_sq))
  # or do this instead to use the background region as the source of noise:
  #signal_to_noise = y_smoothed/math.sqrt(flex.mean(noise_sq[50:190]))
  signal_to_noise = y_smoothed/flex.sqrt(sigma_sq)
  #signal_to_noise.set_selected(x < 50, 0)
  #signal_to_noise.set_selected(x > 375, 0)
  if plot:
    from matplotlib import pyplot
    linewidth=2
    pyplot.plot(x, y_noisy, linewidth=linewidth)
    pyplot.plot(x, y_smoothed, linewidth=linewidth)
    pyplot_label_axes()
    pyplot.show()
    pyplot.plot(x, noise, linewidth=linewidth, label="noise")
    pyplot.plot(x, flex.sqrt(sigma_sq), linewidth=linewidth, label="sigma")
    pyplot_label_axes()
    pyplot.legend(loc=2, prop={'size':20})
    pyplot.show()
    pyplot.plot(x, signal_to_noise, linewidth=linewidth)
    pyplot_label_axes()
    pyplot.show()

  return signal_to_noise


if __name__ == '__main__':
  run(sys.argv[1:])
