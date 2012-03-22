import math
import os
import sys

from matplotlib import pyplot

from libtbx.option_parser import option_parser
from scitbx.array_family import flex
from scitbx import smoothing

def run(args):

  command_line = (option_parser()
                  .option("--savitzky_golay_half_window", "--sg_half_window",
                          type="int",
                          help="The number of values either side of a data point"
                              "to use.")
                  .option("--savitzky_golay_degree", "--sg_degree",
                          type="int",
                          help="The degree of polynomial to fit to data points.")
                  .option("--fourier_filter_cutoff",
                          type="int",
                          help="Cutoff frequency for Fourier filtering.")
                  ).process(args=args)
  args = command_line.args
  savitzky_golay_half_window = command_line.options.savitzky_golay_half_window
  savitzky_golay_degree = command_line.options.savitzky_golay_degree
  fourier_cutoff = command_line.options.fourier_filter_cutoff

  assert (fourier_cutoff is not None or
          [savitzky_golay_degree, savitzky_golay_half_window].count(None) == 0)

  for i, filename in enumerate(args):
    print filename
    f = open(filename, 'rb')
    x, y = zip(*[line.split() for line in f.readlines() if not line.startswith("#")])
    x = flex.double(flex.std_string(x))
    y = flex.double(flex.std_string(y))
    x_orig = x.deep_copy()
    y_orig = y.deep_copy()

    x, y = interpolate(x, y)

    if savitzky_golay_degree is not None:
      x, y_smoothed = smoothing.savitzky_golay_filter(
        x, y, savitzky_golay_half_window, savitzky_golay_degree)

    elif fourier_cutoff is not None:
      x, y_smoothed = fourier_filter(x, y, cutoff_frequency=fourier_cutoff)

    pyplot.plot(x_orig, y_orig, color='black', linestyle='dotted', linewidth=2)
    pyplot.plot(x, y_smoothed, linewidth=2, color='red')
    pyplot.show()
    pyplot.plot(x, y - y_smoothed)
    pyplot.show()

    filename = os.path.join(os.path.dirname(filename), "smoothed_spectrum.txt")
    f = open(filename, "wb")
    print >> f, "\n".join(["%i %f" %(xi, yi)
                           for xi, yi in zip(x, y_smoothed)])
    f.close()
    print "Smoothed spectrum written to %s" %filename

    x_interp_size = x.size()
    for i, x_i in enumerate(reversed(x)):
      if x_i not in x_orig:
        assert x[x_interp_size - i - 1] == x_i
        del x[x_interp_size - i - 1]
        del y[x_interp_size - i - 1]
        del y_smoothed[x_interp_size - i - 1]

    signal_to_noise = estimate_signal_to_noise(x, y, y_smoothed)

    pyplot.plot(x[:375], signal_to_noise[:375])
    pyplot.show()


def interpolate(x, y):
  perm = flex.sort_permutation(x)
  x = x.select(perm)
  y = y.select(perm)
  x_all = flex.double()
  y_all = flex.double()
  for i in range(x.size()):
    x_all.append(x[i])
    y_all.append(y[i])
    if i < x.size()-1 and (x[i+1] - x[i]) > 1:
      window_left = min(10, i)
      window_right = min(10, x.size() - i)
      x_ = x[i-window_left:i+window_right]
      y_ = y[i-window_left:i+window_right]
      from scitbx.math import curve_fitting
      # fit a 2nd order polynomial through the missing points
      polynomial = curve_fitting.univariate_polynomial(1, 1, 1)
      fit = curve_fitting.lbfgs_minimiser([polynomial], x_,y_).functions[0]
      missing_x = flex.double(range(int(x[i]), int(x[i+1])))
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
       Rev. Sci. Instrum. 63, 856 (1992); http://dx.doi.org/10.1063/1.1142627
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
    pyplot.plot(x, y_noisy, linewidth=2)
    pyplot.plot(x, y_smoothed, linewidth=2)
    pyplot.show()
    pyplot.plot(x, noise, linewidth=2)
    pyplot.plot(x, flex.sqrt(sigma_sq), linewidth=2)
    pyplot.show()
    pyplot.plot(x[:375], signal_to_noise[:375])
    pyplot.show()

  return signal_to_noise


if __name__ == '__main__':
  run(sys.argv[1:])
