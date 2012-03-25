from stdlib import math
from stdlib import random

from libtbx.utils import frange
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import scitbx.math
import scitbx.random
from scitbx.math import curve_fitting
from scitbx.smoothing import savitzky_golay_filter, savitzky_golay_coefficients
from scitbx.smoothing import convolve

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)


def exercise_convolve():
  data = flex.double(20, 0)
  for i in range(5,10): data[i] = 1
  response = flex.double([1,1,1,0,0,0,1,1])
  response = flex.double([1,1,1,1,1])
  convolved = convolve(data, response)
  assert approx_equal(convolve(data, response, mode="same"),
                      [0,0,0,1,2,3,4,5,4,3,2,1,0,0,0,0,0,0,0,0])

  try:
    import numpy
  except ImportError:
    print "Skipping numpy compatibility..."
    return

  # convolution of two rectangles gives a triangle
  for data_size in (20,21):
    data = flex.double(data_size, 0)
    for i in range(5,10): data[i] = 1
    for mode in ("full", "same", "valid"):
      for response_size in (5,6):
        response = flex.double([1]*response_size)
        assert approx_equal(
          convolve(data, response, mode=mode),
          numpy.convolve(data, response, mode=mode))


def exercise_savitzky_golay_coefficients():
  coeffs = savitzky_golay_coefficients(5, 5, 4, wraparound=False)
  assert approx_equal(
    coeffs,
    (0.042, -0.105, -0.023, 0.140, 0.280, 0.333, 0.280, 0.140, -0.023, -0.105, 0.042), eps=1e-3)

  coeffs = savitzky_golay_coefficients(4, 4, 4, wraparound=False)
  assert approx_equal(
    coeffs,
    (0.035, -0.128, 0.070, 0.315, 0.417, 0.315, 0.070, -0.128, 0.035), eps=1e-3)

  coeffs = savitzky_golay_coefficients(4, 0, 2, wraparound=False)
  assert approx_equal(
    coeffs,
    (0.086, -0.143, -0.086, 0.257, 0.886), eps=1e-3)


def exercise_savitzky_golay_smoothing():

  plot = False

  def rms(flex_double):
    return math.sqrt(flex.mean(flex.pow2(flex_double)))

  for sigma_frac in (0.005, 0.01, 0.05, 0.1):
    mean = random.randint(-5,5)
    scale = flex.random_double() * 10
    sigma = flex.random_double() * 5 + 1
    gaussian = curve_fitting.gaussian(scale, mean, sigma)

    x = flex.double(frange(-20,20,0.1))
    y = gaussian(x)
    rand_norm = scitbx.random.normal_distribution(
      mean=0, sigma=sigma_frac*flex.max_absolute(y))
    g = scitbx.random.variate(rand_norm)
    noise = g(y.size())
    y_noisy = y + noise
    # according to numerical recipes the best results are obtained where the
    # full window width is between 1 and 2 times the number of points at fwhm
    # for polynomials of degree 4
    half_window = int(round(0.5 * 2.355 * sigma * 10))
    y_filtered = savitzky_golay_filter(x, y_noisy, half_window=half_window, degree=4)[1]
    extracted_noise = y_noisy - y_filtered
    rms_noise = rms(noise)
    rms_extracted_noise = rms(extracted_noise)

    assert abs(rand_norm.sigma - rms_noise)/rand_norm.sigma < 0.15
    assert abs(rand_norm.sigma - rms_extracted_noise)/rand_norm.sigma < 0.15

    diff = y_filtered - y
    assert (rms(diff)/ rand_norm.sigma) < 0.4

    if plot:
      from matplotlib import pyplot
      pyplot.plot(x, y)
      pyplot.plot(x, noise)
      pyplot.scatter(x, y_noisy, marker="x")
      pyplot.plot(x, y_filtered)
      pyplot.show()
      pyplot.plot(x, extracted_noise)
      pyplot.plot(x, noise)
      pyplot.show()

  return


if __name__ == '__main__':
  exercise_convolve()
  exercise_savitzky_golay_smoothing()
  exercise_savitzky_golay_coefficients()
  print "OK"
