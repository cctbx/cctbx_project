from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial as Poly

"""When an energy scan is conducted at LCLS, we acquire FEE spectra of the incident beam with a varying, known, narrow energy band removed -- the "notch". Energy calibration is the process of identifying the notch in each scan and using the known pixel-energy pairs to generate a function (linear fit) returning the energy for any given pixel position on the spectrometer (reported FEE energy in xtc streams). This file supplies functions necessary for automating this process. Functionality custom to LCLS (calling the psana api to read xtc streams) and the end-to-end automated process are provided in cctbx_project/xfel/command_line."""

notch_phil_string = """
fit_half_range = 3
  .type = int
  .help = how many steps (pixels) above and below the notch do we use to fit the notch shape
kernel_size = 30
  .type = int
  .help = how many steps (pixels) do we incorporate into the moving average smoothing kernel
baseline_cutoff = 0.01
  .type = float
  .help = consider this fraction of the beam profile maximum to be the baseline, and when we inflate values to account for the curved shape of the profile, don't inflate (divide) by any smaller fractions than this
reference_spectrum = None
  .type = path
  .help = if we have a notch-free FEE spectrum to use as a reference, subtract this directly instead of a smoothed spectrum to convert the notch to an absolute minimum
per_run_plots = False
  .type = bool
  .help = for each run, additionally plot original pixel data, smoothed or reference data, and scaled flattened data from which notch position is finally determined
"""

def find_notch(data_x, data_y, kernel_size, fit_half_range, baseline_cutoff, ref_spectrum=None, minima=False):
  """Given a series that smoothly increases and then decreases, with a notch cut out somewhere,
     find the position in data_x that we expect (via interpolation) corresponds to the true
     minimum or maximum of the notch in data_y. This is empirically a good way to find the notches of known
     energies in the FEE spectra."""
  if ref_spectrum is not None: # subtract baseline from data_y
    flattened_y = data_y - ref_spectrum
  else: # subtract smoothed data_y from data_y
    # smooth data_y
    kernel = np.ones(kernel_size)/kernel_size
    smoothed_y = np.convolve(data_y, kernel, mode='same')
    if minima:
      # subtract smoothed shape to isolate notch as an absolute minimum
      flattened_y = data_y - smoothed_y
      # scale to account for scan shape
      minimum = max(smoothed_y)*baseline_cutoff
      scale = np.max([smoothed_y, np.ones(len(smoothed_y))*minimum], axis=0)
      flattened_y /= scale
    else:
      flattened_y = smoothed_y
  # isolate a range where the notch can be approximated by a small degree polynomial
  if minima:
    rough_min_idx = np.argmin(flattened_y)
  else:
    rough_max_idx = np.argmax(flattened_y)

  start = rough_max_idx - fit_half_range
  end = rough_max_idx + fit_half_range
  # model it as a parabola
  notch_shape = Poly.fit(data_x[start:end+1], data_y[start:end+1], 2)
  fitted_x = notch_shape.deriv().roots()[0]
  fitted_y = notch_shape.__call__(fitted_x)
  return ((fitted_x, fitted_y), notch_shape, smoothed_y, flattened_y)

def plot_notches(runs, rundata, notches, per_run_plots=False, use_figure=None):
  """Plot the energy scan, optionally one plot per spectrum for troubleshooting misses when automatically identifying notch positions, and always as an overlay of spectra in the scan with notch positions marked."""
  if per_run_plots:
    for run, data_y, notch in zip(runs, rundata, notches):
      data_x = range(len(data_y))
      (notch_x, notch_y), notch_shape, smoothed_y, flattened_y = notch
      plt.title("FEE spectrum for run %d"%run)
      plt.xlabel("FEE spectrometer pixels")
      plt.ylabel("Mean counts")
      plt.plot(data_x, data_y, 'g-', label="raw energy scan")
      plt.plot(data_x, smoothed_y, 'y-', label="smoothed scan")
      plt.plot(data_x, flattened_y*max(smoothed_y), 'b-', label="scaled flattened scan")
      plt.plot([notch_x], [notch_y], 'k+', label=f"notch found at {int(notch_x)} pixels".format())
      plt.legend()
      plt.figure()

  fig = use_figure or plt.figure()
  ax = fig.subplots()
  for run, data, notch in zip(runs, rundata, notches):
    (notch_x, notch_y), notch_shape, smoothed_y, flattened_y = notch
    ax.plot(range(len(data)), data, '-', label=f"run {run}: notch at {int(notch_x)} pixels".format())
    ax.plot([notch_x], [notch_y], 'k+', label="_nolegend_")
  # repeat last one to add legend
  ax.plot([notch_x], [notch_y], 'k+', label="identified notches")

  ax.legend()
  ax.set_title("Energy scan")
  ax.set_xlabel("FEE spectrometer pixels")
  ax.set_ylabel("Mean counts")

def calibrate_energy(notches, energies, return_trendline=False, use_figure=None):
  """Having identified the pixel positions on the FEE spectrometer corresponding to known energy values, get a linear fit of these ordered pairs and report the eV offset and eV per pixel matching the fit."""
  pixels = [n[0][0] for n in notches]
  linear_fit = Poly.fit(pixels, energies, 1).convert()
  eV_offset, eV_per_pixel = linear_fit.coef
  print(f"Calibrated eV offset of {eV_offset} and eV per pixel of {eV_per_pixel}".format())
  fig = use_figure or plt.figure()
  ax = fig.subplots()
  ax.scatter(pixels, energies, color='k', label="known energy positions")
  px_min = int(min(pixels))
  px_max = int(max(pixels))
  px_range = px_max - px_min
  trendline_x = np.arange(px_min-0.1*px_range, px_max+0.1*px_range, int(px_range/10))
  trendline_y = linear_fit.__call__(trendline_x)
  ax.plot(trendline_x, trendline_y, 'b-', label=f"linear fit: y = {eV_per_pixel:.4f}*x + {eV_offset:.4f}".format())
  ax.set_xlabel("FEE spectrometer pixels")
  ax.set_ylabel("Energy (eV)")
  ax.set_title("Energy calibration")
  ax.legend()
  if use_figure is None:
    plt.show()
  if return_trendline:
    return ((eV_offset, eV_per_pixel), (linear_fit, trendline_x, trendline_y))
  else:
    return (eV_offset, eV_per_pixel)

