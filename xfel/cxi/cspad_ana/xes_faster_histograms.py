from __future__ import division
"""Specialized version of xes_histograms.
1) no support for background region of interest
2) photon_counting method only; uses integrated area under 1-photon Gaussian
3) no support for >1 photon
4) no multiprocessing
5) no output of gain map; no input gain correction
6) Fixed constraints for ratio of peak widths 1-photon: 0-photon
7) Fixed constraints for ratio of 1-photon gain: 0-photon sigma
8) 60-fold speed improvement over xes_histograms.py; takes 7 seconds.
"""

import os
import sys
import math

from libtbx import easy_pickle
import iotbx.phil
from scitbx.array_family import flex
import scitbx.math

from xfel.command_line import view_pixel_histograms
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import xes_finalise
from scitbx.lstbx import normal_eqns
from scitbx.lstbx import normal_eqns_solving
from scitbx.math import curve_fitting

from xfel.cxi.cspad_ana.xes_histograms import master_phil_str

master_phil_str = master_phil_str + """
xes {
  fudge_factor {
    gain_to_sigma = 6.75
      .type = float
      .help = On the assumption that one-mean is zero_mean + zero_sigma * gain_to_sigma
      .help = with gain_to_sigma being a constant for the pixel array detector.
      .help = approx 6.75 for LB67 r0100, 6.00 for LG36 r0025
  }
  fit_limits = (20,150)
    .type = ints(size=2)
    .help = x-Limits for histogram fitting, relative to the presumably -50 ADU histogram origin
}
"""

def run(args):
  processed = iotbx.phil.process_command_line(
    args=args, master_string=master_phil_str)
  args = processed.remaining_args
  work_params = processed.work.extract().xes
  processed.work.show()
  assert len(args) == 1
  output_dirname = work_params.output_dirname
  roi = cspad_tbx.getOptROI(work_params.roi)
  bg_roi = cspad_tbx.getOptROI(work_params.bg_roi)
  gain_map_path = work_params.gain_map
  estimated_gain = work_params.estimated_gain

  print output_dirname
  if output_dirname is None:
    output_dirname = os.path.join(os.path.dirname(args[0]), "finalise")
    print output_dirname
  hist_d = easy_pickle.load(args[0])
  if len(hist_d.keys())==2:
    hist_d = hist_d['histogram']
  pixel_histograms = faster_methods_for_pixel_histograms(
    hist_d, work_params)

  result = xes_from_histograms(
    pixel_histograms, output_dirname=output_dirname,
    gain_map_path=gain_map_path, estimated_gain=estimated_gain,
    roi=roi, run=work_params.run)

class xes_from_histograms(object):

  def __init__(self, pixel_histograms, output_dirname=".", gain_map_path=None,
               gain_map=None, estimated_gain=30,roi=None,run=None):

    self.sum_img = flex.double(flex.grid(370,391), 0) # XXX define the image size some other way?
    gain_img = flex.double(self.sum_img.accessor(), 0)

    assert [gain_map, gain_map_path].count(None) > 0
    if gain_map_path is not None:
      d = easy_pickle.load(gain_map_path)
      gain_map = d["DATA"]

    mask = flex.int(self.sum_img.accessor(), 0)

    start_row = 370
    end_row = 0
    print len(pixel_histograms.histograms)

    pixels = list(pixel_histograms.pixels())
    n_pixels = len(pixels)
    if roi is not None:
      for k, (i, j) in enumerate(reversed(pixels)):
        if (   i < roi[2]
            or i > roi[3]
            or j < roi[0]
            or j > roi[1]):
          del pixels[n_pixels-k-1]

    if gain_map is None:
      fixed_func = pixel_histograms.fit_one_histogram
    else:
      def fixed_func(pixel):
        return pixel_histograms.fit_one_histogram(pixel, n_gaussians=1)

    for i, pixel in enumerate(pixels):
      #print i,pixel
      LEG = False
      start_row = min(start_row, pixel[0])
      end_row = max(end_row, pixel[0])
      n_photons = 0

      try:
          if LEG: gaussians, two_photon_flag = pixel_histograms.fit_one_histogram(pixel)
          alt_gaussians = pixel_histograms.fit_one_histogram_two_gaussians(pixel)
      except ZeroDivisionError:
          print "HEY DIVIDE BY ZERO"
          #pixel_histograms.plot_combo(pixel, gaussians)
          mask[pixel] = 1
          continue
      except RuntimeError, e:
          print "Error fitting pixel %s" %str(pixel)
          print str(e)
          mask[pixel] = 1
          continue

      hist = pixel_histograms.histograms[pixel]

      if not LEG:
        gs = alt_gaussians[1].params
        fit_photons = gs[0] * gs[2] * math.sqrt(2.*math.pi)
        n_photons = int(round(fit_photons,0))
        if n_photons< 0:
          print "altrn %d photons from curvefitting"%( n_photons )
          pixel_histograms.plot_combo(pixel, alt_gaussians)
          mask[pixel]=1
          continue
        if False:# or two_photon_flag:
          #pixel_histograms.plot_combo(pixel, gaussians)
          pixel_histograms.plot_combo(pixel, alt_gaussians)
        #pixel_histograms.plot_one_histogram(pixel_histograms.histograms[pixel])
        #pixel_histograms.plot_gaussians(pixel)
        self.sum_img[pixel] = n_photons

    mask.set_selected(self.sum_img == 0, 1)
    unbound_pixel_mask = xes_finalise.cspad_unbound_pixel_mask()
    mask.set_selected(unbound_pixel_mask > 0, 1)
    bad_pixel_mask = xes_finalise.cspad2x2_bad_pixel_mask_cxi_run7()
    mask.set_selected(bad_pixel_mask > 0, 1)

    for row in range(self.sum_img.all()[0]):
      self.sum_img[row:row+1,:].count(0)

    spectrum_focus = self.sum_img[start_row:end_row,:]
    mask_focus = mask[start_row:end_row,:]

    spectrum_focus.set_selected(mask_focus > 0, 0)

    xes_finalise.filter_outlying_pixels(spectrum_focus, mask_focus)

    print "Number of rows: %i" %spectrum_focus.all()[0]
    print "Estimated no. photons counted: %i" %flex.sum(spectrum_focus)
    print "Number of images used: %i" %flex.sum(
      pixel_histograms.histograms.values()[0].slots())

    d = cspad_tbx.dpack(
      address='CxiSc1-0|Cspad2x2-0',
      data=spectrum_focus,
      distance=1,
      ccd_image_saturation=2e8, # XXX
    )
    if run is not None: runstr="_%04d"%run
    else: runstr=""
    cspad_tbx.dwritef(d, output_dirname, 'sum%s_'%runstr)


    plot_x, plot_y = xes_finalise.output_spectrum(
      spectrum_focus.iround(), mask_focus=mask_focus,
      output_dirname=output_dirname, run=run)
    self.spectrum = (plot_x, plot_y)
    self.spectrum_focus = spectrum_focus
    xes_finalise.output_matlab_form(spectrum_focus, "%s/sum%s.m" %(output_dirname,runstr))
    print output_dirname

class faster_methods_for_pixel_histograms(view_pixel_histograms.pixel_histograms):

  def __init__(self,hist_dict,work_params):
    self.work_params = work_params
    super(faster_methods_for_pixel_histograms,self
      ).__init__(hist_dict,work_params.estimated_gain)

  def plot_combo(self, pixel, gaussians,
                         window_title=None, title=None,
                         log_scale=False, normalise=False, save_image=False):
    histogram = self.histograms[pixel]
    from matplotlib import pyplot
    from xfel.command_line.view_pixel_histograms import hist_outline
    slots = histogram.slots().as_double()
    if normalise:
      normalisation = (flex.sum(slots) + histogram.n_out_of_slot_range()) / 1e5
      print "normalising by factor: ", normalisation
      slots /= normalisation
    bins, data = hist_outline(histogram)
    if log_scale:
      data.set_selected(data == 0, 0.1) # otherwise lines don't get drawn when we have some empty bins
      pyplot.yscale("log")
    pyplot.plot(bins, data, '-k', linewidth=2)
    pyplot.plot(bins, data/1000., '-k', linewidth=2)
    pyplot.suptitle(title)
    data_min = min([slot.low_cutoff for slot in histogram.slot_infos() if slot.n > 0])
    data_max = max([slot.low_cutoff for slot in histogram.slot_infos() if slot.n > 0])
    pyplot.xlim(data_min, data_max+histogram.slot_width())
    pyplot.xlim(-50, 100)
    pyplot.ylim(-10,40)
    x = histogram.slot_centers()
    y_calc = flex.double(x.size(), 0)
    for g in gaussians:
      print "Height %7.2f mean %4.1f sigma %3.1f"%(g.params)
      y = g(x)
      y_calc += y
      pyplot.plot(x, y, linewidth=2)

    pyplot.show()

  def fit_one_histogram_two_gaussians(self,pixel):
    histogram = self.histograms[pixel]
    fitted_gaussians = []
    GAIN_TO_SIGMA = self.work_params.fudge_factor.gain_to_sigma
    low_idx = self.work_params.fit_limits[0]
    high_idx = self.work_params.fit_limits[1]

    slot_centers = histogram.slot_centers()
    free_x = slot_centers[low_idx:high_idx]
    #print list(free_x)
    slots = histogram.slots().as_double()
    free_y = slots[low_idx:high_idx]

    # zero_mean = 0. # originally intended mean=0
    maxidx = flex.max_index(free_y) # but if dark subtraction (pedstal correction) is off
    zero_mean = free_x[maxidx] # use this non-zero maximum instead

    zero_amplitude = flex.max(free_y)
    assert 1./zero_amplitude #guard against division by zero
    total_population = flex.sum(free_y)
    zero_sigma = self.estimated_gain / GAIN_TO_SIGMA
    one_amplitude = 0.001
    helper = self.per_frame_helper_factory(initial_estimates =
      (zero_mean, 1.0, zero_sigma, one_amplitude),
      GAIN_TO_SIGMA=GAIN_TO_SIGMA,
      free_x = free_x,
      free_y = free_y/zero_amplitude) # put y values on 0->1 scale for normal eqn solving
    helper.restart()
    iterations = normal_eqns_solving.levenberg_marquardt_iterations(
          non_linear_ls = helper,
          n_max_iterations = 40,
          gradient_threshold = 1.E-5)
    #print "current values after iterations", list(helper.x),

    fitted_gaussians = helper.as_gaussians()
    for item in fitted_gaussians: item.params = (item.params[0] * zero_amplitude,
                                  item.params[1], item.params[2]) # convert back to full scale
    return fitted_gaussians

  @staticmethod
  def per_frame_helper_factory(initial_estimates,GAIN_TO_SIGMA,free_x,free_y):
      SIGMAFAC = 1.15

      class per_frame_helper(normal_eqns.non_linear_ls, normal_eqns.non_linear_ls_mixin):
        def __init__(pfh):
          super(per_frame_helper, pfh).__init__(n_parameters=4)
          pfh.x_0 = flex.double(initial_estimates)
          pfh.restart()

        def restart(pfh):
          pfh.x = pfh.x_0.deep_copy()
          pfh.old_x = None

        def step_forward(pfh):
          pfh.old_x = pfh.x.deep_copy()
          pfh.x += pfh.step()

        def step_backward(pfh):
          assert pfh.old_x is not None
          pfh.x, pfh.old_x = pfh.old_x, None

        def parameter_vector_norm(pfh):
          return pfh.x.norm()

        def build_up(pfh, objective_only=False):
          residuals = pfh.fvec_callable(pfh.x)

          pfh.reset()
          if objective_only:
            pfh.add_residuals(residuals, weights=None)
          else:
            grad_r = pfh.jacobian_callable(pfh.x)
            jacobian = flex.double(
              flex.grid(len(grad_r[0]), pfh.n_parameters))
            for j, der_r in enumerate(grad_r):
              jacobian.matrix_paste_column_in_place(der_r,j)
            pfh.add_equations(residuals, jacobian, weights=None)

        def fvec_callable(pfh,current_values):
          z_mean = current_values[0]
          z_ampl = current_values[1]
          z_sigm = current_values[2]
          o_ampl = current_values[3]
          o_mean = z_mean + z_sigm * GAIN_TO_SIGMA
          o_sigm = z_sigm * SIGMAFAC

          # model minus obs
          # sqrt2pi_inv = 1./math.sqrt(2.*math.pi)
          # the gaussian function is sqrt2pi_inv * exp( - (x-mean)**2 / (2.*(sigma**2)))/sigma
          # take off the coefficient sqrt2pi_inv / sigma, use ampl

          terms = [

          ampl * flex.exp(-flex.pow2(free_x - mean) / (2.*sigm*sigm))

          for mean,ampl,sigm in [(z_mean,z_ampl,z_sigm),(o_mean,o_ampl,o_sigm)]]

          model = terms[0]+terms[1]
          pfh.terms = terms
          return model - free_y

        def jacobian_callable(pfh,current_values):
          z_mean = current_values[0]
          z_ampl = current_values[1]
          z_sigm = current_values[2]
          o_ampl = current_values[3]
          o_mean = z_mean + z_sigm * GAIN_TO_SIGMA
          o_sigm = z_sigm * SIGMAFAC

          udiff = free_x - z_mean
          Afactor = udiff/(z_sigm * z_sigm)
          #Sfactor = (udiff*udiff)*math.pow(z_sigm, -3.)
          Sfactor = Afactor * udiff/z_sigm

          # leaving out small cross terms where one-photon peak influences
          # derivatives with respect to z_mean and z_sigm.

          return (pfh.terms[0]*Afactor,
                  pfh.terms[0]/z_ampl,
                  pfh.terms[0]*Sfactor,
                  pfh.terms[1]/o_ampl)

        def as_gaussians(pfh):
          return [curve_fitting.gaussian( a = pfh.x[1], b = pfh.x[0], c = pfh.x[2] ),
                  curve_fitting.gaussian( a = pfh.x[3], b = pfh.x[0] + pfh.x[2] * GAIN_TO_SIGMA,
                                          c = pfh.x[2] * SIGMAFAC )]

      value = per_frame_helper()
      return value


if __name__ == '__main__':
  run(sys.argv[1:])
