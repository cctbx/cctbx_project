from __future__ import division
from __future__ import print_function
from six.moves import range
import h5py
import math
from matplotlib import pyplot as plt
import numpy as np
import iotbx.phil
from scitbx.array_family import flex
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
      .help = approx 6.75 for LB67 r0100
      .help = manually optimized for LG36: r0025,6.0 / r0080,5.8 (MnCl2) / r0137,5.9 (PSII solution)
  }
  estimated_gain = 30.
    .type = float
    .help = CSPAD gain for photons recorded from inelastic scattering.
    .help = Not critical, only used to calculate the zero photon peak width
  fit_limits = (20,150)
    .type = ints(size=2)
    .help = x-Limits for histogram fitting, relative to the presumably -50 ADU histogram origin
    .help = 20,150 used for LG36
  first_slot_value = -50
    .type = int
    .help = x-Limit for left hand boundary of histogram, presumably -50 ADU
    .help = -50 used for LG36
  gaussian_3 {
    zero_sigma = 5.
      .type = float
    inelastic_gain_to_sigma = 6.5
      .type = float
      .help = On the assumption that one-mean is zero_mean + zero_sigma * gain_to_sigma
      .help = with gain_to_sigma being a constant for the pixel array detector.
      .help = approx 6.75 for LB67 r0100
      .help = manually optimized for LG36: r0025,6.0 / r0080,5.8 (MnCl2) / r0137,5.9 (PSII solution)
    elastic_gain_to_sigma = 10.1
      .type = float
      .help = On the assumption that one-mean is zero_mean + zero_sigma * gain_to_sigma
      .help = with gain_to_sigma being a constant for the pixel array detector.
      .help = approx 6.75 for LB67 r0100
      .help = manually optimized for LG36: r0025,6.0 / r0080,5.8 (MnCl2) / r0137,5.9 (PSII solution)
  }
}
"""
Usage = """cctbx.python xes_faster_histograms.py output_dirname=[outdir] \
   roi=0:388,99:126 [datadir]/hist_r0[run_no].pickle run=[run_no] fudge_factor.gain_to_sigma=5.9
   ...converts histogram.pickle file (described elsewhere) into spectrum by fitting
      0- and 1-photon Gaussians to histograms representing each pixel on the XES spectrometer."""

def get_work_params(args=[],verbose=True):
  processed = iotbx.phil.process_command_line(
    args=args, master_string=master_phil_str)
  args = processed.remaining_args
  work_params = processed.work.extract().xes
  if verbose:
    processed.work.show()

  return work_params

from xfel.cxi.cspad_ana.xes_faster_histograms import faster_methods_for_pixel_histograms
class derived_class(faster_methods_for_pixel_histograms):
  def __init__(self,params):
    self.work_params = params
  def fit_one_histogram_two_gaussians(self,histogram):
    fitted_gaussians = []
    GAIN_TO_SIGMA = self.work_params.fudge_factor.gain_to_sigma
    low_idx = self.work_params.fit_limits[0]
    high_idx = self.work_params.fit_limits[1]

    slot_centers = flex.double(range(self.work_params.first_slot_value,
                                      self.work_params.first_slot_value + len(histogram)))
    free_x = slot_centers[low_idx:high_idx]
    #print list(free_x)

    slots = flex.double(histogram.astype(np.float64))
    free_y = slots[low_idx:high_idx]

    # zero_mean = 0. # originally intended mean=0
    maxidx = flex.max_index(free_y) # but if dark subtraction (pedstal correction) is off
    zero_mean = free_x[maxidx] # use this non-zero maximum instead

    zero_amplitude = flex.max(free_y)

    assert 1./zero_amplitude #guard against division by zero
    total_population = flex.sum(free_y)
    zero_sigma = self.work_params.estimated_gain / GAIN_TO_SIGMA
    one_amplitude = 0.001
    helper = self.per_pixel_helper_factory(initial_estimates =
      (zero_mean, 1.0, zero_sigma, one_amplitude),
      GAIN_TO_SIGMA=GAIN_TO_SIGMA,
      free_x = free_x,
      free_y = free_y/zero_amplitude) # put y values on 0->1 scale for normal eqn solving
    helper.restart()
    iterations = normal_eqns_solving.levenberg_marquardt_iterations(
          non_linear_ls = helper,
          n_max_iterations = 7,
          gradient_threshold = 1.E-3)
    print("current values after iterations", list(helper.x), end=' ')

    fitted_gaussians = helper.as_gaussians()
    for item in fitted_gaussians: item.params = (item.params[0] * zero_amplitude,
                                  item.params[1], item.params[2]) # convert back to full scale
    return fitted_gaussians

  def plot_combo(self, histogram, gaussians,
                         window_title=None, title=None,
                         log_scale=False, normalise=False, save_image=False, interpretation=None):

    from matplotlib import pyplot
    #from xfel.command_line.view_pixel_histograms import hist_outline
    slots = flex.double(histogram.astype(np.float64))
    if normalise:
      normalisation = (flex.sum(slots) + histogram.n_out_of_slot_range()) / 1e5
      print("normalising by factor: ", normalisation)
      slots /= normalisation
    slot_centers = flex.double(range(self.work_params.first_slot_value,
                                      self.work_params.first_slot_value + len(histogram)))
    bins, data = hist_outline(slot_width=1,slots=slots,slot_centers=slot_centers)
    if log_scale:
      data.set_selected(data == 0, 0.1) # otherwise lines don't get drawn when we have some empty bins
      pyplot.yscale("log")
    pyplot.plot(bins, data, '-k', linewidth=2)
    pyplot.plot(bins, data/1000., '-k', linewidth=2)
    pyplot.suptitle(title)
    #data_min = min([slot.low_cutoff for slot in histogram.slot_infos() if slot.n > 0])
    #data_max = max([slot.low_cutoff for slot in histogram.slot_infos() if slot.n > 0])
    #pyplot.xlim(data_min, data_max+histogram.slot_width())
    pyplot.xlim(-50, 150)
    pyplot.ylim(-10, 40)
    x = slot_centers
    for g in gaussians:
      print("Height %7.2f mean %4.1f sigma %3.1f"%(g.params))
      pyplot.plot(x, g(x), linewidth=2)

    if interpretation is not None:
      interpretation.plot_multiphoton_fit(pyplot)
      interpretation.plot_quality(pyplot)
    pyplot.show()

  def multiphoton_and_fit_residual(self,histogram,gaussians):

    class per_pixel_analysis:

      def __init__(OO):

        slots = flex.double(histogram.astype(np.float64))
        slot_centers = flex.double(range(self.work_params.first_slot_value,
                                      self.work_params.first_slot_value + len(histogram)))
        x = slot_centers
        y_calc = flex.double(x.size(), 0)
        for g in gaussians:
          y_calc += g(x)

        #figure a good window for plotting the residual, taken as 5 sigma away from extreme gaussian
        ceilings = [gaussians[n].params[1] + 5.*gaussians[n].params[2] for n in range(len(gaussians))]
        ceiling = max(ceilings)
        floors = [gaussians[n].params[1] - 5.*gaussians[n].params[2] for n in range(len(gaussians))]
        floor = min(floors)
        #print floors
        #print ceilings
        #print "floor",floor,"ceiling",ceiling
        #xfloor = gaussians[1].params[1] + 3.*gaussians[0].params[2]
        selection = (slot_centers<floor).__or__(slot_centers>ceiling)

        OO.fit_xresid = slot_centers.select(selection)
        OO.fit_yresid = slots.select(selection) - y_calc.select(selection)
        OO.xweight = (OO.fit_xresid - gaussians[0].params[1])/(gaussians[1].params[1] - gaussians[0].params[1])
        OO.additional_photons = flex.sum( OO.xweight * OO.fit_yresid )

        #Now the other half of the data; the part supposedly fit by the 0- and 1-photon gaussians
        OO.qual_xresid = slot_centers.select(~selection)
        ysignal = slots.select(~selection)
        OO.qual_yresid = ysignal - y_calc.select(~selection)
        OO.qual_y_fit = y_calc.select(~selection)
        # Not sure how to treat weights for channels with zero observations; default to 1
        _variance = ysignal.deep_copy().set_selected(ysignal==0., 1.)
        OO.weight = 1./_variance
        OO.weighted_numerator = OO.weight * (OO.qual_yresid * OO.qual_yresid)
        OO.sumsq_signal = flex.sum(ysignal * ysignal)
        OO.sumsq_residual = flex.sum(OO.qual_yresid * OO.qual_yresid)

      def get_multiphoton_count(OO):
        # XXX insert a test here as to whether the analysis has been carried out
        #   far enough along x-axis to capture all the high multi-photon signal
        #   if not, raise an exception
        return int(round(OO.additional_photons,0))

      def plot_multiphoton_fit(OO,plotter):
        print("counted %.0f multiphoton photons on this pixel"%OO.additional_photons)
        #plotter.plot(OO.fit_xresid, 10*OO.xweight, "b.")
        plotter.plot(OO.fit_xresid,OO.fit_yresid,"r.")

      def plot_quality(OO,plotter):
        #plotter.plot(OO.qual_xresid,OO.qual_yresid/5.,"m.")
        plotter.plot(OO.qual_xresid,OO.qual_y_fit,"k-")
        plotter.plot(OO.qual_xresid,OO.qual_yresid,"m.")
        print(OO.sumsq_signal,OO.sumsq_residual, OO.quality_factor, math.sqrt(OO.sumsq_signal))

      def chi_squared(OO):
        return flex.sum(OO.weighted_numerator)/len(OO.weighted_numerator)

    E = per_pixel_analysis()
    E.quality_factor = E.sumsq_signal/E.sumsq_residual
    return E

  def fit_3_gaussians(self,histogram):
    fitted_gaussians = []
    low_idx = self.work_params.fit_limits[0]
    high_idx = self.work_params.fit_limits[1]

    slot_centers = flex.double(range(self.work_params.first_slot_value,
                                      self.work_params.first_slot_value + len(histogram)))
    free_x = slot_centers[low_idx:high_idx]

    slots = flex.double(histogram.astype(np.float64))
    free_y = slots[low_idx:high_idx]
    total_population = flex.sum(free_y)

    # zero_mean = 0. # originally intended mean=0
    maxidx = flex.max_index(free_y) # but if dark subtraction (pedstal correction) is off
    zero_mean = free_x[maxidx] # use this non-zero maximum instead
    zero_amplitude = flex.max(free_y)
    assert 1./zero_amplitude #guard against division by zero
    zero_sigma = self.work_params.gaussian_3.zero_sigma

    inelastic_amplitude = 0.001
    elastic_amplitude = 0.001
    elastic_sigma = self.work_params.gaussian_3.zero_sigma
    #elastic_mean = zero_mean+ zero_sigma*self.work_params.gaussian_3.elastic_gain_to_sigma #**

    helper = self.helper_3_gaussian_factory(initial_estimates =
      (zero_mean, 1.0, zero_sigma, inelastic_amplitude, elastic_amplitude, elastic_sigma),
      constants=flex.double((self.work_params.gaussian_3.inelastic_gain_to_sigma,
                             self.work_params.gaussian_3.elastic_gain_to_sigma,
                             )),
      free_x = free_x,
      free_y = free_y/zero_amplitude) # put y values on 0->1 scale for normal eqn solving
    helper.restart()
    iterations = normal_eqns_solving.levenberg_marquardt_iterations(
          non_linear_ls = helper,
          n_max_iterations = 7,
          gradient_threshold = 1.E-3)

    fitted_gaussians = helper.as_gaussians()
    for item in fitted_gaussians: item.params = (item.params[0] * zero_amplitude,
                                  item.params[1], item.params[2]) # convert back to full scale
    return fitted_gaussians

  @staticmethod
  def helper_3_gaussian_factory(initial_estimates,constants,free_x,free_y):

      from xfel.vonHamos import gaussian_3fit_inheriting_from_non_linear_ls

      class helper_3_gaussian(gaussian_3fit_inheriting_from_non_linear_ls, normal_eqns.non_linear_ls_mixin):
        def __init__(pfh):
          super(helper_3_gaussian, pfh).__init__(n_parameters=6)
          pfh.x_0 = flex.double(initial_estimates)
          pfh.restart()
          pfh.set_cpp_data(constants,free_x,free_y)
          pfh.counter = 0

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
          pfh.reset()
          #rely on C++ and go directly for add_equation singular
          pfh.access_cpp_build_up_directly(objective_only, current_values = pfh.x)
          if not objective_only:
            objective = pfh.objective()
            pfh.counter+=1
            #print "iteration",pfh.counter,"Objective %.4f"%objective


        def as_gaussians(pfh):          # a: amplitude, b: mean, c: sigma
          return [curve_fitting.gaussian( a = pfh.x[1], b = pfh.x[0], c = pfh.x[2] ),
                  curve_fitting.gaussian( a = pfh.x[3], b = pfh.x[0] + pfh.x[2] * constants[0],
                                          c = (constants[0]/constants[1])*(pfh.x[5] - pfh.x[2])+pfh.x[2]  ),
                  curve_fitting.gaussian( a = pfh.x[4],
                                          b = pfh.x[0] + pfh.x[2] * constants[1],
                                          c = pfh.x[5])]

      value = helper_3_gaussian()
      return value


def hist_outline(slot_width,slots,slot_centers):

  step_size = slot_width
  half_step_size = 0.5 * step_size
  n_slots = len(slots)

  bins = flex.double(n_slots * 2 + 2, 0)
  data = flex.double(n_slots * 2 + 2, 0)
  for i in range(n_slots):
    bins[2 * i + 1] = slot_centers[i] - half_step_size
    bins[2 * i + 2] = slot_centers[i] + half_step_size
    data[2 * i + 1] = slots[i]
    data[2 * i + 2] = slots[i]

  bins[0] = bins[1] - step_size
  bins[-1] = bins[-2] + step_size
  data[0] = 0
  data[-1] = 0

  return (bins, data)

"""
start here
Suppose we take the following gaussians
0:
free amplitude, free mean, free sigma
max_y, 0., estimated sigma
1: inelastic peak
free amplitude
0.001, fixed mean, fixed sigma
2: elastic peak
free amplitude, free sigma
0.001, fixed mean, estimated_sigma
3) introduce a third gaussian
4) how good is the fit.  what constraints should be used?
"""


def test_fit(histo_data,plot=True):
  work_params = get_work_params(verbose=plot)
  pixel_histograms = derived_class(work_params)
  #alt_gaussians = pixel_histograms.fit_one_histogram_two_gaussians(histogram = histo_data)
  alt_gaussians = pixel_histograms.fit_3_gaussians(histogram = histo_data)
  #pixel_histograms.plot_combo(histo_data,alt_gaussians)
  gs = alt_gaussians[1].params
  fit_photons = gs[0] * gs[2] * math.sqrt(2.*math.pi)
  n_photons = int(round(fit_photons,0))
  gs3 = alt_gaussians[2].params
  fit_photons3 = gs3[0] * gs3[2] * math.sqrt(2.*math.pi)
  n_photons3 = int(round(fit_photons3,0))
  fit_interpretation=pixel_histograms.multiphoton_and_fit_residual(
                     histo_data, alt_gaussians)
  multi_photons = fit_interpretation.get_multiphoton_count()
  total_photons = n_photons + n_photons3
  if plot:
    print("%d single + %d elastic = %d total"%(n_photons,n_photons3, total_photons))
    pixel_histograms.plot_combo(histo_data, alt_gaussians,
                              interpretation=fit_interpretation)
  return n_photons
'''
start here fit doesn't seem to work
the particular one being illustrated here makes no sense i1==129 j1==155
1) the red curve isn't the right height or mean
2) is the objective under the purple curve consistent with the objective calculated by build up?
3) it doesn't seem like the purple curve is being minimized.
4) it seems like red curve mean should be refined
5) XXX maybe I should fit each gaussian on a region of interest +/-5 sigma
6) XXX maybe I should weight the gaussians equally instead of zero peak outweighin
'''

if __name__=="__main__":
  fname = './xppe0314_241.h5'
  f = h5py.File(fname, 'r')

  histograms = f['mydata/histograms']
  cspad = f['mydata/cspad_sum']

  print(histograms)

  image = cspad[:,:,1]
  fig=False
  if fig:
    plt.figure()
    plt.imshow(image,interpolation="none")
    plt.colorbar()
    plt.clim(0, 25000)

  i1 = 111
  j1 = 155
  i2 = 73#113
  j2 = 264#140

  histo1 = histograms[i1*388+j1,:]
  #print list(histo1)
  #exit()

  histo2 = histograms[i2*388+j2,:]
  print(histo1.sum(), histo2.sum())

  x = np.arange(-50,histo1.shape[0]-50)
  if fig:
    plt.figure()
    plt.plot(x,histo1,'r',label='pixel 1')
    plt.plot(x,histo2,'b',label='pixel 2')
    #plt.yscale('log')
    plt.xlim(-20,150)
    plt.show()

  test_fit(histo1)
  for i1 in range(80, 150):
    print("Pair i j ",i1,j1)
    histo1 = histograms[i1*388+j1,:]
    test_fit(histo1)
