from __future__ import division
from __future__ import print_function
from six.moves import range
import math
import sys
from dials.array_family import flex
from six.moves import cPickle as pickle
from rstbx.dials_core.integration_core import show_observations

class reduction(object):
  """Reduced data container.  Merges the following concepts:

  filename : the integration pickle
  experiment : the dxtbx-experimental model
  HKL : the original-index miller set taken from joined database
  """
  def __init__(self, filename, experiment, HKL, i_sigi, measurements, params):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())
    self.stash_type = None
    self.stash_res_filter = None

    from dxtbx.model import DetectorFactory
    self.dummy_detector = DetectorFactory.simple(
      sensor = DetectorFactory.sensor("PAD"),
      distance = 100,
      beam_centre = [1000, 1000],
      fast_direction = "+x",
      slow_direction = "+y",
      pixel_size = [0.2,0.2],
      image_size = [2000,2000],
      )
    # simple view of post-integration, no longer need to know detector

  def get_delta_psi_deg(self):
    from scitbx import matrix
    if self.stash_type is None:
       self.stash_type = self.experiment.crystal.get_space_group().type()
    UC = self.experiment.crystal.get_unit_cell()
    from dials.algorithms.spot_prediction \
      import StillsDeltaPsiReflectionPredictor
    S = StillsDeltaPsiReflectionPredictor(
      self.experiment.beam, \
      self.dummy_detector, \
      self.experiment.crystal.get_A(), \
      UC, \
      self.stash_type, \
      10.0) #dummy value for dmin

    length = len(self.HKL)
    R= flex.reflection_table.empty_standard(length)
    R['miller_index'] = self.HKL
    S.for_reflection_table(R,matrix.sqr(self.experiment.crystal.get_A()))
    degrees = (180./math.pi)*R["delpsical.rad"]
    return degrees

  def get_two_theta_deg(self):
    wavelength=self.experiment.beam.get_wavelength()
    UC = self.experiment.crystal.get_unit_cell()
    two_theta = UC.two_theta(
      miller_indices=self.HKL, wavelength=wavelength, deg=True)
    return two_theta

  def get_imposed_res_filter(self, out):
    if self.stash_res_filter is not None:  return self.stash_res_filter
    if self.params.significance_filter.apply is True: #------------------------------------

      print("Step 5. Frame by frame resolution filter", file=out)
      # Apply an I/sigma filter ... accept resolution bins only if they
      #   have significant signal; tends to screen out higher resolution observations
      #   if the integration model doesn't quite fit
      N_obs_pre_filter = self.i_sigi.size()
      N_bins_small_set = N_obs_pre_filter // self.params.significance_filter.min_ct
      N_bins_large_set = N_obs_pre_filter // self.params.significance_filter.max_ct

      # Ensure there is at least one bin.
      N_bins = max(
        [min([self.params.significance_filter.n_bins,N_bins_small_set]),
         N_bins_large_set, 1]
      )
      print("Total obs %d Choose n bins = %d"%(N_obs_pre_filter,N_bins), file=out)
      bin_results = show_observations(self.measurements, out=sys.stdout, n_bins=N_bins)

      if True: # no fuller kapton -- not implemented here,
               # but code can and should be borrowed from cxi.merge
        acceptable_resolution_bins = [
          bin.mean_I_sigI > self.params.significance_filter.sigma for bin in bin_results]
        acceptable_nested_bin_sequences = [i for i in range(len(acceptable_resolution_bins))
                                           if False not in acceptable_resolution_bins[:i+1]]
        N_acceptable_bins = max(acceptable_nested_bin_sequences) + 1
        imposed_res_filter = float(bin_results[N_acceptable_bins-1].d_range.split()[2])
        print("New resolution filter at %7.2f"%imposed_res_filter,self.filename, file=out)
        print("N acceptable bins",N_acceptable_bins, file=out)
      print("Old n_obs: %d, new n_obs: %d"%(N_obs_pre_filter,self.measurements.size()), file=out)
      # Finished applying the binwise I/sigma filter---------------------------------------
    else:
      imposed_res_filter=None
    self.stash_res_filter = imposed_res_filter
    return imposed_res_filter
