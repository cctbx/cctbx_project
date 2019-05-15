from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import miller
import math
from six.moves import range

class merge_data_handler(object):
  """keep arrays of unmerged data"""
  def __init__(self):
    #extract data
    self.miller_indices_merge = flex.miller_index()
    self.I_merge = flex.double()
    self.sigI_merge = flex.double()
    self.r_meas_div = flex.double()
    self.r_meas_divisor = flex.double()
    self.multiplicities = flex.int()
    self.I_even = flex.double()
    self.I_odd  = flex.double()
    self.I_even_h = flex.double()
    self.I_odd_h = flex.double()
    self.I_even_k = flex.double()
    self.I_odd_k = flex.double()
    self.I_even_l = flex.double()
    self.I_odd_l = flex.double()
    self.miller_array_merge = None

  def extend_data(self, miller_indices_merge, I_merge, sigI_merge, stat_all,
                   I_two_halves_tuple, uc_mean, wavelength_mean):
    self.miller_indices_merge.extend(miller_indices_merge)
    self.I_merge.extend(I_merge)
    self.sigI_merge.extend(sigI_merge)
    r_meas_div, r_meas_divisor, multiplicities = stat_all
    self.r_meas_div.extend(r_meas_div)
    self.r_meas_divisor.extend(r_meas_divisor)
    self.multiplicities.extend(multiplicities)
    I_even, I_odd, I_even_h, I_odd_h, I_even_k, I_odd_k, I_even_l, I_odd_l = I_two_halves_tuple
    self.I_even.extend(I_even)
    self.I_odd.extend(I_odd)
    self.I_even_h.extend(I_even_h)
    self.I_odd_h.extend(I_odd_h)
    self.I_even_k.extend(I_even_k)
    self.I_odd_k.extend(I_odd_k)
    self.I_even_l.extend(I_even_l)
    self.I_odd_l.extend(I_odd_l)
    self.uc_mean = uc_mean
    self.wavelength_mean = wavelength_mean

  def extend(self, mdh):
    self.miller_indices_merge.extend(mdh.miller_indices_merge)
    self.I_merge.extend(mdh.I_merge)
    self.sigI_merge.extend(mdh.sigI_merge)
    self.r_meas_div.extend(mdh.r_meas_div)
    self.r_meas_divisor.extend(mdh.r_meas_divisor)
    self.multiplicities.extend(mdh.multiplicities)
    self.I_even.extend(mdh.I_even)
    self.I_odd.extend(mdh.I_odd)
    self.I_even_h.extend(mdh.I_even_h)
    self.I_odd_h.extend(mdh.I_odd_h)
    self.I_even_k.extend(mdh.I_even_k)
    self.I_odd_k.extend(mdh.I_odd_k)
    self.I_even_l.extend(mdh.I_even_l)
    self.I_odd_l.extend(mdh.I_odd_l)
    self.uc_mean = mdh.uc_mean
    self.wavelength_mean = mdh.wavelength_mean

  def generate_miller_array_from_miller_set(self, miller_set, target_anomalous_flag):
    self.miller_array_merge = miller_set.array() \
        .customized_copy(indices=self.miller_indices_merge, \
            data=self.I_merge, sigmas=self.sigI_merge, anomalous_flag=target_anomalous_flag) \
        .map_to_asu() \
        .set_observation_type_xray_intensity()
    return self.miller_array_merge

  def get_size(self):
    return len(self.I_merge)

  def get_cc_anom(self):
    cc_anom_acentric, nrefl_anom_acentric = (0,0)
    if self.miller_array_merge.anomalous_flag():
      ma_anom_dif_even = self.miller_array_merge.customized_copy(data = self.I_even).anomalous_differences()
      ma_anom_dif_odd = self.miller_array_merge.customized_copy(data = self.I_odd).anomalous_differences()
      i_acentric = (ma_anom_dif_even.centric_flags().data() == False)
      cc_anom_acentric = flex.linear_correlation(ma_anom_dif_even.data().select(i_acentric), ma_anom_dif_odd.data().select(i_acentric)).coefficient()
      nrefl_anom_acentric = len(ma_anom_dif_even.data().select(i_acentric))
    return cc_anom_acentric, nrefl_anom_acentric

  def get_multiplicity(self):
    return flex.sum(self.multiplicities)/self.get_size() if self.get_size() else 0

  def get_r_meas(self):
    return flex.sum(self.r_meas_div)/ flex.sum(self.r_meas_divisor) if flex.sum(self.r_meas_divisor) else 0

  def get_r_split(self):
    try:
      r_split_bin = (1/math.sqrt(2))*(flex.sum(flex.abs(self.I_even - self.I_odd))/(flex.sum(self.I_even + self.I_odd)*0.5))
    except Exception as e:
      print("Warning: R_split calculation failed.")
      print(e)
      r_split_bin =0
    return r_split_bin if self.get_size() else 0

  def get_cc12(self):
    cc12 = flex.linear_correlation(self.I_even, self.I_odd).coefficient()
    return (cc12, self.get_size()) if self.get_size() else (0,0)

  def get_cciso(self, miller_array_iso):
    cciso, n_refl_cciso = (0,0)
    if miller_array_iso:
      matches_iso = miller.match_multi_indices(
          miller_indices_unique=miller_array_iso.indices(),
          miller_indices=self.miller_array_merge.indices())
      I_iso = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_iso.pairs()])
      I_merge_match_iso = flex.double([self.I_merge[pair[1]] for pair in matches_iso.pairs()])
      n_refl_cciso = len(matches_iso.pairs())
      if len(matches_iso.pairs()) > 0 : cciso = flex.linear_correlation(I_merge_match_iso, I_iso).coefficient()
    return cciso, n_refl_cciso

  def get_mean_IoversigI(self):
    return flex.mean(self.I_merge/self.sigI_merge) if self.get_size() else 0

  def get_mean_I(self):
    return flex.mean(self.I_merge) if self.get_size() else 0

  def get_mean_sigI(self):
    return flex.mean(self.sigI_merge) if self.get_size() else 0

  def get_second_moment(self):
    second_moment = 0
    if self.get_size():
      i_acentric = (self.miller_array_merge.centric_flags().data() == False)
      I_acentric = self.I_merge.select(i_acentric)
      second_moment = flex.mean(I_acentric**2)/(flex.mean(I_acentric)**2)
    return second_moment

  def reduce_by_selection(self, selections):
    """from selections (flex.bool), filter all class attributes"""
    self.miller_indices_merge = self.miller_indices_merge.select(selections)
    self.I_merge = self.I_merge.select(selections)
    self.sigI_merge = self.sigI_merge.select(selections)
    self.r_meas_div = self.r_meas_div.select(selections)
    self.r_meas_divisor = self.r_meas_divisor.select(selections)
    self.multiplicities = self.multiplicities.select(selections)
    self.I_even = self.I_even.select(selections)
    self.I_odd = self.I_odd.select(selections)
    self.I_even_h = self.I_even_h.select(selections)
    self.I_odd_h = self.I_odd_h.select(selections)
    self.I_even_k = self.I_even_k.select(selections)
    self.I_odd_k = self.I_odd_k.select(selections)
    self.I_even_l = self.I_even_l.select(selections)
    self.I_odd_l = self.I_odd_l.select(selections)
    if self.miller_array_merge:
      self.miller_array_merge = self.miller_array_merge.select(selections)

  def reduce_by_miller_index(self, miller_indices):
    sequences = list(range(self.get_size()))
    matches = miller.match_multi_indices(
                  miller_indices_unique=miller_indices,
                  miller_indices=self.miller_indices_merge)
    pair_1 = flex.int([pair[1] for pair in matches.pairs()])
    sequences_bin = flex.size_t([sequences[pair_1[j]] for j in range(len(matches.pairs()))])
    self.reduce_by_selection(sequences_bin)

  def reduce_to_cone_on_axis(self, axis_point_2, fraction_percent):
    miller_array_reduced = self.miller_array_merge.remove_cone(fraction_percent, axis_point_2=axis_point_2, negate=True)
    self.reduce_by_miller_index(miller_array_reduced.indices())
