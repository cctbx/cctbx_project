from __future__ import division


def calc_spot_radius(a_star_matrix, miller_indices, wavelength):
  #calculate spot_radius based on rms delta_S for all spots
  delta_S_all = flex.double()
  for miller_index in miller_indices:
    S0 = -1*col((0,0,1./wavelength))
    h = col(miller_index)
    x = a_star_matrix * h
    S = x + S0
    delta_S = S.length() - (1./wavelength)
    delta_S_all.append(delta_S)

  spot_radius = math.sqrt(flex.mean(delta_S_all*delta_S_all))

  return spot_radius

def get_overall_correlation (data_a, data_b) :
  """
  Correlate the averaged intensities to the intensities from the
  reference data set.
  """

  assert len(data_a) == len(data_b)
  corr = 0
  slope = 0
  try:
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
    sum_x  = 0
    sum_y  = 0
    N      = 0
    for i in xrange(len(data_a)):
      I_r       = data_a[i]
      I_o       = data_b[i]
      N      += 1
      sum_xx += I_r**2
      sum_yy += I_o**2
      sum_xy += I_r * I_o
      sum_x  += I_r
      sum_y  += I_o
    slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
    corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
               math.sqrt(N * sum_yy - sum_y**2))
  except:
    pass

  return corr, slope

def compare_quality(miller_array_iso, miller_array_compare):
  target_unit_cell = (90.309,90.309,45.204,90,90,120)
  target_space_group = 'P6'
  target_anomalous_flag = False
  d_min = 1.5

  miller_set = symmetry(
      unit_cell=target_unit_cell,
      space_group_symbol=target_space_group
    ).build_miller_set(
      anomalous_flag=target_anomalous_flag,
      d_min=d_min)

  #filter negative intensities
  i_I_positive = (miller_array_compare.data() > 0)
  miller_indices_positive = miller_array_compare.indices().select(i_I_positive)
  I_positive = miller_array_compare.data().select(i_I_positive)
  sigI_positive = miller_array_compare.sigmas().select(i_I_positive)

  miller_array_compare = miller_array_compare.customized_copy(indices=miller_indices_positive,
      data=I_positive,
      sigmas=sigI_positive,
      anomalous_flag=target_anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )


  #Prepare original, asu, and rev (for polar case)
  miller_array_compare_asu = miller_array_compare.map_to_asu()

  from cctbx import sgtbx
  cb_op = sgtbx.change_of_basis_op('k,h,-l')
  miller_array_compare_rev = miller_array_compare_asu.change_basis(cb_op).map_to_asu()

  matches_asu = miller.match_multi_indices(
                miller_indices_unique=miller_array_iso.indices(),
                miller_indices=miller_array_compare_asu.indices())

  I_ref_asu = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_asu.pairs()])
  I_obs_asu = flex.double([miller_array_compare_asu.data()[pair[1]] for pair in matches_asu.pairs()])

  cc_asu, slope_asu = get_overall_correlation (I_obs_asu, I_ref_asu)
  r_asu = sum(flex.abs(I_ref_asu - (slope_asu*I_obs_asu)))/sum(slope_asu * I_obs_asu)

  matches_rev = miller.match_multi_indices(
                miller_indices_unique=miller_array_iso.indices(),
                miller_indices=miller_array_compare_rev.indices())

  I_ref_rev = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_rev.pairs()])
  I_obs_rev = flex.double([miller_array_compare_rev.data()[pair[1]] for pair in matches_rev.pairs()])

  cc_rev, slope_rev = get_overall_correlation (I_obs_rev, I_ref_rev)
  r_rev = sum(flex.abs(I_ref_rev - (slope_rev*I_obs_rev)))/sum(slope_rev * I_obs_rev)

  cc_iso = cc_asu
  slope_iso = slope_asu
  r_iso = r_asu
  I_ref_iso = I_ref_asu
  I_obs_iso = I_obs_asu
  if cc_rev > cc_iso:
    cc_iso = cc_rev
    slope_iso = slope_rev
    r_iso = r_rev
    I_ref_iso = I_ref_rev
    I_obs_iso = I_obs_rev

  plt.scatter(I_ref_iso, I_obs_iso,s=10, marker='x', c='r')
  plt.title('CC=%6.5f Slope=%6.5f R=%6.5f'%(cc_iso, slope_iso, r_iso))
  plt.xlabel('Reference intensity')
  plt.ylabel('Observed intensity')
  plt.show()


def calc_partiality_mproc(i_sel, observations, crystal_init_orientation):

  h = observations.indices()[i_sel]
  partiality_c_set = flex.double()

  a_star_init = sqr(crystal_init_orientation.reciprocal_matrix())

  eh = energy_handler()
  eh.get_energy_info(file_name_in_img, file_name_in_energy, pickle_filename)


  #calculate spot_radisu (rs) from mean_wavelength
  spot_radius = calc_spot_radius(a_star_init, observations.indices(), eh.mean_wavelength)

  for wavelength, c_weight in zip(eh.wavelength_at_counts, eh.photon_counts_normalized):

    ph = partiality_handler(wavelength, spot_radius)
    p, dummy = ph.calc_partiality(a_star_init, h)

    partiality_c_set.append(p*c_weight)

  #correct intensity usinging single color
  ph = partiality_handler(eh.mean_wavelength, spot_radius)
  p_single, dummy = ph.calc_partiality(a_star_init, h)
  I_p_single = observations.data()[i_sel]/p_single
  sigI_p_single = observations.sigmas()[i_sel]/p_single

  p_multi = sum(partiality_c_set)/sum(eh.photon_counts_normalized)
  I_p_multi = observations.data()[i_sel]/p_multi
  sigI_p_multi = observations.sigmas()[i_sel]/p_multi

  return i_sel, h, p_single, I_p_single, sigI_p_single, p_multi, I_p_multi, sigI_p_multi


if __name__=="__main__":

  from iotbx import reflection_file_reader
  from cctbx.array_family import flex
  from cctbx import miller
  from cctbx import crystal
  from cctbx.crystal import symmetry
  from scitbx.matrix import col, sqr
  from libtbx.easy_mp import pool_map, get_processes
  import os,cPickle as pickle,math
  from mod_partiality import partiality_handler
  from mod_energy import energy_handler
  import matplotlib.pyplot as plt
  import numpy as np

  pickle_filename ='/net/viper/raid1/mu238/L614/results/myo_stills/015/integration/int-data_00198.pickle'

  trial_results = pickle.load(open(pickle_filename,"rb"))
  crystal_init_orientation = trial_results["current_orientation"][0]
  wavelength = trial_results["wavelength"]
  observations = trial_results["observations"][0]
  crystal_pointgroup = trial_results["pointgroup"]
  unit_cell = trial_results["current_orientation"][0].unit_cell()
  target_unit_cell = unit_cell.parameters()

  file_name_in_energy = 'first_energy.csv'
  file_name_in_img = 'serial_img_trial_015.out'



  #for all miller indices calculate the integral of partiality
  i_set = range(0, len(observations.indices()))

  def calc_partiality_mproc_wrapper(arg):
    return calc_partiality_mproc(arg,
        observations,
        crystal_init_orientation)

  calc_partiality_result = pool_map(
          args=i_set,
          func=calc_partiality_mproc_wrapper,
          processes=None)

  I_p_single = flex.double([0]*len(observations.indices()))
  sigI_p_single = flex.double([0]*len(observations.indices()))
  I_p_multi = flex.double([0]*len(observations.indices()))
  sigI_p_multi = flex.double([0]*len(observations.indices()))
  for result in calc_partiality_result:
    if result is not None:
      i_sel, h, p_single, i_p_single, sigi_p_single, p_multi, i_p_multi, sigi_p_multi = result
      I_p_single[i_sel] = i_p_single
      sigI_p_single[i_sel] = sigi_p_single
      I_p_multi[i_sel] = i_p_multi
      sigI_p_multi[i_sel] = sigi_p_multi
      print i_sel, h, 'p single', p_single, 'p multi', p_multi

  observations_p_single = observations.customized_copy(data=I_p_single,
      sigmas=sigI_p_single)

  observations_p_multi = observations.customized_copy(data=I_p_multi,
      sigmas=sigI_p_multi)


  #calculate cc_iso, r_iso for I_p_single and I_p_multi
  file_name_iso_mtz = '3u3e-sf.mtz'
  reflection_file_iso = reflection_file_reader.any_reflection_file(file_name_iso_mtz)
  miller_arrays_iso=reflection_file_iso.as_miller_arrays()
  miller_array_iso = miller_arrays_iso[1].generate_bijvoet_mates()

  compare_quality(miller_array_iso, observations)
  compare_quality(miller_array_iso, observations_p_single)
  compare_quality(miller_array_iso, observations_p_multi)
