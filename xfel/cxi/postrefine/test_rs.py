from __future__ import division
from cctbx.array_family import flex
from scitbx.matrix import sqr, col
import math


def get_observations(dir_name,data_subset):
  file_names = []
  for file_name in os.listdir(dir_name):
    if (file_name.endswith("_00000.pickle")):
      if data_subset==0 or \
        (data_subset==1 and (int(os.path.basename(file_name).split("_00000.pickle")[0][-1])%2==1)) or \
        (data_subset==2 and (int(os.path.basename(file_name).split("_00000.pickle")[0][-1])%2==0)):
        file_names.append(os.path.join(dir_name, file_name))
    elif (file_name.endswith(".pickle")):
      if data_subset==0 or \
        (data_subset==1 and (int(os.path.basename(file_name).split(".pickle")[0][-1])%2==1)) or \
        (data_subset==2 and (int(os.path.basename(file_name).split(".pickle")[0][-1])%2==0)):
        file_names.append(os.path.join(dir_name, file_name))
  print "Number of pickle files found:", len(file_names)
  return file_names

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
  except Exception:
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
  r_asu = sum(flex.abs((I_obs_asu*slope_asu) - I_ref_asu))/sum(slope_asu * I_obs_asu)

  matches_rev = miller.match_multi_indices(
                miller_indices_unique=miller_array_iso.indices(),
                miller_indices=miller_array_compare_rev.indices())

  I_ref_rev = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_rev.pairs()])
  I_obs_rev = flex.double([miller_array_compare_rev.data()[pair[1]] for pair in matches_rev.pairs()])

  cc_rev, slope_rev = get_overall_correlation (I_obs_rev, I_ref_rev)
  r_rev = sum(flex.abs((I_obs_rev*slope_rev) - I_ref_rev))/sum(slope_rev * I_obs_rev)

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

  """
  plt.scatter(I_ref_iso, I_obs_iso,s=10, marker='x', c='r')
  plt.title('CC=%6.5f Slope=%6.5f R=%6.5f'%(cc_iso, slope_iso, r_iso))
  plt.xlabel('Reference intensity')
  plt.ylabel('Observed intensity')
  plt.show()
  """
  return cc_iso, slope_iso, r_iso


def calc_partiality_mproc(calc_partiality_arg, frames):

  calc_partiality_arg_tag = calc_partiality_arg.split(':')
  frame_no = int(calc_partiality_arg_tag[0])
  ry_shift = float(calc_partiality_arg_tag[1])

  pickle_filename = frames[frame_no]

  trial_results = pickle.load(open(pickle_filename,"rb"))
  crystal_init_orientation = trial_results["current_orientation"][0]
  wavelength = trial_results["wavelength"]
  observations = trial_results["observations"][0]
  crystal_pointgroup = trial_results["pointgroup"]
  unit_cell = trial_results["current_orientation"][0].unit_cell()
  target_unit_cell = unit_cell.parameters()

  a_star_init = sqr(crystal_init_orientation.reciprocal_matrix())

  eh = energy_handler()
  eh.get_energy_info(file_name_in_img, file_name_in_energy, pickle_filename)

  #calc spot_radius y and z component
  spot_radius = calc_spot_radius(a_star_init, observations.indices(), wavelength)
  ry = spot_radius + ry_shift
  rz = ry * 0.5

  #calc alpha_angle
  two_theta = observations.two_theta(wavelength=wavelength).data()

  ph = partiality_handler(wavelength, spot_radius)
  I_p_aniso = flex.double()
  sigI_p_aniso = flex.double()
  I_p = flex.double()
  sigI_p = flex.double()
  for i_sel in range(len(observations.indices())):
    h = observations.indices()[i_sel]

    p, dummy = ph.calc_partiality(a_star_init, h)
    I_p.append(observations.data()[i_sel]/ p)
    sigI_p.append(observations.sigmas()[i_sel]/ p)

    alpha_angle = two_theta[i_sel]
    p_aniso, rs_aniso = ph.calc_partiality_anisotropy(a_star_init, h, ry, rz, alpha_angle)
    I_p_aniso.append(observations.data()[i_sel]/ p_aniso)
    sigI_p_aniso.append(observations.sigmas()[i_sel]/ p_aniso)



  observations_p = observations.customized_copy(data=I_p, sigmas=sigI_p)
  observations_p_aniso = observations.customized_copy(data=I_p_aniso, sigmas=sigI_p_aniso)

  return frame_no, observations_p, observations_p_aniso, (ry, rz, spot_radius)


if __name__=="__main__":

  from iotbx import reflection_file_reader
  from cctbx.array_family import flex
  from cctbx import miller
  from cctbx.crystal import symmetry
  from scitbx.matrix import col, sqr
  from libtbx.easy_mp import pool_map
  import os,cPickle as pickle,math
  from mod_partiality import partiality_handler
  from mod_energy import energy_handler
  import matplotlib.pyplot as plt
  import numpy as np



  file_name_in_energy = 'first_energy.csv'
  file_name_in_img = 'serial_img_trial_015.out'
  pickle_dir = '/net/viper/raid1/mu238/L614/results/myo_stills/015/integration/'

  file_name_iso_mtz = '3u3e-sf.mtz'
  reflection_file_iso = reflection_file_reader.any_reflection_file(file_name_iso_mtz)
  miller_arrays_iso=reflection_file_iso.as_miller_arrays()
  miller_array_iso = miller_arrays_iso[1].generate_bijvoet_mates()

  #for each frame calculate the corrected intensities using partiality
  frames = range(10,11)
  frame_files = get_observations(pickle_dir, 0)

  ry_shift = 0.0001
  ry_set = np.arange(-0.005,0.005, ry_shift)

  calc_partiality_arg = []
  for fr in frames:
    for ry in ry_set:
       calc_partiality_arg.append(str(fr)+':'+str(ry))

  def calc_partiality_mproc_wrapper(arg):
    return calc_partiality_mproc(arg,
        frame_files)

  calc_partiality_result = pool_map(
          args=calc_partiality_arg,
          func=calc_partiality_mproc_wrapper,
          processes=None)


  for result in calc_partiality_result:
    if result is not None:
      frame_no = result[0]
      observations_p = result[1]
      observations_p_aniso = result[2]
      ry, rz, spot_radius = result[3]

      d_min = min(observations_p.d_spacings().data())
      d_max = max(observations_p.d_spacings().data())

      #calculate cc_iso, r_iso for I_p_single and I_p_multi
      cc_iso_p, slope_iso_p, r_iso_p = compare_quality(miller_array_iso, observations_p)
      cc_iso_p_aniso, slope_iso_p_aniso, r_iso_p_aniso = compare_quality(miller_array_iso, observations_p_aniso)

      print '%03d %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f' %(frame_no, d_min, d_max, ry, rz, spot_radius, cc_iso_p, cc_iso_p_aniso, r_iso_p, r_iso_p_aniso)
