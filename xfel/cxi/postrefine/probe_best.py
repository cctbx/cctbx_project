from __future__ import division
#probe.py reads in mtz file and build the reciprocal lattices
#from pickle file

def get_observations (dir_name,data_subset):
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


def calc_cc_mproc(frame_no,
      frame_files,
      miller_array_iso):

  pickle_filename = frame_files[frame_no]

  trial_results = pickle.load(open(pickle_filename,"rb"))
  observations = trial_results["observations"][0]

  #for observations, only filter out those refletions with I/sigI < 2.0
  sigma_min = 2.0
  I_over_sigi = observations.data()/ observations.sigmas()
  i_I_obs_sel = (I_over_sigi > sigma_min)

  observations_asu = observations.customized_copy(indices=observations.indices().select(i_I_obs_sel),
      data=observations.data().select(i_I_obs_sel),
      sigmas=observations.sigmas().select(i_I_obs_sel)).map_to_asu()

  if flag_polar:
    from cctbx import sgtbx
    cb_op = sgtbx.change_of_basis_op('k,h,-l')
    observations_rev = observations_asu.change_basis(cb_op).map_to_asu()
  else:
    observations_rev = observations_asu

  matches_asu = miller.match_multi_indices(
                    miller_indices_unique=miller_array_iso.indices(),
                    miller_indices=observations_asu.indices())
  I_iso_asu = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_asu.pairs()])
  I_obs_asu = flex.double([observations_asu.data()[pair[1]] for pair in matches_asu.pairs()])
  cc_asu, slope_asu = get_overall_correlation(I_iso_asu, I_obs_asu)

  matches_rev = miller.match_multi_indices(
                    miller_indices_unique=miller_array_iso.indices(),
                    miller_indices=observations_rev.indices())
  I_iso_rev = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_rev.pairs()])
  I_obs_rev = flex.double([observations_rev.data()[pair[1]] for pair in matches_rev.pairs()])
  cc_rev, slope_rev = get_overall_correlation(I_iso_rev, I_obs_rev)

  cc_all = cc_asu
  observations_sel = observations_asu
  if cc_rev > cc_asu:
    cc_all = cc_rev
    observations_sel = observations_rev


  return frame_no, observations_sel.d_min(), len(observations_sel.indices()), cc_all

if __name__=="__main__":

  from iotbx import reflection_file_reader
  import os,cPickle as pickle,math
  from scitbx.matrix import sqr, col
  import numpy as np
  from libtbx.easy_mp import pool_map, get_processes
  from cctbx.array_family import flex
  from cctbx import miller
  import matplotlib.pyplot as plt
  import matplotlib.mlab as mlab

  flag_polar = True

  #read mtz file and expand to P1
  file_name_iso_mtz = '3u3e-sf.mtz'
  reflection_file_iso = reflection_file_reader.any_reflection_file(file_name_iso_mtz)
  miller_arrays_iso=reflection_file_iso.as_miller_arrays()
  miller_array_iso = miller_arrays_iso[1].expand_to_p1().generate_bijvoet_mates()

  #read pickle file (observations and crystal orientation)
  #pickle_dir = '/net/cci/mu238/XfelProject/t4_base/data/'
  pickle_dir = '/net/viper/raid1/mu238/L614/results/myo_stills/015/integration/'
  frame_files = get_observations(pickle_dir, 0)
  frames = range(0, len(frame_files))

  n_bins = 25

  #for each frame find out true reflections and observed reflections
  def calc_cc_mproc_wrapper(arg):
    return calc_cc_mproc(arg,
        frame_files,
        miller_array_iso)

  probe_complete_result = pool_map(
          args=frames,
          func=calc_cc_mproc_wrapper,
          processes=None)

  frame_no_all = flex.int()
  d_min_obs_all = flex.double()
  sum_refl_obs_all = flex.int()
  cc_all = flex.double()
  sort_weight_all = flex.double()
  for result_row in probe_complete_result:
    if result_row is not None:
      frame_no = result_row[0]
      d_min_obs = result_row[1]
      sum_refl_obs = result_row[2]
      cc = result_row[3]

      frame_no_all.append(frame_no)
      d_min_obs_all.append(d_min_obs)
      sum_refl_obs_all.append(sum_refl_obs)
      cc_all.append(cc)
      sort_weight_all.append(cc*sum_refl_obs)

      #print frame_no, d_min_iso, d_min_obs, sum_refl_iso, sum_refl_obs, cc
      #print txt_out


  perm = flex.sort_permutation(sort_weight_all, reverse = True)
  frame_no_all_sort = frame_no_all.select(perm)
  d_min_obs_all_sort = d_min_obs_all.select(perm)
  sum_refl_obs_all_sort = sum_refl_obs_all.select(perm)
  cc_all_sort = cc_all.select(perm)

  n_show_top = 10
  cn_n_show = 0
  for frame_no, d_min_obs, sum_obs, cc in zip(frame_no_all_sort,
        d_min_obs_all_sort, sum_refl_obs_all_sort, cc_all_sort):
    print frame_no, d_min_obs, sum_obs, cc
    cn_n_show += 1
    if cn_n_show == n_show_top:
      break


  print list(frame_no_all_sort[0:100])
