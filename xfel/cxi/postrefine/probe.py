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

def compute_detector_d_min(d_distance, dsx, dsy, wavelength):
  half_diag = (dsx**2 + dsy**2)**0.5 / 2
  import math
  theta_rad = math.atan2(half_diag, d_distance) / 2
  assert theta_rad != 0
  denom = 2 * math.sin(theta_rad)
  assert denom != 0
  return round(wavelength / denom, 2)

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
  
def determine_reflection_visible(miller_array_raw, 
      wavelength, 
      crystal_init_orientation, 
      I_over_sigI_min=0,
      x_len_lim = 99,
      ewald_proximity=0.001):
  
  #filter reflections by given I_over_sigI_min
  I_over_sigI = miller_array_raw.data()/ miller_array_raw.sigmas()
  i_sel = (I_over_sigI > I_over_sigI_min)
  miller_array_raw = miller_array_raw.customized_copy(indices=miller_array_raw.indices().select(i_sel),
    data=miller_array_raw.data().select(i_sel),
    sigmas=miller_array_raw.data().select(i_sel))
    
  #determine rh (offset from Ewald sphere)
  S0 = -1*col((0,0,1/wavelength))
  a_star_init = sqr(crystal_init_orientation.reciprocal_matrix())
  ewald_offset_list = flex.double()
  
  x_len_list = flex.double()
  miller_indices_visible = flex.miller_index()
  I_visible = flex.double()
  sigI_visible = flex.double()
  for i in range(len(miller_array_raw.indices())): 
    miller_index = miller_array_raw.indices()[i]
    h = col(miller_index)
    x = a_star_init * h
    S = x + S0
    delta_S = S.length() - (1/wavelength)
    
    x_len_list.append(x.length())
    ewald_offset_list.append(delta_S)
    #a reflection is on this frame when
    #1. Length of x is less than certain limit (says 0.1)
    #2. delta_S is smaller than ewald_proximity
    if x.length() < x_len_lim and abs(delta_S) < ewald_proximity:
      miller_indices_visible.append(miller_index)
      I_visible.append(miller_array_raw.data()[i])
      sigI_visible.append(miller_array_raw.sigmas()[i])
      
  
  """
  x = x_len_list.as_numpy_array()
  mu = np.mean(x) 
  med = np.median(x)
  sigma = np.std(x)
  num_bins = 50
  n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
  y = mlab.normpdf(bins, mu, sigma)
  plt.plot(bins, y, 'r--')
  plt.ylabel('Frequencies')
  plt.title('X length mean=%3.2f median=%3.2f sigma=%3.2f'%(mu,med,sigma))
  plt.show()
  """
  max_x_len = max(x_len_list)  
  max_abs_ewald_offset = max(flex.abs(ewald_offset_list))
  
  miller_array_visible = miller_array_raw.customized_copy(indices=miller_indices_visible,
    data=I_visible,
    sigmas=sigI_visible)
  
  return miller_array_visible, max_x_len, max_abs_ewald_offset
  
def probe_complete_mproc(frame_no,
      frame_files, 
      miller_array_iso):
  print 'iso_raw', len(miller_array_iso.indices())
  
  pickle_filename = frame_files[frame_no]
  
  trial_results = pickle.load(open(pickle_filename,"rb"))
  observations = trial_results["observations"][0]
  wavelength = trial_results["wavelength"]
  crystal_init_orientation = trial_results["current_orientation"][0]
  
  d_min_visible = compute_detector_d_min(detector_distance, detector_size_x, detector_size_y, wavelength)
  print 'From detector setup and lambda, d_min_visible=', d_min_visible
  
  miller_array_iso = miller_array_iso.resolution_filter(d_min=d_min_visible, d_max=45)
  print 'iso filter_res', len(miller_array_iso.indices())
  
  #use observations to determine max_x_len and max_ewald_offset
  observations_visible, max_x_len, max_abs_ewald_offset = determine_reflection_visible(observations, wavelength, crystal_init_orientation)
  print 'From observations, max_x_len=', max_x_len, 'max_abs_ewald_offset=', max_abs_ewald_offset
  
  #for iso set, determine visible observations
  miller_array_iso_visible, dummy, dummy = determine_reflection_visible(miller_array_iso, wavelength, crystal_init_orientation, 
      I_over_sigI_min=0, x_len_lim = 0.9, ewald_proximity=0.005)
  print 'iso visible', len(miller_array_iso_visible.indices())
  
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
                    miller_indices_unique=miller_array_iso_visible.indices(),
                    miller_indices=observations_asu.indices())
  I_iso_asu = flex.double([miller_array_iso_visible.data()[pair[0]] for pair in matches_asu.pairs()])
  I_obs_asu = flex.double([observations_asu.data()[pair[1]] for pair in matches_asu.pairs()])
  cc_asu, slope_asu = get_overall_correlation(I_iso_asu, I_obs_asu)
  
  matches_rev = miller.match_multi_indices(
                    miller_indices_unique=miller_array_iso_visible.indices(),
                    miller_indices=observations_rev.indices())
  I_iso_rev = flex.double([miller_array_iso_visible.data()[pair[0]] for pair in matches_rev.pairs()])
  I_obs_rev = flex.double([observations_rev.data()[pair[1]] for pair in matches_rev.pairs()])
  cc_rev, slope_rev = get_overall_correlation(I_iso_rev, I_obs_rev)
  
  cc_all = cc_asu
  observations_sel = observations_asu
  if cc_rev > cc_asu:
    cc_all = cc_rev
    observations_sel = observations_rev
  print len(matches_asu.pairs()), len(matches_rev.pairs())
  #bin the reflection set and report no. of reflections by resolutions
  binner_merge = miller_array_iso_visible.setup_binner(n_bins=n_bins)
  binner_merge_indices = binner_merge.bin_indices()
  sum_refl_iso = 0
  sum_refl_obs = 0
  txt_out = '#Bin Res   #exp   #obs  <I_iso> <I_obs> CCiso\n'
  for i in range(1,n_bins+1):
    i_binner = (binner_merge_indices == i)
    miller_indices_iso_bin = miller_array_iso_visible.indices().select(i_binner)
    matches_bin = miller.match_multi_indices(
                    miller_indices_unique=miller_indices_iso_bin,
                    miller_indices=observations_sel.indices())
    I_iso_bin = flex.double([miller_array_iso_visible.data().select(i_binner)[pair[0]] for pair in matches_bin.pairs()])
    I_obs_bin = flex.double([observations_sel.data()[pair[1]] for pair in matches_bin.pairs()])
    cc_bin, slope_bin = get_overall_correlation(I_iso_bin, I_obs_bin)
    mean_I_iso_bin = 0
    mean_I_obs_bin = 0
    if len(matches_bin.pairs()) > 0:
      mean_I_iso_bin = np.mean(I_iso_bin)
      mean_I_obs_bin = np.mean(I_obs_bin)
      
    sum_refl_iso += len(miller_indices_iso_bin)
    sum_refl_obs += len(matches_bin.pairs())
    
    txt_out += '%02d %5.2f - %5.2f %5.0f %5.0f %7.2f %7.2f %0.4f' \
          %(i, binner_merge.bin_d_range(i)[0], binner_merge.bin_d_range(i)[1], 
          len(miller_indices_iso_bin), len(matches_bin.pairs()), mean_I_iso_bin, mean_I_obs_bin, cc_bin)
    txt_out += '\n'
  
  return frame_no, miller_array_iso_visible.d_min(), observations_sel.d_min(), sum_refl_iso, sum_refl_obs, cc_all, txt_out
  
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
  
  #from detector distance and size, determine at maximum resolution observable
  detector_distance = 116
  detector_size_x = 325.001216
  detector_size_y = 325.001216
  
  #read pickle file (observations and crystal orientation)
  #pickle_dir = '/net/cci/mu238/XfelProject/t4_base/data/'
  pickle_dir = '/net/viper/raid1/mu238/L614/results/myo_stills/015/integration/'
  frame_files = get_observations(pickle_dir, 0)
  #frames = range(0, len(frame_files))
  frames = range(1,2)
  
  n_bins = 25
  
  #for each frame find out true reflections and observed reflections
  def probe_complete_mproc_wrapper(arg):
    return probe_complete_mproc(arg, 
        frame_files,
        miller_array_iso)
    
  probe_complete_result = pool_map(
          args=frames,
          func=probe_complete_mproc_wrapper,
          processes=None)
  
  frame_no_all = flex.int()
  d_min_iso_all = flex.double()
  d_min_obs_all = flex.double()
  sum_refl_iso_all = flex.int()
  sum_refl_obs_all = flex.int()
  cc_all = flex.double()     
  sort_weight_all = flex.double()   
  for result_row in probe_complete_result:
    if result_row is not None:
      frame_no = result_row[0]
      d_min_iso = result_row[1]
      d_min_obs = result_row[2]
      sum_refl_iso = result_row[3]
      sum_refl_obs = result_row[4]
      cc = result_row[5]
      txt_out = result_row[6]
      
      frame_no_all.append(frame_no)
      d_min_iso_all.append(d_min_iso)
      d_min_obs_all.append(d_min_obs)
      sum_refl_iso_all.append(sum_refl_iso)
      sum_refl_obs_all.append(sum_refl_obs)
      cc_all.append(cc)
      sort_weight_all.append(cc*sum_refl_obs) 
     
      #print frame_no, d_min_iso, d_min_obs, sum_refl_iso, sum_refl_obs, cc
      #print txt_out
  
  
  perm = flex.sort_permutation(sort_weight_all, reverse = True)
  frame_no_all_sort = frame_no_all.select(perm)
  d_min_iso_all_sort = d_min_iso_all.select(perm)
  d_min_obs_all_sort = d_min_obs_all.select(perm)
  sum_refl_iso_all_sort = sum_refl_iso_all.select(perm)
  sum_refl_obs_all_sort = sum_refl_obs_all.select(perm)
  cc_all_sort = cc_all.select(perm)
  
  n_show_top = 10
  cn_n_show = 0
  for frame_no, d_min_iso, d_min_obs, sum_iso, sum_obs, cc in zip(frame_no_all_sort, 
        d_min_iso_all_sort, d_min_obs_all_sort, sum_refl_iso_all_sort, sum_refl_obs_all_sort, cc_all_sort):
    print frame_no, d_min_iso, d_min_obs, sum_iso, sum_obs, cc
    cn_n_show += 1
    if cn_n_show == n_show_top:
      break
  
      
      
      
      
      
