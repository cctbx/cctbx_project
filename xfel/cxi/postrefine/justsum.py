"""
Read pickle files output from integration process to determine the polarity and
refine the rotation matrix
Usage:
hklin=file_ref.mtz 
pickle_dir=/path/to/integration/folder/ 
frame_start=0 
frame_end=1 
flag_polar=1
target_unit_cell = a,b,c,alpha,beta,gamma
target_space_group = 'P1'
target_anomalous_flag = False
"""
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
  except ZeroDivisionError:
    pass
    
  return corr, slope

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

def organize_input(observations, anomalous_flag=False):
  """Given the supplied observations and Miller array, return a tuple of
  the original observations, the ASU-mapped observations, and the
  observations reindexed using (k, h, -l).  Merge Bijvoet mates unless
  anomalous_flag is True.
  """
  
  observations = observations.resolution_filter(d_min=d_min, d_max=d_max)
  
  miller_set = symmetry(
      unit_cell=target_unit_cell,
      space_group_symbol=target_space_group
    ).build_miller_set(
      anomalous_flag=anomalous_flag,
      d_min=d_min)
  
  i_I_positive = (observations.data() > 0)
  miller_indices_positive = observations.indices().select(i_I_positive)
  I_positive = observations.data().select(i_I_positive)
  sigI_positive = observations.sigmas().select(i_I_positive)
  
  observations = observations.customized_copy(indices=miller_indices_positive,
      data=I_positive,
      sigmas=sigI_positive,
      anomalous_flag=anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  
  I_obs_scaled = (observations.data() - np.mean(observations.data()))/np.std(observations.data())
  i_I_obs_sel = (flex.abs(I_obs_scaled) < sigma_max)
  
  observations = observations.customized_copy(indices=observations.indices().select(i_I_obs_sel), 
      data=observations.data().select(i_I_obs_sel), 
      sigmas=observations.sigmas().select(i_I_obs_sel),
      anomalous_flag=anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  
  observations_asu = observations.customized_copy(
      anomalous_flag=anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      ).map_to_asu()
  
  observations_original = observations.customized_copy(
      anomalous_flag=anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  
  if flag_polar:
    from cctbx import sgtbx
    cb_op = sgtbx.change_of_basis_op('k,h,-l')
    observations_rev = observations_asu.change_basis(cb_op).map_to_asu()
  else:
    observations_rev = observations_asu
  
  return observations_original, observations_asu, observations_rev
  
  
def determine_polar(observations_original, observations_asu, observations_rev, miller_array_main):
  
  """
  determine polarity based on input data
  """
  
  matches_asu = miller.match_multi_indices(
              miller_indices_unique=miller_array_main.indices(),
              miller_indices=observations_asu.indices())
  
  I_ref_asu = flex.double([miller_array_main.data()[pair[0]] for pair in matches_asu.pairs()])
  I_obs_asu = flex.double([observations_asu.data()[pair[1]] for pair in matches_asu.pairs()])
  miller_indices_ori_asu = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_asu.pairs()))
  
  corr_raw_asu, slope_raw_asu = get_overall_correlation(I_ref_asu, I_obs_asu)
  
  matches_rev = miller.match_multi_indices(
              miller_indices_unique=miller_array_main.indices(),
              miller_indices=observations_rev.indices())
  
  I_ref_rev = flex.double([miller_array_main.data()[pair[0]] for pair in matches_rev.pairs()])
  I_obs_rev = flex.double([observations_rev.data()[pair[1]] for pair in matches_rev.pairs()])
  miller_indices_ori_rev = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_rev.pairs()))
  
  corr_raw_rev, slope_raw_rev = get_overall_correlation(I_ref_rev, I_obs_rev)
  
  if corr_raw_asu >= corr_raw_rev:
    polar_hkl = 'h,k,l'
    I_ref = I_ref_asu
    I_obs = I_obs_asu
    miller_indices_ori = miller_indices_ori_asu
  else:
    polar_hkl = 'k,h,-l'
    I_ref = I_ref_rev
    I_obs = I_obs_rev
    miller_indices_ori = miller_indices_ori_rev
  
  if len(I_ref) != len(miller_indices_ori):
    print "Error"
    exit()        
    
  return polar_hkl, I_ref, I_obs, miller_indices_ori, corr_raw_asu, corr_raw_rev
  
    
def justsum_mproc(frame_no, miller_array_main):
  """
  Main module to support multi-processor
  """
  
  #1. Organize input
  pickle_filename = frame_files[frame_no]
  trial_results = pickle.load(open(pickle_filename,"rb"))
  observations = trial_results["observations"][0]
  wavelength = trial_results["wavelength"]
        
  observations_original, observations_asu, observations_rev = organize_input(observations, anomalous_flag=target_anomalous_flag)
  
  #2. Determine polarity
  polar_hkl, I_ref, I_obs, miller_indices_ori, corr_raw_asu, corr_raw_rev = determine_polar(observations_original,
            observations_asu, observations_rev, miller_array_main)
  #reset polarity if flag_polar is off
  if flag_polar==False:
    polar_hkl = 'h,k,l'
  
  #3 output mtz file
  if polar_hkl == 'h,k,l':
    observations_sel = observations_asu
  elif polar_hkl == 'k,h,-l':
    observations_sel = observations_rev
  
  return frame_no, (corr_raw_asu, corr_raw_rev), observations_sel.indices(), observations_sel.data(), observations_sel.sigmas(), observations_sel, wavelength
  
  
if __name__=="__main__":

  from iotbx import reflection_file_reader
  from cctbx.array_family import flex
  from cctbx import miller
  from cctbx import crystal
  from cctbx.crystal import symmetry
  from scitbx.matrix import sqr, col
  
  from libtbx.easy_mp import pool_map, get_processes
  from cctbx.crystal_orientation import crystal_orientation
  from scitbx.lstbx import normal_eqns
  
  import matplotlib
  import matplotlib.cm as cm
  import matplotlib.mlab as mlab
  import matplotlib.pyplot as plt
  
  import sys
  import numpy as np
  import os,cPickle as pickle,math
  from datetime import date, datetime, time, timedelta
  
  from mod_polar import polar_manager
  from mod_partiality import partiality_handler
  from mod_util import intensities_scaler
  
  
  #capture starting time
  time_global_start=datetime.now()
  
  #Extracting input
  frame_start = 0
  frame_end = 0
  flag_plot = False
  file_name_in_mtz=''
  pickle_dir = '' 
  d_min = 1.5
  d_max = 45.0
  flag_polar = False
  
  n_bins = 25
  sigma_max = 15
  
  for i in range(len(sys.argv)):
    pair=sys.argv[i].split('=')
    if pair[0]=='frame_start':
      frame_start=int(pair[1])
    elif pair[0]=='frame_end':
      frame_end=int(pair[1])
    elif pair[0]=='flag_plot':
      if int(pair[1])==1:
        flag_plot=True
    elif pair[0]=='flag_polar':
      if int(pair[1])==1:
        flag_polar=True
    elif pair[0]=='hklin':
      file_name_in_mtz=pair[1]
    elif pair[0]=='pickle_dir':
      pickle_dir=pair[1]
  
  
  if file_name_in_mtz == '':
    print "Please input mtz file with reference intensity"
    exit() 
  
  reflection_file = reflection_file_reader.any_reflection_file(file_name_in_mtz)
  miller_arrays=reflection_file.as_miller_arrays()
    
  miller_array_main = miller_arrays[0].as_intensity_array()
  
  data_subset = 0
  frame_files = get_observations(pickle_dir,data_subset)
  
  
  frames = range(frame_start, frame_end)
    
  def justsum_mproc_wrapper(arg):
    return justsum_mproc(arg, miller_array_main)
  
  justsum_result = pool_map(
        args=frames,
        func=justsum_mproc_wrapper,
        processes=None)
  
  
  cc_correct_polar = flex.double()
  cc_wrong_polar = flex.double()
  miller_indices_all = flex.miller_index()
  I_refined_all = flex.double()
  sigI_refined_all = flex.double()
  cn_accept_frame = 0
  observations_all = []
  wavelength_all = flex.double()
  for frame_no, cc_polar, miller_indices_refined, I_obs_refined, sig_I_obs_refined, observations_sel, wavelength  in justsum_result:
    print frame_no, cc_polar[0], cc_polar[1]
    cn_accept_frame += 1
    observations_all.append(observations_sel)
    wavelength_all.append(wavelength)
    if cc_polar[0] > cc_polar[1]:
      cc_correct_polar.append(cc_polar[0])
      cc_wrong_polar.append(cc_polar[1])
    else:
      cc_correct_polar.append(cc_polar[1])
      cc_wrong_polar.append(cc_polar[0])
    
    for index_new, I_new, sigI_new in zip(miller_indices_refined, I_obs_refined, sig_I_obs_refined):
      miller_indices_all.append(index_new)
      I_refined_all.append(I_new)
      sigI_refined_all.append(sigI_new)
        
  
  miller_array_all = miller_array_main.customized_copy(indices=miller_indices_all,
        data=I_refined_all,
        sigmas=sigI_refined_all)
  
  
  
  #scale and output to mtz files
  inten_scaler = intensities_scaler()
  
  miller_arrays_merge, cc_merge_sets, slope_merge_sets, multiplicities = inten_scaler.output_mtz_files(target_unit_cell, 
        target_space_group, 
        target_anomalous_flag,
        miller_indices_all, 
        I_refined_all, 
        sigI_refined_all,
        flex.double([1]*len(I_refined_all)),
        miller_array_main,
        'justsum_757')
  
   
  miller_array_merge = miller_arrays_merge[0]
  binner_merge = miller_array_merge.setup_binner(n_bins=n_bins)
  binner_merge_indices = binner_merge.bin_indices()
  completeness_merge = miller_array_merge.completeness(use_binning=True)
  
  print 'Completeness (%) Multiplicities (mean, median, std., max, min) and correlation (mean, median, max)'
  bin_multiplicities = []
  for i in range(1,n_bins+1):
    i_binner = (binner_merge_indices == i)
    multiplicities_sel = multiplicities.select(i_binner)
    bin_multiplicities.append((np.mean(multiplicities_sel), np.median(multiplicities_sel), np.std(multiplicities_sel), flex.max(multiplicities_sel), flex.min(multiplicities_sel)))
    
    #calculate C.C.iso per bin
    matches_bin = miller.match_multi_indices(
              miller_indices_unique=miller_array_main.indices(),
              miller_indices=miller_arrays_merge[0].indices().select(i_binner))
  
    I_ref_bin = flex.double([miller_array_main.data()[pair[0]] for pair in matches_bin.pairs()])
    I_obs_bin_mean = flex.double([miller_arrays_merge[0].data().select(i_binner)[pair[1]] for pair in matches_bin.pairs()])
    I_obs_bin_median = flex.double([miller_arrays_merge[1].data().select(i_binner)[pair[1]] for pair in matches_bin.pairs()])
    I_obs_bin_max = flex.double([miller_arrays_merge[2].data().select(i_binner)[pair[1]] for pair in matches_bin.pairs()])
    
    corr_bin_mean, slope_bin_mean = get_overall_correlation(I_obs_bin_mean, I_ref_bin)
    corr_bin_median, slope_bin_median = get_overall_correlation(I_obs_bin_median, I_ref_bin)
    corr_bin_max, slope_bin_max = get_overall_correlation(I_obs_bin_max, I_ref_bin)
    
    completeness = completeness_merge.data[i]
          
    print 'bin ', i, ': %2.4f - %2.4f' %binner_merge.bin_d_range(i), \
      '%3.2f'%(completeness*100), '%4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f'%(np.mean(multiplicities_sel), np.median(multiplicities_sel), np.std(multiplicities_sel), flex.max(multiplicities_sel), flex.min(multiplicities_sel), corr_bin_mean, corr_bin_median, corr_bin_max)
    
  print 'Summary'  
  print 'C.C. correct polar mean=%6.5f median=%6.5f sigma=%6.5f' %(np.mean(cc_correct_polar), np.median(cc_correct_polar), np.std(cc_correct_polar))
  print 'C.C. wrong polar mean=%6.5f median=%6.5f sigma=%6.5f' %(np.mean(cc_wrong_polar), np.median(cc_wrong_polar), np.std(cc_wrong_polar))
  print 'No. of accepted frames:', cn_accept_frame
  print 'No. of all reflections', len(I_refined_all)
  print 'No. of reflections after merge', len(miller_arrays_merge[0].indices())
  print 'Final-merge mean (iso) C.C=%6.5f Slope=%6.5f' %(cc_merge_sets[0], slope_merge_sets[0])
  print 'Final-merge median (iso) C.C=%6.5f Slope=%6.5f' %(cc_merge_sets[1], slope_merge_sets[1])
  print 'Final-merge max (iso) C.C=%6.5f Slope=%6.5f' %(cc_merge_sets[2], slope_merge_sets[2])
      
  """
  plt.scatter(I_ref_all,I_obs_all,s=[10]*len(I_ref_all), marker='o', c='b')
  plt.title('CC Final='+str(cc_all))
  plt.xlabel('Reference intensity')
  plt.ylabel('Observed intensity')
  plt.show()
  """
