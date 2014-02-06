"""
Read pickle files output from integration process to determine the polarity and
refine the rotation matrix
Usage:
libtbx.python postrefine.py input=input.inp
Best protocal:
LP correction CC postref=0.4 CC merge=0.5
"""
class ScoringContainer:
  """
  Compare two rotation matrix and report angular error (Nick's legacy code)
  """
  pass
  
  def angular_rotation(self):
    """Compare the two rotation operators self.model to self.reference and
    return the angle (in radians) between them.
    """

    #print "model"
    model_unit = self.direct_matrix_as_unit_vectors(self.model)
    #print "reference"
    reference_unit = self.direct_matrix_as_unit_vectors(self.reference) 
    rotation = model_unit*reference_unit.inverse()
    UQ = rotation.r3_rotation_matrix_as_unit_quaternion()
    UQ = UQ.normalize() # bugfix; without this many evaluations come to 0.00000 degrees
    angle, axis = UQ.unit_quaternion_as_axis_and_angle()
    #print "axis length",axis.length()
    #print "axis %7.4f %7.4f %7.4f"%(axis[0],axis[1],axis[2])
    return angle
  
  def direct_matrix_as_unit_vectors(self, ori):
    """Return ori's direct_matrix(), with normalized rows.
    """
    direct = sqr(ori.direct_matrix())
    A = col((direct[0],direct[1],direct[2])).normalize()
    B = col((direct[3],direct[4],direct[5])).normalize()
    C = col((direct[6],direct[7],direct[8])).normalize()
    #print "gamma deg",math.acos(A.dot(B))*180./math.pi
    #print "alpha deg",math.acos(B.dot(C))*180./math.pi
    #print " beta deg",math.acos(C.dot(A))*180./math.pi
    direct_as_unit = sqr((A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2]))
    return direct_as_unit

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

def organize_input(observations, 
      d_min, 
      d_max, 
      sigma_max, 
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      flag_polar,
      wavelength,
      miller_array_iso=None):
  """Given the supplied observations and Miller array, return a tuple of
  the original observations, the ASU-mapped observations, and the
  observations reindexed using (k, h, -l).  Merge Bijvoet mates unless
  anomalous_flag is True.
  """
  
  #Lorentz-polarization correction
  two_theta = observations.two_theta(wavelength=wavelength).data()
  one_over_LP = (2 * flex.sin(two_theta))/(1 + (flex.cos(two_theta)**2))
  one_over_P = 2/(1 + (flex.cos(two_theta)**2))
  observations = observations.customized_copy(data=observations.data()*one_over_P)
  
  observations = observations.resolution_filter(d_min=d_min, d_max=d_max)
  
  miller_set = symmetry(
      unit_cell=target_unit_cell,
      space_group_symbol=target_space_group
    ).build_miller_set(
      anomalous_flag=target_anomalous_flag,
      d_min=d_min)
  
  #Filter negative intensities
  i_I_positive = (observations.data() > 0)
  miller_indices_positive = observations.indices().select(i_I_positive)
  I_positive = observations.data().select(i_I_positive)
  sigI_positive = observations.sigmas().select(i_I_positive)
  
  observations = observations.customized_copy(indices=miller_indices_positive,
      data=I_positive,
      sigmas=sigI_positive,
      anomalous_flag=target_anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  
  #Filter possible outliers (use keyword sigma_max=value to control the outliers)
  I_over_sigi = observations.data()/ observations.sigmas()
  i_I_obs_sel = (I_over_sigi > sigma_max)
  
  observations = observations.customized_copy(indices=observations.indices().select(i_I_obs_sel), 
      data=observations.data().select(i_I_obs_sel), 
      sigmas=observations.sigmas().select(i_I_obs_sel),
      anomalous_flag=target_anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  observations = observations.resolution_filter(d_min=d_min, d_max=d_max)
  
  #Prepare original, asu, and rev (for polar case)
  observations_asu = observations.customized_copy(
      anomalous_flag=target_anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      ).map_to_asu()
  
  observations_original = observations.customized_copy(
      anomalous_flag=target_anomalous_flag,
      crystal_symmetry=miller_set.crystal_symmetry()
      )
  
  if flag_polar:
    from cctbx import sgtbx
    cb_op = sgtbx.change_of_basis_op('k,h,-l')
    observations_rev = observations_asu.change_basis(cb_op).map_to_asu()
  else:
    observations_rev = observations_asu
  
  return observations_original, observations_asu, observations_rev
  
  
def determine_polar(observations_original, observations_asu, observations_rev, miller_array_iso, flag_polar):
  
  """
  Determine polarity based on input data.
  The function still needs isomorphous reference so, if flag_polar is True,
  miller_array_iso must be supplied in input file. 
  """
  if flag_polar:
    matches_asu = miller.match_multi_indices(
                miller_indices_unique=miller_array_iso.indices(),
                miller_indices=observations_asu.indices())
    
    I_ref_asu = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_asu.pairs()])
    miller_indices_ref_asu = flex.miller_index((miller_array_iso.indices()[pair[0]] for pair in matches_asu.pairs()))
    I_obs_asu = flex.double([observations_asu.data()[pair[1]] for pair in matches_asu.pairs()])
    sigI_obs_asu = flex.double([observations_asu.sigmas()[pair[1]] for pair in matches_asu.pairs()])
    miller_indices_ori_asu = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_asu.pairs()))
    
    corr_raw_asu, slope_raw_asu = get_overall_correlation(I_obs_asu, I_ref_asu)
    
    matches_rev = miller.match_multi_indices(
                miller_indices_unique=miller_array_iso.indices(),
                miller_indices=observations_rev.indices())
    
    I_ref_rev = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_rev.pairs()])
    miller_indices_ref_rev = flex.miller_index((miller_array_iso.indices()[pair[0]] for pair in matches_rev.pairs()))
    I_obs_rev = flex.double([observations_rev.data()[pair[1]] for pair in matches_rev.pairs()])
    sigI_obs_rev = flex.double([observations_rev.sigmas()[pair[1]] for pair in matches_rev.pairs()])
    miller_indices_ori_rev = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_rev.pairs()))
    
    corr_raw_rev, slope_raw_rev = get_overall_correlation(I_obs_rev, I_ref_rev)
  
  
    if corr_raw_asu >= corr_raw_rev:
      polar_hkl = 'h,k,l'
    else:
      polar_hkl = 'k,h,-l'
      
  else:
    polar_hkl = 'h,k,l'
    corr_raw_asu = 0
    corr_raw_rev = 0
  
  
  return polar_hkl, corr_raw_asu, corr_raw_rev
  
  
def get_pickle_info(pickle_filename):
  """
  #in case there is a reference matrix, calculate the angular error.
  """
  
  trial_results = pickle.load(open(pickle_filename,"rb"))
  observations = trial_results["observations"][0]
  current_hexagonal_ori = trial_results["current_orientation"][0]
  wavelength = trial_results["wavelength"]
  
  current_cb_op_to_primitive = trial_results["current_cb_op_to_primitive"][0]
  current_triclinic_ori = current_hexagonal_ori.change_basis(current_cb_op_to_primitive)
  
  crystal_pointgroup = trial_results["pointgroup"]
   
  CONTAINER_SZ=1000
  mosflm = sqr((1,0,0,0,1,0,0,0,1)) # mosflm lab vectors in their own frame
  labelit = sqr((0,0,-1,-1,0,0,0,1,0)) # mosflm basis vectors in labelit frame
  LM = mosflm * labelit.inverse() # converts labelit frame coords to mosflm frame
  SWAPXY = sqr((0,1,0,1,0,0,0,0,1)) # in labelit frame
  SWAPZ  = sqr((1,0,0,0,1,0,0,0,-1)) # in labelit frame
  R90    = sqr((0,-1,0,1,0,0,0,0,1)) # in labelit frame, rotation 90 on beam axis
  
  i_found_item = pickle_filename.find('/int-data_')
  item_no = int(pickle_filename[i_found_item+10:i_found_item+15].strip())
  container_no = item_no//CONTAINER_SZ
    
    
  filename = "/net/viper/raid1/sauter/fake/holton/mosflm_matrix"
  filename = os.path.join( filename, "%02d"%container_no, "%05d.mat"%item_no)
    
    
  lines = open(filename).readlines()

  A0 = lines[0].strip().split()
  A1 = lines[1].strip().split()
  A2 = lines[2].strip().split()
  A = sqr((float(A0[0]), float(A0[1]), float(A0[2]),
             float(A1[0]), float(A1[1]), float(A1[2]),
             float(A2[0]), float(A2[1]), float(A2[2])))
    
            
  A = A/wavelength
  Holton_hexagonal_ori = crystal_orientation(SWAPZ*SWAPXY*R90*LM*A,True)
  Holton_triclinic_ori = Holton_hexagonal_ori.change_basis(current_cb_op_to_primitive)

  c_inv_r_best = Holton_triclinic_ori.best_similarity_transformation(
                 other=current_triclinic_ori,
                 fractional_length_tolerance=50.,
                 unimodular_generator_range=1)
  c_inv_r_int = tuple([int(round(ij,0)) for ij in c_inv_r_best])
  from cctbx import sgtbx
  c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r_int))
  cb_op = sgtbx.change_of_basis_op(c_inv)
  comparison_triclinic = Holton_triclinic_ori.change_basis(cb_op)
  comparison_hexagonal = comparison_triclinic.change_basis(current_cb_op_to_primitive.inverse())
  
  return current_hexagonal_ori, comparison_hexagonal, wavelength, observations, crystal_pointgroup  

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
  

def get_refined_observations(observations,
      refined_parameters,
      crystal_init_orientation,
      crystal_compare_orientation,
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      wavelength,
      d_min_merge,
      d_min,
      d_max,
      sigma_max,
      flag_polar,
      polar_hkl,
      partiality_threshold_merge,
      miller_array_ref,
      miller_array_iso):
  
  
  observations_original, observations_asu, observations_rev = organize_input(observations, 
      d_min_merge, 
      d_max, 
      sigma_max, 
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      flag_polar,
      wavelength)
  
  observations_sel = observations_asu.deep_copy()
  if flag_polar:
    if polar_hkl == 'k,h,-l':
      observations_sel = observations_rev.deep_copy()
  
  #for B_factor calculation
  observations_sel_two_theta = observations_sel.two_theta(wavelength=wavelength)
  two_theta = observations_sel_two_theta.data()
  sin_theta_over_lambda_sq = observations_sel_two_theta.sin_theta_over_lambda_sq().data()
  
  #grab refined parameters
  G_best = refined_parameters[0]
  B_factor_best = refined_parameters[1]
  rotx_best = refined_parameters[2]
  roty_best = refined_parameters[3]
  ry_best = refined_parameters[4]
  rz_best = refined_parameters[5]
  
  #check how much the rotation matrix has changed
  effective_orientation = crystal_init_orientation.rotate_thru((1,0,0), rotx_best
             ).rotate_thru((0,1,0), roty_best)
  SC = ScoringContainer()
  SC.model = crystal_init_orientation
  SC.reference = crystal_compare_orientation
  lstsqr_angle_error_before = SC.angular_rotation()
  SC.model = effective_orientation
  SC.reference = crystal_compare_orientation
  lstsqr_angle_error_after = SC.angular_rotation()
  
  effective_a_star = sqr(effective_orientation.reciprocal_matrix())
  ph = partiality_handler(wavelength, 0)
  partiality_final = ph.calc_partiality_anisotropy_set(effective_a_star, observations_original.indices(), ry_best, rz_best, two_theta)
  
  I_obs_refined = (G_best * (flex.exp(-2*B_factor_best*sin_theta_over_lambda_sq)) * observations_sel.data())/partiality_final
  sigI_obs_refined = (G_best * (flex.exp(-2*B_factor_best*sin_theta_over_lambda_sq)) * observations_sel.sigmas())/partiality_final
  
  observations_weight = (partiality_final**2)/((G_best**2)*(sigI_obs_refined**2))
  
          
  observations_refined = observations_sel.customized_copy(data=I_obs_refined,
              sigmas=sigI_obs_refined)
              
  #calculate cc_ref, r_ref
  matches = miller.match_multi_indices(
                miller_indices_unique=miller_array_ref.indices(),
                miller_indices=observations_refined.indices())
  
  I_ref = flex.double([miller_array_ref.data()[pair[0]] for pair in matches.pairs()])
  I_obs = flex.double([observations_sel.data()[pair[1]] for pair in matches.pairs()])
  I_obs_refined_match = flex.double([I_obs_refined[pair[1]] for pair in matches.pairs()])
  
  cc_ref_init, slope_ref_init = get_overall_correlation(I_obs, I_ref)
  cc_ref_final, dummy = get_overall_correlation(I_obs_refined_match, I_ref)
  r_ref_init = sum(flex.abs(I_ref - (I_obs*slope_ref_init)))/sum(I_obs*slope_ref_init)
  r_ref_final = sum(flex.abs(I_ref - (I_obs_refined_match)))/sum(I_obs_refined_match)
  
  #calculate cc_iso, r_iso
  cc_iso_init = 0
  cc_iso_final = 0
  r_iso_init = 0
  r_iso_final = 0
  if miller_array_iso is not None:
    matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_iso.indices(),
                  miller_indices=observations_refined.indices())
    I_ref = flex.double([miller_array_iso.data()[pair[0]] for pair in matches.pairs()])
    I_obs = flex.double([observations_sel.data()[pair[1]] for pair in matches.pairs()])
    I_obs_refined_match = flex.double([I_obs_refined[pair[1]] for pair in matches.pairs()])
    
    cc_iso_init, slope_iso_init = get_overall_correlation(I_obs, I_ref)
    cc_iso_final, slope_iso_final = get_overall_correlation(I_obs_refined_match, I_ref)
    r_iso_init = sum(flex.abs(I_ref - (I_obs*slope_iso_init)))/sum(I_obs*slope_iso_init)
    r_iso_final = sum(flex.abs(I_ref - (I_obs_refined_match)))/sum(I_obs_refined_match)
  
  cc_r_results = (cc_ref_init, cc_ref_final, cc_iso_init, cc_iso_final, 
      r_ref_init, r_ref_final, r_iso_init, r_iso_final)
  """
  plt.subplot(211)
  plt.scatter(I_ref, I_obs,s=10, marker='x', c='r')
  plt.title('CC=%6.5f Slope=%6.5f R=%6.5f'%(cc_iso_init, slope_iso_init, r_iso_init))
  plt.xlabel('Reference intensity')
  plt.ylabel('Observed intensity')
  plt.subplot(212)
  plt.scatter(I_ref, I_obs_refined_match,s=10, marker='x', c='r')
  plt.title('CC=%6.5f Slope=%6.5f R=%6.5f'%(cc_iso_final, slope_iso_final, r_iso_final))
  plt.xlabel('Reference intensity')
  plt.ylabel('Observed intensity')
  plt.show()
  """
  return observations_refined, lstsqr_angle_error_before, lstsqr_angle_error_after, cc_r_results, partiality_final, effective_orientation, observations_weight   
  
  
def postrefine_mproc(postref_arg, 
      miller_array_ref, 
      d_min, 
      d_max, 
      sigma_max,
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      flag_polar,
      flag_plot,
      file_name_in_img,
      file_name_in_energy,
      flag_reference_a_matrix,
      crystal_orientations_info,
      miller_array_iso=None):
  """
  Main module to support multi-processor
  """
  postref_arg_tag = postref_arg.split(':')
  frame_no = int(postref_arg_tag[0])
  wavelength_shift = float(postref_arg_tag[1])
  weight_i_ev = 1
  
  #1. Organize input
  pickle_filename = frame_files[frame_no]
  
  if flag_reference_a_matrix:
    crystal_init_orientation, crystal_compare_orientation, wavelength, observations, crystal_pointgroup = get_pickle_info(pickle_filename)
  else:
    trial_results = pickle.load(open(pickle_filename,"rb"))
    crystal_init_orientation = trial_results["current_orientation"][0]
    crystal_compare_orientation = trial_results["current_orientation"][0]
    wavelength = trial_results["wavelength"]
    observations = trial_results["observations"][0]
    crystal_pointgroup = trial_results["pointgroup"]
    unit_cell = trial_results["current_orientation"][0].unit_cell()
    target_unit_cell = unit_cell.parameters()
  
    
  #if crystal_orientations_info is given (next cycle), use the given orientation
  if crystal_orientations_info is not None:
    coi_found = False
    for coi_tuple in crystal_orientations_info:
      if coi_tuple[0] == frame_no and coi_tuple[1] == wavelength_shift:
        crystal_init_orientation = coi_tuple[2]
        coi_found = True
        break
    
    if coi_found == False:
      print 'Error - not found given crystal orientation'
      return None
      
  #handle energy - if energy file is given, use the shifted (sigma) to 
  #determine which wavelength to be used for post-refinement
  #otherwise, use the estimated sigma.
  eh = energy_handler()
  if file_name_in_energy != '' and file_name_in_img !='':
    eh.get_energy_info(file_name_in_img, file_name_in_energy, pickle_filename)  
    wavelength = eh.wavelength_at_counts_list[int(wavelength_shift)]
    weight_i_ev = eh.ev_weight_list[int(wavelength_shift)]
  else:
    wavelength = wavelength + (wavelength_shift*0.009)
  
  imgname = pickle_filename
  if file_name_in_img != '':
    fh = file_handler()
    imgname = fh.get_imgname_from_pickle_filename(iph.file_name_in_img, pickle_filename)
  
  
  observations_original, observations_asu, observations_rev = organize_input(observations, 
      d_min, 
      d_max, 
      sigma_max, 
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      flag_polar,
      wavelength,
      miller_array_iso=miller_array_iso)
  
  #2. Determine polarity - always do this even if flag_polar = False 
  #the function will take care of it.
  polar_hkl, cc_iso_raw_asu, cc_iso_raw_rev = determine_polar(observations_original, 
        observations_asu, observations_rev, miller_array_iso, flag_polar)
  
  #3. Prepare data for post-refinement
  if polar_hkl == 'h,k,l':
    matches_asu = miller.match_multi_indices(
                miller_indices_unique=miller_array_ref.indices(),
                miller_indices=observations_asu.indices())
    
    I_ref_asu = flex.double([miller_array_ref.data()[pair[0]] for pair in matches_asu.pairs()])
    miller_indices_ref_asu = flex.miller_index((miller_array_ref.indices()[pair[0]] for pair in matches_asu.pairs()))
    I_obs_asu = flex.double([observations_asu.data()[pair[1]] for pair in matches_asu.pairs()])
    sigI_obs_asu = flex.double([observations_asu.sigmas()[pair[1]] for pair in matches_asu.pairs()])
    miller_indices_ori_asu = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_asu.pairs()))
    
    references_sel = miller_array_ref.customized_copy(data=I_ref_asu, indices=miller_indices_ref_asu)
    observations_original_sel = observations_original.customized_copy(data=I_obs_asu, 
        sigmas=sigI_obs_asu, 
        indices=miller_indices_ori_asu,
        anomalous_flag=target_anomalous_flag)
  elif polar_hkl == 'k,h,-l':
    matches_rev = miller.match_multi_indices(
                miller_indices_unique=miller_array_ref.indices(),
                miller_indices=observations_rev.indices())
    
    I_ref_rev = flex.double([miller_array_ref.data()[pair[0]] for pair in matches_rev.pairs()])
    miller_indices_ref_rev = flex.miller_index((miller_array_ref.indices()[pair[0]] for pair in matches_rev.pairs()))
    I_obs_rev = flex.double([observations_rev.data()[pair[1]] for pair in matches_rev.pairs()])
    sigI_obs_rev = flex.double([observations_rev.sigmas()[pair[1]] for pair in matches_rev.pairs()])
    miller_indices_ori_rev = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_rev.pairs()))
    
    references_sel = miller_array_ref.customized_copy(data=I_ref_rev, indices=miller_indices_ref_rev)
    observations_original_sel = observations_original.customized_copy(data=I_obs_rev, 
        sigmas=sigI_obs_rev, 
        indices=miller_indices_ori_rev,
        anomalous_flag=target_anomalous_flag)
    
  
  #3. Refine rotation matrix
  a_star_init = sqr(crystal_init_orientation.reciprocal_matrix())
  
  I_ref = references_sel.data()
  I_obs = observations_original_sel.data()
  I_obs = I_obs/ weight_i_ev
  sigI_obs = observations_original_sel.sigmas()
  
  #for B-factor weighting, precalculate (sin(2theta)/lambda)^2
  observations_original_sel_two_theta = observations_original_sel.two_theta(wavelength=wavelength)
  two_theta = observations_original_sel_two_theta.data()
  sin_theta_over_lambda_sq = observations_original_sel_two_theta.sin_theta_over_lambda_sq().data()
  
  #use grid search to approximate intitial value for ry and rz in
  #III.4 (Winkler et al. 1979 A35 p.901)
  spot_radius = calc_spot_radius(a_star_init, observations_original_sel.indices(), wavelength)
  r_shift = 0.0001
  ry_set = np.arange(-0.001,0.001, r_shift)
  r_ref_set = flex.double()
  ry_add_set = flex.double()
  rz_add_set = flex.double()
  
  cc_ref_init, G_init = get_overall_correlation(I_obs, I_ref)
  ph = partiality_handler(wavelength, 0)
  for ry_shift in ry_set:
    ry = spot_radius + ry_shift
    rz = ry * 0.5
    partiality_aniso = ph.calc_partiality_anisotropy_set(a_star_init, observations_original_sel.indices(), ry, rz, two_theta)
      
    r_ref = sum(flex.abs((G_init * (I_obs/partiality_aniso)) - I_ref))/sum(G_init * (I_obs/partiality_aniso))
    r_ref_set.append(r_ref)
    ry_add_set.append(ry)
    rz_add_set.append(rz)
    
  #select ry and rz with the smallest r_ref
  perm = flex.sort_permutation(r_ref_set)
  r_ref_sort = r_ref_set.select(perm)
  ry_add_sort = ry_add_set.select(perm)
  rz_add_sort = rz_add_set.select(perm)
  
  ry_init = ry_add_sort[0]
  rz_init = rz_add_sort[0]
  
  partiality_init = ph.calc_partiality_anisotropy_set(a_star_init, observations_original_sel.indices(), ry, rz, two_theta)
  
  #start solving non-linear least-squares      
  cc_ref_init, G_init = get_overall_correlation(I_obs/partiality_init, I_ref)
  B_factor_init = 0
  rotx_init = 0
  roty_init = 0
  
  lstsqr_before = sum((((G_init * (flex.exp(-2*B_factor_init*sin_theta_over_lambda_sq)) * I_obs)/partiality_init) - I_ref)**2)
  
  if (cc_ref_init < corr_threshold_postref) or (frame_no in frame_excludes):
    G_best = G_init
    B_factor_best = B_factor_init
    rotx_best = rotx_init
    roty_best = roty_init
    ry_best = ry_init
    rz_best = rz_init
  else:
    """
    J_history = []
    grad_history = []
    grad_thres = 1e-5
    parameters = (G_init, rotx_init, roty_init, ry_init, rz_init) 
    const_params = (ry_init, rz_init)
    
    neh = normal_eqns_handler()
    helper = neh.per_frame_helper_factory(I_ref, observations_original_sel, 
                wavelength, crystal_init_orientation,  
                parameters, J_history, grad_history)
    helper.restart()
    iterations = normal_eqns_solving.naive_iterations(
                 non_linear_ls = helper,
                 gradient_threshold = grad_thres)
    lstsqr_results =  helper.x
    G_best = lstsqr_results[0]
    B_factor_best = 0
    rotx_best = lstsqr_results[1]
    roty_best = lstsqr_results[2]
    ry_best = lstsqr_results[3]
    rz_best = lstsqr_results[4]
    """
    n_iters = 10
    parameters = flex.double([rotx_init, roty_init, ry_init, rz_init]) 
    neh = normal_eqns_handler()
    lstsqr_results = neh.postref_with_separable_scale_factor(I_ref, observations_original_sel, 
              wavelength, crystal_init_orientation,  
              parameters, n_iters)
    
    B_factor_best = 0
    rotx_best = lstsqr_results[0]
    roty_best = lstsqr_results[1]
    ry_best = lstsqr_results[2]
    rz_best = lstsqr_results[3]
    
    #find G_bset
    effective_orientation = crystal_init_orientation.rotate_thru((1,0,0), rotx_best
             ).rotate_thru((0,1,0), roty_best)
    effective_a_star = sqr(effective_orientation.reciprocal_matrix())
    ph = partiality_handler(wavelength, 0)
    partiality_final = ph.calc_partiality_anisotropy_set(effective_a_star, observations_original.indices(), ry_best, rz_best, two_theta)
    
    I_obs_refined = (flex.exp(-2*B_factor_best*sin_theta_over_lambda_sq) * observations_original_sel.data())/partiality_final
    dummy, G_best = get_overall_correlation(I_obs_refined, I_ref)
  
  
  lstsqr_results_best = (G_best, B_factor_best, rotx_best, roty_best, ry_best, rz_best)
  observations_refined, lstsqr_angle_error_before, lstsqr_angle_error_after, cc_r_results, partiality_final, effective_orientation, observations_weight = get_refined_observations(observations,
      lstsqr_results_best,
      crystal_init_orientation,
      crystal_compare_orientation,
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      wavelength,
      d_min_merge,
      d_min,
      d_max,
      sigma_max,
      flag_polar,
      polar_hkl,
      partiality_threshold_merge,
      miller_array_ref,
      miller_array_iso)
  
    
  lstsqr_after = 0
  
  cc_ref_init, cc_ref_final, cc_iso_init, cc_iso_final, r_ref_init, r_ref_final, r_iso_init, r_iso_final = cc_r_results
  
  print 'frame:%04d %3.2f %2.2f %7.0f %7.0f %2.4f %2.4f %2.4f %2.4f'%(frame_no, wavelength, weight_i_ev, len(observations.indices()), len(observations_refined.indices()), cc_ref_init, cc_ref_final, r_ref_init, r_ref_final)
  
  return frame_no, (len(observations.indices()), len(observations_refined.indices())), \
    cc_r_results, \
    polar_hkl, pickle_filename, lstsqr_results_best, (lstsqr_angle_error_before, lstsqr_angle_error_after), \
    (lstsqr_before, lstsqr_after), \
    eh, wavelength, observations_refined, partiality_init, partiality_final, wavelength_shift, effective_orientation, observations_weight
  

def justsum_mproc(frame_no, 
      d_min_merge, 
      d_max, 
      sigma_max,
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      flag_polar,
      flag_plot,
      miller_array_iso=None):
  """
  generate justavg.mtz file as a reference set
  """
  pickle_filename = frame_files[frame_no]
  
  trial_results = pickle.load(open(pickle_filename,"rb"))
  observations = trial_results["observations"][0]
  wavelength = trial_results["wavelength"]
  pointgroup = trial_results["pointgroup"]
  
  observations_original, observations_asu, observations_rev = organize_input(observations, 
      d_min_merge, 
      d_max, 
      sigma_max, 
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      flag_polar,
      wavelength)
  polar_hkl, cc_iso_raw_asu, cc_iso_raw_rev = determine_polar(observations_original, 
        observations_asu, observations_rev, miller_array_iso, flag_polar)
        
  if polar_hkl == 'h,k,l':
    observations_sel = observations_asu
  elif polar_hkl == 'k,h,-l':
    observations_sel = observations_rev
  
  
  return frame_no, observations_sel

def scale_justsum_mproc(frame_no_new_order,
      observations_all, mean_observations, frames_no_real=None):
  
  this_observations = observations_all[frame_no_new_order]
  this_observations_G = mean_observations/ np.mean(this_observations.data())
  this_observations_scaled = this_observations.customized_copy(data = this_observations_G * this_observations.data(),
      sigmas = this_observations_G * this_observations.sigmas())
      
  if frames_no_real is None:
    frame_no = frame_no_new_order
  else:
    frame_no = frames_no_real[frame_no_new_order]
  
  return this_observations_scaled, this_observations_G, frame_no
   
#if __name__=="__main__":
def start(args):
  from iotbx import reflection_file_reader
  from cctbx.array_family import flex
  from cctbx import miller
  from cctbx import crystal
  from cctbx.crystal import symmetry
  from scitbx.matrix import sqr, col
  
  from libtbx.easy_mp import pool_map, get_processes
  from cctbx.crystal_orientation import crystal_orientation
  from scitbx.lstbx import normal_eqns
  from scitbx.lstbx import normal_eqns_solving
  
  import matplotlib
  import matplotlib.cm as cm
  import matplotlib.mlab as mlab
  import matplotlib.pyplot as plt
  
  import sys
  import numpy as np
  import random
  import os,cPickle as pickle,math
  from datetime import date, datetime, time, timedelta
  
  from mod_polar import polar_manager
  from mod_partiality import partiality_handler
  from mod_energy import energy_handler
  from mod_util import intensities_scaler, file_handler, input_handler
  from mod_normal_eqns import normal_eqns_handler
  
  
  #capture starting time
  time_global_start=datetime.now()
  
  #non-customizable parameters for program
  corr_threshold_postref = 0.25
  corr_threshold_merge = 0.25
  partiality_threshold_postref = 0.0
  partiality_threshold_merge = 0.0
  n_bins = 25
  d_min_merge = 1.5
  frame_excludes = (9999,9999)
  
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='input':
      file_name_input = pair[1]
  
  if file_name_input == '':
    print "Please provide input-parameters file (usage: input=yourinput.inp)"
    exit() 
  
  iph = input_handler()
  iph.read_input(file_name_input)
  
  #make run_no folder
  if not os.path.exists(iph.run_no):
    os.makedirs(iph.run_no)
  
  #fetch isomorphous structure
  miller_array_iso = None
  if iph.file_name_iso_mtz != '':
    reflection_file_iso = reflection_file_reader.any_reflection_file(iph.file_name_iso_mtz)
    miller_arrays_iso=reflection_file_iso.as_miller_arrays()
    miller_array_iso = miller_arrays_iso[1].generate_bijvoet_mates()
    
    #special for t4_base (M. Rossmann - A5.mtz 3 miller arrays and intesity is the last one)
    #miller_array_iso = miller_arrays_iso[2].generate_bijvoet_mates()
    
    
  frame_files = get_observations(iph.pickle_dir, 0)
  
  frames = range(iph.frame_start, iph.frame_end)
  frames_rand_200 = [347, 3, 526, 726, 172, 692, 436, 670, 701, 426, 279, 663, 137, 225, 390, 363, 110, 266, 150, 485, 508, 582, 600, 129, 407, 349, 56, 194, 441, 4, 562, 11, 496, 629, 381, 323, 410, 377, 333, 133, 34, 156, 59, 525, 199, 254, 527, 164, 309, 189, 606, 31, 229, 166, 592, 68, 265, 191, 248, 669, 675, 651, 424, 307, 78, 661, 736, 596, 601, 152, 487, 41, 70, 158, 745, 107, 329, 169, 126, 117, 210, 171, 585, 432, 359, 188, 438, 306, 462, 178, 258, 239, 51, 687, 623, 691, 27, 589, 184, 212, 529, 547, 520, 690, 135, 443, 626, 397, 396, 356, 47, 299, 593, 574, 476, 345, 69, 519, 219, 157, 186, 657, 18, 567, 537, 557, 565, 654, 76, 221, 83, 8, 754, 170, 46, 319, 656, 492, 388, 200, 636, 447, 569, 278, 5, 676, 723, 82, 461, 403, 486, 247, 314, 311, 634, 357, 44, 344, 452, 737, 700, 584, 560, 287, 732, 594, 375, 48, 506, 729, 650, 206, 142, 187, 748, 120, 538, 330, 342, 62, 331, 246, 662, 42, 412, 680, 98, 89, 23, 553, 118, 431, 433, 613, 36, 6, 75, 459, 481, 37]
  frames_rand_100 = [212, 380, 396, 445, 71, 455, 182, 553, 492, 732, 612, 218, 421, 260, 65, 77, 412, 626, 700, 590, 89, 441, 57, 171, 125, 654, 354, 448, 487, 64, 216, 68, 72, 728, 243, 78, 4, 280, 321, 390, 564, 210, 646, 620, 417, 1, 431, 495, 381, 499, 283, 551, 289, 184, 488, 504, 82, 9, 144, 652, 531, 284, 336, 708, 484, 247, 369, 598, 98, 534, 599, 497, 438, 587, 521, 165, 116, 115, 674, 465, 563, 738, 242, 450, 391, 663, 337, 503, 415, 404, 45, 332, 343, 684, 588, 399, 146, 265, 407, 481]
  frames_rand_best_110 = [198, 192, 208, 172, 41, 221, 684, 23, 246, 81, 382, 546, 446, 402, 308, 223, 53, 637, 343, 185, 140, 173, 211, 358, 160, 725, 346, 16, 169, 204, 293, 360, 170, 302, 516, 368, 189, 430, 432, 466, 183, 344, 148, 456, 476, 299, 244, 250, 231, 270, 119, 202, 251, 188, 200, 333, 574, 667, 59, 46, 245, 452, 284, 547, 92, 187, 328, 597, 696, 95, 209, 528, 711, 310, 156, 589, 641, 217, 233, 549, 721, 273, 353, 619, 48, 714, 228, 118, 491, 212, 214, 424, 289, 461, 659, 612, 268, 174, 523, 676, 735, 435, 332, 743, 577, 506, 414, 138, 184, 559]
  n_frame = 100
  #frames = frames_rand_best_110[0:n_frame]
  
  #fetch reference set
  txt_out_ref_set = ''
  if iph.file_name_ref_mtz != '':
    reflection_file_ref = reflection_file_reader.any_reflection_file(iph.file_name_ref_mtz)
    miller_arrays_ref = reflection_file_ref.as_miller_arrays()
    miller_array_ref = miller_arrays_ref[1].generate_bijvoet_mates()
    
    #special for t4_base (M. Rossmann - A5.mtz 3 miller arrays and intesity is the last one)
    #miller_array_ref = miller_arrays_ref[2].generate_bijvoet_mates()
    
  else:
    print 'Generate reference set by averaging'
    txt_out_ref_set += 'Generate reference set by averaging\n'
    
    #grab all observations for <I> overall distributions.
    def justsum_mproc_wrapper(arg):
      return justsum_mproc(arg, 
        d_min_merge, 
        iph.d_max, 
        iph.sigma_max,
        iph.target_unit_cell,
        iph.target_space_group,
        iph.target_anomalous_flag,
        iph.flag_polar,
        iph.flag_plot,
        miller_array_iso)
    
    justsum_result = pool_map(
          args=frames,
          func=justsum_mproc_wrapper,
          processes=None)
    
    observations_justsum = []
    miller_indices_all = flex.miller_index()
    I_refined_all = flex.double()
    sigI_refined_all = flex.double()
    mean_I_isnot_nan = flex.double()
    for result_row in justsum_result:
      if result_row is not None:
        frame_no = result_row[0]
        observations_sel = result_row[1]
        
        observations_justsum.append(observations_sel)
        frame_mean_i = np.mean(observations_sel.data())
        if math.isnan(frame_mean_i) == False:
          mean_I_isnot_nan.append(frame_mean_i)
      
    
    #scale each frame
    frames_new_order = range(len(observations_justsum))
    def scale_justsum_mproc_wrapper(arg):
      return scale_justsum_mproc(arg, 
        observations_justsum,
        np.mean(mean_I_isnot_nan))
    
    scale_justsum_result = pool_map(
          args=frames_new_order,
          func=scale_justsum_mproc_wrapper,
          processes=None)
       
    miller_indices_all = flex.miller_index()
    I_refined_all = flex.double()
    sigI_refined_all = flex.double()
    frame_no_all = flex.double()
    for observations_scaled, mean_G, frame_no_new_order  in scale_justsum_result:
      for miller_index_scaled, I_scaled, sigI_scaled in zip(observations_scaled.indices(), 
            observations_scaled.data(), observations_scaled.sigmas()):
        miller_indices_all.append(miller_index_scaled)
        I_refined_all.append(I_scaled)
        sigI_refined_all.append(sigI_scaled)
        frame_no_all.append(frame_no_new_order)
    
    
        
    inten_scaler = intensities_scaler()
    miller_array_merge_justsumG, cc_merge, slope_merge,txt_output_mtz = inten_scaler.output_mtz_files(iph.target_unit_cell, 
        iph.target_space_group, 
        iph.target_anomalous_flag,
        miller_indices_all, 
        I_refined_all, 
        sigI_refined_all,
        frame_no_all,
        miller_array_iso,
        iph.run_no+'/justavgG',
        n_bins=n_bins,
        flag_on_screen_output=True,
        sigI_thres=99)
    
    txt_out_ref_set += txt_output_mtz
    #assign the mean-intensity-scaled as the reference set
    miller_array_ref = miller_array_merge_justsumG.generate_bijvoet_mates()
    
    
  wavelength_shift = 1
  wavelength_set = np.arange(0,1, wavelength_shift)
  
  postref_arg = []
  for fr in frames:
    for wl in wavelength_set:  
       postref_arg.append(str(fr)+':'+str(wl))
       
  #do post-refinement
  crystal_orientations_info = None
  txt_out_postref = ''
  for i_postref_cycle in range(iph.n_postref_cycle):
    print 'Running post-refinement cycle#', i_postref_cycle
    txt_out_postref += 'Running post-refinement cycle#'+str(i_postref_cycle+1)+'\n'
    def postrefine_mproc_wrapper(arg):
      return postrefine_mproc(arg, miller_array_ref, 
        iph.d_min, 
        iph.d_max, 
        iph.sigma_max,
        iph.target_unit_cell,
        iph.target_space_group,
        iph.target_anomalous_flag,
        iph.flag_polar,
        iph.flag_plot,
        iph.file_name_in_img,
        iph.file_name_in_energy,
        iph.flag_reference_a_matrix,
        crystal_orientations_info,
        miller_array_iso)
    
    post_refine_result = pool_map(
          args=postref_arg,
          func=postrefine_mproc_wrapper,
          processes=None)
    
  
    ae_ini_all = flex.double()
    ae_fin_all = flex.double()
    G_all = flex.double()
    std_e_all = flex.double()
    partiality_init_mean_all = flex.double()
    partiality_refined_mean_all = flex.double()
    partiality_init_median_all = flex.double()
    partiality_refined_median_all = flex.double()
    
    observations_postref_set = []
    frames_no_set = flex.double()
    mean_I_isnot_nan = flex.double()
    crystal_orientations_info = [] # list of (frame_no, wavelength_shift, crystal_orientation)
    cn_indices = 0
    
    miller_indices_all = flex.miller_index()
    I_refined_all = flex.double()
    sigI_refined_all = flex.double()
    I_weight_all = flex.double()
    frames_no_all = flex.double()
    for result_row  in post_refine_result:
      if result_row is not None:
        
        frame_no = result_row[0]
        cc_r_results = result_row[2]
        pickle_filename = result_row[4]
        lstsqr_results = result_row[5]
        lstsqr_angle_error = result_row[6]
        lstsqr = result_row[7]
        eh = result_row[8]
        wavelength = result_row[9]
        observations_postref = result_row[10].deep_copy()
        partiality_init = result_row[11]
        partiality_refined = result_row[12]
        wavelength_shift = result_row[13]
        crystal_orientation_refined = result_row[14]
        observations_weight = result_row[15]
        
        crystal_orientations_info.append((frame_no, wavelength_shift, crystal_orientation_refined))
        
        lstsqr_angle_error_before = lstsqr_angle_error[0]
        lstsqr_angle_error_after = lstsqr_angle_error[1]
        lstsqr_before = lstsqr[0]
        lstsqr_after = lstsqr[1]
        
        partiality_init_mean_all.append(flex.mean(partiality_init))
        partiality_refined_mean_all.append(flex.mean(partiality_refined))
        
        partiality_init_median_all.append(np.median(partiality_init))
        partiality_refined_median_all.append(np.median(partiality_refined))
        
        imgname = pickle_filename
        if iph.file_name_in_img != '':
          fh = file_handler()
          imgname = fh.get_imgname_from_pickle_filename(iph.file_name_in_img, pickle_filename)
        
        
        cc_ref_init, cc_ref_final, cc_iso_init, cc_iso_final, r_ref_init, r_ref_final, r_iso_init, r_iso_final = cc_r_results
        
        img_exclude = 'XX'  
        
        if (cc_ref_final > corr_threshold_merge) and \
          (imgname.find(img_exclude) == -1):
          
          ae_ini_all.append(lstsqr_angle_error_before)
          ae_fin_all.append(lstsqr_angle_error_after)
          G_all.append(lstsqr_results[0])
          std_e_all.append(eh.std_energy)
          
          observations_postref_set.append(observations_postref)
          frames_no_set.append(frame_no)
          
          for miller_index_scaled, I_scaled, sigI_scaled, I_weight in zip(observations_postref.indices(), 
              observations_postref.data(), observations_postref.sigmas(), observations_weight):
            miller_indices_all.append(miller_index_scaled)
            I_refined_all.append(I_scaled)
            sigI_refined_all.append(sigI_scaled)
            frames_no_all.append(frame_no)
            I_weight_all.append(I_weight)
          
          frame_mean_i = np.mean(observations_postref.data())
          if math.isnan(frame_mean_i) == False:
            mean_I_isnot_nan.append(frame_mean_i)
          
        else:
          #output image name not merged
          print imgname, ' - not merged'
        
        
        txt_out_postref += '%04d %2.4f %7.0f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f'%(frame_no, wavelength,
        len(observations_postref.data()), cc_ref_init, cc_ref_final, cc_iso_init, cc_iso_final,
        r_ref_init, r_ref_final, r_iso_init, r_iso_final)
        txt_out_postref += '\n'  
    
    
    """
    #rescale each frame (after post-refinement)
    frames_new_order = range(len(observations_postref_set))
    def scale_postref_mproc_wrapper(arg):
      return scale_justsum_mproc(arg, 
        observations_postref_set,
        np.mean(mean_I_isnot_nan),
        frames_no_real=frames_no_set)
      
    scale_postref_result = pool_map(
            args=frames_new_order,
            func=scale_postref_mproc_wrapper,
            processes=None)
         
    miller_indices_all = flex.miller_index()
    I_refined_all = flex.double()
    sigI_refined_all = flex.double()
    frames_no_all = flex.double()
    I_weight_all = flex.double()
    for observations_scaled, mean_G, frame_no  in scale_postref_result:
      for miller_index_scaled, I_scaled, sigI_scaled in zip(observations_scaled.indices(), 
              observations_scaled.data(), observations_scaled.sigmas()):
        miller_indices_all.append(miller_index_scaled)
        I_refined_all.append(I_scaled)
        sigI_refined_all.append(sigI_scaled)
        frames_no_all.append(frame_no)
        I_weight_all.append(1)
    """   
    #collect caculating time
    time_global_end=datetime.now()
    time_global_spent=time_global_end-time_global_start
    
    file_mtz_prefix = 'postrefine'
    if iph.file_name_ref_mtz != '':
      file_mtz_prefix = 'postrefine_iso'
      
    if len(I_refined_all) > 0:
      #scale and output to mtz files
      inten_scaler = intensities_scaler()
      miller_array_merge, cc_merge_mean, slope_merge_mean, txt_output_mtz = inten_scaler.output_mtz_files(iph.target_unit_cell, 
            iph.target_space_group, 
            iph.target_anomalous_flag,
            miller_indices_all, 
            I_refined_all, 
            sigI_refined_all,
            frames_no_all,
            miller_array_iso,
            iph.run_no+'/'+file_mtz_prefix+'_cycle_'+str(i_postref_cycle+1),
            I_weight_all=I_weight_all,
            n_bins=n_bins,
            flag_on_screen_output=True,
            sigI_thres=99)
      txt_out_postref += txt_output_mtz
    
    
    #count unique frame
    from collections import Counter
    c = Counter(frames_no_all)
    cn_accept_frame =  len(c.items())
    
    txt_out_postref_sum = ''
    txt_out_postref_sum += 'Post-refinement summary for cycle# %02d' %(i_postref_cycle+1)
    txt_out_postref_sum += '\n'  
    txt_out_postref_sum += 'Angular error (initial) mean=%6.5f med=%6.5f' %(np.mean(ae_ini_all)*180/math.pi, np.median(ae_ini_all)*180/math.pi)
    txt_out_postref_sum += '\n'
    txt_out_postref_sum += 'Angular error (final) lstsqr mean=%6.5f med=%6.5f' %(np.mean(ae_fin_all)*180/math.pi, np.median(ae_fin_all)*180/math.pi)  
    txt_out_postref_sum += '\n'
    txt_out_postref_sum += 'G mean=%6.5f med=%6.5f' %(np.mean(G_all), np.median(G_all)) 
    txt_out_postref_sum += '\n'
    txt_out_postref_sum += 'Partiality (initial) mean=%6.5f med=%6.5f std=%6.5f' %(np.mean(partiality_init_mean_all), np.median(partiality_init_mean_all), np.std(partiality_init_mean_all))
    txt_out_postref_sum += '\n'
    txt_out_postref_sum += 'Partiality (final) mean=%6.5f med=%6.5f std=%6.5f' %(np.mean(partiality_refined_mean_all), np.median(partiality_refined_mean_all), np.std(partiality_refined_mean_all))
    txt_out_postref_sum += '\n'
    txt_out_postref_sum += 'Total calculation time: '+'{0:.2f}'.format(time_global_spent.seconds)+' seconds'
    txt_out_postref_sum += '\n'
    txt_out_postref_sum += 'No. of accepted frames:'+str(cn_accept_frame)
    txt_out_postref_sum += '\n'
    print txt_out_postref_sum
    
    txt_out_postref += txt_out_postref_sum
    #set the reference set to the refine set
    miller_array_ref = miller_array_merge.deep_copy()
  
  #output log file
  file_name_log = 'log.txt'
  if iph.file_name_ref_mtz != '':
    file_name_log = 'log_iso.txt'
  
  
  f = open(iph.run_no+'/'+file_name_log, 'w')
  f.write(iph.txt_out+txt_out_ref_set+txt_out_postref)
  f.close()
  
