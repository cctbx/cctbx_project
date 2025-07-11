"""Helper tools for auto-sharpening"""
from __future__ import absolute_import, division, print_function

import sys,os
from libtbx.utils import Sorry
from cctbx.array_family import flex
from copy import deepcopy
from libtbx import group_args
from libtbx.utils import null_out

import scitbx.lbfgs
import math
from cctbx.maptbx.segment_and_split_map import map_and_b_object
from six.moves import range
from six.moves import zip
from scitbx import matrix
from cctbx import adptbx

def amplitude_quasi_normalisations(ma, d_star_power=1, set_to_minimum=None,
    pseudo_likelihood=False):
    """Calculations for normalizing miller arrays
    for pseudo-likelihood calculation"""
    epsilons = ma.epsilons().data().as_double()
    mean_f_sq_over_epsilon = flex.double()
    for i_bin in ma.binner().range_used():
      sel = ma.binner().selection(i_bin)
      if pseudo_likelihood:
        sel_f_sq = flex.pow2(ma.data().select(sel)) # original method used
      else: # usual
        sel_f_sq = ma.data().select(sel)
      if (sel_f_sq.size() > 0):
        sel_epsilons = epsilons.select(sel)
        sel_f_sq_over_epsilon = sel_f_sq / sel_epsilons
        mean_f_sq_over_epsilon.append(flex.mean(sel_f_sq_over_epsilon))
      else:
        mean_f_sq_over_epsilon.append(0)
    mean_f_sq_over_epsilon_interp = ma.binner().interpolate(
      mean_f_sq_over_epsilon, d_star_power)
    if set_to_minimum and not mean_f_sq_over_epsilon_interp.all_gt(0):
      # HACK NO REASON THIS SHOULD WORK BUT IT GETS BY THE FAILURE
      sel = (mean_f_sq_over_epsilon_interp <= set_to_minimum)
      mean_f_sq_over_epsilon_interp.set_selected(sel,-mean_f_sq_over_epsilon_interp)
      sel = (mean_f_sq_over_epsilon_interp <= set_to_minimum)
      mean_f_sq_over_epsilon_interp.set_selected(sel,set_to_minimum)
    assert mean_f_sq_over_epsilon_interp.all_gt(0)
    from cctbx.miller import array
    return array(ma, flex.sqrt(mean_f_sq_over_epsilon_interp))
    # XXX was below before 2017-10-25
    # return array(ma, mean_f_sq_over_epsilon_interp)

def quasi_normalize_structure_factors(ma, d_star_power=1, set_to_minimum=None,
     pseudo_likelihood=False):
    """Normalize miller arrays for likelihood calculation"""
    normalisations = amplitude_quasi_normalisations(ma, d_star_power,
       set_to_minimum=set_to_minimum,pseudo_likelihood=pseudo_likelihood)
    if pseudo_likelihood:
      print("Norms:")
      for n,d in zip(normalisations[:100],ma.data()[:100]): print(n,d)

    q = ma.data() / normalisations.data()
    from cctbx.miller import array
    return array(ma, q)

def get_array(file_name=None,labels=None):
  """Read from a reflection file. Used in testing these tools"""
  print("Reading from %s" %(file_name))
  from iotbx import reflection_file_reader
  reflection_file = reflection_file_reader.any_reflection_file(
       file_name=file_name)
  array_to_use=None
  if labels:
    for array in reflection_file.as_miller_arrays():
      if ",".join(array.info().labels)==labels:
        array_to_use=array
        break
  else:
    for array in reflection_file.as_miller_arrays():
      if array.is_complex_array() or array.is_xray_amplitude_array() or\
          array.is_xray_intensity_array():
        array_to_use=array
        break
  if not array_to_use:
    text=""
    for array in reflection_file.as_miller_arrays():
      text+=" %s " %(",".join(array.info().labels))

    raise Sorry("Cannot identify array to use...possibilities: %s" %(text))

  print("Using the array %s" %(",".join(array_to_use.info().labels)))
  return array_to_use


def get_amplitudes(args):
  """Read in map coefficients or amplitudes and sharpen, used in testing"""
  if not args or 'help' in args or '--help' in args:
    print("\nsharpen.py")
    print("Read in map coefficients or amplitudes and sharpen")
    return

  new_args=[]
  file_name=None
  for arg in args:
    if os.path.isfile(arg) and arg.endswith(".mtz"):
      file_name=arg
    else:
      new_args.append(arg)
  args=new_args
  labels=None

  array_list=[]

  array_list.append(get_array(file_name=file_name,labels=labels))
  array=array_list[-1]
  phases=None
  assert array.is_complex_array()
  return array


def get_effective_b_values(d_min_ratio=None,resolution_dependent_b=None,
    resolution=None):
  """Return effective b values at sthol2_1 2 and 3
  see adjust_amplitudes_linear below """

  d_min=resolution*d_min_ratio
  sthol2_2=0.25/resolution**2
  sthol2_1=sthol2_2*0.5
  sthol2_3=0.25/d_min**2
  b1=resolution_dependent_b[0]
  b2=resolution_dependent_b[1]
  b3=resolution_dependent_b[2]

  b3_use=b3+b2

  res_1=(0.25/sthol2_1)**0.5
  res_2=(0.25/sthol2_2)**0.5
  res_3=(0.25/sthol2_3)**0.5

  #  Scale factor is exp(b3_use) at res_3 for example
  #  f=exp(-b sthol2)
  #  b= - ln(f)/sthol2
  b1=-b1/sthol2_1
  b2=-b2/sthol2_2
  b3_use=-b3_use/sthol2_3


  return [res_1,res_2,res_3],[b1,b2,b3_use]

def adjust_amplitudes_linear(f_array,b1,b2,b3,resolution=None,
    d_min_ratio=None):
  """do something to the amplitudes.
  b1=delta_b at midway between d=inf and d=resolution,b2 at resolution,
  b3 at d_min (added to b2)
  pseudo-B at position of b1= -b1/sthol2_2= -b1*4*resolution**2
  or...b1=-pseudo_b1/(4*resolution**2)
  typical values of say b1=1 at 3 A -> pseudo_b1=-4*9=-36 """

  data_array=f_array.data()
  sthol2_array=f_array.sin_theta_over_lambda_sq()
  scale_array=flex.double()
  import math
  #d_min=f_array.d_min()
  #if resolution is None: resolution=d_min
  d_min=d_min_ratio*resolution

  sthol2_2=0.25/resolution**2
  sthol2_1=sthol2_2*0.5
  sthol2_3=0.25/d_min**2
  b0=0.0
  d_spacings=f_array.d_spacings()
  b3_use=b3+b2
  for x,(ind,sthol2),(ind1,d) in zip(data_array,sthol2_array,d_spacings):
      if sthol2 > sthol2_2:
        value=b2+(sthol2-sthol2_2)*(b3_use-b2)/(sthol2_3-sthol2_2)
      elif sthol2 > sthol2_1:
        value=b1+(sthol2-sthol2_1)*(b2-b1)/(sthol2_2-sthol2_1)
      else:
        value=b0+(sthol2-0.)*(b1-b0)/(sthol2_1-0.)
      scale_array.append(math.exp(value))
  data_array=data_array*scale_array
  return f_array.customized_copy(data=data_array)

def get_model_map_coeffs_normalized(pdb_hierarchy=None,
   si=None,
   f_array=None,
   overall_b=None,
   resolution=None,
   n_bins=None,
   target_b_iso_model_scale=0,
   target_b_iso_ratio = 5.9,  # empirical, see params for segment_and_split_map
   out=sys.stdout):
  """define Wilson B for the model, generate model map coeffs,
      and normalize the model-map coefficients"""
  if not pdb_hierarchy: return None
  if not si:
    from cctbx.maptbx.segment_and_split_map import sharpening_info
    si=sharpening_info(resolution=resolution,
     target_b_iso_model_scale=0,
     target_b_iso_ratio = target_b_iso_ratio,
     n_bins=n_bins)
  if overall_b is None:
    if si.resolution:
      overall_b=si.get_target_b_iso()*si.target_b_iso_model_scale
    else:
      overall_b=0
    print("Setting Wilson B = %5.1f A" %(overall_b), file=out)

  # create model map using same coeffs
  from cctbx.maptbx.segment_and_split_map import get_f_phases_from_model
  try:
    model_map_coeffs=get_f_phases_from_model(
     pdb_hierarchy=pdb_hierarchy,
     f_array=f_array,
     overall_b=overall_b,
     k_sol=si.k_sol,
     b_sol=si.b_sol,
     out=out)
  except Exception as e:
    print ("Failed to get model map coeffs...going on",file=out)
    return None


  from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi,get_b_iso
  model_f_array,model_phases=map_coeffs_as_fp_phi(model_map_coeffs)
  (d_max,d_min)=f_array.d_max_min(d_max_is_highest_defined_if_infinite=True)
  model_f_array.setup_binner(n_bins=si.n_bins,d_max=d_max,d_min=d_min)

  # Set overall_b....
  final_b_iso=get_b_iso(model_f_array,d_min=resolution)
  print("Effective b_iso of "+\
     "adjusted model map:  %6.1f A**2" %(final_b_iso), file=out)
  model_map_coeffs_normalized=model_f_array.phase_transfer(
     phase_source=model_phases,deg=True)
  return model_map_coeffs_normalized

def get_b_eff(si=None,out=sys.stdout):
  """Get value of effective B factor from si object"""
  if si.rmsd is None:
    b_eff=None
  else:
    b_eff=8*3.14159*si.rmsd**2
    print("Setting b_eff for fall-off at %5.1f A**2 based on model uncertainty of %5.1f A" \
       %( b_eff,si.rmsd), file=out)
  return b_eff

def cc_fit(s_value_list=None,scale=None,value_zero=None,baseline=None,
     scale_using_last=None):
  """Fit CC values to simple exponential.
    For scale_using_last, fix final value at zero"""
  fit=flex.double()
  s_zero=s_value_list[0]
  for s in s_value_list:
    fit.append(value_zero*math.exp(-scale*(s-s_zero)))
  if scale_using_last:
    fit=fit-fit[-1]
  return fit

def get_baseline(scale=None,scale_using_last=None,max_cc_for_rescale=None):
  """Estimate baseline for exponential fit"""
  if not scale_using_last:
    return 0
  else:
    baseline=min(0.99,max(0.,scale[-scale_using_last:].min_max_mean().mean))
    if baseline > max_cc_for_rescale:
       return None
    else:
       return baseline

def fit_cc(cc_list=None,s_value_list=None,
    scale_min=None,scale_max=None,n_tries=None,scale_using_last=None):
  """find value of scale in range scale_min,scale_max that minimizes rms diff
  between cc_list and cc_list[0]*exp(-scale*(s_value-s_value_list[0]))
  for scale_using_last, require it to go to zero at end """

  best_scale=None
  best_rms=None
  for i in range(n_tries):
    scale=scale_min+(scale_max-scale_min)*i/n_tries
    fit=cc_fit(s_value_list=s_value_list,scale=scale,value_zero=cc_list[0],
       scale_using_last=scale_using_last)
    fit=fit-cc_list
    rms=fit.rms()
    if best_rms is None or rms<best_rms:
      best_rms=rms
      best_scale=scale
  return cc_fit(s_value_list=s_value_list,scale=best_scale,value_zero=cc_list[0],
      scale_using_last=scale_using_last)

def get_fitted_cc(cc_list=None,s_value_list=None, cc_cut=None,
   scale_using_last=None,keep_cutoff_point=False,force_scale_using_last=False,
   cutoff_after_last_high_point = False):
  """Alternative fit to CC values.
  only do this if there is some value of s where cc is at least 2*cc_cut or
  (1-c_cut/2), whichever is smaller. """

  min_cc=min(2*cc_cut,1-0.5*cc_cut)
  if cc_list.min_max_mean().max < min_cc and (not force_scale_using_last):
    return cc_list
  # 2020-10-08 instead last point where cc >=cc_cut cutoff_after_last_high_point
  # find first point after point where cc>=min_cc that cc<=cc_cut
  #   then back off by 1 point  # 2019-10-12 don't back off if keep_cutoff_point
  found_high=False
  s_cut=None
  i_cut=0
  if (not cutoff_after_last_high_point):
    for s,cc in zip(s_value_list,cc_list):
      if cc > min_cc:
        found_high=True
      if found_high and cc < cc_cut:
        s_cut=s
        break
      i_cut+=1
  else:
    ii=0
    for s,cc in zip(s_value_list,cc_list):
      if cc > min_cc:
        found_high=True
      if found_high and cc >= cc_cut:
        s_cut = s
        i_cut = ii
      ii += 1

  if force_scale_using_last:
    scale_using_last=True
    s_cut=s_value_list[0]
    i_cut=1
  if s_cut is None or i_cut==0:
    return cc_list

  if keep_cutoff_point:
    i_cut=max(1,i_cut-1)

  #Fit remainder
  s_value_remainder_list=s_value_list[i_cut:]
  cc_remainder_list=cc_list[i_cut:]
  n=cc_remainder_list.size()
  scale_min=10 # soft
  scale_max=500 # hard
  n_tries=200

  fitted_cc_remainder_list=fit_cc(
     cc_list=cc_remainder_list,s_value_list=s_value_remainder_list,
     scale_min=scale_min,scale_max=scale_max,n_tries=n_tries,
     scale_using_last=scale_using_last)

  new_cc_list=cc_list[:i_cut]
  new_cc_list.extend(fitted_cc_remainder_list)
  return new_cc_list

def estimate_cc_star(cc_list=None,s_value_list=None, cc_cut=None,
    scale_using_last=None,
    keep_cutoff_point=False):
  """Estimate value of CC* (True map correlation)
  cc ~ sqrt(2*half_dataset_cc/(1+half_dataset_cc))
  however for small cc the errors are very big and we think cc decreases
  rapidly towards zero once cc is small
  So find value of s_value_zero that gives cc about cc_cut...for
  s_value >s_value_zero use fit of cc_zero * exp(-falloff*(s_value-s_value_zero))
  for scale_using_last set: subtract off final values so it goes to zero."""

  fitted_cc=get_fitted_cc(
    cc_list=cc_list,s_value_list=s_value_list,cc_cut=cc_cut,
    scale_using_last=scale_using_last,
    keep_cutoff_point=keep_cutoff_point)

  cc_star_list=flex.double()
  for cc in fitted_cc:
     cc=max(cc,0.)
     cc_star=(2.*cc/(1.+cc))**0.5
     cc_star_list.append(cc_star)
  return cc_star_list


def rescale_cc_list(cc_list=None,scale_using_last=None,
    max_cc_for_rescale=None):
  """Rescale CC values by replacing cc with (cc-baseline)/(1-baseline)"""
  baseline=get_baseline(scale=cc_list,
      scale_using_last=scale_using_last,
      max_cc_for_rescale=max_cc_for_rescale)
  if baseline is None:
     return cc_list,baseline

  scaled_cc_list=flex.double()
  for cc in cc_list:
    scaled_cc_list.append((cc-baseline)/(1-baseline))
  return scaled_cc_list,baseline

def get_calculated_scale_factors(
      s_value_list=None,
      effective_b=None,
      b_zero = None,
      cc_list = None,
      dv = None,
      uc = None):
    """Calculate expected scale factors from s_value_list and cc_list"""
    recip_space_vectors = flex.vec3_double()
    scale_values = flex.double()
    original_scale_values = flex.double()
    indices = flex.miller_index()
    if effective_b is None and not cc_list:
      #no info
      return group_args(
        indices = indices,
        scale_values = scale_values,
        original_scale_values = original_scale_values)

    for s,cc in zip(s_value_list,cc_list):
      d = 1/s
      sthol2 = 0.25/d**2
      recip_space_vectors.append(matrix.col(dv) * s)
      if effective_b is not None:
        scale_values.append(b_zero*math.exp (max(-20.,min(20.,
         -effective_b * sthol2))))
      else:
        scale_values.append(cc)

      original_scale_values.append(cc)
    indices = flex.miller_index(tuple(get_nearest_lattice_points(
        uc,recip_space_vectors)))
    return group_args(
      indices = indices,
      scale_values = scale_values,
      original_scale_values = original_scale_values)

def get_rms_fo_values(fo = None, direction_vectors = None, d_avg = None,
      apply_aniso_dict_by_dv = None,
       aniso_scale_factor_array = None,
     i_bin = None,
     first_bin = None):
    '''Estimate rms Fo values along direction vectors.
       For first bin everything is just the average'''

    abs_fo_data = flex.abs(fo.data())
    sqr_fo = fo.customized_copy(data=flex.pow2(abs_fo_data))
    rms_fo= sqr_fo.data().min_max_mean().mean**0.5

    if first_bin:
      return flex.double(direction_vectors.size(), rms_fo)

    if aniso_scale_factor_array:
      abs_fo = fo.customized_copy(
        data = abs_fo_data * aniso_scale_factor_array.data())
    else:
      abs_fo = fo.customized_copy(data = abs_fo_data)
    sar = shell_aniso_refinery(abs_fo)
    sar.run()
    # Now apply the values to a dummy array with indices for our dv vectors
    recip_space_vectors = flex.vec3_double()
    apply_aniso_values = flex.double()
    i = -1
    for dv in direction_vectors:
      i+=1
      s = 1/d_avg
      recip_space_vectors.append(matrix.col(dv) * s)
      if  aniso_scale_factor_array:
        apply_aniso_values.append(apply_aniso_dict_by_dv[i][i_bin-1])
    from iotbx.map_model_manager import create_fine_spacing_array
    fo_to_get_values = create_fine_spacing_array( abs_fo.crystal_symmetry().unit_cell())
    indices = flex.miller_index(tuple(get_nearest_lattice_points(
        fo_to_get_values.crystal_symmetry().unit_cell(),recip_space_vectors)))

    fo_to_get_values = fo_to_get_values.customized_copy(indices = indices,
      data = flex.double(indices.size(),1))
    values_array = sar.get_calc_array(
      sar.x,
      array_with_indices=fo_to_get_values)
    # Return estimates of rms F in each direction
    values = values_array.data() * (1./.785)**0.5  # expect: <F>**2/<F**2> = 0.785
    if aniso_scale_factor_array:
       values *= apply_aniso_values
    return values

def calculate_fsc(**kw):
  '''
    Calculate FSC of 2 maps and estimate scale factors
    If direction_vector or direction_vectors supplied, calculate
     CC values weighted by abs() component along direction vector
    If list of direction vectors, return group_args with si objects and
     also include an overall si object
  '''
  si = kw.get('si',None)
  f_array = kw.get('f_array',None)
  map_coeffs = kw.get('map_coeffs',None)
  model_map_coeffs = kw.get('model_map_coeffs',None)
  external_map_coeffs = kw.get('external_map_coeffs',None)
  first_half_map_coeffs = kw.get('first_half_map_coeffs',None)
  second_half_map_coeffs = kw.get('second_half_map_coeffs',None)
  resolution = kw.get('resolution',None)
  n_bins_use = kw.get('n_bins_use',None)
  fraction_complete = kw.get('fraction_complete',None)
  min_fraction_complete = kw.get('min_fraction_complete', 0.05)
  is_model_based = kw.get('is_model_based',None)
  cc_cut = kw.get('cc_cut',None)
  scale_using_last = kw.get('scale_using_last',None)
  max_cc_for_rescale = kw.get('max_cc_for_rescale',None)
  equalize_power = kw.get('equalize_power',False)
  verbose = kw.get('verbose',None)
  maximum_scale_factor = kw.get('maximum_scale_factor',None)
  maximum_ratio = kw.get('maximum_ratio', 25)
  use_dv_weighting = kw.get('use_dv_weighting',None)
  run_analyze_anisotropy = kw.get('run_analyze_anisotropy',None)
  pseudo_likelihood = kw.get('pseudo_likelihood',False)
  skip_scale_factor = kw.get('skip_scale_factor',False)
  optimize_b_eff = kw.get('optimize_b_eff',None)
  direction_vector = kw.get('direction_vector',None)
  direction_vectors = kw.get('direction_vectors',None)
  smooth_fsc = kw.get('smooth_fsc',None)
  cutoff_after_last_high_point = kw.get('cutoff_after_last_high_point',None)
  expected_rms_fc_list = kw.get('expected_rms_fc_list',None)
  expected_ssqr_list = kw.get('expected_ssqr_list',None)
  rmsd_resolution_factor = kw.get('rmsd_resolution_factor',0.25)
  low_res_bins = kw.get('low_res_bins',3)
  low_res_bins = kw.get('low_res_bins',3)
  remove_anisotropy_before_analysis = kw.get(
      'remove_anisotropy_before_analysis',False)
  aniso_scale_factor_array = kw.get( 'aniso_scale_factor_array',None)
  aniso_scale_factor_as_u_cart = kw.get( 'aniso_scale_factor_as_u_cart',None)
  out = kw.get('out',sys.stdout)


  if direction_vectors and direction_vectors != [None]:
     kw['direction_vectors'] = None
     print("Getting overall analysis first",file = out)
     overall_si = calculate_fsc(**kw)
     print("Done with getting overall analysis ",file = out)
  else:
     overall_si = None

  # calculate anticipated fall-off of model data with resolution
  if si.rmsd is None and is_model_based:
    if not rmsd_resolution_factor:
      rmsd_resolution_factor = si.rmsd_resolution_factor
    if not resolution:
      resolution = si.resolution
    si.rmsd=resolution*rmsd_resolution_factor
    print("Setting rmsd to %5.1f A based on resolution of %5.1f A" %(
       si.rmsd,resolution), file=out)
  elif is_model_based:
    if not resolution:
      resolution = si.resolution
    print("RMSD is %5.1f A and resolution is %5.1f A" %(
       si.rmsd,resolution), file=out)

  # get f and model_f vs resolution and FSC vs resolution and apply

  # If external_map_coeffs then simply scale f to external_map_coeffs

  # scale to f_array and return sharpened map
  dsd = f_array.d_spacings().data()
  from cctbx.maptbx.segment_and_split_map import map_coeffs_to_fp

  if is_model_based:
    mc1=map_coeffs
    mc2=model_map_coeffs
    fo_map=map_coeffs # scale map_coeffs to model_map_coeffs*FSC
    fc_map=model_map_coeffs
    b_eff=get_b_eff(si=si,out=out)
  elif external_map_coeffs:
    mc1=map_coeffs
    mc2=external_map_coeffs
    fo_map=map_coeffs # scale map_coeffs to external_map_coeffs
    fc_map=external_map_coeffs
    b_eff=None

  else: # half_dataset
    mc1=first_half_map_coeffs
    mc2=second_half_map_coeffs
    fo_map=map_coeffs # scale map_coeffs to cc*
    fc_map=model_map_coeffs
    b_eff=None

  if not mc1 or not mc2:  # nothing to do
    si.target_scale_factors = None
    return si


  ratio_list=flex.double()
  target_sthol2=flex.double()
  s_value_list=flex.double()
  d_min_list=flex.double()
  rms_fo_list=flex.double()
  rms_fc_list=flex.double()
  max_possible_cc=None
  n_list = flex.double()

  if direction_vectors:
    pass # already ok
  elif direction_vector:
    direction_vectors = [direction_vector]
  else:
    direction_vectors = [None]

  cc_dict_by_dv = {}
  rms_fo_dict_by_dv = {}
  rms_fc_dict_by_dv = {}
  ratio_dict_by_dv = {}
  i = 0
  for dv in direction_vectors:
    cc_dict_by_dv [i] = flex.double()
    rms_fo_dict_by_dv [i] = flex.double()
    rms_fc_dict_by_dv [i] = flex.double()
    ratio_dict_by_dv [i] = flex.double()
    i += 1

  first_bin =True
  weights_para_list = [] # NOTE: this makes N * f_array.size() arrays!!
  for dv in direction_vectors:
    if dv:
      weights_para_list.append(
        get_normalized_weights_para(f_array,direction_vectors, dv,
          include_all_in_lowest_bin = True))
    else:
      weights_para_list.append(None)

  # Set up to work with anisotropy-removed data

  if remove_anisotropy_before_analysis and not aniso_scale_factor_as_u_cart:
   # Get aniso_scale_factor_as_u_cart
    aniso_scale_factor_as_u_cart = get_aniso_scale_info(
      fo_map, resolution = resolution)

  if remove_anisotropy_before_analysis and aniso_scale_factor_as_u_cart and (
     not aniso_scale_factor_array):
    # Get aniso_scale_factor_array

    # Array of scale factors to remove anisotropy
    aniso_scale_factor_array = get_aniso_scale_factor_array(fo_map,
      aniso_scale_factor_as_u_cart)
  if aniso_scale_factor_as_u_cart:
    apply_aniso_dict_by_dv = get_apply_aniso_dict_by_dv(direction_vectors,
       f_array, aniso_scale_factor_as_u_cart)
  else:
    apply_aniso_dict_by_dv = None



  if n_bins_use is None:
    n_bins = len(list(f_array.binner().range_used()))
    n_bins_use = min(n_bins,max(3,n_bins//3))
    set_n_bins_use = True
  else:
    set_n_bins_use = False

  ratio_fo_to_fc_zero = None
  for i_bin in f_array.binner().range_used():
    sel       = f_array.binner().selection(i_bin)
    d         = dsd.select(sel)
    if d.size()<1:
      raise Sorry("Please reduce number of bins (no data in bin "+
        "%s) from current value of %s" %(i_bin,f_array.binner().n_bins_used()))
    d_min     = flex.min(d)
    d_max     = flex.max(d)
    d_avg     = flex.mean(d)

    if set_n_bins_use and i_bin-1 > n_bins_use and (
          (not resolution) or (d_avg >= resolution)):
      n_bins_use = i_bin - 1

    n         = d.size()
    m1        = mc1.select(sel)
    m2        = mc2.select(sel)

    cc = None
    i = 0
    if fc_map:
      fc        = fc_map.select(sel)
    else:
      fc = None
    if fo_map:
          fo        = fo_map.select(sel)
    else:
      fo = None

    if remove_anisotropy_before_analysis and aniso_scale_factor_as_u_cart:
      rms_fo_values = len(direction_vectors) * [None]
    elif (direction_vectors and direction_vectors[0] != None) and (
        aniso_scale_factor_array is not None):
      # Let's fit rms fo in just this shell
      rms_fo_values = get_rms_fo_values(
       fo = fo,
       apply_aniso_dict_by_dv = apply_aniso_dict_by_dv,
       aniso_scale_factor_array = aniso_scale_factor_array.select(sel),
       direction_vectors = direction_vectors,
       d_avg = d_avg,
       i_bin = i_bin,
       first_bin = first_bin)
    else:
      rms_fo_values = len(direction_vectors) * [None]

    for dv, weights_para, rms_fo in zip(direction_vectors, weights_para_list,
         rms_fo_values):
      if dv:
        weights_para_sel = weights_para.select(sel)
        weights_para_sel_sqrt = flex.sqrt(weights_para_sel)
        if remove_anisotropy_before_analysis and aniso_scale_factor_as_u_cart\
            and aniso_scale_factor_array and aniso_scale_factor_array.data():
          scale_sel = aniso_scale_factor_array.data().select(sel)
        else:
          scale_sel = None

        if scale_sel:
          m1a=m1.customized_copy(
             data = m1.data() * weights_para_sel_sqrt * scale_sel)
          m2a=m2.customized_copy(
             data = m2.data() * weights_para_sel_sqrt * scale_sel)
        else: # usual
          m1a=m1.customized_copy(data = m1.data() * weights_para_sel_sqrt)
          m2a=m2.customized_copy(data = m2.data() * weights_para_sel_sqrt)
        cca        = m1a.map_correlation(other = m2a)

        if external_map_coeffs: # only for no direction vectors
          cc=1.
        if cca is None:
          cca=0.
        cc_dict_by_dv[i].append(cca)
        normalization = 1./max(1.e-10,weights_para_sel.rms())
        if scale_sel:
          fo_a = fo.customized_copy(data=fo.data()*weights_para_sel*scale_sel)
          f_array_fo=map_coeffs_to_fp(fo_a)
          rms_fo=normalization * f_array_fo.data().rms() \
              * apply_aniso_dict_by_dv[i][i_bin-1]
        elif (not rms_fo) and fo_map:
          fo_a = fo.customized_copy(data=fo.data()*weights_para_sel)
          f_array_fo=map_coeffs_to_fp(fo_a)
          rms_fo=normalization * f_array_fo.data().rms()
        elif (not rms_fo):
          rms_fo=1.

        if expected_rms_fc_list:
          rms_fc = expected_rms_fc_list[i_bin-1]
        elif fc_map:
          fc_a  = fc.customized_copy(data=fc.data()*weights_para_sel)
          f_array_fc=map_coeffs_to_fp(fc_a)
          rms_fc=normalization *f_array_fc.data().rms()
        else:
          rms_fc=1.
        # Normalize rms_fc to make rms_fc[0] == rms_fo[0]
        if rms_fo and rms_fc and not ratio_fo_to_fc_zero:
          ratio_fo_to_fc_zero = rms_fo/rms_fc
        rms_fc *= (1 if ratio_fo_to_fc_zero is None else ratio_fo_to_fc_zero)


        rms_fo_dict_by_dv[i].append(rms_fo)
        rms_fc_dict_by_dv[i].append(rms_fc)
        ratio_dict_by_dv[i].append(max(1.e-10,rms_fc)/max(1.e-10,rms_fo))

        i += 1
        if (cca is not None) and (cc is None):
          cc = cca # save first one
      else:
        cc        = m1.map_correlation(other = m2)
        if external_map_coeffs: # only for no direction vectors
          cc=1.
        if cc is None:
          cc= 0
        cc_dict_by_dv[i].append(cc)
        if fo_map:
          f_array_fo=map_coeffs_to_fp(fo)
          rms_fo=f_array_fo.data().rms()
        else:
          rms_fo=1.

        if expected_rms_fc_list:
          rms_fc = expected_rms_fc_list[i_bin-1]
        elif fc_map:
          f_array_fc=map_coeffs_to_fp(fc)
          rms_fc=f_array_fc.data().rms()
        else:
          rms_fc=1.4

        # Normalize rms_fc to make rms_fc[0] == rms_fo[0]
        if rms_fo and rms_fc and not ratio_fo_to_fc_zero:
          ratio_fo_to_fc_zero = rms_fo/rms_fc
        rms_fc *= (1 if ratio_fo_to_fc_zero is None else ratio_fo_to_fc_zero)

        rms_fo_dict_by_dv[i].append(rms_fo)
        rms_fc_dict_by_dv[i].append(rms_fc)
        ratio_dict_by_dv[i].append(max(1.e-10,rms_fc)/max(1.e-10,rms_fo))


    sthol2=0.25/d_avg**2 # note this is 0.25 * s_value**2
    target_sthol2.append(sthol2)
    s_value_list.append(1/d_avg)
    d_min_list.append(d_min)
    n_list.append(m1.size())


    if b_eff is not None:
      max_cc_estimate=cc* math.exp(min(20.,sthol2*b_eff))
    else:
      max_cc_estimate=cc
    max_cc_estimate=max(0.,min(1.,max_cc_estimate))
    if max_possible_cc is None or (
        max_cc_estimate > 0 and max_cc_estimate > max_possible_cc):
      max_possible_cc=max_cc_estimate
    if verbose:
      print("d_min: %5.1f  FC: %7.1f  FOBS: %7.1f   CC: %5.2f" %(
      d_avg,rms_fc,rms_fo,cc), file=out)
    first_bin = False

  input_info = group_args(
     f_array = f_array,
     n_list = n_list,
     target_sthol2 = target_sthol2,
     d_min_list = d_min_list,
     pseudo_likelihood = pseudo_likelihood,
     equalize_power = equalize_power,
     is_model_based = is_model_based,
     skip_scale_factor = skip_scale_factor,
     maximum_scale_factor = maximum_scale_factor,
     out = out)

  # Now apply analyses on each cc_list (if more than one)
  si_list = []
  for i in range(len(direction_vectors)):
    ratio_list = remove_values_if_necessary(ratio_dict_by_dv[i])
    rms_fo_list = remove_values_if_necessary(rms_fo_dict_by_dv[i])
    rms_fc_list = remove_values_if_necessary(rms_fc_dict_by_dv[i])
    cc_list = smooth_values(cc_dict_by_dv[i],
         overall_values=getattr(overall_si,'cc_list',None),
         smooth=smooth_fsc)
    if len(direction_vectors) > 1:
      working_si = deepcopy(si)
      dv = direction_vectors[i]
    else:
      dv = None
      working_si = si  # so we can modify it in place

    working_si = complete_cc_analysis(
       dv,
       cc_list,
       rms_fc_list,
       rms_fo_list,
       ratio_list,
       scale_using_last,
       max_cc_for_rescale,
       optimize_b_eff,
       is_model_based,
       s_value_list,
       cc_cut,
       max_possible_cc,
       fraction_complete,
       min_fraction_complete,
       low_res_bins,
       working_si,
       b_eff,
       input_info,
       cutoff_after_last_high_point,
       expected_rms_fc_list,
       expected_ssqr_list,
       overall_si,
       out)
    si_list.append(working_si)
  if direction_vectors == [None]:
    return si_list[0]

  # Results so far
  scale_factor_info = group_args(
     group_args_type = 'scaling_info objects, one set per direction_vector',
     direction_vectors = direction_vectors,
     scaling_info_list = si_list,
     overall_si = overall_si,
     )

  # Analyze anisotropy

  if run_analyze_anisotropy:
    aniso_info = analyze_anisotropy(
        mc1,
        mc2,
        f_array,
        overall_si,
        si_list,
        s_value_list,
        direction_vectors,
        weights_para_list,
        resolution,
        n_bins_use,
        use_dv_weighting,
        expected_ssqr_list,
        maximum_ratio = maximum_ratio,
        out = out)
  else:
    aniso_info = group_args(
     a_zero_values = None,
     fo_b_cart = None,
     fo_b_cart_as_u_cart = None,
     aa_b_cart = None,
     aa_b_cart_as_u_cart = None,
     bb_b_cart = None,
     bb_b_cart_as_u_cart = None,
     ss_b_cart = None,
     ss_b_cart_as_u_cart = None,
     uu_b_cart = None,
     uu_b_cart_as_u_cart = None,
     )

  # Merge in aniso info to scale_factor_info
  scale_factor_info.merge(aniso_info) # aniso_info overwrites scale_factor_info
  return scale_factor_info

def get_apply_aniso_dict_by_dv(direction_vectors,
    f_array, aniso_scale_factor_as_u_cart):
  """Calculate anisotropic values to apply for each direction vector"""
  s_value_list = flex.double()
  dsd = f_array.d_spacings().data()
  for i_bin in f_array.binner().range_used():
    sel       = f_array.binner().selection(i_bin)
    d         = dsd.select(sel)
    d_avg     = flex.mean(d)
    s_value_list.append(1/d_avg)

  apply_aniso_dict_by_dv = {}
  i = -1
  for dv in direction_vectors:
    i+=1

    if not dv:
      apply_aniso_dict_by_dv[i] = flex.double(s_value_list.size(),1.)
    else: # usual

      from iotbx.map_model_manager import create_fine_spacing_array
      fine_array = create_fine_spacing_array(
         f_array.crystal_symmetry().unit_cell())

      indices = get_calculated_scale_factors( # get the indices
          s_value_list=s_value_list,
          cc_list = flex.double(s_value_list.size(),1),
          dv = dv,
          uc = fine_array.crystal_symmetry().unit_cell(),
           ).indices
      fine_array=fine_array.customized_copy(indices = indices,
        data = flex.double(indices.size(), 1.))
      scale_array = get_aniso_scale_factor_array(fine_array,
        tuple(-flex.double(aniso_scale_factor_as_u_cart)))
      apply_aniso_dict_by_dv[i] = scale_array.data()
  return apply_aniso_dict_by_dv


def get_aniso_scale_factor_array(fo_map,
   aniso_scale_factor_as_u_cart):
  """Apply anisotropic scale factor"""
  scale_values_array = fo_map.customized_copy(
       data = flex.double(fo_map.size(),1.))
  u_star= adptbx.u_cart_as_u_star(
     scale_values_array.unit_cell(),
     tuple(matrix.col(aniso_scale_factor_as_u_cart)))
  from mmtbx.scaling import absolute_scaling
  return absolute_scaling.anisotropic_correction(
     scale_values_array,0.0, u_star ,must_be_greater_than=-0.0001)

def get_aniso_scale_info(fo_map, resolution = None):
  ''' Get overall anisotropic scale for fo_map'''
  if not fo_map:
    return
  assert resolution is not None
  from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
  f_local,phases_local=map_coeffs_as_fp_phi(fo_map)
  aniso_obj=analyze_aniso_object()
  aniso_obj.set_up_aniso_correction(f_array=f_local,d_min=resolution)
  if (not aniso_obj) or (not aniso_obj.b_cart):
    return

  return adptbx.b_as_u(aniso_obj.b_cart)



def analyze_anisotropy(
  half_map_coeffs_1,
  half_map_coeffs_2,
  f_array,
  overall_si,
  si_list,
  s_value_list,
  direction_vectors,
  weights_para_list,
  resolution,  # nominal resolution
  n_bins_use,
  use_dv_weighting,
  expected_ssqr_list,
  n_display = 6,
  weight_by_variance = True,
  minimum_sd = 0.1, # ratio to average
  maximum_ratio = None, # maximum ratio of target_scale_factor to overall
  update_scale_values_to_match_s_matrix = False,
  update_scale_values_to_match_qq_values= True,
  out = sys.stdout):
  """Analyze anisotropy in a map.

  Define:
  F1a_obs = rmsFc(|s|) * Ao(|s|) * (A(s) * F1a + B(s) * s1)
  ssqr = rms(s1)**2 = ssqr(|s|)
  At(|s|)  = rmsFc(|s|)

  Then we can ssqr(s) as:
  ssqr(|s|), ssqr(s)  = (1/CC_half(s)) - 1  # s**2, overall and  by direction
  rmsFo(|s|), rmsFo(s)  # rms F obs (mean of half maps) overall and by direction

  Define D(s), D(|s|):
  D(s) = 1./sqrt(1 +  0.5* ssqr(s))

  These can be used to calculate the anisotropy values A(s) and B(s):
  A(s) = (rmsFo(s)/rmsFo(|s|))*sqrt((1 +  0.5* ssqr(|s|))/(1 +  0.5* ssqr(s)))
  B(s) = A(s) * sqrt (ssqr(s) /ssqr(|s|))

  Ao(|s|)  = rmsFo(|s|) / (rmsFc(|s|) * sqrt(1 +  0.5* ssqr(|s|)))
    # overall true fall-off

  The scale factor to apply to Fobs is:
  Q(s) = (rmsFc(|s|)/rmsFo(s))  * sqrt(1 +  0.5 * ssqr(s))/(1 + ssqr(s))

  Normalize these to Q(s=0).
  """
  print ("\n",79*"=","\nAnalyzing anisotropy","\n",79*"=", file = out)

  cc_b_cart = get_overall_anisotropy(
     # Note: sets a_zero_values in overall_si
     overall_si = overall_si,
     n_bins_use = n_bins_use,
     out = out)

  aniso_info = get_aniso_info(
    f_array = f_array,
    overall_si = overall_si,
    si_list = si_list,
    expected_ssqr_list = expected_ssqr_list,
    s_value_list = s_value_list,
    n_bins_use = n_bins_use,
    direction_vectors = direction_vectors,
    weights_para_list = weights_para_list,
    resolution = resolution,
    weight_by_variance = weight_by_variance,
    minimum_sd = minimum_sd,
    maximum_ratio = maximum_ratio,
    use_dv_weighting = use_dv_weighting,
    cc_b_cart = cc_b_cart,
    out = out)


  aniso_info = estimate_s_matrix(aniso_info)

  # Total correction
  print("\nD matrix (Uncertainty-based scaling correction factor)\n"+
       "(Negative means uncertainties increase more in this direction)\n " +
      "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f) " %(
     tuple(aniso_info.dd_b_cart)), file = out)

  print("\nS matrix (anisotropy correction to scale factor)\n"+
       "(Negative means amplitudes fall off more in this direction)\n " +
      "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f) " %(
     tuple(aniso_info.ss_b_cart)), file = out)

  print("\nU matrix (Overall anisotropic fall-off relative to ideal)\n"+
       "(Negative means amplitudes fall off more in this direction)\n " +
      "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f) " %(
     tuple(aniso_info.uu_b_cart)), file = out)

  print("\nS scale values by direction vector", file = out)
  print("  D-min  A-zero   E**2   ", file = out, end = "")

  display_scale_values(
    aniso_info = aniso_info,
    n_display=n_display,
    values_by_dv = aniso_info.ss_values_by_dv,
    out = out)

  print("\nS scale values (calculated) by direction vector", file = out)
  print("  D-min  A-zero   E**2   ", file = out, end = "")

  display_scale_values(
      aniso_info = aniso_info,
      n_display=n_display,
      values_by_dv = aniso_info.ss_calc_values_by_dv,
      out = out)

  # Get summary information

  # Use qq_values instead of target_scale_factors (very similar)
  if update_scale_values_to_match_qq_values:
    print("\nUpdating scale values with qq_values (new version)\n",file = out)
    aniso_info = update_scale_values_with_qq_values(
      aniso_info = aniso_info,
      out = out)

  # Update scale factors from overall values if desired
  if update_scale_values_to_match_s_matrix:
    print("\nRecalculating scale values from S matrix and overall scale",
       file = out)
    # Recalculate target_scale_factors for each direction vector and save in
    #   corresponding si
    aniso_info = update_scale_values(  # update with S matrix
      aniso_info = aniso_info,
      out = out)


  scale_factor_info = group_args(
    group_args_type = """
     scale_factor_info:
  overall_si:  si overall
  scaling_info_list: si (scaling_info) objects, one for each direction vector
    each si:  si.target_scale_factors   # scale factors vs sthol2
              si.target_sthol2 # sthol2 values  d = 1/s = 0.25/sthol2**0.5
  a_zero_values:  isotropic fall-off of the data
  fo_b_cart_as_u_cart: anisotropic fall-off of data
  aa_b_cart_as_u_cart: anisotropy of the data
  bb_b_cart_as_u_cart: anisotropy of the uncertainties
  ss_b_cart_as_u_cart: anisotropy of the scale factors
  uu_b_cart_as_u_cart: anisotropic fall-off of data relative to ideal
  overall_scale: radial part of overall correction factor


    """,
    scaling_info_list = aniso_info.si_list,
    overall_si = overall_si,
    a_zero_values = overall_si.a_zero_values,
    fo_b_cart_as_u_cart = adptbx.b_as_u(aniso_info.fo_b_cart),
    aa_b_cart_as_u_cart = adptbx.b_as_u(aniso_info.aa_b_cart),
    bb_b_cart_as_u_cart = adptbx.b_as_u(aniso_info.bb_b_cart),
    ss_b_cart_as_u_cart = adptbx.b_as_u(aniso_info.ss_b_cart),
    uu_b_cart_as_u_cart = adptbx.b_as_u(aniso_info.uu_b_cart),
    overall_scale = overall_si.target_scale_factors,
   )
  return scale_factor_info

def estimate_s_matrix(aniso_info):

  """Estimate errors and calculate A,B,D S and Q matrices"""

  aniso_info = get_starting_sd_info(
    aniso_info = aniso_info,)

  aniso_info.ss_b_cart = get_aniso_from_scale_values(
    aniso_info,
    aniso_info.ss_scale_values,
    aniso_info.ss_sd_values,
  )

  aniso_info.ss_calc_values_by_dv = get_calc_values(
      aniso_info = aniso_info,
      b_cart = aniso_info.ss_b_cart,
      apply_b = True)

  # Update errors and run again:

  aniso_info = update_sd_values(aniso_info)

  # Now real thing for S matrix
  aniso_info.ss_b_cart = get_aniso_from_scale_values(
    aniso_info,
    aniso_info.ss_scale_values,
    aniso_info.ss_sd_values,
    b_cart = aniso_info.ss_b_cart,
  )
  aniso_info.ss_calc_values_by_dv = get_calc_values(
      aniso_info = aniso_info,
      b_cart = aniso_info.ss_b_cart,
      apply_b = True)

  # Repeat for A matrix
  aniso_info.aa_b_cart = get_aniso_from_scale_values(
    aniso_info,
    aniso_info.aa_scale_values,
    aniso_info.ss_sd_values,  # using sigmas for ss_b_cart
  )

  # Repeat for B matrix
  aniso_info.bb_b_cart = get_aniso_from_scale_values(
    aniso_info,
    aniso_info.bb_scale_values,
    aniso_info.ss_sd_values,  # using sigmas for ss_b_cart
  )

  # Repeat for D matrix
  aniso_info.dd_b_cart = get_aniso_from_scale_values(
    aniso_info,
    aniso_info.dd_scale_values,
    aniso_info.ss_sd_values,  # using sigmas for ss_b_cart
  )

  # Repeat for F (anisotropy of data)
  aniso_info.fo_b_cart = get_aniso_from_scale_values(
    aniso_info,
    aniso_info.fo_scale_values,
    aniso_info.ss_sd_values,
  )

  # Repeat for U (anisotropy of data relative to ideal)
  aniso_info.uu_b_cart = get_aniso_from_scale_values(
    aniso_info,
    aniso_info.uu_scale_values,
    aniso_info.ss_sd_values,
  )


  return aniso_info

def get_overall_anisotropy(overall_si,
  n_bins_use,
  out = sys.stdout):
  """Estimate overall fall-off with resolution of data
    and corrections for it including uncertainties"""

  print("\nEstimating overall fall-off with resolution of data"+
    " and correction including uncertainties", file = out)

  overall_si.cc_list.set_selected(overall_si.cc_list<1.e-4,1.e-4)
  overall_si.rms_fo_list.set_selected(overall_si.rms_fo_list <= 1.e-10, 1.e-10)

  overall_si.ssqr_values = 2. * ( -1 + 1./overall_si.cc_list)  # cc* is cc_list

  correction_factor = 1./flex.sqrt(1. + 0.5*overall_si.ssqr_values)
       # == sqrt(overall_si.cc_list)

  # Ao(|s|) = ( rmsFobs(|s|)/rmsFc(|s|) ) / sqrt (1 + 0.5*ssqr(|s|)) == C(|s|)
  # D(|s|) = 1./sqrt(1 +  0.5* ssqr(|s|))

  # So Ao(|s|) = ( rmsFobs(|s|)/rmsFc(|s|) )  * D(|s|)

  overall_si.dd_values = correction_factor

  overall_si.a_zero_values = correction_factor *(
      overall_si.rms_fo_list/ overall_si.rms_fc_list)

  info = get_effective_b(
    values = overall_si.a_zero_values,
    sthol2_values = overall_si.target_sthol2,
    n_bins_use = n_bins_use)

  print("Overall approximate B-value: %.2f A**2 (Scale = %.4f  rms = %.4f)" %(
    info.effective_b,
    info.b_zero,
    info.rms,), file = out)

  cc_b_cart = (
    -info.effective_b,-info.effective_b,-info.effective_b,0,0,0)

  print (" D-min   Ao(|s|)   Calc", file = out)
  for i in range(info.sthol2_values.size()):
    dd = 0.5/info.sthol2_values[i]**0.5
    print ("%6.2f  %7.4f   %6.4f  " %(
      dd, overall_si.a_zero_values[i], info.calc_values[i]),file = out)
  return cc_b_cart

def get_aniso_info(
    f_array = None,
    overall_si = None,
    si_list = None,
    expected_ssqr_list = None,
    s_value_list = None,
    n_bins_use = None,
    direction_vectors = None,
    weights_para_list = None,
    resolution = None,  # nominal resolution
    weight_by_variance = None,
    minimum_sd = None,
    maximum_ratio = None,
    small_ssqr_for_maximum_ratio = 5.,
    use_dv_weighting= None,
    cc_b_cart = None,
    out = sys.stdout):

  """
  Calculate anisotropic information for an array.

  ssqr(s) = (1/CC_half(s) - 1)  ...along any direction s
  Ao(|s|)**2 = rmsFobs(s*)**2 * (1 + 0.5*ssqr(s*))    ...along s*

  A(s) = (rmsFo(s)/rmsFc(|s|))*sqrt((1 +  0.5* ssqr(|s|))/(1 +  0.5* ssqr(s)))
  B(s) = A(s) * sqrt (ssqr(s) /ssqr(|s|))
  Target scale factors:
  Q(s) = (rmsFc(|s|)/rmsFo(s)) * sqrt(1 + 0.5* ssqr(s))/ (1 + ssqr(s))

  S(s) =  Q(s)/Q(|s|) # anisotropy of target scale factors

  """

  print("\nUsing overall average fall-off as baseline", file = out)

  overall_si.ssqr_values = 2. * ( -1 + 1./overall_si.cc_list)  # cc* is in cc_list

  fo_scale_values = flex.double()
  fo_values_by_dv = []
  aa_scale_values = flex.double()
  aa_values_by_dv = []
  bb_scale_values = flex.double()
  bb_values_by_dv = []
  dd_scale_values = flex.double()
  dd_values_by_dv = []
  ss_scale_values = flex.double()
  ss_values_by_dv = []
  uu_scale_values = flex.double()
  uu_values_by_dv = []

  indices = flex.miller_index()
  indices_by_dv = []
  # Limit range of scale factors relative to overall if errors are not small
  overall_ssqr_values_are_large = (overall_si.ssqr_values > small_ssqr_for_maximum_ratio)

  # We are going to recalculate target_scale_factors here with slightly
  #   different formula and call them qq_values
  #  target_scale_factors: T = (rmsFc(s)/rmsFo(s)) * 1/sqrt(1+0.5*ssqr(s))
  #  qq_values = (rmsFc(|s|)/rmsFo(s)) * sqrt(1 + 0.5* ssqr(s))/ (1 + ssqr(s))
  #  Normalized these:  qq_values = qq_values/qq_values[0]

  #  These differ by the ratio:
  #   (1 +  0.5 * ssqr(s))**2/(1 + ssqr(s))
  #  ... which has a value of 1 for ssqr=0 and 1.04 for ssqr=5...

  #  dd_values = 1/sqrt(1+0.5*ssqr)

  overall_si.qq_values = (
     overall_si.rms_fc_list *
     flex.sqrt(1. + 0.5 * overall_si.ssqr_values)) / (
     overall_si.rms_fo_list *
     (1. + overall_si.ssqr_values))

  overall_normalization = overall_si.target_scale_factors.min_max_mean().mean/ \
         max(1.e-10,overall_si.qq_values.min_max_mean().mean)
  overall_si.target_scale_factors *= overall_normalization

  overall_si.target_scale_factors.set_selected(
    overall_si.target_scale_factors < 1.e-10, 1.e-10)
  overall_si.qq_values.set_selected(
    overall_si.qq_values< 1.e-10, 1.e-10)

  overall_si.aa_values = flex.double(overall_si.qq_values.size(),1)
  overall_si.bb_values = flex.double(overall_si.qq_values.size(),1)

  for dv,si in zip(direction_vectors, si_list):
    si.cc_list.set_selected(si.cc_list < 1.e-10,1.e-10)
    si.dd_values = flex.sqrt(si.cc_list)  # 1/(1. + 0.5 * ssqr)**0.5
    si.ssqr_values = 2. * ( -1 + 1./si.cc_list)  # cc* is in cc_list
    si.qq_values = (
     si.rms_fc_list *
     flex.sqrt(1. + 0.5 * si.ssqr_values)) / (
     si.rms_fo_list *
     (1. + si.ssqr_values))

    # Normalize to overall_si.target_scale_factors[0] if possible
    local_normalization=\
        overall_si.qq_values.min_max_mean().mean/max(1.e-10,
       si.qq_values.min_max_mean().mean)
    si.qq_values *= local_normalization

    if overall_si.qq_values[0] > 1.e-10 and \
       si.qq_values[0] > 1.e-10 and \
       si.qq_values[0]/overall_si.qq_values[0] > 0.1 and\
       si.qq_values[0]/overall_si.qq_values[0] < 10:
      local_normalization = \
         overall_si.qq_values[0]/si.qq_values[0]
      si.qq_values *= local_normalization

    # Anisotropy of the data: A
    si.aa_values = ( si.rms_fo_list * si.dd_values ) / (
      overall_si.rms_fo_list * overall_si.dd_values )

    # Anisotropy of the errors: B
    si.bb_values = si.aa_values * flex.sqrt(
       si.ssqr_values/overall_si.ssqr_values)

    # Anisotropy of scale factors S:
    #     S = ratio of target scale factors to overall values
    si.ss_values = si.qq_values/overall_si.qq_values
    si.ss_values.set_selected(
       (overall_ssqr_values_are_large) & (si.ss_values > maximum_ratio),
           maximum_ratio)
    si.ss_values.set_selected(
       (overall_ssqr_values_are_large) & (si.ss_values < 1/maximum_ratio),
           1/maximum_ratio)

    # Estimate of true fall-off of data relative to ideal U
    # U(s) = rmsFo(s) / (rmsFc(|s|) sqrt(1 +  0.5* ssqr(s)))
    si.uu_values = si.rms_fo_list * si.dd_values / overall_si.rms_fc_list

    info = get_calculated_scale_factors( # get the indices
          s_value_list=s_value_list,
          cc_list = flex.double(s_value_list.size(),1),
          dv = dv,
          uc = f_array.unit_cell(),
           )
    indices_by_dv.append(info.indices)
    indices.extend(info.indices)
    fo_values_by_dv.append(si.rms_fo_list)
    fo_scale_values.extend(si.rms_fo_list)
    aa_values_by_dv.append(si.aa_values)
    aa_scale_values.extend(si.aa_values)
    bb_values_by_dv.append(si.bb_values)
    bb_scale_values.extend(si.bb_values)
    dd_values_by_dv.append(si.dd_values)
    dd_scale_values.extend(si.dd_values)
    ss_values_by_dv.append(si.ss_values)
    ss_scale_values.extend(si.ss_values)
    uu_values_by_dv.append(si.uu_values)
    uu_scale_values.extend(si.uu_values)

  return group_args(
    f_array = f_array,
    direction_vectors = direction_vectors,
    weights_para_list = weights_para_list,
    resolution = resolution,
    ssqr_values = overall_si.ssqr_values,
    s_value_list = s_value_list,
    minimum_sd = minimum_sd,
    maximum_ratio = maximum_ratio,
    weight_by_variance = weight_by_variance,
    use_dv_weighting = use_dv_weighting,
    n_bins = overall_si.target_sthol2.size(),
    n_dv = direction_vectors.size(),
    si_list = si_list,
    overall_si = overall_si,
    n_bins_use = n_bins_use,
    cc_b_cart = cc_b_cart,
    indices = indices,
    indices_by_dv = indices_by_dv,
    fo_scale_values = fo_scale_values,
    fo_values_by_dv = fo_values_by_dv,
    aa_scale_values = aa_scale_values,
    aa_values_by_dv = aa_values_by_dv,
    bb_scale_values = bb_scale_values,
    bb_values_by_dv = bb_values_by_dv,
    dd_scale_values = dd_scale_values,
    dd_values_by_dv = dd_values_by_dv,
    ss_scale_values = ss_scale_values,
    ss_values_by_dv = ss_values_by_dv,
    uu_scale_values = uu_scale_values,
    uu_values_by_dv = uu_values_by_dv,
  )

def get_starting_sd_info(aniso_info = None):

  """Try to get variances of aa, bb values within each resolution bin"""
  sd_ss = flex.double()
  for i in range(aniso_info.n_bins):
    ss_in_bin = flex.double()
    for k in range(aniso_info.n_dv):
      ss_in_bin.append(aniso_info.ss_values_by_dv[k][i])
    sd_value = ss_in_bin.sample_standard_deviation()
    # Try to catch cases where the values are stopped by bounds
    if ss_in_bin.min_max_mean().mean >= 0.9* aniso_info.maximum_ratio:
      sd_value = aniso_info.maximum_ratio
    elif ss_in_bin.min_max_mean().mean <= 1.1/aniso_info.maximum_ratio:
      sd_value = 1
    sd_ss.append(sd_value)
  ss_sd_values = flex.double()
  for k in range(aniso_info.n_dv):
    ss_sd_values.extend(sd_ss)
  average_sd = ss_sd_values.min_max_mean().mean
  minimum_sd_use = max(1.e-5,aniso_info.minimum_sd * average_sd)
  ss_sd_values.set_selected(ss_sd_values < minimum_sd_use, minimum_sd_use)
  aniso_info.ss_sd_values = ss_sd_values
  return aniso_info

def get_calc_values_with_dv_weighting(
   aniso_info,
   b_cart,
   direction_vector_k_list = None,
   apply_b = None):
  """Calculate anisotropic scale factors in each direction dv
      that removes anisotropic b_cart from data, applying weighting
      depending on the direction dv"""

  if direction_vector_k_list is None:
    direction_vector_k_list = list(range(aniso_info.n_dv))

  # Calculate values from matrices but weight as in dv weighting
  calc_values= get_scale_from_aniso_b_cart(
      f_array = aniso_info.f_array,
      indices = aniso_info.f_array.indices(),
      b_cart = b_cart,
      apply_b = apply_b)
  calc_values_by_dv=[]

  if direction_vector_k_list is None:
    direction_vector_k_list = list(range(aniso_info.n_dv))

  for k in range(aniso_info.n_dv):
    calc_values_by_dv.append(flex.double())

  for i_bin in aniso_info.f_array.binner().range_used():
    sel       = aniso_info.f_array.binner().selection(i_bin)
    for k in direction_vector_k_list:
      weights_para = aniso_info.weights_para_list[k]
      weights_para_sel = weights_para.select(sel)
      mean_weight = max(1.e-10,weights_para_sel.min_max_mean().mean)

      calc_values_by_dv[k].append( (calc_values.select(sel) *
         weights_para_sel).min_max_mean().mean/mean_weight)
  return calc_values_by_dv

def get_calc_values(
   aniso_info,
   b_cart,
   use_dv_weighting = None,
   direction_vector_k_list = None,
   apply_b = None):
  """Calculate anisotropic scale factors in each direction dv
      that removes anisotropic b_cart from data"""
  if use_dv_weighting is None:
    use_dv_weighting = aniso_info.use_dv_weighting
  if use_dv_weighting:
   values = get_calc_values_with_dv_weighting(
     aniso_info,
     b_cart,
     direction_vector_k_list = direction_vector_k_list,
     apply_b = apply_b,
    )
   return values

  if direction_vector_k_list is None:
    direction_vector_k_list = list(range(aniso_info.n_dv))

  # Calculate values from matrices
  calc_values_by_dv=[]
  for k in range(aniso_info.n_dv):
    if k in direction_vector_k_list:
      calc_values_by_dv.append(get_scale_from_aniso_b_cart(
        f_array = aniso_info.f_array,
        indices = aniso_info.indices_by_dv[k],
        b_cart = b_cart,
        apply_b = apply_b))
    else:
      calc_values_by_dv.append(flex.double())  # empty
  return calc_values_by_dv

def update_sd_values(aniso_info):

  """Update error estimates using calculated values as reference
  Try to get variances of aa, bb values within each resolution bin"""

  delta_ss = flex.double(aniso_info.n_bins,0)
  for k in range(aniso_info.n_dv):
    delta_ss += flex.pow2(
       aniso_info.ss_values_by_dv[k] - aniso_info.ss_calc_values_by_dv[k])
  scale = max(1,aniso_info.n_dv - 1)
  ss_sd_vs_resolution= flex.sqrt(delta_ss/scale) # sd vs resolution

  # And create sd vector
  ss_sd_values = flex.double()
  for k in range(aniso_info.n_dv):
    ss_sd_values.extend(ss_sd_vs_resolution)

  aniso_info.ss_sd_values = ss_sd_values
  return aniso_info

def get_scale_from_aniso_b_cart(f_array = None,
    indices = None,
    b_cart = None,
    apply_b = False):
  """Calculate anisotropic scale factor that removes anisotropic b_cart from
  data"""
  scale_values_array = f_array.customized_copy(
    data = flex.double(indices.size(),1),
    indices = indices)
  from mmtbx.scaling import absolute_scaling
  scale_values_array.set_observation_type_xray_amplitude()

  # NOTE: this removes b_cart from data. To apply it,  use apply_b=True
  if apply_b:
    overall_u_cart_to_apply = adptbx.b_as_u(tuple(-flex.double(tuple(b_cart))))
  else:
    overall_u_cart_to_apply = adptbx.b_as_u(tuple( flex.double(tuple(b_cart))))
  u_star= adptbx.u_cart_as_u_star(
    scale_values_array.unit_cell(),
     tuple(matrix.col(overall_u_cart_to_apply)))
  scaled_f_array = absolute_scaling.anisotropic_correction(
          scale_values_array,0.0, u_star ,must_be_greater_than=-0.0001)
  return scaled_f_array.data()

def update_scale_values_with_qq_values(
    aniso_info = None,
    out = sys.stdout):
  '''
   Replace target_scale_factors with qq_values
  '''

  aniso_info.overall_si.target_scale_factors = aniso_info.overall_si.qq_values

  for k in range(aniso_info.direction_vectors.size()):
    si = aniso_info.si_list[k]
    original_target_scale_factors = si.target_scale_factors
    si.target_scale_factors = si.qq_values

    print("Target scale factor replacement with qq_values for direction_vector",
       k,file=out)
    for sthol2,o,t in zip(
      si.target_sthol2,original_target_scale_factors,si.target_scale_factors):
      dd = 0.5/sthol2**0.5
      print(" %.2f %.3f %.3f" %(dd,o,t), file = out)

  return aniso_info

def update_scale_values(
    aniso_info = None,
    out = sys.stdout):
  '''
   Update the values of target_scale_factors using the information in
   aniso_info and s_matrix and target_scale_factors from overall_si (overall)
  '''
  for k in range(aniso_info.direction_vectors.size()):
    si = aniso_info.si_list[k]
    original_target_scale_factors = si.target_scale_factors
    normalized_scale_factors = get_any_list_from_any(
      aniso_info = aniso_info,
      any_matrix = aniso_info.ss_b_cart,
      use_dv_weighting = True, # REQUIRED
      target_sthol2_list = si.target_sthol2,
      direction_vector_k = k)
    si.target_scale_factors = normalized_scale_factors * \
       aniso_info.overall_si.target_scale_factors

    print("Target scale factor replacement with S matrix for direction_vector ",
       k,file=out)
    for sthol2,o,t in zip(
      si.target_sthol2,original_target_scale_factors,si.target_scale_factors):
      dd = 0.5/sthol2**0.5
      print(" %.2f %.3f %.3f" %(dd,o,t), file = out)

  return aniso_info

def get_any_list_from_any(
   aniso_info = None,
   any_matrix = None,
   target_sthol2_list = None,
   use_dv_weighting = None,
   direction_vector_k = None):
  '''
  Calculate estimated amplitude fall-off along this direction vector vs sthol2
  Value of any_matrix = A(s) (any anisotropic tensor)
  '''

  direction_vector_k_list = [direction_vector_k]

  calc_values_by_dv = get_calc_values(
    aniso_info,
    any_matrix,
    use_dv_weighting,
    direction_vector_k_list,
    apply_b = True)

  return calc_values_by_dv[direction_vector_k]

def display_scale_values(
    aniso_info = None,
    n_display=None,
    values_by_dv = None,
    out = sys.stdout):
  """Display scale values as function of direction vectors"""
  for k in range(min(aniso_info.n_dv,n_display)):
    print("  %4s " %(k+1), file = out, end = "")
  print("", file = out)
  for i in range(aniso_info.n_bins):
    dd = 0.5/aniso_info.overall_si.target_sthol2[i]**0.5
    print ("%6.2f  %6.3f %8.3f   " %(dd,
      aniso_info.overall_si.a_zero_values[i],
         aniso_info.overall_si.ssqr_values[i]),
        file = out, end = "")
    for k in range(min(aniso_info.n_dv,n_display)):
      print (" %5.2f " %(values_by_dv[k][i]), file = out, end= "")
    print("", file=out)

def get_aniso_from_scale_values(
   aniso_info = None,
   scale_values = None,
   sd_values = None,
   b_cart = None,
   ):
  """Estimate anistropy from scale values"""
  # If nothing present yet,
  #   Get a first cut for aniso_obj.b_cart with analyze_aniso
  if not b_cart:
    scale_values_array = aniso_info.f_array.customized_copy(
        data = scale_values,
        indices = aniso_info.indices)
    scaled_array,aniso_obj=analyze_aniso(
        b_iso=0,
        f_array=scale_values_array,resolution=aniso_info.resolution,
        remove_aniso=True,out=null_out())
    if not aniso_obj:
      return None
    elif aniso_obj.b_cart:
      b_cart = aniso_obj.b_cart
    else:
      b_cart = (0,0,0,0,0,0)

  # Optimize this (with weighting if weight_by_variance)

  ar = aniso_refinery(
    aniso_info,
    b_cart,
    scale_values,
    sd_values,
    eps = .01)
  ar.run()
  b_cart = ar.get_b()

  # If we want overall scale factor, it is here:
  # overall_scale = ar.get_overall_scale()
  # resid=ar.residual(b_cart)

  return b_cart

class aniso_refinery:
  """Refine anisotropic parameters to match scale factors along direction
   vectors dv. Specific for tools in this file"""
  def __init__(self,
    aniso_info,
    b_cart,
    scale_values,
    sd_values,
    overall_scale = 1.0,
    eps=0.01,
    tol=1.e-6,
    max_iterations=20,
    start_with_grid_search = True,
    grid_delta = 10,
    grid_n = 2,
    get_overall_scale_factor = True,
    ):

    self.aniso_info = aniso_info
    self.b_cart=b_cart
    self.overall_scale=overall_scale
    self.scale_values=scale_values
    self.sd_values=sd_values

    self.tol=tol
    self.eps=eps
    self.max_iterations=max_iterations
    self.get_overall_scale_factor=get_overall_scale_factor

    self.x = flex.double(b_cart)
    self.start_with_grid_search = start_with_grid_search
    self.grid_delta = grid_delta
    self.grid_n = grid_n



  def run(self):

    if self.start_with_grid_search:
      self.grid_search()
      b = self.get_b()
      resid = self.residual(b)
      if resid < self.tol: # done
        return

    else:
      best_b = self.get_b(self.x)
      resid=self.residual(self.x)
      scitbx.lbfgs.run(target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test_eps=self.tol,
                     max_iterations=self.max_iterations,
       ))

  def grid_search(self):
    best_b = self.get_b()
    working_b = self.get_b()
    best_resid = self.residual(best_b)
    for i in range(-self.grid_n,self.grid_n+1):
      for j in range(-self.grid_n,self.grid_n+1):
        for k in range(-self.grid_n,self.grid_n+1):
          b = list(
          matrix.col(working_b) +
          matrix.col((i*self.grid_delta,
           j*self.grid_delta,k*self.grid_delta,0,0,0)
           ))
          resid = self.residual(b)
          if resid < best_resid:
            best_resid = resid
            best_b = b
    self.x = flex.double(best_b)

  def show_result(self,out=sys.stdout):

    b=self.get_b()
    value = self.residual(b)
    return value

  def compute_functional_and_gradients(self):
    b = self.get_b()
    f = self.residual(b)
    g = self.gradients(b)
    return f, g

  def calculate_overall_scale(self,calc_values):
     if self.get_overall_scale_factor:
       return max(0, self.scale_values.min_max_mean().mean/max(1.e-10,
         calc_values.min_max_mean().mean))
     else:
       return 1

  def residual(self,b):
    calc_values_by_dv = get_calc_values(
      aniso_info = self.aniso_info,
      b_cart = b,
      apply_b = True)

    calc_values = flex.double()
    for c in calc_values_by_dv:
      calc_values.extend(c)

    # Get overall scale
    self.overall_scale = self.calculate_overall_scale(calc_values)

    diffs = calc_values*self.overall_scale - self.scale_values
    if self.aniso_info.weight_by_variance:
      diffs = diffs/self.sd_values
    diffs_dc=diffs.deep_copy()
    sd_dc= self.sd_values.deep_copy()
    if self.aniso_info.n_bins_use:  # just take first n_bins_use of each group
      i_pos = 0
      n_bins = calc_values_by_dv[0].size()
      diffs_use = flex.double()
      for k in range(len(calc_values_by_dv)):
        diffs_use.extend(diffs[i_pos:i_pos+self.aniso_info.n_bins_use])
        i_pos += n_bins
      diffs = diffs_use
    resid = diffs.rms()
    return resid

  def gradients(self,b):

    result = flex.double()
    for i in range(len(list(b))):
      rs = []
      for signed_eps in [self.eps, -self.eps]:
        params_eps = deepcopy(b)
        params_eps[i] += signed_eps
        rs.append(self.residual(params_eps))
      result.append((rs[0]-rs[1])/(2*self.eps))
    return result

  def get_overall_scale(self):
    return self.overall_scale

  def get_b(self):
    return list(self.x)

  def callback_after_step(self, minimizer):
    pass # can do anything here

class shell_aniso_refinery(aniso_refinery):
  def __init__(self,
    f_array,
    power = 1.,
    b_cart = None,
    overall_scale = None,
    eps=.1,
    tol=1.e-6,
    max_iterations=20,
    ):

    if not b_cart:
      b_cart = (0,0,0,0,0,0)
    if not overall_scale:
      overall_scale = f_array.data().min_max_mean().mean

    self.f_array = f_array.deep_copy()

    self.power=power
    self.tol=tol
    self.eps=eps
    self.max_iterations=max_iterations

    self.x = self.get_x(b_cart,overall_scale)

    self.start_with_grid_search = False

  def get_b(self, x = None):
    if x is None:
      x = self.x
    return [x[0],x[1],-1.*(x[0]+x[1]),x[2],x[3],x[4]]

  def get_overall_scale(self, x = None):
    if x is None:
      x = self.x
    return x[5]


  def get_x(self, b_cart, overall_scale):
    return flex.double((b_cart[0],b_cart[1],b_cart[3],b_cart[4],b_cart[5],
      overall_scale))

  def compute_functional_and_gradients(self):
    x = self.x
    f = self.residual(x)
    g = self.gradients(x)
    return f, g

  def gradients(self,x):

    result = flex.double()
    for i in range(len(list(x))):
      rs = []
      for signed_eps in [self.eps, -self.eps]:
        params_eps = deepcopy(x)
        params_eps[i] += signed_eps
        rs.append(self.residual(params_eps))
      result.append((rs[0]-rs[1])/(2*self.eps))
    return result


  def get_calc_array(self, x,  array_with_indices = None):
    b = self.get_b(x)
    overall_scale = self.get_overall_scale(x)
    from mmtbx.scaling import absolute_scaling
    overall_u_cart_to_apply = adptbx.b_as_u(tuple(
       -1*self.power*flex.double(tuple(b))))

    if not array_with_indices:
      array_with_indices = self.f_array

    u_star= adptbx.u_cart_as_u_star(
      array_with_indices.unit_cell(),
       tuple(matrix.col(overall_u_cart_to_apply)))
    fit_array = array_with_indices.customized_copy(
       data=flex.double(array_with_indices.size(),overall_scale))
    fit_array.set_observation_type_xray_amplitude()
    calc_array = absolute_scaling.anisotropic_correction(
          fit_array,0.0, u_star ,must_be_greater_than=-0.0001)

    return calc_array

  def residual(self,x):
    calc_array = self.get_calc_array(x)
    if not calc_array or calc_array.size()==0:
       return 1.e+30
    diffs = calc_array.data()-self.f_array.data()
    residual = flex.pow2(diffs).min_max_mean().mean
    rms = residual**0.5
    return rms

class iso_refinery(aniso_refinery):
  def __init__(self,
    values,
    sthol2_values,
    n_bins_use = None,
    b_iso = None,
    eps=0.01,
    tol=1.e-6,
    max_iterations=20,
    start_with_grid_search = False,
    grid_delta = 50,
    grid_n = 10,
    ):

    if not b_iso:
      b_iso = 0

    self.values = values
    self.sthol2_values = sthol2_values
    self.n_bins_use = n_bins_use
    self.tol=tol
    self.eps=eps
    self.max_iterations=max_iterations

    self.x = flex.double((b_iso,))

    self.start_with_grid_search = start_with_grid_search
    self.grid_delta = grid_delta
    self.grid_n = grid_n

  def get_info(self):
    b = self.get_b()
    b_value = b[0]
    info=get_b_calc(b_value,
       self.sthol2_values,
       self.values,
       n_bins_use = self.n_bins_use)
    return info

  def grid_search(self):
    best_b = self.get_b()
    working_b = self.get_b()
    grid_delta = self.grid_delta
    for cycle in range(self.grid_n):
      best_resid = self.residual(best_b)
      for i in range(-self.grid_n,self.grid_n+1):
        b = [working_b[0] + i*grid_delta]
        resid = self.residual(b)
        if resid < best_resid:
          best_resid = resid
          best_b = b
      self.x = flex.double(best_b)
      grid_delta = grid_delta/(self.grid_n*2)
      if grid_delta < self.eps: break

  def residual(self,b):
    b_value = b[0]
    info=get_b_calc(b_value,
       self.sthol2_values,
       self.values,
       n_bins_use = self.n_bins_use)
    return info.rms

def calculate_kurtosis(ma,phases,b,resolution,n_real=None,
    d_min_ratio=None):
  """Calculate the kurtosis of a map"""
  map_data=get_sharpened_map(ma,phases,b,resolution,n_real=n_real,
  d_min_ratio=d_min_ratio)
  return get_kurtosis(map_data.as_1d())



def get_aniso_obj_from_direction_vectors(
        f_array = None,
        resolution = None,
        direction_vectors = None,
        si_list = None,
        s_value_list = None,
        expected_rms_fc_list = None,
        invert_in_scaling = False,
        ):
      """Calculate anisotropic information from direction vectors and
         a list of CC values or expected rms FC values """
      scale_values = flex.double()
      indices = flex.miller_index()
      extra_b = -15.  # Just so ml scaling gives about the right answer for
                # constant number that are not really reflection data
      if not si_list:
        si_list = []
        for dv in direction_vectors:
          si_list.append(group_args(
            effective_b = None,
            effective_b_f_obs = None,
            b_zero = 0,
            cc_list =expected_rms_fc_list,
           ))

      for dv,si in zip(direction_vectors,si_list):
        info = get_calculated_scale_factors(
          s_value_list=s_value_list,
          effective_b=(None if (si.effective_b is None) else (si.effective_b + extra_b)),  # so xtriage gives about right
          b_zero = si.b_zero,
          cc_list = si.cc_list,
          dv = dv,
          uc = f_array.unit_cell(),
           )
        indices.extend(info.indices)
        scale_values.extend(info.scale_values)

      if invert_in_scaling:
        scale_values = 1/scale_values
      scale_values_array = f_array.customized_copy(
        data = scale_values,
        indices = indices)

      scaled_array,aniso_obj=analyze_aniso(
        b_iso=0,
        invert = invert_in_scaling,
        f_array=scale_values_array,resolution=resolution,
        remove_aniso=True,out=null_out())
      return aniso_obj

def get_resolution_for_aniso(s_value_list=None, si_list=None,
    minimum_ratio = 0.05):
  """Estimate a reasonable resolution to use for anisotropic calculations"""
  highest_d_min = None
  for si in si_list:
    first_rms_fo = None
    for rms_fo,s_value in  zip(si.rms_fo_list,s_value_list):
      d = 1/s_value
      if first_rms_fo is None:
        first_rms_fo = rms_fo
      else:
        ratio = rms_fo/max(1.e-10, first_rms_fo)
        if ratio < minimum_ratio and (
            highest_d_min is None or d > highest_d_min):
          highest_d_min = d
  return highest_d_min

def remove_values_if_necessary(f_values,max_ratio=100, min_ratio=0.01):
  """Make sure values are within a factor of 100 of low_res ones...if they are
  not it was probably something like zero values or near-zero values"""
  f_values=flex.double(f_values)
  low_res = f_values[:3].min_max_mean().mean
  new_values=flex.double()
  last_value = low_res
  for x in f_values:
    if x > max_ratio * low_res or x < min_ratio * low_res:
      new_values.append(last_value)
    else:
      new_values.append(x)
      last_value = x
  return new_values

def smooth_values(cc_values, max_relative_rms=10, n_smooth = None,
    skip_first_frac = 0.1,
    overall_values = None,
    smooth = True,
    max_ratio = 2.,
    min_ratio = 0.,
      ): # normally do not smooth the very first ones
  '''  If overall_values are supplied, make sure all values are within
      min_ratio to max_ratio of those values'''


  if overall_values and (max_ratio is not None or min_ratio is not None):
    new_cc_values = flex.double()
    for cc,cc_overall in zip(cc_values,overall_values):
      new_cc_values.append(max(
        min_ratio*cc_overall, min(max_ratio*cc_overall, cc)))
    cc_values = new_cc_values

  if not smooth:
    return cc_values

  skip_first = max (1, int(0.5+skip_first_frac*cc_values.size()))

  if n_smooth:
    # smooth with window of n_smooth
    new_cc_values = flex.double()
    for i in range(cc_values.size()):
      if i < skip_first:
        new_cc_values.append(cc_values[i])
      else:
        sum=0.
        sum_n=0.
        for j in range(-n_smooth, n_smooth+1):
          weight = 1/(1+(abs(j)/n_smooth)) # just less as we go out
          k = i+j
          if k < 0 or k >= cc_values.size(): continue
          sum += cc_values[k] * weight
          sum_n += weight
        new_cc_values.append(sum/max(1.e-10,sum_n))
    return new_cc_values

  # Smooth values in cc_values  max_relative_rms is avg rms / avg delta
  if relative_rms(cc_values) <= max_relative_rms:
    return cc_values
  for i in range(1,cc_values.size()//2):
    smoothed_cc_values=smooth_values(cc_values,n_smooth=i)
    if relative_rms(smoothed_cc_values) <= max_relative_rms:
      return smoothed_cc_values
  smoothed_cc_values=smooth_values(cc_values,n_smooth=cc_values.size()//2)
  return smoothed_cc_values


def relative_rms(cc_values):
  """Calculate relative rms of CC values"""
  diffs = cc_values[:-1] - cc_values[1:]
  avg_delta = abs(diffs.min_max_mean().mean)
  rms = diffs.sample_standard_deviation()
  return rms/max(1.e-10,avg_delta)

def complete_cc_analysis(
       direction_vector,
       cc_list,
       rms_fc_list,
       rms_fo_list,
       ratio_list,
       scale_using_last,
       max_cc_for_rescale,
       optimize_b_eff,
       is_model_based,
       s_value_list,
       cc_cut,
       max_possible_cc,
       fraction_complete,
       min_fraction_complete,
       low_res_bins,
       si,
       b_eff,
       input_info,
       cutoff_after_last_high_point,
       expected_rms_fc_list,
       expected_ssqr_list,
       overall_si,
       out):
  """Analyze CC values as function of direction vectors. Return
    an si object with this analysis"""
  if scale_using_last: # rescale to give final value average==0
    cc_list,baseline=rescale_cc_list(
       cc_list=cc_list,scale_using_last=scale_using_last,
       max_cc_for_rescale=max_cc_for_rescale)
    if baseline is None: # don't use it
      scale_using_last=None

  if expected_ssqr_list and overall_si:
    # Replace with expected scaled by avg fo/fo; if
    #  this direction has small fo then errors are bigger
    directional_ssqr_list = expected_ssqr_list * flex.pow2(
       overall_si.rms_fo_list/rms_fo_list)
    cc_list = 1/(max(1.e-10,0.5*directional_ssqr_list + 1))

  original_cc_list=deepcopy(cc_list)
  if is_model_based: # jut smooth cc if nec
    fitted_cc=get_fitted_cc(
      cc_list=cc_list,s_value_list=s_value_list,cc_cut=cc_cut,
      scale_using_last=scale_using_last,
      cutoff_after_last_high_point = cutoff_after_last_high_point,)
    cc_list=fitted_cc
    text=" FIT "
  else:
    cc_list=estimate_cc_star(cc_list=cc_list,s_value_list=s_value_list,
      cc_cut=cc_cut,scale_using_last=scale_using_last)
    text=" CC* "


  if not max_possible_cc:
    max_possible_cc=0.01
  if si.target_scale_factors: # not using these
    max_possible_cc=1.
    fraction_complete=1.
  elif (not is_model_based):
    max_possible_cc=1.
    fraction_complete=1.
  else:
    # Define overall CC based on model completeness (CC=sqrt(fraction_complete))

    if fraction_complete is None:
      fraction_complete=max_possible_cc**2

      print(
     "Estimated fraction complete is %5.2f based on low_res CC of %5.2f" %(
          fraction_complete,max_possible_cc), file=out)
    else:
      print(
      "Using fraction complete value of %5.2f "  %(fraction_complete), file=out)
      max_possible_cc=fraction_complete**0.5

  if optimize_b_eff and is_model_based:
    ''' Find b_eff that maximizes expected map-model-cc to model with B=0'''
    best_b_eff = b_eff
    best_weighted_cc = get_target_scale_factors(
      cc_list=cc_list,
      rms_fo_list=rms_fo_list,
      ratio_list=ratio_list,
      b_eff=b_eff,
      max_possible_cc=max_possible_cc,
      **input_info()).weighted_cc
    for i in range(20):
      b_eff_working = 0.1 * i * b_eff
      weighted_cc=get_target_scale_factors(
         cc_list=cc_list,
         rms_fo_list=rms_fo_list,
         ratio_list=ratio_list,
         b_eff = b_eff_working,
         max_possible_cc=max_possible_cc,
         **input_info()).weighted_cc
      if weighted_cc > best_weighted_cc:
        best_b_eff = b_eff_working
        best_weighted_cc = weighted_cc
    print("Optimized effective B value: %.3f A**2 " %(best_b_eff),file=out)
    b_eff = best_b_eff

  info = get_target_scale_factors(
      cc_list=cc_list,
      rms_fo_list=rms_fo_list,
      ratio_list=ratio_list,
      b_eff=b_eff,
      max_possible_cc=max_possible_cc,
      **input_info())
  target_scale_factors = info.target_scale_factors

  if direction_vector:
    print ("\n Analysis for direction vector (%5.2f, %5.2f, %5.2f): "% (
      direction_vector), file = out)

  if info.effective_b:
    print("\nEffective B value for CC*: %.3f A**2 " %(
        info.effective_b),file=out)

  if fraction_complete < min_fraction_complete:
    print("\nFraction complete (%5.2f) is less than minimum (%5.2f)..." %(
      fraction_complete,min_fraction_complete) + "\nSkipping scaling", file=out)
    target_scale_factors=flex.double(target_scale_factors.size()*(1.0,))
  print ("\nAverage CC: %.3f" %(cc_list.min_max_mean().mean),file=out)
  print("\nScale factors vs resolution:", file=out)
  print("Note 1: CC* estimated from sqrt(2*CC/(1+CC))", file=out)
  print("Note 2: CC estimated by fitting (smoothing) for values < %s" %(cc_cut), file=out)
  print("Note 3: Scale = A  CC*  rmsFc/rmsFo (A is normalization)", file=out)
  print("  d_min     rmsFo       rmsFc    CC      %s  Scale " %(
      text), file=out)

  for sthol2,scale,rms_fo,cc,rms_fc,orig_cc in zip(
     input_info.target_sthol2,target_scale_factors,rms_fo_list,
      cc_list,rms_fc_list,
      original_cc_list):
     print("%7.2f  %9.1f  %9.1f %7.3f  %7.3f  %5.2f " %(
       0.5/sthol2**0.5,rms_fo,rms_fc,orig_cc,cc,scale),
        file=out)

  si.target_scale_factors=target_scale_factors
  si.target_sthol2=input_info.target_sthol2
  si.d_min_list=input_info.d_min_list
  si.original_cc_list=original_cc_list # this is CC(half-map1, half_map2)
  si.cc_list=cc_list
  si.rms_fo_list = rms_fo_list
  si.rms_fc_list = rms_fc_list
  si.low_res_cc = cc_list[:low_res_bins].min_max_mean().mean # low-res average
  si.effective_b = info.effective_b
  si.effective_b_f_obs = info.effective_b_f_obs
  si.b_zero = info.b_zero
  si.rms = info.rms
  si.expected_rms_fc_list = expected_rms_fc_list

  return si

def get_sel_para(f_array, direction_vector, minimum_dot = 0.70):
    """get selections based on |dot(normalized_indices, direction_vector)|"""
    u = f_array.unit_cell()
    rcvs = u.reciprocal_space_vector(f_array.indices())
    norms = rcvs.norms()
    norms.set_selected((norms == 0),1)
    index_directions = rcvs/norms
    sel = (flex.abs(index_directions.dot(direction_vector)) > minimum_dot)
    return sel

def get_normalized_weights_para(f_array,direction_vectors, dv,
    include_all_in_lowest_bin = None):
    """Get normalized weights parallel to direction vectors dv"""

    sum_weights = flex.double(f_array.size(),0)
    current_weights = None
    for direction_vector in direction_vectors:
      weights = get_weights_para(f_array, direction_vector,
        include_all_in_lowest_bin = include_all_in_lowest_bin)
      if direction_vector == dv:
        current_weights = weights
      sum_weights += weights
    sum_weights.set_selected((sum_weights <= 1.e-10), 1.e-10)
    return current_weights * (1/sum_weights)

def get_weights_para(f_array, direction_vector,
       weight_by_cos = True,
       min_dot = 0.7,
       very_high_dot = 0.9,
       pre_factor_scale= 10,
       include_all_in_lowest_bin = None):
    """Get weights parallel to direction vectors dv"""
    u = f_array.unit_cell()
    rcvs = u.reciprocal_space_vector(f_array.indices())
    norms = rcvs.norms()
    norms.set_selected((norms == 0),1)
    index_directions = rcvs/norms
    if weight_by_cos:

      weights = flex.abs(index_directions.dot(direction_vector))
      sel = (weights < min_dot)
      weights.set_selected(sel,0)

      weights += (1-very_high_dot)  # move very_high to 1.0
      sel = (weights > 1)
      weights.set_selected(sel,1) # now from (min_dot+(1-very_high_dot) to 1)

      weights = (weights - 1 ) * pre_factor_scale
      sel = (weights > -20)  &  (weights < 20)
      weights.set_selected(sel, flex.exp(weights.select(sel)))
      weights.set_selected(~sel,0)


    else:
      weights = flex.double(index_directions.size(),0)
      sel = (flex.abs(index_directions.dot(direction_vector)) > min_dot)
      weights.set_selected(sel,1)
    if include_all_in_lowest_bin:
      i_bin = 1
      sel       = f_array.binner().selection(i_bin)
      weights.set_selected(sel, 1.0)  # full weights on low-res in all directions
    return weights

def get_nearest_lattice_points(unit_cell, reciprocal_space_vectors):
  """Return lattice points nearest each of a set of reciprocal-space vectors"""
  lattice_points=flex.vec3_double()
  v = matrix.sqr(unit_cell.fractionalization_matrix()).inverse()
  for x in reciprocal_space_vectors:
    lattice_points.append(v*x)
  lattice_points=lattice_points.iround()
  return lattice_points

def get_target_scale_factors(
     f_array = None,
     ratio_list = None,
     rms_fo_list = None,
     cc_list = None,
     n_list = None,
     target_sthol2 = None,
     d_min_list = None,
     max_possible_cc = None,
     pseudo_likelihood = None,
     equalize_power = None,
     is_model_based = None,
     skip_scale_factor = None,
     maximum_scale_factor = None,
     b_eff = None,
     out = sys.stdout):
  """Estimate target scale factors and return a group_args object
   containing them"""

  weighted_cc = 0
  weighted = 0

  target_scale_factors=flex.double()
  sum_w=0.
  sum_w_scale=0.
  for i_bin in f_array.binner().range_used():
    index=i_bin-1
    ratio=ratio_list[index]
    cc=cc_list[index]
    sthol2=target_sthol2[index]
    d_min=d_min_list[index]


    corrected_cc=max(0.00001,min(1.,cc/max(1.e-10,max_possible_cc)))

    if (not is_model_based): # NOTE: cc is already converted to cc*
      scale_on_fo=ratio * corrected_cc  # CC* * rmsFc/rmsFo
    elif b_eff is not None:
      if pseudo_likelihood:
        scale_on_fo=(cc/max(0.001,1-cc**2))
      else: # usual
        scale_on_fo=ratio * min(1.,
          max(0.00001,corrected_cc) * math.exp(min(20.,sthol2*b_eff)) )
    else:
      scale_on_fo=ratio * min(1.,max(0.00001,corrected_cc))

    w = n_list[index]*(rms_fo_list[index])**2
    sum_w += w
    sum_w_scale += w * scale_on_fo**2
    target_scale_factors.append(scale_on_fo)
    weighted_cc += n_list[index]*rms_fo_list[index] * scale_on_fo * corrected_cc
    weighted += n_list[index]*rms_fo_list[index] * scale_on_fo

  weighted_cc = weighted_cc/max(1.e-10,weighted)
  if not pseudo_likelihood and not skip_scale_factor: # normalize
    avg_scale_on_fo = (sum_w_scale/max(1.e-10,sum_w))**0.5
    if equalize_power and avg_scale_on_fo>1.e-10:
      # XXX do not do this if only 1 bin has values > 0
      scale_factor = 1/avg_scale_on_fo
    else: # usual
      scale_factor=1./target_scale_factors.min_max_mean().max
    target_scale_factors=\
      target_scale_factors*scale_factor
  if maximum_scale_factor and \
     target_scale_factors.min_max_mean().max > maximum_scale_factor:
    truncated_scale_factors = flex.double()
    for x in target_scale_factors:
      truncated_scale_factors.append(min(maximum_scale_factor,x ))
    target_scale_factors = truncated_scale_factors

  # Get effective B for cc_list and target_sthol2
  info = get_effective_b(values = cc_list,
    sthol2_values = target_sthol2)
  effective_b = info.effective_b
  b_zero= info.b_zero
  rms= info.rms

  # Also get effective_b for amplitudes
  amplitude_info = get_effective_b(values = rms_fo_list/max(
     1.e-10,rms_fo_list[0]),
    sthol2_values = target_sthol2)
  effective_b_f_obs = amplitude_info.effective_b

  return group_args(
    target_scale_factors = target_scale_factors,
    weighted_cc = weighted_cc,
    effective_b = effective_b,
    effective_b_f_obs = effective_b_f_obs,
    b_zero = b_zero,
    rms = rms
  )

def get_effective_b(values = None,
      sthol2_values = None,
       n_bins_use = None):
  """Estimate effective B corresponding to falloff of a set of values at
   sin**2(theta)/lambda**2 values"""
  ir = iso_refinery(
    values,
    sthol2_values,
    n_bins_use,
    start_with_grid_search=True)
  ir.run()
  return ir.get_info()

  """
  info object looks like:
    group_args_type =
       'get_b_iso values calc_values = b_zero * exp(-b_value* sthol2)',
    effective_b = refined_b_value,
    b_zero =  scale_factor_b_zero,
    sthol2_values = sthol2_values,
    values = values,
    calc_values = calc_values,
    n_bins_use = n_bins_use,
    rms = rms
  """


def get_b_calc( b_value, sthol2_values, values, n_bins_use = None):
  '''Calculate values using b_value and sthol2_values, normalize to
   first element of values.
   If n_bins_use is set, just use that many for rms value'''
  import math
  sum= 0.
  sumx= 0.
  sumy= 0.
  sum2= 0.
  sumn=0.
  calc_values = flex.double()
  for sthol2, value in zip (sthol2_values, values):
    calc_values.append(math.exp(max(-20.,min(20.,
      - b_value* sthol2))))
  b_zero = values[0]/calc_values[0]
  calc_values *= b_zero
  from libtbx.test_utils import approx_equal
  assert approx_equal(values[0],calc_values[0])
  if n_bins_use is None:
    n_bins_use = values.size()
  delta = values[:n_bins_use]-calc_values[:n_bins_use]
  rms = delta.rms()
  return group_args(
    group_args_type = \
      'get_b_iso values calc_values = b_zero * exp(-b_value* sthol2)',
    effective_b = b_value,
    b_zero=b_zero,
    sthol2_values = sthol2_values,
    values=values,
    calc_values=calc_values,
    n_bins_use = n_bins_use,
    rms=rms)


def analyze_aniso(f_array=None,map_coeffs=None,b_iso=None,resolution=None,
     get_remove_aniso_object=True,
     invert = False,
     remove_aniso=None, aniso_obj=None, out=sys.stdout):
  """Analyze anisotropy in f_array or map_coeffs.
     Optionally remove anisotropy and set all directions to mean value.
  return array and analyze_aniso_object
  resolution can be None, b_iso can be None
  if remove_aniso is None, just analyze and return original array"""

  if map_coeffs:  # convert to f and apply
    from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
    f_local,phases_local=map_coeffs_as_fp_phi(map_coeffs)
    f_local,f_local_aa=analyze_aniso(f_array=f_local,
       aniso_obj=aniso_obj,
       get_remove_aniso_object=get_remove_aniso_object,
       remove_aniso=remove_aniso, resolution=resolution,out=out)
    return f_local.phase_transfer(phase_source=phases_local,deg=True),f_local_aa

  elif not get_remove_aniso_object:
    return f_array,aniso_obj # don't do anything

  else:  # have f_array and resolution
    if not aniso_obj:
      aniso_obj=analyze_aniso_object()
      aniso_obj.set_up_aniso_correction(f_array=f_array,d_min=resolution,
        b_iso=b_iso, invert = invert)

    if remove_aniso and aniso_obj and aniso_obj.b_cart:
      f_array=aniso_obj.apply_aniso_correction(f_array=f_array)
      print("Removing anisotropy with b_cart=(%7.2f,%7.2f,%7.2f)\n" %(
        aniso_obj.b_cart[:3]), file=out)
    return f_array,aniso_obj

def scale_amplitudes(model_map_coeffs=None,
    map_coeffs=None,
    external_map_coeffs=None,
    first_half_map_coeffs=None,
    second_half_map_coeffs=None,
    si=None,resolution=None,overall_b=None,
    fraction_complete=None,
    min_fraction_complete=0.05,
    map_calculation=True,
    verbose=False,
    out=sys.stdout):
  """Figure out resolution_dependent sharpening to optimally
  match map and model. Then apply it as usual.
  if second_half_map_coeffs instead of model,
  use second_half_map_coeffs same as
  normalized model map_coeffs, except that the target fall-off should be
  skipped (could use fall-off based on a dummy model...)"""

  if model_map_coeffs and (
      not first_half_map_coeffs or not second_half_map_coeffs):
    is_model_based=True
  elif si.target_scale_factors or (
       first_half_map_coeffs and second_half_map_coeffs) or (
        external_map_coeffs):
    is_model_based=False
  else:
    assert map_coeffs
    if si.is_model_sharpening():
      is_model_based=True
    else:
      is_model_based=False

  if si.verbose and not verbose:
    verbose=True

  # if si.target_scale_factors is set, just use those scale factors

  from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi,get_b_iso

  f_array,phases=map_coeffs_as_fp_phi(map_coeffs)

  (d_max,d_min)=f_array.d_max_min(d_max_is_highest_defined_if_infinite=True)
  n_bins_use = si.n_bins
  while not f_array.binner():
    try:
      f_array.setup_binner(n_bins=si.n_bins,d_max=d_max,d_min=d_min)
      f_array.binner().require_all_bins_have_data(min_counts=1,
        error_string="Please use a lower value of n_bins")
    except Exception as e:
      if n_bins_use > 1:
        n_bins_use -= 1
      else:
        raise Exception(e)

  if resolution is None:
    resolution=si.resolution
  if resolution is None:
    raise Sorry("Need resolution for model sharpening")

  obs_b_iso=get_b_iso(f_array,d_min=resolution)
  print("\nEffective b_iso of observed data: %6.1f A**2" %(obs_b_iso), file=out)

  if not si.target_scale_factors: # get scale factors if don't already have them
    si=calculate_fsc(si=si,
      f_array=f_array,  # just used for binner
      map_coeffs=map_coeffs,
      model_map_coeffs=model_map_coeffs,
      first_half_map_coeffs=first_half_map_coeffs,
      second_half_map_coeffs=second_half_map_coeffs,
      external_map_coeffs=external_map_coeffs,
      resolution=resolution,
      fraction_complete=fraction_complete,
      min_fraction_complete=min_fraction_complete,
      is_model_based=is_model_based,
      cc_cut=si.cc_cut,
      scale_using_last=si.scale_using_last,
      max_cc_for_rescale=si.max_cc_for_rescale,
      pseudo_likelihood=si.pseudo_likelihood,
      verbose=verbose,
      out=out)
    # now si.target_scale_factors array are the scale factors

  # Now create resolution-dependent coefficients from the scale factors

  if not si.target_scale_factors: # nothing to do
    print("\nNo scaling applied", file=out)
    map_data=calculate_map(map_coeffs=map_coeffs,n_real=si.n_real)
    return map_and_b_object(map_data=map_data)
  elif not map_calculation:
    return map_and_b_object()
  else:  # apply scaling
    if si.pseudo_likelihood:
      print("Normalizing structure factors", file=out)
      f_array=quasi_normalize_structure_factors(f_array,set_to_minimum=0.01,
        pseudo_likelihood=si.pseudo_likelihood)
      f_array.setup_binner(n_bins=si.n_bins,d_max=d_max,d_min=d_min)
    map_and_b=apply_target_scale_factors(
      f_array=f_array,phases=phases,resolution=resolution,
      target_scale_factors=si.target_scale_factors,
      n_real=si.n_real,
      out=out)
    return map_and_b

def apply_target_scale_factors(f_array=None,phases=None,
   resolution=None,target_scale_factors=None,
   n_real=None,
   return_map_coeffs=None,out=sys.stdout):
    """Apply target_scale_factors to f_array.  Returns map_and_b object
    with map_data and starting and final b_iso values. """
    from cctbx.maptbx.segment_and_split_map import get_b_iso
    f_array_b_iso=get_b_iso(f_array,d_min=resolution)
    scale_array=f_array.binner().interpolate(
      target_scale_factors, 1) # d_star_power=1
    scaled_f_array=f_array.customized_copy(data=f_array.data()*scale_array)
    scaled_f_array_b_iso=get_b_iso(scaled_f_array,d_min=resolution)
    print("\nInitial b_iso for "+\
      "map: %5.1f A**2     After applying scaling: %5.1f A**2" %(
      f_array_b_iso,scaled_f_array_b_iso), file=out)
    new_map_coeffs=scaled_f_array.phase_transfer(phase_source=phases,deg=True)
    assert new_map_coeffs.size() == f_array.size()
    if return_map_coeffs:
      return new_map_coeffs

    map_data=calculate_map(map_coeffs=new_map_coeffs,n_real=n_real)
    return map_and_b_object(map_data=map_data,starting_b_iso=f_array_b_iso,
      final_b_iso=scaled_f_array_b_iso)
def calculate_map(map_coeffs=None,crystal_symmetry=None,n_real=None):
  """Calculate a map from map_coeffs and crystal_symmetry"""
  if crystal_symmetry is None: crystal_symmetry=map_coeffs.crystal_symmetry()
  from cctbx.development.create_models_or_maps import get_map_from_map_coeffs
  map_data=get_map_from_map_coeffs(
     map_coeffs=map_coeffs,crystal_symmetry=crystal_symmetry, n_real=n_real)
  return map_data

def get_sharpened_map(ma=None,phases=None,b=None,resolution=None,
    n_real=None,d_min_ratio=None):
  """Calculate a sharpened map from amplitudes in ma, the values in b,
   and phases"""
  assert n_real is not None
  sharpened_ma=adjust_amplitudes_linear(ma,b[0],b[1],b[2],resolution=resolution,
     d_min_ratio=d_min_ratio)
  new_map_coeffs=sharpened_ma.phase_transfer(phase_source=phases,deg=True)
  map_data=calculate_map(map_coeffs=new_map_coeffs,n_real=n_real)
  return map_data

def calculate_match(target_sthol2=None,target_scale_factors=None,b=None,resolution=None,d_min_ratio=None,rmsd=None,fraction_complete=None):
  """Calculate residual between target_scale_factors and those calculated
   using effective b values sthol2_1 2 and 3"""

  if fraction_complete is None:
    pass # XXX not implemented for fraction_complete

  if rmsd is None:
    rmsd=resolution/3.
    print("Setting rmsd to %5.1f A based on resolution of %5.1f A" %(
       rmsd,resolution), file=out)

  if rmsd is None:
    b_eff=None
  else:
    b_eff=8*3.14159*rmsd**2

  d_min=d_min_ratio*resolution
  sthol2_2=0.25/resolution**2
  sthol2_1=sthol2_2*0.5
  sthol2_3=0.25/d_min**2
  b0=0.0
  b1=b[0]
  b2=b[1]
  b3=b[2]
  b3_use=b3+b2

  resid=0.
  import math
  value_list=flex.double()
  scale_factor_list=flex.double()

  for sthol2,scale_factor in zip(target_sthol2,target_scale_factors):
    if sthol2 > sthol2_2:
      value=b2+(sthol2-sthol2_2)*(b3_use-b2)/(sthol2_3-sthol2_2)
    elif sthol2 > sthol2_1:
      value=b1+(sthol2-sthol2_1)*(b2-b1)/(sthol2_2-sthol2_1)
    else:
      value=b0+(sthol2-0.)*(b1-b0)/(sthol2_1-0.)

    value=math.exp(value)
    if b_eff is not None:
      value=value*math.exp(-sthol2*b_eff)
    value_list.append(value)
    scale_factor_list.append(scale_factor)
  mean_value=value_list.min_max_mean().mean
  mean_scale_factor=scale_factor_list.min_max_mean().mean
  ratio=mean_scale_factor/mean_value
  value_list=value_list*ratio
  delta_list=value_list-scale_factor_list
  delta_sq_list=delta_list*delta_list
  resid=delta_sq_list.min_max_mean().mean
  return resid

def calculate_adjusted_sa(ma,phases,b,
    resolution=None,
    d_min_ratio=None,
    solvent_fraction=None,
    region_weight=None,
    max_regions_to_test=None,
    sa_percent=None,
    fraction_occupied=None,
    wrapping=None,
    n_real=None):
  """Calculate adjusted surface area for a map"""
  map_data=get_sharpened_map(ma,phases,b,resolution,n_real=n_real,
    d_min_ratio=d_min_ratio)
  from cctbx.maptbx.segment_and_split_map import score_map

  si=score_map(
    map_data=map_data,
    solvent_fraction=solvent_fraction,
    fraction_occupied=fraction_occupied,
    wrapping=wrapping,
    sa_percent=sa_percent,
    region_weight=region_weight,
    max_regions_to_test=max_regions_to_test,
    out=null_out())
  return si.adjusted_sa

def get_kurtosis(data=None):
  """Calculate kurtosis of a list of values in a flex.double array"""
  mean=data.min_max_mean().mean
  sd=data.sample_standard_deviation()
  x=data-mean
  return (x**4).min_max_mean().mean/sd**4

class analyze_aniso_object:
  """ Object to hold information about anisotropy in structure factors"""
  def __init__(self):

    self.b_cart=None
    self.b_cart_aniso_removed=None
    self.b_iso=None

  def set_up_aniso_correction(self,f_array=None,b_iso=None,d_min=None,
     b_cart_to_remove = None, invert = False):

    assert f_array is not None
    if not d_min:
      (d_max,d_min)=f_array.d_max_min(d_max_is_highest_defined_if_infinite=True)

    if b_cart_to_remove and b_iso:
      self.b_cart=b_cart_to_remove
      self.b_cart_aniso_removed = [ -b_iso, -b_iso, -b_iso, 0, 0, 0] # change
      self.b_iso=b_iso
    else:
      from cctbx.maptbx.segment_and_split_map import get_b_iso
      b_mean,aniso_scale_and_b=get_b_iso(f_array,d_min=d_min,
        return_aniso_scale_and_b=True)
      if not aniso_scale_and_b or not aniso_scale_and_b.b_cart:
        return # failed

      if b_iso is None:
        b_iso=b_mean  # use mean
      self.b_iso=b_iso

      self.b_cart=aniso_scale_and_b.b_cart  # current
      self.b_cart_aniso_removed = [ -b_iso, -b_iso, -b_iso, 0, 0, 0] # change
      if invert:
        self.b_cart = tuple([-x for x in self.b_cart])

      # ready to apply

  def apply_aniso_correction(self,f_array=None):

    if self.b_cart is None or self.b_cart_aniso_removed is None:
      return f_array  # nothing to do

    from mmtbx.scaling import absolute_scaling

    u_star= adptbx.u_cart_as_u_star(
      f_array.unit_cell(), adptbx.b_as_u( self.b_cart) )

    f_array.set_observation_type_xray_amplitude()

    u_star_aniso_removed = adptbx.u_cart_as_u_star(
      f_array.unit_cell(), adptbx.b_as_u( self.b_cart_aniso_removed  ) )

    no_aniso_array = absolute_scaling.anisotropic_correction(
      f_array,0.0, u_star ,must_be_greater_than=-0.0001)

    no_aniso_array.set_observation_type_xray_amplitude()
    no_aniso_array = absolute_scaling.anisotropic_correction(
      no_aniso_array,0.0,u_star_aniso_removed,must_be_greater_than=-0.0001)

    no_aniso_array=no_aniso_array.set_observation_type( f_array)
    return no_aniso_array


class refinery:
  """Refine machinery specific for auto-sharpening"""
  def __init__(self,ma,phases,b,resolution,
    residual_target=None,
    solvent_fraction=None,
    region_weight=None,
    max_regions_to_test=None,
    sa_percent=None,
    fraction_occupied=None,
    wrapping=None,
    eps=0.01,
    tol=0.01,
    max_iterations=20,
    n_real=None,
    target_sthol2=None,
    target_scale_factors=None,
    d_min_ratio=None,
    rmsd=None,
    fraction_complete=None,
    dummy_run=False):

    self.ma=ma
    self.n_real=n_real
    self.phases=phases
    self.resolution=resolution
    self.d_min_ratio=d_min_ratio
    self.rmsd=rmsd
    self.fraction_complete=fraction_complete

    self.target_sthol2=target_sthol2
    self.target_scale_factors=target_scale_factors

    self.tol=tol
    self.eps=eps
    self.max_iterations=max_iterations

    self.solvent_fraction=solvent_fraction
    self.region_weight=region_weight
    self.max_regions_to_test=max_regions_to_test
    self.residual_target=residual_target
    self.sa_percent=sa_percent
    self.fraction_occupied=fraction_occupied
    self.wrapping=wrapping

    self.x = flex.double(b)

  def run(self):

    scitbx.lbfgs.run(target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test_eps=self.tol,
                     max_iterations=self.max_iterations,
       ))

  def show_result(self,out=sys.stdout):

    b=self.get_b()
    value = -1.*self.residual(b)
    print("Result: b1 %7.2f b2 %7.2f b3 %7.2f resolution %7.2f %s: %7.3f" %(
     b[0],b[1],b[2],self.resolution,self.residual_target,value), file=out)

    if self.ma:
      self.sharpened_ma=adjust_amplitudes_linear(
         self.ma,b[0],b[1],b[2],resolution=self.resolution,
         d_min_ratio=self.d_min_ratio)
    else:
      self.sharpened_ma=None
    return value

  def compute_functional_and_gradients(self):
    b = self.get_b()
    f = self.residual(b)
    g = self.gradients(b)
    return f, g

  def residual(self,b,restraint_weight=100.):

    if self.residual_target=='kurtosis':
      resid=-1.*calculate_kurtosis(self.ma,self.phases,b,self.resolution,
         n_real=self.n_real,d_min_ratio=self.d_min_ratio)

    elif self.residual_target=='adjusted_sa':
      resid=-1.*calculate_adjusted_sa(self.ma,self.phases,b,
        resolution=self.resolution,
        d_min_ratio=self.d_min_ratio,
        solvent_fraction=self.solvent_fraction,
        region_weight=self.region_weight,
        max_regions_to_test=self.max_regions_to_test,
        sa_percent=self.sa_percent,
        fraction_occupied=self.fraction_occupied,
        wrapping=self.wrapping,n_real=self.n_real)

    elif self.residual_target=='model':
      resid=calculate_match(target_sthol2=self.target_sthol2,
        target_scale_factors=self.target_scale_factors,
        b=b,
        resolution=self.resolution,
        d_min_ratio=self.d_min_ratio,
        emsd=self.rmsd,
        fraction_complete=self.complete)

    else:
      raise Sorry("residual_target must be kurtosis or adjusted_sa or match_target")

    # put in restraint so b[1] is not bigger than b[0]
    if b[1]>b[0]:  resid+=(b[1]-b[0])*restraint_weight
    # put in restraint so b[2] <=0
    if b[2]>0:  resid+=b[2]*restraint_weight
    return resid

  def gradients(self,b):

    result = flex.double()
    for i in range(len(list(b))):
      rs = []
      for signed_eps in [self.eps, -self.eps]:
        params_eps = deepcopy(b)
        params_eps[i] += signed_eps
        rs.append(self.residual(params_eps))
      result.append((rs[0]-rs[1])/(2*self.eps))
    return result

  def get_b(self):
    return list(self.x)

  def callback_after_step(self, minimizer):
    pass # can do anything here

def calculate_kurtosis(ma,phases,b,resolution,n_real=None,
    d_min_ratio=None):
  """Calculate kurtosis from amplitudes and phases for a map"""
  map_data=get_sharpened_map(ma,phases,b,resolution,n_real=n_real,
  d_min_ratio=d_min_ratio)
  return get_kurtosis(map_data.as_1d())

def run(map_coeffs=None,
  b=[0,0,0],
  sharpening_info_obj=None,
  resolution=None,
  residual_target=None,
  solvent_fraction=None,
  region_weight=None,
  max_regions_to_test=None,
  sa_percent=None,
  fraction_occupied=None,
  n_bins=None,
  eps=None,
  wrapping=False,
  n_real=False,
  target_sthol2=None,
  target_scale_factors=None,
  d_min_ratio=None,
  rmsd=None,
  fraction_complete=None,
  normalize_amplitudes_in_resdep=None,
  out=sys.stdout):

  """Run map auto-sharpening"""
  if sharpening_info_obj:
    solvent_fraction=sharpening_info_obj.solvent_fraction
    wrapping=sharpening_info_obj.wrapping
    n_real=sharpening_info_obj.n_real
    fraction_occupied=sharpening_info_obj.fraction_occupied
    sa_percent=sharpening_info_obj.sa_percent
    region_weight=sharpening_info_obj.region_weight
    max_regions_to_test=sharpening_info_obj.max_regions_to_test
    residual_target=sharpening_info_obj.residual_target
    resolution=sharpening_info_obj.resolution
    d_min_ratio=sharpening_info_obj.d_min_ratio
    rmsd=sharpening_info_obj.rmsd
    fraction_complete=sharpening_info_obj.fraction_complete
    eps=sharpening_info_obj.eps
    n_bins=sharpening_info_obj.n_bins
    normalize_amplitudes_in_resdep= \
    sharpening_info_obj.normalize_amplitudes_in_resdep
  else:
    from cctbx.maptbx.segment_and_split_map import sharpening_info
    sharpening_info_obj=sharpening_info()


  if map_coeffs:
    phases=map_coeffs.phases(deg=True)
    ma=map_coeffs.as_amplitude_array()
  else:
    phases=None
    ma=None

  # set some defaults
  if residual_target is None: residual_target='kurtosis'

  assert (solvent_fraction is not None ) or residual_target=='kurtosis'
  assert resolution is not None

  if residual_target=='adjusted_sa' and solvent_fraction is None:
    raise Sorry("Solvent fraction is required for residual_target=adjusted_sa")

  if eps is None and residual_target=='kurtosis':
    eps=0.01
  elif eps is None:
    eps=0.5

  if fraction_complete is None:
    pass # XXX not implemented

  if rmsd is None:
    rmsd=resolution/3.
    print("Setting rmsd to %5.1f A based on resolution of %5.1f A" %(
       rmsd,resolution), file=out)

  if fraction_occupied is None: fraction_occupied=0.20
  if region_weight is None: region_weight=20.
  if sa_percent is None: sa_percent=30.
  if n_bins is None: n_bins=20
  if max_regions_to_test is None: max_regions_to_test=30


  # Get initial value

  best_b=b
  print("Getting starting value ...",residual_target, file=out)
  refined = refinery(ma,phases,b,resolution,
    residual_target=residual_target,
    solvent_fraction=solvent_fraction,
    region_weight=region_weight,
    sa_percent=sa_percent,
    max_regions_to_test=max_regions_to_test,
    fraction_occupied=fraction_occupied,
    wrapping=wrapping,
    n_real=n_real,
    target_sthol2=target_sthol2,
    target_scale_factors=target_scale_factors,
    d_min_ratio=d_min_ratio,
    rmsd=rmsd,
    fraction_complete=fraction_complete,
    eps=eps)


  starting_result=refined.show_result(out=out)
  print("Starting value: %7.2f" %(starting_result), file=out)

  if ma:
    (d_max,d_min)=ma.d_max_min(d_max_is_highest_defined_if_infinite=True)
    ma.setup_binner(n_bins=n_bins,d_max=d_max,d_min=d_min)
    if normalize_amplitudes_in_resdep:
      print("Normalizing structure factors...", file=out)
      ma=quasi_normalize_structure_factors(ma,set_to_minimum=0.01)
  else:
    assert resolution is not None

  refined = refinery(ma,phases,b,resolution,
    residual_target=residual_target,
    solvent_fraction=solvent_fraction,
    region_weight=region_weight,
    max_regions_to_test=max_regions_to_test,
    sa_percent=sa_percent,
    fraction_occupied=fraction_occupied,
    wrapping=wrapping,
    n_real=n_real,
    target_sthol2=target_sthol2,
    target_scale_factors=target_scale_factors,
    d_min_ratio=d_min_ratio,
    rmsd=rmsd,
    fraction_complete=fraction_complete,
    eps=eps)

  starting_normalized_result=refined.show_result(out=out)
  print("Starting value after normalization: %7.2f" %(
     starting_normalized_result), file=out)
  best_sharpened_ma=ma
  best_result=starting_normalized_result
  best_b=refined.get_b()

  refined.run()

  final_result=refined.show_result(out=out)
  print("Final value: %7.2f" %(
     final_result), file=out)

  if final_result>best_result:
    best_sharpened_ma=refined.sharpened_ma
    best_result=final_result
    best_b=refined.get_b()
  print("Best overall result: %7.2f: " %(best_result), file=out)

  sharpening_info_obj.resolution_dependent_b=best_b
  return best_sharpened_ma,phases


if (__name__ == "__main__"):
  args=sys.argv[1:]
  residual_target='kurtosis'
  if 'adjusted_sa' in args:
    residual_target='adjusted_sa'
  resolution=2.9 # need to set this as nominal resolution
  # get data
  map_coeffs=get_amplitudes(args)

  new_map_coeffs=run(map_coeffs=map_coeffs,
    resolution=resolution,
    residual_target=residual_target)
  mtz_dataset=new_map_coeffs.as_mtz_dataset(column_root_label="FWT")
  mtz_dataset.mtz_object().write(file_name='sharpened.mtz')
