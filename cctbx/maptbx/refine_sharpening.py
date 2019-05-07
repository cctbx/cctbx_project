from __future__ import division, print_function

from six.moves import range
import sys,os
from libtbx.utils import Sorry
from cctbx.array_family import flex
from copy import deepcopy

from cctbx.array_family import flex
import scitbx.lbfgs
import math
from cctbx.maptbx.segment_and_split_map import map_and_b_object

def write_mtz(ma=None,phases=None,file_name=None):
  mtz_dataset=ma.as_mtz_dataset(column_root_label="FWT")
  mtz_dataset.add_miller_array(miller_array=phases,column_types="P", column_root_label="PHWT")
  mtz_dataset.mtz_object().write(file_name=file_name)

# XXX copied these from cctbx.miller; made small change to catch weird case
#  where normalizations are negative.  Just multiply these *-1 and it seems to
#  close to what we want. Figure this out later...
# XXX Also set means=1 not mean square = 1

def amplitude_quasi_normalisations(ma, d_star_power=1, set_to_minimum=None,
    pseudo_likelihood=False):  # Used for pseudo-likelihood calculation
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
    normalisations = amplitude_quasi_normalisations(ma, d_star_power,
       set_to_minimum=set_to_minimum,pseudo_likelihood=pseudo_likelihood)
    if pseudo_likelihood:
      print("Norms:")
      for n,d in zip(normalisations[:100],ma.data()[:100]): print(n,d)

    q = ma.data() / normalisations.data()
    from cctbx.miller import array
    return array(ma, q)

def get_array(file_name=None,labels=None):

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
  # Return effective b values at sthol2_1 2 and 3
  # see adjust_amplitudes_linear below

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
  # do something to the amplitudes.
  #   b1=delta_b at midway between d=inf and d=resolution,b2 at resolution,
  #   b3 at d_min (added to b2)
  # pseudo-B at position of b1= -b1/sthol2_2= -b1*4*resolution**2
  #  or...b1=-pseudo_b1/(4*resolution**2)
  #  typical values of say b1=1 at 3 A -> pseudo_b1=-4*9=-36

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

def get_scale_factors(f_array,target_scale_factors=None):
  scale_array=flex.double(f_array.data().size()*(-1.,))
  assert len(target_scale_factors)==len(list(f_array.binner().range_used()))
  for i_bin in f_array.binner().range_used():
    sel       = f_array.binner().selection(i_bin)
    scale_array.set_selected(sel,target_scale_factors[i_bin-1])
  assert scale_array.count(-1.)==0
  return scale_array

def get_model_map_coeffs_normalized(pdb_inp=None,
   si=None,
   f_array=None,
   overall_b=None,
   resolution=None,
   out=sys.stdout):
  if not pdb_inp: return None

  # define Wilson B for the model
  if overall_b is None:
    if si.resolution:
      overall_b=si.get_target_b_iso()*si.target_b_iso_model_scale
    else:
      overall_b=0
    print("Setting Wilson B = %5.1f A" %(overall_b), file=out)

  # create model map using same coeffs
  from cctbx.maptbx.segment_and_split_map import get_f_phases_from_model
  model_map_coeffs=get_f_phases_from_model(
     pdb_inp=pdb_inp,
     f_array=f_array,
     overall_b=overall_b,
     k_sol=si.k_sol,
     b_sol=si.b_sol,
     out=out)

  from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi,get_b_iso
  model_f_array,model_phases=map_coeffs_as_fp_phi(model_map_coeffs)
  (d_max,d_min)=f_array.d_max_min()
  model_f_array.setup_binner(n_bins=si.n_bins,d_max=d_max,d_min=d_min)

  # Set overall_b....
  final_b_iso=get_b_iso(model_f_array,d_min=resolution)
  print("Effective b_iso of "+\
     "adjusted model map:  %6.1f A**2" %(final_b_iso), file=out)
  model_map_coeffs_normalized=model_f_array.phase_transfer(
     phase_source=model_phases,deg=True)
  return model_map_coeffs_normalized

def get_b_eff(si=None,out=sys.stdout):
  if si.rmsd is None:
    b_eff=None
  else:
    b_eff=8*3.14159*si.rmsd**2
    print("Setting b_eff for fall-off at %5.1f A**2 based on model error of %5.1f A" \
       %( b_eff,si.rmsd), file=out)
  return b_eff

def cc_fit(sthol_list=None,scale=None,value_zero=None,baseline=None,
     scale_using_last=None):
  # for scale_using_last, fix final value at zero
  fit=flex.double()
  s_zero=sthol_list[0]
  for s in sthol_list:
    fit.append(value_zero*math.exp(-scale*(s-s_zero)))
  if scale_using_last:
    fit=fit-fit[-1]
  return fit

def get_baseline(scale=None,scale_using_last=None,max_cc_for_rescale=None):
  if not scale_using_last:
    return 0
  else:
    baseline=min(0.99,max(0.,scale[-scale_using_last:].min_max_mean().mean))
    if baseline > max_cc_for_rescale:
       return None
    else:
       return baseline

def fit_cc(cc_list=None,sthol_list=None,
    scale_min=None,scale_max=None,n_tries=None,scale_using_last=None):
  # find value of scale in range scale_min,scale_max that minimizes rms diff
  # between cc_list and cc_list[0]*exp(-scale*(sthol-sthol_list[0]))
  # for scale_using_last, require it to go to zero at end

  best_scale=None
  best_rms=None
  for i in range(n_tries):
    scale=scale_min+(scale_max-scale_min)*i/n_tries
    fit=cc_fit(sthol_list=sthol_list,scale=scale,value_zero=cc_list[0],
       scale_using_last=scale_using_last)
    fit=fit-cc_list
    rms=fit.norm()
    if best_rms is None or rms<best_rms:
      best_rms=rms
      best_scale=scale
  return cc_fit(sthol_list=sthol_list,scale=best_scale,value_zero=cc_list[0],
      scale_using_last=scale_using_last)

def get_fitted_cc(cc_list=None,sthol_list=None, cc_cut=None,
   scale_using_last=None):
  # only do this if there is some value of s where cc is at least 2*cc_cut or
  #  (1-c_cut/2), whichever is smaller
  min_cc=min(2*cc_cut,1-0.5*cc_cut)
  if cc_list.min_max_mean().max < min_cc:
    return cc_list
  # find first point after point where cc>=min_cc that cc<=cc_cut
  #   then back off by 1 point
  found_high=False
  s_cut=None
  i_cut=0
  for s,cc in zip(sthol_list,cc_list):
    if cc > min_cc:
      found_high=True
    if found_high and cc < cc_cut:
      s_cut=s
      break
    i_cut+=1
  if s_cut is None or i_cut==0:
    return cc_list

  #Fit remainder
  sthol_remainder_list=sthol_list[i_cut:]
  cc_remainder_list=cc_list[i_cut:]
  n=cc_remainder_list.size()
  scale_min=10 # soft
  scale_max=500 # hard
  n_tries=200

  fitted_cc_remainder_list=fit_cc(
     cc_list=cc_remainder_list,sthol_list=sthol_remainder_list,
     scale_min=scale_min,scale_max=scale_max,n_tries=n_tries,
     scale_using_last=scale_using_last)

  new_cc_list=cc_list[:i_cut]
  new_cc_list.extend(fitted_cc_remainder_list)
  return new_cc_list

def estimate_cc_star(cc_list=None,sthol_list=None, cc_cut=None,
    scale_using_last=None):
  # cc ~ sqrt(2*half_dataset_cc/(1+half_dataset_cc))
  # however for small cc the errors are very big and we think cc decreases
  #  rapidly towards zero once cc is small
  # So find value of sthol_zero that gives cc about cc_cut...for
  #   sthol >sthol_zero use fit of cc_zero * exp(-falloff*(sthol-sthol_zero))
  #  for scale_using_last set: subtract off final values so it goes to zero.

  fitted_cc=get_fitted_cc(
    cc_list=cc_list,sthol_list=sthol_list,cc_cut=cc_cut,
    scale_using_last=scale_using_last)

  cc_star_list=flex.double()
  for cc in fitted_cc:
     cc=max(cc,0.)
     cc_star=(2.*cc/(1.+cc))**0.5
     cc_star_list.append(cc_star)
  return cc_star_list


def rescale_cc_list(cc_list=None,scale_using_last=None,
    max_cc_for_rescale=None):
  baseline=get_baseline(scale=cc_list,
      scale_using_last=scale_using_last,
      max_cc_for_rescale=max_cc_for_rescale)
  if baseline is None:
     return cc_list,baseline

  # replace cc with (cc-baseline)/(1-baseline)
  scaled_cc_list=flex.double()
  for cc in cc_list:
    scaled_cc_list.append((cc-baseline)/(1-baseline))
  return scaled_cc_list,baseline


def calculate_fsc(si=None,
     f_array=None,  # just used for binner
     map_coeffs=None,
     model_map_coeffs=None,
     first_half_map_coeffs=None,
     second_half_map_coeffs=None,
     resolution=None,
     fraction_complete=None,
     min_fraction_complete=None,
     is_model_based=None,
     cc_cut=None,
     scale_using_last=None,
     max_cc_for_rescale=None,
     pseudo_likelihood=False,
     skip_scale_factor=False,
     verbose=None,
     out=sys.stdout):

  # calculate anticipated fall-off of model data with resolution
  if si.rmsd is None and is_model_based:
    si.rmsd=resolution*si.rmsd_resolution_factor
    print("Setting rmsd to %5.1f A based on resolution of %5.1f A" %(
       si.rmsd,resolution), file=out)


  # get f and model_f vs resolution and FSC vs resolution and apply
  # scale to f_array and return sharpened map
  dsd = f_array.d_spacings().data()
  from cctbx.maptbx.segment_and_split_map import map_coeffs_to_fp

  if is_model_based:
    mc1=map_coeffs
    mc2=model_map_coeffs
    fo_map=map_coeffs # scale map_coeffs to model_map_coeffs*FSC
    fc_map=model_map_coeffs
    b_eff=get_b_eff(si=si,out=out)
  else: # half_dataset
    mc1=first_half_map_coeffs
    mc2=second_half_map_coeffs
    fo_map=map_coeffs # scale map_coeffs to cc*
    fc_map=model_map_coeffs
    b_eff=None


  ratio_list=flex.double()
  target_sthol2=flex.double()
  cc_list=flex.double()
  sthol_list=flex.double()
  d_min_list=flex.double()
  rms_fo_list=flex.double()
  rms_fc_list=flex.double()
  max_possible_cc=None
  for i_bin in f_array.binner().range_used():
    sel       = f_array.binner().selection(i_bin)
    d         = dsd.select(sel)
    if d.size()<1:
      raise Sorry("Please reduce number of bins (no data in bin "+
        "%s) from current value of %s" %(i_bin,f_array.binner().n_bins_used()))
    d_min     = flex.min(d)
    d_max     = flex.max(d)
    d_avg     = flex.mean(d)
    n         = d.size()
    m1        = mc1.select(sel)
    m2        = mc2.select(sel)
    cc        = m1.map_correlation(other = m2)


    if fo_map:
      fo        = fo_map.select(sel)
      f_array_fo=map_coeffs_to_fp(fo)
      rms_fo=f_array_fo.data().norm()
    else:
      rms_fo=1.

    if fc_map:
      fc        = fc_map.select(sel)
      f_array_fc=map_coeffs_to_fp(fc)
      rms_fc=f_array_fc.data().norm()
    else:
      rms_fc=1.

    sthol2=0.25/d_avg**2
    ratio_list.append(max(1.e-10,rms_fc)/max(1.e-10,rms_fo))
    target_sthol2.append(sthol2)
    if cc is None: cc=0.
    cc_list.append(cc)
    sthol_list.append(1/d_avg)
    d_min_list.append(d_min)
    rms_fo_list.append(rms_fo)
    rms_fc_list.append(rms_fc)

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

  if scale_using_last: # rescale to give final value average==0
    cc_list,baseline=rescale_cc_list(
       cc_list=cc_list,scale_using_last=scale_using_last,
       max_cc_for_rescale=max_cc_for_rescale)
    if baseline is None: # don't use it
      scale_using_last=None


  original_cc_list=deepcopy(cc_list)
  if not is_model_based:  # calculate cc* for half-dataset cc
    cc_list=estimate_cc_star(cc_list=cc_list,sthol_list=sthol_list,
      cc_cut=cc_cut,scale_using_last=scale_using_last)

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

      print("Estimated fraction complete is %5.2f based on low_res CC of %5.2f" %(
          fraction_complete,max_possible_cc), file=out)
    else:
      print("Using fraction complete value of %5.2f "  %(fraction_complete), file=out)
      max_possible_cc=fraction_complete**0.5

  target_scale_factors=flex.double()
  for i_bin in f_array.binner().range_used():
    index=i_bin-1
    ratio=ratio_list[index]
    cc=cc_list[index]
    sthol2=target_sthol2[index]
    d_min=d_min_list[index]

    corrected_cc=max(0.00001,min(1.,cc/max_possible_cc))

    if (not is_model_based): # cc is already cc*
      scale_on_fo=ratio * corrected_cc
    elif b_eff is not None:
      if pseudo_likelihood:
        scale_on_fo=(cc/max(0.001,1-cc**2))
      else: # usual
        scale_on_fo=ratio * min(1.,
          max(0.00001,corrected_cc) * math.exp(min(20.,sthol2*b_eff)) )
    else:
      scale_on_fo=ratio * min(1.,max(0.00001,corrected_cc))


    target_scale_factors.append(scale_on_fo)

  if not pseudo_likelihood and not skip_scale_factor: # normalize
    scale_factor=1./target_scale_factors.min_max_mean().max
    target_scale_factors=\
      target_scale_factors*scale_factor
    print("Scale factor A: %.5f" %(scale_factor), file=out)
  if fraction_complete < min_fraction_complete:
    print("\nFraction complete (%5.2f) is less than minimum (%5.2f)..." %(
      fraction_complete,min_fraction_complete) + "\nSkipping scaling", file=out)
    target_scale_factors=flex.double(target_scale_factors.size()*(1.0,))

  print("\nScale factors vs resolution:", file=out)
  print("Note 1: CC* estimated from sqrt(2*CC/(1+CC))", file=out)
  print("Note 2: CC estimated by fitting (smoothing) for values < %s" %(cc_cut), file=out)
  print("Note 3: Scale = A  CC*  rmsFc/rmsFo (A is normalization)", file=out)
  print("  d_min     rmsFo       rmsFc    CC       CC*   Scale", file=out)

  for sthol2,scale,rms_fo,cc,rms_fc,orig_cc in zip(
     target_sthol2,target_scale_factors,rms_fo_list,cc_list,rms_fc_list,
      original_cc_list):
     print("%7.1f  %9.1f  %9.1f %7.3f  %7.3f  %5.2f" %(
       0.5/sthol2**0.5,rms_fo,rms_fc,orig_cc,cc,scale), file=out)

  si.target_scale_factors=target_scale_factors
  si.target_sthol2=target_sthol2
  si.d_min_list=d_min_list
  si.cc_list=cc_list

  return si

def analyze_aniso(f_array=None,map_coeffs=None,b_iso=None,resolution=None,
     get_remove_aniso_object=True,
     remove_aniso=None, aniso_obj=None, out=sys.stdout):
  # optionally remove anisotropy and set all directions to mean value
  #  return array and analyze_aniso_object
  #  resolution can be None, b_iso can be None
  #  if remove_aniso is None, just analyze and return original array

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
      aniso_obj.set_up_aniso_correction(f_array=f_array,d_min=resolution)

    if remove_aniso and aniso_obj and aniso_obj.b_cart:
      f_array=aniso_obj.apply_aniso_correction(f_array=f_array)
      print("Removing anisotropy with b_cart=(%7.2f,%7.2f,%7.2f)\n" %(
        aniso_obj.b_cart[:3]), file=out)
    return f_array,aniso_obj

def scale_amplitudes(model_map_coeffs=None,
    map_coeffs=None,
    first_half_map_coeffs=None,
    second_half_map_coeffs=None,
    si=None,resolution=None,overall_b=None,
    fraction_complete=None,
    min_fraction_complete=0.05,
    map_calculation=True,
    verbose=False,
    out=sys.stdout):
  # Figure out resolution_dependent sharpening to optimally
  #  match map and model. Then apply it as usual.
  #  if second_half_map_coeffs instead of model, use second_half_map_coeffs same as
  #    normalized model map_coeffs, except that the target fall-off should be
  #    skipped (could use fall-off based on a dummy model...)

  if model_map_coeffs and (
      not first_half_map_coeffs or not second_half_map_coeffs):
    is_model_based=True
  else:
    assert si.target_scale_factors or (
       first_half_map_coeffs and second_half_map_coeffs)
    is_model_based=False

  if si.verbose and not verbose:
    verbose=True

  # if si.target_scale_factors is set, just use those scale factors

  from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi,get_b_iso

  f_array,phases=map_coeffs_as_fp_phi(map_coeffs)

  (d_max,d_min)=f_array.d_max_min()
  if not f_array.binner():
    f_array.setup_binner(n_bins=si.n_bins,d_max=d_max,d_min=d_min)
    f_array.binner().require_all_bins_have_data(min_counts=1,
      error_string="Please use a lower value of n_bins")

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
    from cctbx.maptbx.segment_and_split_map import get_b_iso
    f_array_b_iso=get_b_iso(f_array,d_min=resolution)
    scale_array=get_scale_factors(f_array,
        target_scale_factors=target_scale_factors)
    scaled_f_array=f_array.customized_copy(data=f_array.data()*scale_array)
    scaled_f_array_b_iso=get_b_iso(scaled_f_array,d_min=resolution)
    print("\nInitial b_iso for "+\
      "map: %5.1f A**2     After applying scaling: %5.1f A**2" %(
      f_array_b_iso,scaled_f_array_b_iso), file=out)
    new_map_coeffs=scaled_f_array.phase_transfer(phase_source=phases,deg=True)
    if return_map_coeffs:
      return new_map_coeffs

    map_data=calculate_map(map_coeffs=new_map_coeffs,n_real=n_real)
    return map_and_b_object(map_data=map_data,starting_b_iso=f_array_b_iso,
      final_b_iso=scaled_f_array_b_iso)
def calculate_map(map_coeffs=None,crystal_symmetry=None,n_real=None):

  if crystal_symmetry is None: crystal_symmetry=map_coeffs.crystal_symmetry()
  from cctbx.maptbx.segment_and_split_map import get_map_from_map_coeffs
  map_data=get_map_from_map_coeffs(
     map_coeffs=map_coeffs,crystal_symmetry=crystal_symmetry, n_real=n_real)
  return map_data

def get_sharpened_map(ma=None,phases=None,b=None,resolution=None,
    n_real=None,d_min_ratio=None):
  assert n_real is not None
  sharpened_ma=adjust_amplitudes_linear(ma,b[0],b[1],b[2],resolution=resolution,
     d_min_ratio=d_min_ratio)
  new_map_coeffs=sharpened_ma.phase_transfer(phase_source=phases,deg=True)
  map_data=calculate_map(map_coeffs=new_map_coeffs,n_real=n_real)
  return map_data

def calculate_match(target_sthol2=None,target_scale_factors=None,b=None,resolution=None,d_min_ratio=None,rmsd=None,fraction_complete=None):

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

  map_data=get_sharpened_map(ma,phases,b,resolution,n_real=n_real,
    d_min_ratio=d_min_ratio)
  from cctbx.maptbx.segment_and_split_map import score_map

  from libtbx.utils import null_out
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
  mean=data.min_max_mean().mean
  sd=data.standard_deviation_of_the_sample()
  x=data-mean
  return (x**4).min_max_mean().mean/sd**4

class analyze_aniso_object:
  def __init__(self):

    self.b_iso=None # target b_iso, default is mean of existing
    self.b_cart=None
    self.b_cart_aniso_removed=None

  def set_up_aniso_correction(self,f_array=None,b_iso=None,d_min=None):

    assert f_array is not None
    if not d_min:
      (d_max,d_min)=f_array.d_max_min()

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

    # ready to apply

  def apply_aniso_correction(self,f_array=None):

    if self.b_cart is None or self.b_cart_aniso_removed is None:
      return f_array  # nothing to do

    from mmtbx.scaling import absolute_scaling
    from cctbx import adptbx

    u_star= adptbx.u_cart_as_u_star(
      f_array.unit_cell(), adptbx.b_as_u( self.b_cart) )

    u_star_aniso_removed = adptbx.u_cart_as_u_star(
      f_array.unit_cell(), adptbx.b_as_u( self.b_cart_aniso_removed  ) )

    no_aniso_array = absolute_scaling.anisotropic_correction(
      f_array,0.0, u_star ,must_be_greater_than=-0.0001)

    no_aniso_array = absolute_scaling.anisotropic_correction(
      no_aniso_array,0.0,u_star_aniso_removed,must_be_greater_than=-0.0001)

    no_aniso_array=no_aniso_array.set_observation_type( f_array)
    return no_aniso_array


class refinery:
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
    (d_max,d_min)=ma.d_max_min()
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
