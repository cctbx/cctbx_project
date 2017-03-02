from __future__ import division

import sys,os
from libtbx.utils import Sorry
from cctbx.array_family import flex
from copy import deepcopy

from cctbx.array_family import flex
import scitbx.lbfgs

def write_mtz(ma=None,phases=None,file_name=None):
  mtz_dataset=ma.as_mtz_dataset(column_root_label="FWT")
  mtz_dataset.add_miller_array(miller_array=phases,column_types="P", column_root_label="PHWT")
  mtz_dataset.mtz_object().write(file_name=file_name)

# XXX copied these from cctbx.miller; made small change to catch weird case
#  where normalizations are negative.  Just multiply these *-1 and it seems to
#  close to what we want. Figure this out later...
# XXX Also set means=1 not mean square = 1

def amplitude_quasi_normalisations(ma, d_star_power=1, set_to_minimum=None):
    epsilons = ma.epsilons().data().as_double()
    mean_f_sq_over_epsilon = flex.double()
    for i_bin in ma.binner().range_used():
      sel = ma.binner().selection(i_bin)
      #sel_f_sq = flex.pow2(ma.data().select(sel))
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
      # HACK NO REASON THIS SHOULD WORK
      sel = (mean_f_sq_over_epsilon_interp <= set_to_minimum)
      mean_f_sq_over_epsilon_interp.set_selected(sel,-mean_f_sq_over_epsilon_interp)
      sel = (mean_f_sq_over_epsilon_interp <= set_to_minimum)
      mean_f_sq_over_epsilon_interp.set_selected(sel,set_to_minimum)
    assert mean_f_sq_over_epsilon_interp.all_gt(0)
    from cctbx.miller import array
    #return array(ma, flex.sqrt(mean_f_sq_over_epsilon_interp))
    return array(ma, mean_f_sq_over_epsilon_interp)

def quasi_normalize_structure_factors(ma, d_star_power=1, set_to_minimum=None):
    normalisations = amplitude_quasi_normalisations(ma, d_star_power,
       set_to_minimum=set_to_minimum)
    q = ma.data() / normalisations.data()
    from cctbx.miller import array
    return array(ma, q)

def get_array(file_name=None,labels=None):

  print "Reading from %s" %(file_name)
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

  print "Using the array %s" %(",".join(array_to_use.info().labels))
  return array_to_use


def get_amplitudes(args):
  if not args or 'help' in args or '--help' in args:
    print "\nsharpen.py"
    print "Read in map coefficients or amplitudes and sharpen"
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

def scale_amplitudes(pdb_inp=None,map_coeffs=None,
    si=None,resolution=None,overall_b=None,
    fraction_complete=None,
    verbose=False,
    out=sys.stdout):

  # Figure out resolution_dependent sharpening to optimally
  #  match map and model. Then apply it as usual.
  from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi,get_b_iso
  from cctbx.maptbx.segment_and_split_map import get_f_phases_from_model
  from cctbx.maptbx.segment_and_split_map import map_coeffs_to_fp
  import math

  f_array,phases=map_coeffs_as_fp_phi(map_coeffs)

  from cctbx.maptbx.refine_sharpening import get_sharpened_map,\
     quasi_normalize_structure_factors
  (d_max,d_min)=f_array.d_max_min()
  if not f_array.binner():
    f_array.setup_binner(n_bins=si.n_bins,d_max=d_max,d_min=d_min)
  f_array_normalized=quasi_normalize_structure_factors(
        f_array,set_to_minimum=0.01)

  if resolution is None:
    resolution=si.resolution
  if resolution is None:
    raise Sorry("Need resolution for model sharpening")
  # define Wilson B for the model
  if overall_b is None:
    overall_b=10*resolution
    print >>out,"Setting Wilson B = %5.1f A based on resolution of %5.1f A" %(
      overall_b,resolution)

  # Define expected fall-off of CC due to model error
  if si.rmsd is None:
    b_eff=None
  else:
    b_eff=8*3.14159*si.rmsd**2
    print >>out,\
    "Setting b_eff for fall-off at %5.1f A**2 based on model error of %5.1f A" \
       %( b_eff,si.rmsd)

  # create model map using same coeffs
  model_map_coeffs=get_f_phases_from_model(
     pdb_inp=pdb_inp,
     f_array=f_array,
     overall_b=overall_b,
     out=out)

  model_f_array,model_phases=map_coeffs_as_fp_phi(model_map_coeffs)
  model_f_array.setup_binner(n_bins=si.n_bins,d_max=d_max,d_min=d_min)

  # Set overall_b....
  starting_b_iso=get_b_iso(model_f_array,d_min=resolution)
  model_f_array=\
   model_f_array.apply_debye_waller_factors(
      b_iso=overall_b-starting_b_iso)
  final_b_iso=get_b_iso(model_f_array,d_min=resolution)
  print "Effective b_iso of initial and "+\
     "adjusted model map: %6.1f A**2  %6.1f A**2" %(starting_b_iso,final_b_iso)
  model_map_coeffs_normalized=model_f_array.phase_transfer(
     phase_source=model_phases,deg=True)

  # get f and model_f vs resolution and FSC vs resolution and apply
  # scale to f_array and return sharpened map
  dsd = f_array.d_spacings().data()

  ratio_list=flex.double()
  target_sthol2=flex.double()
  cc_list=flex.double()
  d_min_list=flex.double()
  rms_fo_list=flex.double()
  max_possible_cc=None
  for i_bin in f_array.binner().range_used():
    sel       = f_array.binner().selection(i_bin)
    d         = dsd.select(sel)
    d_min     = flex.min(d)
    d_max     = flex.max(d)
    d_avg     = flex.mean(d)
    n         = d.size()
    fo        = map_coeffs.select(sel)
    fc        = model_map_coeffs_normalized.select(sel)
    cc        = fc.map_correlation(other = fo)
    f_array_fc=map_coeffs_to_fp(fc)
    f_array_fo=map_coeffs_to_fp(fo)

    rms_fc=f_array_fc.data().norm()
    rms_fo=f_array_fo.data().norm()
    sthol2=0.25/d_avg**2
    ratio_list.append(rms_fc/rms_fo)
    target_sthol2.append(sthol2)
    cc_list.append(cc)
    d_min_list.append(d_min)
    rms_fo_list.append(rms_fo)
    if b_eff is not None:
      max_cc_estimate=cc* math.exp(min(20.,sthol2*b_eff))
    else:
      max_cc_estimate=cc
    max_cc_estimate=max(0.,min(1.,max_cc_estimate))
    if max_possible_cc is None or (max_cc_estimate > max_possible_cc):
      max_possible_cc=max_cc_estimate
    if verbose:
      print >>out,"d_min: %5.1f  FC: %7.1f  FOBS: %7.1f   CC: %5.2f" %(
      d_avg,rms_fc,rms_fo,cc)

  # Define overall CC based on model completeness (CC=sqrt(fraction_complete))
  if fraction_complete is None:
    fraction_complete=max_possible_cc**2

    print >>out,\
      "Estimated fraction complete is is %5.2f based on low_res CC of %5.2f" %(
        fraction_complete,max_possible_cc)
  else:
    print >>out,"Using fraction complete value of %5.2f "  %(fraction_complete)
    max_possible_cc=fraction_complete**0.5


  target_scale_factors=flex.double()
  max_scale=None
  for i_bin in f_array.binner().range_used():
    index=i_bin-1
    ratio=ratio_list[index]
    cc=cc_list[index]
    sthol2=target_sthol2[index]
    d_min=d_min_list[index]
    corrected_cc=max(0.00001,min(1.,cc/max_possible_cc))
    if b_eff is not None:
      scale_on_fo=ratio * min(1.,
        max(0.00001,corrected_cc) * math.exp(min(20.,sthol2*b_eff)) )
    else:
      scale_on_fo=ratio * min(1.,max(0.00001,corrected_cc))

    if max_scale is None or d_min >= resolution:
      max_scale=scale_on_fo
    else: # at less than resolution don't allow any bigger scale than last one
      scale_on_fo=min(scale_on_fo,max_scale)

    target_scale_factors.append(scale_on_fo)


  target_scale_factors=\
     target_scale_factors/target_scale_factors.min_max_mean().max
  print >>out,"\nScale factors vs resolution:"
  for sthol2,scale,rms_fo in zip(
     target_sthol2,target_scale_factors,rms_fo_list):
     print >>out,"d_min: %5.1f   scale:  %5.2f  rms F:  %7.1f" %(
       0.5/math.sqrt(sthol2),scale,rms_fo*scale)

  # Now create resolution-dependent coefficients from the scale factors
  from cctbx.maptbx.refine_sharpening import run as refine_sharpening
  si.sharpening_method='model_sharpening'
  si.residual_target='model'
  si.b_sharpen=0
  si.b_iso=None

  si.target_scale_factors=target_scale_factors
  si.target_sthol2=target_sthol2

  scale_array=flex.double(f_array.data().size()*(-1.,))
  for i_bin in f_array.binner().range_used():
    sel       = f_array.binner().selection(i_bin)
    scale_array.set_selected(sel,target_scale_factors[i_bin-1])
  assert scale_array.count(-1.)==0

  f_array_b_iso=get_b_iso(f_array,d_min=resolution)

  scaled_f_array=f_array.customized_copy(data=f_array.data()*scale_array)
  scaled_f_array_b_iso=get_b_iso(scaled_f_array,d_min=resolution)
  print >>out,"\nInitial b_iso for "+\
      "map: %5.1f A**2     After adjustment: %5.1f A**2" %(
      f_array_b_iso,scaled_f_array_b_iso)

  assert si.n_real is not None
  new_map_coeffs=scaled_f_array.phase_transfer(phase_source=phases,deg=True)
  return calculate_map(map_coeffs=new_map_coeffs,n_real=si.n_real)

def calculate_map(map_coeffs=None,crystal_symmetry=None,n_real=None):

  if crystal_symmetry is None: crystal_symmetry=map_coeffs.crystal_symmetry()
  from cctbx.maptbx.segment_and_split_map import get_map_from_map_coeffs
  return get_map_from_map_coeffs(
     map_coeffs=map_coeffs,crystal_symmetry=crystal_symmetry, n_real=n_real)

def get_sharpened_map(ma=None,phases=None,b=None,resolution=None,
    n_real=None,d_min_ratio=None):
  assert n_real is not None
  sharpened_ma=adjust_amplitudes_linear(ma,b[0],b[1],b[2],resolution=resolution,
     d_min_ratio=d_min_ratio)
  new_map_coeffs=sharpened_ma.phase_transfer(phase_source=phases,deg=True)
  return calculate_map(map_coeffs=new_map_coeffs,n_real=n_real)

def calculate_match(target_sthol2=None,target_scale_factors=None,b=None,resolution=None,d_min_ratio=None,rmsd=None,fraction_complete=None):

  if fraction_complete is None:
    pass # XXX not implemented for fraction_complete

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
    print >>out,"Result: b1 %7.2f b2 %7.2f b3 %7.2f resolution %7.2f %s: %7.3f" %(
     b[0],b[1],b[2],self.resolution,self.residual_target,value)

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
        rmsd=self.rmsd,
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
    for i in xrange(len(list(b))):
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
  else:
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
    print >>out,"Setting rmsd to %5.1f A based on resolution of %5.1f A" %(
       rmsd,resolution)

  if fraction_occupied is None: fraction_occupied=0.20
  if region_weight is None: region_weight=20.
  if sa_percent is None: sa_percent=30.
  if n_bins is None: n_bins=20
  if max_regions_to_test is None: max_regions_to_test=30


  # Get initial value

  best_b=b
  print >>out,"Getting starting value ...",residual_target
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
  print >>out,"Starting value: %7.2f" %(starting_result)

  if ma:
    print >>out,"Normalizing structure factors..." # XXX why here?
    (d_max,d_min)=ma.d_max_min()
    ma.setup_binner(n_bins=n_bins,d_max=d_max,d_min=d_min)
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
  print >>out,"Starting value after normalization: %7.2f" %(
     starting_normalized_result)
  best_sharpened_ma=ma
  best_result=starting_normalized_result
  best_b=refined.get_b()

  refined.run()

  final_result=refined.show_result(out=out)
  print >>out,"Final value: %7.2f" %(
     final_result)

  if final_result>best_result:
    best_sharpened_ma=refined.sharpened_ma
    best_result=final_result
    best_b=refined.get_b()
  print >>out,"Best overall result: %7.2f: " %(best_result)

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
