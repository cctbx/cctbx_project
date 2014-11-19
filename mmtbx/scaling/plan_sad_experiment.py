
from __future__ import division
import mmtbx.scaling
from iotbx.bioinformatics import any_sequence_format
from libtbx.utils import Sorry, null_out
from libtbx import table_utils
import math
import sys,os

def get_vol_per_residue(chain_type='PROTEIN'):
  if chain_type=='PROTEIN':
    vol_per_residue=135.3
  else:
    vol_per_residue=495.
  return vol_per_residue

def get_atoms_per_residue(chain_type='PROTEIN'):
  if chain_type=='PROTEIN':
    atoms_per_residue=7
  else:
    atoms_per_residue=26
  return atoms_per_residue

def get_nrefl(residues=200,dmin=2.5,solvent_fraction=0.5,chain_type='PROTEIN'):
  nrefl=2.*3.14159*0.3333*residues*get_vol_per_residue(chain_type=chain_type)/ \
        (1.-solvent_fraction)*(1./dmin**3)
  return nrefl

def get_sigf(nrefl,nsites,natoms,z,fpp,target_s_ano=15.,ntries=1000,
     min_cc_ano=None,
     fa2=None,fb2=None,disorder_parameter=None,
     fo_list=None,fo_number_list=None,occupancy=None,include_zero=True,
     ratio_for_failure=0.95,
     resolution=None,
     cc_ano_estimators=None,
     signal_estimators=None,
     max_i_over_sigma=None):
  closest_sigf=None
  closest_dist2=None
  start_value=1
  if include_zero: start_value=0
  for i in xrange(start_value,ntries):
    sigf=i*1./float(ntries)
    if max_i_over_sigma and get_i_over_sigma_from_sigf(sigf)>max_i_over_sigma:
      continue
    s_ano,s_ano_sig,cc_ano,cc_ano_sig,cc_half,cc_half_sig,fpp_weak,\
        cc_ano_weak,cc_half_weak,i_over_sigma=\
      get_values_from_sigf(nrefl,nsites,natoms,z,fpp,sigf,
          fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
          fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy,
          get_fpp_weak=False,resolution=resolution,
          cc_ano_estimators=cc_ano_estimators,
          signal_estimators=signal_estimators)

    if min_cc_ano is not None and cc_ano_weak < min_cc_ano:
      continue

    dist2=(s_ano-target_s_ano)**2
    if closest_dist2 is None or dist2<closest_dist2:
     closest_dist2=dist2
     closest_sigf=sigf
  if closest_sigf==0:  # try with target of 90% of maximum available
    s_ano,s_ano_sig,cc_ano,cc_ano_sig,cc_half,cc_half_sig,\
       fpp_weak,cc_ano_weak,\
       cc_half_weak,i_over_sigma=\
       get_values_from_sigf(nrefl,nsites,natoms,z,fpp,closest_sigf,
          fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
          fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy,
          get_fpp_weak=False,resolution=resolution,
          cc_ano_estimators=cc_ano_estimators,
          signal_estimators=signal_estimators)
    closest_sigf=get_sigf(
     nrefl,nsites,natoms,z,fpp,target_s_ano=s_ano*ratio_for_failure,
     ntries=ntries,
     min_cc_ano=min_cc_ano,
     fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
     fo_list=fo_list,fo_number_list=fo_number_list,
     occupancy=occupancy,include_zero=False,resolution=resolution,
     cc_ano_estimators=cc_ano_estimators,
     signal_estimators=signal_estimators,
     max_i_over_sigma=max_i_over_sigma)
  return closest_sigf

def get_sano(nrefl,nsites,natoms,z,fpp,sigf,
       fa2=None,fb2=None,disorder_parameter=None,
       fo_list=None,fo_number_list=None,occupancy=None):

  # nrefl = number of anomalous data
  # nsites = number of anomalously scattering atoms in asymmetric unit
  # natoms = number of non-anomalously scattering atoms in asymmetric unit
  # z = typical Z of non-anomalously scattering atoms
  # sigf = rms(sigF)/rms(F) overall
  # delta_b_ano=B for anom atoms minus B for rest of structure.
  #  Ignoring this effect here
  #sano2=(4./5.)*(nrefl/nsites)/(1.+(natoms/nsites)*(z*sigf/fpp)**2)

  # more detailed calculation with fa2 and fb2 and disorder_parameter:
  # sano=cc_ano*sqrt((4/5)*nrefl/nsites)
  cc_ano=get_cc_ano(nrefl,nsites,natoms,z,fpp,sigf,
      fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
      fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)
  sano_new=cc_ano*math.sqrt((4./5.)*nrefl/nsites)

  # note this is very close to:
  #  sano2=(4./5.)*nrefl*fpp**2/(natoms*z**2*sigf**2)
  #  sano ~ sqrt(4./5.)*(fpp/z*sigf)*sqrt(nrefl/natoms)
  #     which is proportional to fpp and 1/sigf
  #sano=math.sqrt(sano2)

  return sano_new

def get_cc_ano(nrefl,nsites,natoms,z,fpp,sigf,
      fa2=None,fb2=None,disorder_parameter=None,
      fo_list=None,fo_number_list=None,occupancy=None):
  "cc_ano**2=1/[1 + <sigF**2>/(n/N)(Za/Z)**2<F**2>]"
  #ccano2=1./(1.+sigf**2*natoms*z**2/(nsites*fpp**2))

  # recalculate fa2 and fb2 based on fpp value:
  target_fpp_list=[fpp]
  target_fpp_number_list=[nsites]

  local_fa2=get_normalized_scattering(
      fpp_list=target_fpp_list,
      fpp_number_list=target_fpp_number_list,
      fo_list=fo_list,
      fo_number_list=fo_number_list,
      occupancy=occupancy)

  #  more detailed: cc_ano**2=1/(1+Q+(fb2/fa2)+(sigf2/fa2))
  ccano2_new=1./(1.+disorder_parameter+(fb2/local_fa2)+(sigf**2/local_fa2))

  return math.sqrt(ccano2_new)

def get_cc_half(cc_ano,disorder_parameter=None):
  "cc_half=cc_ano_perf**2/(2-cc_ano_perf**2)"

  #cc_half=cc_ano*cc_ano/(2.-cc_ano*cc_ano)

  # with disorder_parameter, cc**_ano (useful cc_ano)=cc_ano/sqrt(1+Q)
  #   or cc_ano=cc**_ano*sqrt(1+Q), where cc**_ano is cc_ano_perf, useful cc_ano
  ccano2=min(1.,cc_ano*cc_ano*(1+disorder_parameter))
  cc_half_new=ccano2/(2.-ccano2)
  return cc_half_new

def estimate_fpp_weak(nrefl,nsites,natoms,z,fpp,sigf,
      fa2=None,fb2=None,disorder_parameter=None,
      fo_list=None,fo_number_list=None,occupancy=None):
  # estimate fpp that would give half the sano value
  sano=get_sano(nrefl,nsites,natoms,z,fpp,sigf,
      fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
      fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)

  target=sano/2.
  best_scale=1.0
  for i in xrange(1,101):
    scale=float(i)*0.01
    ss=get_sano(nrefl,nsites,natoms,z,fpp*scale,sigf,
      fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
      fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)
    if ss >=target and scale<best_scale:
      best_scale=scale
      break


  return fpp*best_scale

def get_i_over_sigma_from_sigf(sigf):
  # <I/sigI> = 0.88 * rms(F)/rms(sigF) r**2=0.82
  if sigf>0:
    return 0.88/sigf
  else:
    return 999.

def get_sigf_from_i_over_sigma(i_over_sigma):
  # just inverse of get_i_over_sigma_from_sigf
  #  ios=.88/sigf  -> sigf=0.88/ios
  if i_over_sigma > 0:
     sigf=0.88/i_over_sigma
  else:
     sigf=0.001
  return sigf

def get_i_over_sigma(i_obs):
  data=i_obs.data().select(i_obs.sigmas() > 0)
  sigmas=i_obs.sigmas().select(i_obs.sigmas() > 0)
  iover=data/sigmas
  return iover.min_max_mean().mean

def get_values_from_sigf(nrefl,nsites,natoms,z,fpp,sigf,
       fa2=None,fb2=None,disorder_parameter=None,
       fo_list=None,fo_number_list=None,occupancy=None,
       get_fpp_weak=True,resolution=None,
       cc_ano_estimators=None,signal_estimators=None):
    # 2014-11-04 add capability for Bayesian estimation (update) of
    #   signal and cc_ano
    s_ano=get_sano(nrefl,nsites,natoms,z,fpp,sigf,
      fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
      fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)
    cc_ano=get_cc_ano(nrefl,nsites,natoms,z,fpp,sigf,
      fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
      fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)

    if cc_ano_estimators:
      cc_ano,cc_ano_sig=cc_ano_estimators.apply_estimators(
         value_list=[cc_ano],data_items=['pred_cc_perfect'],
         resolution=resolution)
    else:
      cc_ano_sig=None
    if signal_estimators:
      s_ano,s_ano_sig=signal_estimators.apply_estimators(
         value_list=[s_ano],data_items=['pred_signal'],
         resolution=resolution)
    else:
      s_ano_sig=None

    cc_half=get_cc_half(cc_ano,disorder_parameter=disorder_parameter)

    if cc_ano_estimators:
      cc_half_high=get_cc_half(cc_ano+cc_ano_sig,
         disorder_parameter=disorder_parameter)
      cc_half_low=get_cc_half(max(0,cc_ano-cc_ano_sig),
         disorder_parameter=disorder_parameter)
      cc_half_sig=0.5*(cc_half_high-cc_half_low)
    else:
      cc_half_sig=None

    if get_fpp_weak:
      fpp_weak=estimate_fpp_weak(nrefl,nsites,natoms,z,fpp,sigf,
        fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
        fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)
      cc_ano_weak=get_cc_ano(nrefl,nsites,natoms,z,fpp_weak,sigf,
        fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
        fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)

      if cc_ano_estimators:
        cc_ano_weak,cc_ano_weak_sig=cc_ano_estimators.apply_estimators(
           value_list=[cc_ano_weak],data_items=['pred_cc_perfect'],
           resolution=resolution)

      cc_half_weak=get_cc_half(cc_ano_weak,
          disorder_parameter=disorder_parameter)

      i_over_sigma=get_i_over_sigma_from_sigf(sigf)
      return s_ano,s_ano_sig,cc_ano,cc_ano_sig,cc_half,cc_half_sig,fpp_weak,\
          cc_ano_weak,cc_half_weak,i_over_sigma
    else:
       return s_ano,s_ano_sig,cc_ano,cc_ano_sig,cc_half,cc_half_sig,\
         fpp,cc_ano,cc_half,0.0

def get_fp_fdp(atom_type=None,wavelength=None,out=sys.stdout):
  if not atom_type or not wavelength:
    raise Sorry("Please specify either f_double_prime or " +
      "atom_type and wavelength")
  from cctbx.eltbx import sasaki
  try:
    table = sasaki.table(atom_type)
    fp_fdp = table.at_angstrom(wavelength)
  except ValueError :
    raise Sorry("Unable to get scattering factors for %s" %(atom_type))
  return fp_fdp

def get_fo(atom_type=None,wavelength=None,out=sys.stdout):
  if not atom_type or not wavelength:
    raise Sorry("Please specify either f_double_prime or " +
      "atom_type and wavelength")
  from cctbx.eltbx import sasaki
  try:
    table = sasaki.table(atom_type)
    fo=table.atomic_number()
  except ValueError :
    raise Sorry("Unable to get scattering factors for %s" %(atom_type))
  return fo

def get_aa_and_met(sequence):
  sequence=sequence.upper().replace(" ","")
  n_aa=len(sequence)
  n_met=sequence.count("M")
  n_cys=sequence.count("C")
  return n_aa,n_met,n_cys

def get_number_of_sites(atom_type=None,n_met=0,n_cys=0,
      n_aa=0,ncs_copies=1,out=sys.stdout):
    # guess number of sites:
    number_of_sites_lowres=None
    if atom_type is None:
      print >>out, "No heavy atom type set, so no sites estimated"
      number_of_sites=None
    elif atom_type.lower() in ['se']:
      number_of_sites=max(1,ncs_copies*n_met)
    elif atom_type.lower() in ['s']:
      number_of_sites=max(1,ncs_copies*(n_met+n_cys))
      number_of_sites_lowres=ncs_copies*(n_met+int(float(n_cys)/2.))
    elif atom_type.lower() in ['br','i']: # 1 per 20 up to 10 per chain
      number_of_sites=ncs_copies*max(1,min(10,1+int(float(n_aa)/20.)))
    else: # general ha  1 per 100 up to 2 per chain
      number_of_sites=ncs_copies*max(1,min(2,1+int(float(n_aa)/100.)))
    if number_of_sites:
      print >>out,"\nBest guess of number of %s sites: %d" %(
        atom_type.upper(),number_of_sites)
    if number_of_sites_lowres is None: number_of_sites_lowres=number_of_sites
    return number_of_sites,number_of_sites_lowres

def get_residues_and_ha(seq_file,atom_type=None,
      chain_type=None,out=sys.stdout):
  if not seq_file or not os.path.isfile(seq_file):
    raise Sorry("Please supply number of residues or a sequence file")
  objects, non_compliant = any_sequence_format(seq_file)
  if non_compliant:
    raise Sorry("Sorry, unable to read the sequence file %s" %(seq_file))
  n_aa, n_met, n_cys = 0, 0, 0
  for seq_obj in objects :
    n_aa_,n_met_,n_cys_ = get_aa_and_met(sequence=seq_obj.sequence)
    n_aa += n_aa_
    n_met += n_met_
    n_cys += n_cys_
  number_of_s=n_met+n_cys
  number_of_sites,number_of_sites_lowres=get_number_of_sites(
      atom_type=atom_type,n_met=n_met,n_cys=n_cys,
      n_aa=n_aa,ncs_copies=1,out=null_out())
  return n_aa,number_of_sites,number_of_s

def get_disorder_parameter(ideal_cc_anom=None):
  if not ideal_cc_anom:
    return 0.  #
  else:
    # E=1+Q; ideal_cc_anom=1/sqrt(E) -> Q=1/ideal_cc_anom**2 - 1
    cc2=ideal_cc_anom**2
    return (1.-cc2)/cc2

def include_intrinsic_scatterers(
   intrinsic_scatterers_list,target_fpp,wavelength,
    minimum_ratio=0.25):
  all_weak=True # noise if all are weak
  for x in intrinsic_scatterers_list:
    fpp=get_fp_fdp(atom_type=x,wavelength=wavelength).fdp()
    if fpp >= minimum_ratio*target_fpp:
      all_weak=False
  if all_weak:
    return True
  else:
    return False


def get_fo_list(chain_type="PROTEIN",wavelength=1.0,
    residues=1,
    target_fpp=None,target_atom_type=None,target_n=None,
    intrinsic_scatterers_as_noise=None,
    include_weak_anomalous_scattering=None,
    number_of_s=0.,out=sys.stdout):

  if chain_type=="PROTEIN":
     intrinsic_scatterers_list=['S']
     if number_of_s is None:
       number_of_s=0.044*residues # typical s content
     intrinsic_scatterers_number_list=[number_of_s/residues]
     atoms_list=['C','N','O']
     atoms_number_list=[5.15,1.37,1.58] # average (for a2u-globulin)
  else:
     intrinsic_scatterers_list=['P']
     intrinsic_scatterers_number_list=[1]
     atoms_list=['C','N','O']
     atoms_number_list=[8.75,3.75,6.]  # average of DNA
     if chain_type=="RNA":
       atoms_number_list[2]+=1 # an extra O in RNA

  # decide if we include the intrinsic scatterers
  if intrinsic_scatterers_as_noise or (
      intrinsic_scatterers_as_noise is None and include_intrinsic_scatterers(
          intrinsic_scatterers_list,target_fpp,wavelength) ):
       atoms_list+=intrinsic_scatterers_list
       atoms_number_list+=intrinsic_scatterers_number_list

  fo_list=[]
  fo_number_list=[]
  fpp_list=[]
  fpp_number_list=[]
  noise_table_rows = []
  n_atoms=0
  for x,n in zip(atoms_list,atoms_number_list):
    fpp=get_fp_fdp(atom_type=x,wavelength=wavelength).fdp()
    fo=get_fo(atom_type=x,wavelength=wavelength)
    fo_list.append(fo)
    fo_number_list.append(n*residues)
    n_atoms+=n*residues
    if include_weak_anomalous_scattering:
      fpp_list.append(fpp)
      fpp_number_list.append(n*residues)
      contribution=fpp*math.sqrt(n*residues)
      noise_table_rows.append(
         (x,fpp,n*residues,contribution)
      )
  return fo_list,fo_number_list,fpp_list,fpp_number_list,\
      n_atoms,noise_table_rows

def get_normalized_scattering(
      fpp_list=[],
      fpp_number_list=[],
      fo_list=[],
      fo_number_list=[],
      occupancy=1.):
  # estimate ratio of anomalous to real scattering for these atoms
  # sum(Z_b_i**2 N_b_i)/sum(Z_i**2 N_i)
  sum_real=0.
  for fo,fo_n in zip(fo_list,fo_number_list):
    sum_real+=fo**2*fo_n
  sum_anom=0.
  for fpp,fpp_n in zip(fpp_list,fpp_number_list):
    sum_anom+=(occupancy*fpp)**2*fpp_n
  if sum_real <=0:
    sum_real=1.

  return sum_anom/sum_real

def get_local_file_name(estimator_type):
  if estimator_type=='cc_star':
    local_file_name='cc_ano_data.dat'
  elif estimator_type=='signal':
    local_file_name='signal_from_est_signal.dat'
  elif estimator_type=='cc_ano':
    local_file_name='cc_anom_from_cc_anom_star.dat'
  elif estimator_type=='solved':
    local_file_name='percent_solved_vs_signal.dat'
  else:
    raise Sorry("No estimator type %s" %(estimator_type))
  return local_file_name

class interpolator:
  # basic interpolator. Call with list of pairs of (predictor,target)
  # list can have extra data (ignored)
  # then apply with (predictor) and get an estimate of target
  # If off the end in either direction, use last available value (do not
  #   extrapolate)
  def __init__(self,target_predictor_list,extrapolate=False,
      require_monotonic_increase=False):
    assert not extrapolate # (not programmed)
    # make a dict
    self.target_dict={}
    self.keys=[]
    for values in target_predictor_list:
       predictor=values[0]
       target=values[1]
       self.target_dict[predictor]=target
       self.keys.append(predictor)
    self.keys.sort()
    if require_monotonic_increase:
      self.set_monotonic_increase()
    from copy import deepcopy
    self.reverse_keys=deepcopy(self.keys)
    self.reverse_keys.reverse()  # so we can go down too

  def set_monotonic_increase(self):
    # require never to decrease value with increasing key
    # start at middle and make sure that value never goes above
    #  mean of remaining or below previous
    # verify that overall means for each end do not change
    n=len(self.keys)
    if n<2: return
    i_middle=n//2
    from copy import deepcopy
    value_list=[]
    for key in self.keys:
      value_list.append(self.target_dict[key])
    new_values=deepcopy(value_list)
    for i in xrange(i_middle+1,n):
      prev_value=new_values[i-1]
      remainder=value_list[i:]
      value=value_list[i]
      mean_remainder=self.get_mean(remainder)
      new_values[i]=max(prev_value,min(mean_remainder,value))
    for i in xrange(i_middle-1,-1,-1):
      prev_value=new_values[i+1]
      remainder=value_list[:i+1]
      value=value_list[i]
      mean_remainder=self.get_mean(remainder)
      new_values[i]=min(prev_value,max(mean_remainder,value))
    from cctbx.array_family import flex
    a=flex.double()
    b=flex.double()
    delta_list=[]
    for key,value in zip(self.keys,new_values):
      delta_list.append(value-self.target_dict[key])
    # adjust delta mean to zero in each region
    offset_low=self.get_mean(delta_list[:i_middle])
    offset_high=self.get_mean(delta_list[i_middle+1:])
    for i in xrange(i_middle+1,n):
      new_values[i]=new_values[i]-offset_high
    for i in xrange(i_middle-1,-1,-1):
      new_values[i]=new_values[i]-offset_low
    for key,value in zip(self.keys,new_values):
      self.target_dict[key]=value

  def get_mean(self,remainder):
    n=len(remainder)
    if n<1: return 0
    a=0.
    for x in remainder: a+=x
    return a/n



  def interpolate(self,predictor):
    lower_pred=None
    dist_from_lower=None
    higher_pred=None
    dist_from_higher=None
    if predictor <= self.keys[0]:
      return self.target_dict[self.keys[0]]
    elif predictor >= self.keys[-1]:
      return self.target_dict[self.keys[-1]]
    else:  # between two keys. Find highest key < predictor and lowest > pred
      for key in self.reverse_keys: # high to low values of keys
        if predictor > key:
          lower_pred=self.target_dict[key]
          dist_from_lower=predictor-key
          break
      for key in self.keys: # low to high values of keys
        if predictor < key:
          higher_pred=self.target_dict[key]
          dist_from_higher=key-predictor
          break
      if lower_pred is None or higher_pred is None:
        raise Sorry("Unable to interpolate...")
      dist=(dist_from_lower+dist_from_higher)
      value=lower_pred+(higher_pred-lower_pred)*dist_from_lower/dist
      return value


def get_interpolator(estimator_type='solved',predictor_variable='PredSignal',
       require_monotonic_increase=None,out=sys.stdout):
  local_file_name=get_local_file_name(estimator_type)
  print >>out,"\nSetting up interpolator for %s" %(estimator_type)
  import libtbx.load_env
  file_name=libtbx.env.find_in_repositories(
    relative_path=os.path.join("mmtbx","scaling",local_file_name),
      test=os.path.isfile)
  # data looks like:
  """
Signal  %Solved   N  PredSignal  %Solved  N   BayesEstSignal  %Solved  N
   1.0       2    44        1.0       1    69        1.0       0     6
   3.0       0   123        3.0       0    95        3.0       9   104
  """
  # Pick up 2 columns starting with predictor variable

  from mmtbx.scaling.bayesian_estimator import get_table_as_list
  prediction_values_as_list,info_values_as_list,\
     dummy_target_variable,dummy_data_items,dummy_info_items=\
    get_table_as_list(file_name=file_name, select_only_complete=False,
    start_column_header=predictor_variable,out=out)
  return interpolator(prediction_values_as_list,
       require_monotonic_increase=require_monotonic_increase)


def get_estimators(estimator_type='signal',
    resolution_cutoffs=None,out=sys.stdout):
  # get estimators for signal from est_signal or cc_ano from cc*_ano
  local_file_name=get_local_file_name(estimator_type)

  import libtbx.load_env

  print >>out,"\nSetting up estimator for %s" %(estimator_type)
  file_name=libtbx.env.find_in_repositories(
    relative_path=os.path.join("mmtbx","scaling",local_file_name),
      test=os.path.isfile)

  from mmtbx.scaling.bayesian_estimator import estimator_group
  estimators=estimator_group(
    resolution_cutoffs=resolution_cutoffs,out=out)
  estimators.set_up_estimators(
    file_name=file_name,select_only_complete=False)
  estimators.show_summary()
  return estimators

class estimate_necessary_i_sigi (mmtbx.scaling.xtriage_analysis) :
  def __init__ (self,
      chain_type='PROTEIN',
      residues=250,
      number_of_s=0,
      solvent_fraction=0.50,
      nsites=5,
      wavelength=1.0,
      atom_type=None,
      fpp=3.8,
      target_s_ano=30,
      min_cc_ano=0.15,
      dmin=None,
      occupancy=1.,
      ideal_cc_anom=0.76,
      bayesian_estimates=True,
      include_weak_anomalous_scattering=True,
      intrinsic_scatterers_as_noise=None,
      ratio_for_failure=0.95,
      i_over_sigma=None,
      max_i_over_sigma=None,
      quiet=False) :
    self.chain_type = chain_type
    self.residues = residues
    self.solvent_fraction = solvent_fraction
    self.nsites = nsites
    self.fpp = fpp
    self.target_s_ano = target_s_ano
    self.min_cc_ano = min_cc_ano
    self.wavelength = wavelength
    self.max_i_over_sigma = max_i_over_sigma
    self.bayesian_estimates = bayesian_estimates
    if atom_type is None:
      self.atom_type='-'
    else:
      self.atom_type=atom_type
    self.occupancy=occupancy
    self.ratio_for_failure=ratio_for_failure

    vol_per_residue = get_vol_per_residue(chain_type=chain_type)

    # q (disorder_parameter) = normalized mean square anom diff not due
    #   to target atoms, E=1+Q
    self.disorder_parameter=get_disorder_parameter(ideal_cc_anom=ideal_cc_anom)

    # atoms in structure and their numbers, normal and anomalous scattering
    fo_list,fo_number_list,fpp_list,fpp_number_list,\
       self.natoms,self.noise_table_rows=get_fo_list(
        residues=residues,
        chain_type=chain_type,
        wavelength=self.wavelength,
        target_fpp=self.fpp,
        target_atom_type=atom_type,
        target_n=nsites,
        intrinsic_scatterers_as_noise=intrinsic_scatterers_as_noise,
        include_weak_anomalous_scattering=include_weak_anomalous_scattering,
        number_of_s=number_of_s)


    # scattering from target anomalous atoms:
    # fa2=Occ**2*sum(Z_a_i**2 N_a_i)/sum(Z_i**2 N_i)
    target_fpp_list=[self.fpp]
    target_fpp_number_list=[self.nsites]

    self.fa2=get_normalized_scattering(
      fpp_list=target_fpp_list,
      fpp_number_list=target_fpp_number_list,
      fo_list=fo_list,
      fo_number_list=fo_number_list,
      occupancy=occupancy)

    # scattering from all other anomalous atoms
    # fb2=sum(Z_b_i**2 N_b_i)/sum(Z_i**2 N_i)

    self.fb2=get_normalized_scattering(
      fpp_list=fpp_list,
      fpp_number_list=fpp_number_list,
      fo_list=fo_list,
      fo_number_list=fo_number_list)


    z=7 # try 5.6 here ZZZ
    if dmin:
      self.dmin_ranges=[dmin]
    else:
      self.dmin_ranges=[6,5,3,2.5,2,1.5]
    self.table_rows = []
    self.representative_values = None
    self.skipped_resolutions = []
    self.missed_target_resolutions = []
    self.input_i_over_sigma=i_over_sigma
    self.used_max_i_over_sigma=True

    if self.bayesian_estimates:

       # set up estimators of cc_ano from cc_*_ano and signal from est_signal

       # estimate cc_ano from cc_star_ano
       self.cc_ano_estimators=get_estimators(estimator_type='cc_ano',
         resolution_cutoffs=self.dmin_ranges,out=null_out())

       # estimate true signal from estimated signal
       self.signal_estimators=get_estimators(estimator_type='signal',
         resolution_cutoffs=self.dmin_ranges,out=null_out())

       # estimate fom of phasing from cc_ano
       self.fom_estimators=get_estimators(estimator_type='cc_star',
         resolution_cutoffs=self.dmin_ranges,out=null_out())

       # solved (probability of hyss finding >=50% of sites) is just a table
       #  so interpolate the probability.
       self.solved_interpolator=get_interpolator(estimator_type='solved',
          require_monotonic_increase=True,
          predictor_variable='PredSignal',out=null_out())
    else:
      self.cc_ano_estimators=None
      self.signal_estimators=None
      self.fom_estimators=None

    for dmin in self.dmin_ranges:
      # Guess reflections from residues, dmin, solvent fraction

      nrefl=get_nrefl(
         residues=residues,dmin=dmin,solvent_fraction=solvent_fraction)

      # identify rms(sigF)/rms(F) necessary to get target_s_ano with this
      # many reflections, sites, atoms, f" value

      if not i_over_sigma:
        sigf=get_sigf(nrefl,nsites,self.natoms,z,fpp,target_s_ano=target_s_ano,
          min_cc_ano=min_cc_ano,
          fa2=self.fa2,fb2=self.fb2,disorder_parameter=self.disorder_parameter,
          fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy,
          ratio_for_failure=self.ratio_for_failure,resolution=dmin,
          cc_ano_estimators=None,
          signal_estimators=None,
          max_i_over_sigma=self.max_i_over_sigma)
        if sigf is None: # failed so far...
          self.used_max_i_over_sigma=True
          sigf=get_sigf_from_i_over_sigma(self.max_i_over_sigma)
      else:  # input i_over_sigma...estimate sigf
        sigf=get_sigf_from_i_over_sigma(i_over_sigma)
      # what are expected signal, useful cc_ano, cc_half-dataset, <I>/<sigI>


      s_ano,s_ano_sig,cc_ano,cc_ano_sig,cc_half,cc_half_sig,\
          fpp_weak,cc_ano_weak,\
          cc_half_weak,local_i_over_sigma=\
        get_values_from_sigf(nrefl,nsites,self.natoms,z,fpp,sigf,
          fa2=self.fa2,fb2=self.fb2,disorder_parameter=self.disorder_parameter,
          fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy,
          get_fpp_weak=True,
          resolution=dmin,
          cc_ano_estimators=self.cc_ano_estimators,
          signal_estimators=self.signal_estimators)


      if local_i_over_sigma>=999:
        self.skipped_resolutions.append(dmin)
        continue  # hopeless
      if s_ano<0.95*target_s_ano:  # must not be able to get target s_ano
        self.missed_target_resolutions.append(dmin)

      if self.bayesian_estimates:  # use range from s_ano+/-s_ano_sig etc
        solved=self.solved_interpolator.interpolate(s_ano)
        if solved is None: solved=0.0
        if s_ano_sig is None: s_ano_sig=0.0
        s_ano_weak=max(0,s_ano-s_ano_sig)
        cc_half_weak=max(0,cc_half-cc_half_sig)
        cc_ano_weak=max(0.,cc_ano-cc_ano_sig)
        cc_st_square=0.04 #cc_ano**2
        fom,s_fom=self.fom_estimators.apply_estimators(
         value_list=[cc_st_square,None,None],
         data_items=['cc_st_square','skew','e'],
         resolution=dmin)
      else:
        s_ano_weak=s_ano/2
        solved=0.0
        fom=0.

      self.table_rows.append([
        "%5.2f" % dmin,
        "%7d" % nrefl,
        "%6.0f" % local_i_over_sigma,
        "%7.1f" % (100.*sigf),
        "%5.2f" % (cc_half),
        "%5.2f" % ( cc_ano),
        "%3.0f" % ( s_ano),
        "%3.0f" % (solved),
        "%3.2f" % (fom),
      ])
      if ((self.representative_values is None) or
          len(self.dmin_ranges) < 2 or dmin == self.dmin_ranges[-2]) :
        self.representative_values = [dmin,nsites,nrefl,fpp,local_i_over_sigma,
           sigf,cc_half_weak,cc_half,cc_ano_weak,cc_ano,s_ano,solved,fom]


  def representative_dmin(self):
    return self.representative_values[0]

  def representative_nsites(self):
    return self.representative_values[1]

  def representative_nrefl(self):
    return self.representative_values[2]

  def representative_fpp(self):
    return self.representative_values[3]

  def representative_i_over_sigma(self):
    return self.representative_values[4]

  def representative_sigf(self):
    return self.representative_values[5]

  def representative_cc_half_weak(self):
    return self.representative_values[6]

  def representative_cc_half(self):
    return self.representative_values[7]

  def representative_cc_ano_weak(self):
    return self.representative_values[8]

  def representative_cc_ano(self):
    return self.representative_values[9]

  def representative_s_ano(self):
    return self.representative_values[10]

  def representative_solved(self):
    return self.representative_values[11]

  def representative_fom(self):
    return self.representative_values[12]

  def show_summary(self):
    if self.is_solvable():
      print """
I/sigI: %7.1f
Dmin:      %5.2f
cc_half:   %5.2f - %5.2f
cc*_anom:  %5.2f - %5.2f
Signal:   %5.1f - %5.1f
p(Substr):%3d %%
FOM:      %3.2f
""" %(
 self.representative_i_over_sigma(),
 self.representative_dmin(),
 self.representative_cc_half_weak(),
 self.representative_cc_half(),
 self.representative_cc_ano_weak(),
 self.representative_cc_ano(),
 self.representative_s_ano()/2,
 self.representative_s_ano(),
 int(0.5+self.representative_solved()),
 self.representative_fom(),
 )

  def is_solvable (self) :
    return (self.representative_values is not None)

  def _show_impl (self, out) :
    out.show_header("SAD experiment planning")
    out.show_sub_header(
      "Dataset overall I/sigma required to solve a structure")

    out.show_paragraph_header(
      "\nDataset characteristics:")

    out.show_preformatted_text("""\
  Target anomalous signal: %(target_s_ano)7.1f
  Residues: %(residues)d
  Chain-type: %(chain_type)s
  Solvent_fraction: %(solvent_fraction)7.2f
  Atoms: %(natoms)d
  Anomalously-scattering atom: %(atom_type)s
  Wavelength: %(wavelength)7.4f A
  Sites: %(nsites)d
  f-double-prime: %(fpp)7.2f
""" % self.__dict__)

    if self.atom_type:
       t=self.atom_type
    else:
       t='-'
    contribution=self.fpp*math.sqrt(self.nsites)
    out.show_preformatted_text("""\
Target anomalous scatterer:
  Atom: %2s  f": %4.2f  n:%5.0f   rmsF:%7.1f""" %(
         t,self.fpp,self.nsites,contribution))
    if self.noise_table_rows:
      out.show_preformatted_text("""\

Other anomalous scatterers in the structure:""")
      for row in self.noise_table_rows :
        out.show_preformatted_text(
    '  Atom: %2s  f": %4.2f  n:%5.0f   rmsF:%7.1f' %tuple(row))

      fa=100.*math.sqrt(self.fa2)
      fb=100.*math.sqrt(self.fb2)
      fab=math.sqrt(self.fa2/(self.fa2+self.fb2))
      out.show_preformatted_text("""\

Normalized anomalous scattering:
  From target anomalous atoms rms(x**2)/rms(F**2):  %7.2f
  From other anomalous atoms rms(e**2)/rms(F**2):   %7.2f
  Correlation of useful to total anomalous scattering: %4.2f
""" % (fa,fb,fab))

    out.show_sub_header(
      "Dataset <I>/<sigI> needed for anomalous signal of 15-30")
    out.show_preformatted_text("""
-------Targets for entire dataset-------  ----------Likely outcome-----------""")

    if (len(self.table_rows) == 0) :
      out.show_text("SAD solution unlikely with the given parameters.")
      return
    if (not out.gui_output) :
      out.show_preformatted_text("""
                              Anomalous    Useful    Useful
                            Half-dataset  Anom CC   Anomalous
 Dmin   N     I/sigI sigF/F     CC       (cc*_anom)  Signal   P(Substr)   FOM
                      (%)                                        (%)
""")
      for row in self.table_rows :
        out.show_preformatted_text(
        "%s%s%s%s     %s       %s      %s        %s       %s" %
          tuple(row))
    else :
      table = table_utils.simple_table(
        table_rows=self.table_rows,
        column_headers=["d_min", "N", "I/sigI", "sigF/F (%)",
          "Half-dataset CC_ano", "CC*_ano", "Anom. signal","P(Substr)","FOM"])
      out.show_table(table)
    (dmin,nsites,nrefl,fpp,i_over_sigma,sigf,cc_half_weak,cc_half,cc_ano_weak,
      cc_ano,s_ano,solved,fom) = tuple(self.representative_values)

    if self.missed_target_resolutions:
      self.missed_target_resolutions.sort()
      extra_note=""
      if self.used_max_i_over_sigma:
        extra_note="I/sigma shown is value \nof max_i_over_sigma."
      elif not self.input_i_over_sigma:
        extra_note="I/sigma shown achieves about %3.0f%% of \nmaximum anomalous signal." %(self.ratio_for_failure*100.)
      out.show_text("""
Note: Target anomalous signal not achievable with tested I/sigma (up to %d )
for resolutions of %5.2f A and lower. %s
""" % (int(self.max_i_over_sigma),self.missed_target_resolutions[0],extra_note))

    if self.skipped_resolutions:
      self.skipped_resolutions.sort()
      out.show_text("""
Note: No plausible values of I/sigma found for  resolutions of %5.2f A
and lower.
""" % (self.skipped_resolutions[0]))


    out.show_text("""
This table says that if you collect your data to a resolution of %5.1f A with
an overall <I>/<sigma> of about %3.0f then the half-dataset anomalous
correlation should be about %5.2f (typically within a factor of 2).  This
should lead to a correlation of your anomalous data to true anomalous
differences (CC*_ano) of about %5.2f, and a useful anomalous signal around
%3.0f (again within a factor of about two). With this value of estimated
anomalous signal the probability of finding the anomalous substructure is
about %3d%% (based on estimated anomalous signal and actual outcomes for
real structures.)  """ % (dmin, i_over_sigma,  cc_half,  cc_ano,
        s_ano, int(solved)))
    out.show_text("""
The value of rms(sigF)/rms(F) is approximately the inverse of I/sigma. The
calculations are based on rms(sigF)/rms(F).

Note that these values assume data measured with little radiation damage or at
least with anomalous pairs measured close in time. The values also assume that
the anomalously-scattering atoms are nearly as well-ordered as other atoms.
If your crystal does not fit these assumptions it may be necessary to collect
data with even higher I/sigma than indicated here.

Note also that anomalous signal is roughly proportional to the anomalous
structure factors at a given resolution. That means that if you have 50%
occupancy of your anomalous atoms, the signal will be 50% of what it otherwise
would be.  Also it means that if your anomalously scattering atoms only
contribute to 5 A, you should only consider data to 5 A in this analysis.
""")
    out.show_paragraph_header("""What to do next:""")
    out.show_text("""
1. Collect your data, trying to obtain a value of I/sigma for the whole dataset
   at least as high as your target.""")
    out.show_text("""\
2. Scale and analyze your unmerged data with phenix.scale_and_merge to get
   accurate scaled and merged data as well as two half-dataset data files
   that can be used to estimate the quality of your data.""")
    out.show_text("""\
3. Analyze your anomalous data (the scaled merged data and the two half-datdaset
   data files) with phenix.anomalous_signal to estimate the anomalous signal
   in your data. This tool will again guess the fraction of the substructure
   that can be obtained with your data, this time with knowledge of the
   actual anomalous signal.  It will also estimate the figure of merit of
   phasing that you can obtain once you solve the substruture. """)
    out.show_text("""\
4. Compare the anomalous signal in your measured data with the
   estimated values in the table above. If they are lower than expected
   you may need to collect more data to obtain the target anomalous signal.""")
