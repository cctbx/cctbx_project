
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
     ratio_for_failure=0.95):
  closest_sigf=None
  closest_dist2=None
  start_value=1
  if include_zero: start_value=0
  for i in xrange(start_value,ntries):
    sigf=i*1./float(ntries)
    s_ano,cc_ano,cc_half,fpp_weak,cc_ano_weak,cc_half_weak,i_over_sigma=\
      get_values_from_sigf(nrefl,nsites,natoms,z,fpp,sigf,
          fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
          fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy,
          get_fpp_weak=False)
    if min_cc_ano is not None and cc_ano_weak < min_cc_ano:
      continue

    dist2=(s_ano-target_s_ano)**2
    if closest_dist2 is None or dist2<closest_dist2:
     closest_dist2=dist2
     closest_sigf=sigf
  if closest_sigf==0:  # try with target of 90% of maximum available
    s_ano,cc_ano,cc_half,fpp_weak,cc_ano_weak,cc_half_weak,i_over_sigma=\
       get_values_from_sigf(nrefl,nsites,natoms,z,fpp,closest_sigf,
          fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
          fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy,
          get_fpp_weak=False)
    closest_sigf=get_sigf(
     nrefl,nsites,natoms,z,fpp,target_s_ano=s_ano*ratio_for_failure,
     ntries=ntries,
     min_cc_ano=min_cc_ano,
     fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
     fo_list=fo_list,fo_number_list=fo_number_list,
     occupancy=occupancy,include_zero=False)
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
       get_fpp_weak=True):
    s_ano=get_sano(nrefl,nsites,natoms,z,fpp,sigf,
      fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
      fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)
    cc_ano=get_cc_ano(nrefl,nsites,natoms,z,fpp,sigf,
      fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
      fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)
    cc_half=get_cc_half(cc_ano,disorder_parameter=disorder_parameter)

    if get_fpp_weak:
      fpp_weak=estimate_fpp_weak(nrefl,nsites,natoms,z,fpp,sigf,
        fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
        fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)

      cc_ano_weak=get_cc_ano(nrefl,nsites,natoms,z,fpp_weak,sigf,
        fa2=fa2,fb2=fb2,disorder_parameter=disorder_parameter,
        fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy)
      cc_half_weak=get_cc_half(cc_ano_weak,
          disorder_parameter=disorder_parameter)
      i_over_sigma=get_i_over_sigma_from_sigf(sigf)
      return s_ano,cc_ano,cc_half,fpp_weak,cc_ano_weak,cc_half_weak,i_over_sigma
    else:
       return s_ano,cc_ano,cc_half,fpp,cc_ano,cc_half,0.0

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
      include_weak_anomalous_scattering=True,
      intrinsic_scatterers_as_noise=None,
      ratio_for_failure=0.95,
      i_over_sigma=None,
      quiet=False) :
    self.chain_type = chain_type
    self.residues = residues
    self.solvent_fraction = solvent_fraction
    self.nsites = nsites
    self.fpp = fpp
    self.target_s_ano = target_s_ano
    self.min_cc_ano = min_cc_ano
    self.wavelength = wavelength
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


    z=7
    if dmin:
      self.dmin_ranges=[dmin]
    else:
      self.dmin_ranges=[6,5,3,2.5,2,1.5]
    self.table_rows = []
    self.representative_values = None
    self.skipped_resolutions = []
    self.missed_target_resolutions = []
    self.input_i_over_sigma=i_over_sigma
    for dmin in self.dmin_ranges:
      # Guess reflections from residues, dmin, solvent fraction

      nrefl=get_nrefl(
         residues=residues,dmin=dmin,solvent_fraction=solvent_fraction)

      # identify rms(sigF)/rms(F) necessary to get target_s_ano with this
      # many reflections, sites, atoms, f" value

      if i_over_sigma is None:
        sigf=get_sigf(nrefl,nsites,self.natoms,z,fpp,target_s_ano=target_s_ano,
          min_cc_ano=min_cc_ano,
          fa2=self.fa2,fb2=self.fb2,disorder_parameter=self.disorder_parameter,
          fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy,
          ratio_for_failure=self.ratio_for_failure)
      else:  # input i_over_sigma...estimate sigf
        sigf=get_sigf_from_i_over_sigma(i_over_sigma)
      if sigf is None: continue  # hopeless
      # what are expected signal, useful cc_ano, cc_half-dataset, <I>/<sigI>

      s_ano,cc_ano,cc_half,fpp_weak,cc_ano_weak,cc_half_weak,\
         local_i_over_sigma=\
        get_values_from_sigf(nrefl,nsites,self.natoms,z,fpp,sigf,
          fa2=self.fa2,fb2=self.fb2,disorder_parameter=self.disorder_parameter,
          fo_list=fo_list,fo_number_list=fo_number_list,occupancy=occupancy,
          get_fpp_weak=True)

      if local_i_over_sigma>=999:
        self.skipped_resolutions.append(dmin)
        continue  # hopeless
      if s_ano<0.95*target_s_ano:  # must not be able to get target s_ano
        self.missed_target_resolutions.append(dmin)

      self.table_rows.append([
        "%5.2f" % dmin,
        "%7d" % nrefl,
        "%6.0f" % local_i_over_sigma,
        "%7.1f" % (100.*sigf),
        "%5.2f - %5.2f" % (cc_half_weak, cc_half),
        "%5.2f - %5.2f" % (cc_ano_weak, cc_ano),
        "%3.0f - %3.0f" % (s_ano/2, s_ano),
      ])
      if ((self.representative_values is None) or
          len(self.dmin_ranges) < 2 or dmin == self.dmin_ranges[-2]) :
        self.representative_values = [dmin,nsites,nrefl,fpp,local_i_over_sigma,
           sigf,cc_half_weak,cc_half,cc_ano_weak,cc_ano,s_ano]

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

  def show_summary(self):
    if self.is_solvable():
      print """
I/sigI:  %7.1f
Dmin:     %5.2f
cc_half:  %5.2f - %5.2f
cc*_anom: %5.2f - %5.2f
Signal:   %5.1f - %5.1f
""" %(
 self.representative_i_over_sigma(),
 self.representative_dmin(),
 self.representative_cc_half_weak(),
 self.representative_cc_half(),
 self.representative_cc_ano_weak(),
 self.representative_cc_ano(),
 self.representative_s_ano()/2,
 self.representative_s_ano())

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
    out.show_sub_header("(Targets for entire dataset)")
    if (len(self.table_rows) == 0) :
      out.show_text("SAD solution unlikely with the given parameters.")
      return
    if (not out.gui_output) :
      out.show_preformatted_text("""
                                         Anomalous      Useful        Useful
                                        Half-dataset    Anom CC      Anomalous
 Dmin   Nrefl   I/sigI rms(sigF)/           CC         (cc*_anom)     Signal
                         rms(F) (%) """)
      for row in self.table_rows :
        out.show_preformatted_text("%s %s %s   %s        %s  %s   %s" %
          tuple(row))
    else :
      table = table_utils.simple_table(
        table_rows=self.table_rows,
        column_headers=["d_min", "Nrefl", "I/sigI", "rms(sigF)/rms(F) (%)",
          "Half-dataset CC_ano", "CC*_ano", "Anom. signal"])
      out.show_table(table)
    (dmin,nsites,nrefl,fpp,i_over_sigma,sigf,cc_half_weak,cc_half,cc_ano_weak,
      cc_ano,s_ano) = tuple(self.representative_values)

    if self.missed_target_resolutions:
      self.missed_target_resolutions.sort()
      extra_note=""
      if not self.input_i_over_sigma:
        extra_note="I/sigma shown achieves about %3.0f%% of maximum anomalous signal." %(self.ratio_for_failure*100.)
      out.show_text("""
Note: Target anomalous signal not achievable with tested I/sigma for
resolutions of %5.2f A and lower. %s
""" % (self.missed_target_resolutions[0],extra_note))

    if self.skipped_resolutions:
      self.skipped_resolutions.sort()
      out.show_text("""
Note: No plausible values of I/sigma found for  resolutions of %5.2f A
and lower.
""" % (self.skipped_resolutions[0]))


    out.show_text("""
This table says that if you collect your data to a resolution of %5.1f A with
an overall <I>/<sigma> of about %3.0f then the half-dataset anomalous
correlation should be in the range of %5.2f - %5.2f.  This should lead to a
correlation of your anomalous data to the useful anomalous differences
(CC*_ano) in the range of %5.2f - %5.2f and a useful anomalous signal around
%3.0f - %3.0f, where an anomalous correlation of about %5.2f and an
anomalous signal of 10-15 should be sufficient to find the substructure and
calculate phases.
""" % (dmin, i_over_sigma, cc_half_weak, cc_half, cc_ano_weak, cc_ano,
       s_ano/2., s_ano, self.min_cc_ano))
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
2. Scale and analyze your unmerged data with phenix.scale_and_merge to get an
   accurate estimate of your half-dataset anomalous correlation and your
   estimated useful anomalous correlation cc*_anom.""")
    out.show_text("""\
3. Compare the half-dataset anomalous correlation and cc*_anom with the
   estimated values in the table above. If they are lower than expected
   you may need to collect more data to obtain the target half-dataset
   correlation and target anomalous signal.
""")
