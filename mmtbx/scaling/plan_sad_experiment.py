
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
     min_cc_ano=None):
  closest_sigf=None
  closest_dist2=None
  for i in xrange(ntries):
    sigf=i*1./float(ntries)
    s_ano,cc_ano,cc_half,fpp_weak,cc_ano_weak,cc_half_weak,i_over_sigma=\
      get_values_from_sigf(nrefl,nsites,natoms,z,fpp,sigf)
    if min_cc_ano is not None and cc_ano_weak < min_cc_ano: continue

    dist2=(s_ano-target_s_ano)**2
    if closest_dist2 is None or dist2<closest_dist2:
     closest_dist2=dist2
     closest_sigf=sigf
  return closest_sigf

def get_sano(nrefl,nsites,natoms,z,fpp,sigf):
  # nrefl = number of anomalous data
  # nsites = number of anomalously scattering atoms in asymmetric unit
  # natoms = number of non-anomalously scattering atoms in asymmetric unit
  # z = typical Z of non-anomalously scattering atoms
  # sigf = rms(sigF)/rms(F) overall
  # delta_b_ano=B for anom atoms minus B for rest of structure.
  #  Ignoring this effect here
  sano2=(4./5.)*(nrefl/nsites)/(1.+(natoms/nsites)*(z*sigf/fpp)**2)
  # note this is very close to:
  #  sano2=(4./5.)*nrefl*fpp**2/(natoms*z**2*sigf**2)
  #  sano ~ sqrt(4./5.)*(fpp/z*sigf)*sqrt(nrefl/natoms) which is proportional to fpp and 1/sigf

  return math.sqrt(sano2)

def get_cc_ano(nrefl,nsites,natoms,z,fpp,sigf):
  "cc_ano**2=1/[1 + <sigF**2>/(n/N)(Za/Z)**2<F**2>]"
  ccano2=1./(1.+sigf**2*natoms*z**2/(nsites*fpp**2))
  return math.sqrt(ccano2)

def get_cc_half(cc_ano):
  "cc_half=cc_ano_perf**2/(2-cc_ano_perf**2)"
  return cc_ano*cc_ano/(2.-cc_ano*cc_ano)

def estimate_fpp_weak(nrefl,nsites,natoms,z,fpp,sigf):
  # estimate fpp that would give half the sano value
  sano=get_sano(nrefl,nsites,natoms,z,fpp,sigf)
  target=sano/2.
  for i in xrange(10,101):
    scale=i*0.01
    ss=get_sano(nrefl,nsites,natoms,z,fpp*scale,sigf)
    if ss >=target:
      return fpp*scale
  return fpp

def get_i_over_sigma(sigf):
  # <I/sigI> = 0.88 * rms(F)/rms(sigF) r**2=0.82
  if sigf>0:
    return 0.88/sigf
  else:
    return 999.

def get_values_from_sigf(nrefl,nsites,natoms,z,fpp,sigf):
    s_ano=get_sano(nrefl,nsites,natoms,z,fpp,sigf)
    cc_ano=get_cc_ano(nrefl,nsites,natoms,z,fpp,sigf)
    cc_half=get_cc_half(cc_ano)

    fpp_weak=estimate_fpp_weak(nrefl,nsites,natoms,z,fpp,sigf)
    cc_ano_weak=get_cc_ano(nrefl,nsites,natoms,z,fpp_weak,sigf)
    cc_half_weak=get_cc_half(cc_ano_weak)
    i_over_sigma=get_i_over_sigma(sigf)

    return s_ano,cc_ano,cc_half,fpp_weak,cc_ano_weak,cc_half_weak,i_over_sigma

def get_fpp(atom_type=None,wavelength=None,out=sys.stdout):
  if not atom_type or not wavelength:
    raise Sorry("Please specify either f_double_prime or " +
      "atom_type and wavelength")
  from cctbx.eltbx import sasaki
  try:
    table = sasaki.table(atom_type)
    fp_fdp = table.at_angstrom(wavelength)
  except ValueError :
    raise Sorry("Unable to get scattering factors for %s" %(atom_type))
  return fp_fdp.fdp()

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
  number_of_sites,number_of_sites_lowres=get_number_of_sites(
      atom_type=atom_type,n_met=n_met,n_cys=n_cys,
      n_aa=n_aa,ncs_copies=1,out=null_out())
  return n_aa,number_of_sites

class estimate_necessary_i_sigi (mmtbx.scaling.xtriage_analysis) :
  def __init__ (self,
      chain_type='PROTEIN',
      residues=250,
      solvent_fraction=0.50,
      nsites=5,
      fpp=3.8,
      target_s_ano=30,
      min_cc_ano=0.15,
      dmin=None,
      quiet=False) :
    self.chain_type = chain_type
    self.residues = residues
    self.solvent_fraction = solvent_fraction
    self.nsites = nsites
    self.fpp = fpp
    self.target_s_ano = target_s_ano
    self.min_cc_ano = min_cc_ano
    self.natoms = residues*get_atoms_per_residue(chain_type=chain_type)
    vol_per_residue = get_vol_per_residue(chain_type=chain_type)
    z=7
    if dmin:
      self.dmin_ranges=[dmin]
    else:
      self.dmin_ranges=[6,5,3,2.5,2,1.5]
    self.table_rows = []
    self.representative_values = None
    self.skipped_resolutions = []
    for dmin in self.dmin_ranges:
      # Guess reflections from residues, dmin, solvent fraction
      nrefl=get_nrefl(
         residues=residues,dmin=dmin,solvent_fraction=solvent_fraction)

      # identify rms(sigF)/rms(F) necessary to get target_s_ano with this
      # many reflections, sites, atoms, f" value
      sigf=get_sigf(nrefl,nsites,self.natoms,z,fpp,target_s_ano=target_s_ano,
        min_cc_ano=min_cc_ano)

      # what are expected signal, true cc_ano, cc_half-dataset, <I>/<sigI>
      s_ano,cc_ano,cc_half,fpp_weak,cc_ano_weak,cc_half_weak,i_over_sigma=\
        get_values_from_sigf(nrefl,nsites,self.natoms,z,fpp,sigf)
      if i_over_sigma>=999:
        self.skipped_resolutions.append(dmin)
        continue  # hopeless
      self.table_rows.append([
        "%5.2f" % dmin,
        "%7d" % nrefl,
        "%6.0f" % i_over_sigma,
        "%7.3f" % sigf,
        "%5.2f - %5.2f" % (cc_half_weak, cc_half),
        "%5.2f - %5.2f" % (cc_ano_weak, cc_ano),
        "%3.0f - %3.0f" % (s_ano/2, s_ano),
      ])
      if ((self.representative_values is None) or
          len(self.dmin_ranges) < 2 or dmin == self.dmin_ranges[-2]) :
        self.representative_values = [dmin,nsites,nrefl,fpp,i_over_sigma,sigf,
           cc_half_weak,cc_half,cc_ano_weak,cc_ano,s_ano]

  def is_solvable (self) :
    return (self.representative_values is not None)

  def _show_impl (self, out) :
    out.show_header("SAD experiment planning")
    out.show_paragraph_header(
      "Dataset overall I/sigma required to solve a structure")
    out.show_preformatted_text("""\
  Target anomalous signal: %(target_s_ano)7.1f
  Residues: %(residues)d
  Solvent_fraction: %(solvent_fraction)7.2f
  Atoms: %(natoms)d
  Sites: %(nsites)d
  f-double-prime: %(fpp)7.2f
""" % self.__dict__)
    out.show_paragraph_header(
      "Dataset <I>/<sigI> needed to obtain anomalous signal of 15-30")
    out.show_sub_header("Targets for entire dataset")
    if (len(self.table_rows) == 0) :
      out.show_text("SAD solution unlikely with the given parameters.")
      return
    if (not out.gui_output) :
      out.show_preformatted_text("""
                                        Half-dataset     True        Anomalous
Dmin   Nrefl    I/sigI rms(sigF)/rms(F)   CC_ano        CC*_ano       Signal""")
      for row in self.table_rows :
        out.show_preformatted_text("%s %s %s   %s        %s  %s   %s" %
          tuple(row))
    else :
      table = table_utils.simple_table(
        table_rows=self.table_rows,
        column_headers=["d_min", "Nrefl", "I/sigI", "rms(sigF)/rms(F)",
          "Half-dataset CC_ano", "True CC*_ano", "Anom. signal"])
      out.show_table(table)
    (dmin,nsites,nrefl,fpp,i_over_sigma,sigf,cc_half_weak,cc_half,cc_ano_weak,
      cc_ano,s_ano) = tuple(self.representative_values)
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
true correlation of your anomalous data to the true anomalous differences
(CC*_ano) in the range of %5.2f - %5.2f and an anomalous signal in the range
of %3.0f - %3.0f, where a true anomalous correlation of about %5.2f and an
anomalous signal of 10-15 nshould be sufficient to find the substructure and
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
   accurate estimate of your half-dataset anomalous correlation.""")
    out.show_text("""\
3. Compare the half-dataset anomalous correlation with the estimated values in
   the table above. If they are lower than expected you may need to collect
   more data to obtain the target half-dataset correlation and target anomalous
   signal.
""")
