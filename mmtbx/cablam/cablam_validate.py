from __future__ import absolute_import, division, print_function
# (jEdit options) :folding=explicit:collapseFolds=1:
#
#cablam_validate
#Author: Christopher Williams, contact christopher.j.williams@duke.edu
#
#cablam_validate is part of the cablam system for use in Phenix and cctbx,
#  specifically, cablam_validate uses contour information derived from
#  cablam_training to validate the backbone of protein models.
#cablam_validate accepts a protein model (as pdb or Phenix hierarchy) and
#  calculates backbone geometry for each residue.  The backbone geometry is
#  compared against expected behavior and outliers are marked.  Each outlier is
#  then checked for similarity to regular secondary structure types.  Numerical
#  values are returned to indicate which seconday structure type each outlier is
#  likely to be.
#cablam_validate can return this data as machine-friendly, comma-separated text.
#  It can also print validation markup in kinemage format, which can be appended
#  to an existing protein kinemage.
#Within the Phenix environment, cablam_validate can return a list of custom
#  objects (see cablam_validation class below) which includes relevant data and
#  a link to the hierarchy rg object for each outlier residue.  Be sure to pass
#  a hierarchy object into cablam_validate.run() if you want to access this
#  functionality.
#Notes on contour levels: For the CA contours, the choice of 0.5% for "allowed"
#  reflects a level around which a significant change occurs in the behavior of
#  these contours (new minor clusters appear and bridges form between existing
#  clusters) and which provideds good covereage of the secondary structure motif
#  contours.  The choice of 2% for "favored" is an intentional parallel to Rama
#  validation.
#
#2012-08-03: Initial upload
#  A "run whole protein to find probable structure" function might be nice
#2012-09-05:
#  General cleanup. usage() help message and interpretation() guide to output.
#  analyze_pdb() is now the easiest way to run the script from within phenix.
#2012-10-09:
#  The "outliers" object returned by analyze_pdb is now a dict instead of a list
#  The keys are the same as those used by cablam_res. cablam_measures now
#  handles all the calculations setup() needs.
#2013-02-01:
#  Added oneline output. Supports dir of files, as well as single files. Added
#  checks for CA-only contours. Added find_partial_sec_struc and
#  find_whole_sec_struc to piece individual residue assignments together into
#  complete secondary structure elements. Added multicrit kinemage printing.
#2013-09-17:
#  Added support for threeten helices throughout.
#2014-02-07:
#  Added test validation for cis vs tran proline.
#  Added stand-alone output for HELIX/SHEET records.
#2014-03-07: cleaned up commandline control for outputs, now all "output="
#Next: iotbx.file_reader incorporated to control input, cleaned up oneline output
#Next: give_text now in MolProbity format (colons and cnit ids)

import os, sys
import libtbx.phil.command_line #argument parsing
from iotbx import pdb  #contains hierarchy data structure
from iotbx import file_reader
from mmtbx.cablam import cablam_res #contains a data structure derived from
#  hierarchy, but more suited to cablam's needs - specifically it can hold
#  geometric and probe measures and can look forward and backward in sequence
from mmtbx.cablam import cablam_math #contains geometric measure calculators
from mmtbx.rotamer.n_dim_table import NDimTable #handles contours
from libtbx import easy_pickle #NDimTables are stored as pickle files
from libtbx import easy_run
import libtbx.load_env

#{{{ phil
#-------------------------------------------------------------------------------
master_phil = libtbx.phil.parse("""
cablam_validate {
  pdb_infile = None
    .type = path
    .help = '''input PDB file or dirpath'''
  output = *text kin ca_kin full_kin markup_no_ribbons points records records_plus_pdb oneline
    .type = choice
    .help = '''choose output type,
    =text for comma-separated outlier summary
    =kin print markup kinemage to screen
    =ca_kin print CA geo markup kinemage to screen
    =full_kin open in King a multi-crit-type kinemage with markup
    =markup_no_ribbons print allowed, outlier, and ca geom outlier markup kin
    =points print cablam-space points to sys.stdout, in kinemage dotlist format
    =records print HELIX and SHEET-style records based on CaBLAM analysis
    =records_plus_pdb print HELIX and SHHET-style records as part of a pdb file
    =oneline print one-line validation to sys.stdout: percent of cablam outliers
    '''
  outlier_cutoff = 0.05
    .type = float
    .help = '''sets the contour level for detecting outliers, e.g. outlier_cutoff=0.01 for accepting the top 99%, defaults to 0.05'''
  check_prolines = False
    .type = bool
    .help = '''test of a function to differentiate cis vs trans prolines'''
  help = False
    .type = bool
    .help = '''help and data interpretation messages'''
}
""", process_includes=True)
#-------------------------------------------------------------------------------
#}}}

#{{{ classes: cablam_validation, motif_guess, and motif_chunk
#-------------------------------------------------------------------------------
#This object holds information on one residue for purposes of cablam_validate
#It's mostly just a package for data passing and access
class cablam_validation(object):
  def __init__(self, residue=None,outlier_level=None,ca_geom_outlier_level=None):
    self.residue = residue #cablam_res.linked_residue object
    self.outlier_level = outlier_level #percentile level in "peptide_expectations" contours
    self.ca_geom_outlier_level = ca_geom_outlier_level
    self.loose_alpha  = None #percentile level in motif contour
    self.regular_alpha= None #percentile level in motif contour
    self.loose_beta   = None #percentile level in motif contour
    self.regular_beta = None #percentile level in motif contour
    self.loose_threeten = None
    self.regular_threeten = None
    self.suggestion = ""
    if self.residue:
      self.rg = self.residue.rg #rg object from a source hierarchy
    else:
      self.rg = None

#Class for holding percentile contour level values in residue object
class motif_value(object):
  def __init__(self):
    self.loose_alpha = None
    self.regular_alpha = None
    self.loose_beta = None
    self.regular_beta = None
    self.loose_threeten = None
    self.regular_threeten = None

class motif_guess(object):
  def __init__(self):
    self.loose_alpha = False
    self.regular_alpha = False
    self.loose_beta = False
    self.regular_beta = False
    self.loose_threeten = False
    self.regular_threeten = False

class motif_chunk(object):
  def print_record(self,writeto=sys.stdout):
    if self.motif_type == 'helix':
      self.print_helix_record(writeto=writeto)
    elif self.motif_type == 'threeten':
      self.print_threeten_record(writeto=writeto)
    elif self.motif_type == 'sheet':
      self.print_sheet_record(writeto=writeto)
    else:
      sys.stderr.write(
        '\ntried to print a record that was not \'helix\' or \'sheet\'\n')
      sys.exit()

  def print_helix_record(self,helix_num=1,writeto=sys.stdout):
    motif_start_id = self.motif_start.id_with_resname()
    motif_end_id = self.motif_end.id_with_resname()
    motif_len = self.motif_end.resnum - self.motif_start.resnum + 1
    if motif_len >= 3:
      #Note: the following code makes me hate HELIX record formatting
      writeto.write('HELIX  '+ '%3i' %helix_num +' '+ '%3i' %helix_num +' '+
        motif_start_id[:5]+' '+motif_start_id[5:]+' '+  motif_end_id[:5] +' '+
        motif_end_id[5:] +' 1                               '+'%5i'%motif_len+
        '\n')

  def print_threeten_record(self,helix_num=1,writeto=sys.stdout):
    motif_start_id = self.motif_start.id_with_resname()
    motif_end_id = self.motif_end.id_with_resname()
    motif_len = self.motif_end.resnum - self.motif_start.resnum + 1
    if motif_len >= 1:
      #Note: the following code makes me hate HELIX record formatting
      writeto.write('HELIX  '+ '%3i' %helix_num +' '+ '%3i' %helix_num +' '+
        motif_start_id[:5]+' '+motif_start_id[5:]+' '+  motif_end_id[:5] +' '+
        motif_end_id[5:] +' 5                               '+'%5i'%motif_len+
        '\n')

  def print_sheet_record(self,sheet_id='U',strand_num=1,strand_count=1,
    strand_sense=0,writeto=sys.stdout):
    motif_start_id = self.motif_start.id_with_resname()
    motif_end_id = self.motif_end.id_with_resname()
    writeto.write('SHEET  '+'%3i' %strand_num +' '+'%3s'%sheet_id +
      '%2i'%strand_count +' '+ motif_start_id+' '+motif_end_id +
      '%2i'%strand_sense+'\n')

  def __init__(self):
    self.motif_start = None #will hold residue object
    self.motif_end = None #will hold residue object
    self.motif_type = None #'helix' or 'sheet' for the moment
#-------------------------------------------------------------------------------
#}}}

#{{{ usage
#-------------------------------------------------------------------------------
def usage():
  sys.stderr.write("""
phenix.cablam_validate file.pdb [options ...]

Options:

  pdb_infile=filename      input PDB file
                             some output formats support a dir of files in this
                             field
  outlier_cutoff=0.05      sets the contour level for detecting outliers,
                             e.g. outlier_cutoff=0.01 for accepting the top 99%,
                             defaults to 0.05
  output=
                           kin : prints markup kinemage for a single outlier
                             cutoff to screen
                           text : prints machine-readable columnated and comma-
                             separated data to screen
                             (default output)
                           points : prints cablam-space points to screen,
                             in kinemage dotlist format
                           records : prints HELIX and SHEET-style records to
                             screen, based on CaBLAM secondary structure
                             analysis
                           records_plus_pdb : prints HELIX and SHEET records and
                             the associated pdb file to screen
                           ca_kin : prints a kinemage with ca-contour outliers
                             marked
                           full_kin : opens a phenix.king window with multicrit
                             style markup for the submitted structure. Also
                             saves a .pdb with HELIX and SHEET records and a
                             .kin with other markup to working dir
                           markup_no_ribbons : prints 'allowed', 'outlier', and
                             'ca geom outlier' kinemage markup to screen
                           oneline : prints oneline-style output to screen.
                             Supports printing for multiple files

  help=False               prints this usage text, plus notes on data
                             interpretation to screen

Example:

phenix.cablam_validate file.pdb output=full_kin
--------------------------------------------------------------------------------
""")
#-------------------------------------------------------------------------------
#}}}

#{{{ interpretation
#-------------------------------------------------------------------------------
#This function holds a print-to-screen explanation of what the numbers mean
def interpretation():
  sys.stderr.write("""
cablam_validate data interpretation:

Text:
  Text output is provided in comma-separated format:
  residue,contour_level,loose_alpha,regular_alpha,loose_beta,regular_beta

  'residue' is a residue identifier formatted as pdb columns 18 to 27 (1 index)

  'contour_level' is the 3D outlier contour at which this residue was found.
    Smaller values are worse outliers.  Residues with values <=0.05 are worth
    concern.

  'loose_alpha' is a 2D percentile contour for this residue's similarity to
    alpha helix.  A high value indicates some amount of helix character, but not
    necessarily the presence of a regular helix.
    >0.001 is meaningful, >0.01 is likely

  'regular_alpha' is a 2D percentile contour for this residue's similarity to
    regular alpha helix.  The contours for regular helix are very tight.
    >0.00001 is meaningful, >0.0001 is likely

  'loose_beta' is a 2D percentile contour for this residue's similarity to
    beta strand.  A high value indicates some amount of strand character, but
    not necessarily the presence of a regular strand.
    >0.01 is meaningful, >0.1 is likely

  'regular_beta' is a 2D percentile contour for this residue's similarity to
    regular beta sheet.  Beta structure contours are generally more permissive
    than helix contours, and the contour values should not be judged on the same
    scale.
    >0.001 is meaningful, >0.01 is likely

    ('meaningful' and 'likely' scores are subject to change as the system is
    refined)

Kin:
  Kin output is provided in a kinemage format and should be appended to an open
    kinemage of the structure of interest.  This provides validation markup of
    the structure.  Markup takes the form of purple lines, drawn to follow the
    peptide plane diherdrals used to determine outliers.  Clicking on the
    vertices or on a point in the middle of the center line will bring up
    validation numbers.

  The first value corresponds to 'contour_level' in the text output.

  The values preceded by 'a' and 'rega' correspond to 'loose_alpha' and
    'regular_alpha', respectively.

  The values preceded by 'b' and 'regb' correspond to 'loose_beta' and
    'regular_beta', respectively.

Points:
  Points output is provided in kinemage dotlist format.  Each point corresponds
    to one outlier residue and is plotted in the 3D cablam space used to
    determine outliers (dimensions are CA_pseudodihedral_in,
    CA_pseudodihedral_out, and peptide_plane_psedudodihedral).  These points may
    be appended to contour kinemages to visualize outliers in that space.
  This output is intended primarily for developer and exploratory use, rather
    than for validation per se.

""")
#-------------------------------------------------------------------------------
#}}}

#{{{ setup
#-------------------------------------------------------------------------------
#Wrapper for the construction of the resdata object needed for validation
#  All you need is a hierarchy object
def setup(hierarchy,pdbid='pdbid'):
  resdata=cablam_res.construct_linked_residues(hierarchy,
    targetatoms=['CA','C','N','O'],pdbid=pdbid)
  cablam_math.cablam_measures(resdata)
  cablam_math.omegacalc(resdata)
  return resdata
#-------------------------------------------------------------------------------
#}}}

#{{{ fetch_expectations
#-------------------------------------------------------------------------------
#This function finds, unpickles, and returns N-Dim Tables of expected residue
#  behavior for use in determining outliers
#The return object is a dict keyed by residue type: 'general','gly','pro'
#One set of contours defines peptide behavior (CA_d_in,CA_d_out,peptide):
def fetch_peptide_expectations():
  categories = ['general','gly','transpro','cispro']
  unpickled = {}
  for category in categories:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.expected."+category+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for category "+
        category+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[category] = ndt
  return unpickled

#One set of contours defines CA trace (CA_d_in,CA_d_out,CA_a):
def fetch_ca_expectations():
  categories = ['general','gly','transpro','cispro']
  unpickled = {}
  for category in categories:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.expected."+category+"_CA.pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for category "+
        category+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[category] = ndt
  return unpickled
#-------------------------------------------------------------------------------
#}}}

#{{{ fetch_motifs
#-------------------------------------------------------------------------------
#This function finds, unpickles, and returns N-Dim Tables of secondary structure
#  behavior for use in determining what an outlier residue might really be
#The return object is a dict keyed by secondary structure type:
#  'loose_beta','regular_beta','loose_alpha','regular_alpha'
def fetch_motifs():
  motifs = ['loose_beta','regular_beta','loose_alpha','regular_alpha',
  'loose_threeten','regular_threeten']
  unpickled = {}
  for motif in motifs:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.motif."+motif+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for motif "+
        motif+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[motif] = ndt
  return unpickled
#-------------------------------------------------------------------------------
#}}}

#{{{ fetch cis trans pro
#-------------------------------------------------------------------------------
#This function finds, unpickles, and returns N-Dim Tables of proline behavior
#  for use in determining cis versus trans peptide confirmation
#The return object is a dict keyed by: 'cis' and 'trans'
def fetch_cis_trans_proline():
  confs = ['cis','trans']
  unpickled = {}
  for conf in confs:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.proline."+conf+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for "+
        conf+" proline in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[conf] = ndt
  return unpickled
#-------------------------------------------------------------------------------
#}}}

#{{{ find_outliers
#-------------------------------------------------------------------------------
#These are used by the helix_or_sheet function
#For each residue, check the relevant measures against 3D contours and record
#  outliers. Outliercutoffs can be set on use, but default to generally useful
#  values. cutoff=1.0 will consider all residues outliers and is a quick way to
#  get contour data for all residues.
def find_peptide_outliers(resdata,expectations,cutoff=0.05):
  outliers = {}
  reskeys = resdata.keys()
  reskeys.sort()
  for resid in reskeys:
    residue = resdata[resid]
    if 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures and 'CO_d_in' in residue.measures:
      cablam_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CO_d_in']]
      resname = residue.alts[residue.firstalt('CA')]['resname']
      if resname.upper() == 'GLY':
        percentile = expectations['gly'].valueAt(cablam_point)
      elif resname.upper() == 'PRO':
        #Splitting cis vs trans at 90, rather than having a twisted category
        #  This way, no residues are left out of cablam analysis
        #  Use omegalyze for nontrans peptide validation
        if ('omega' not in residue.measures) or (residue.measures['omega'] >= 90) or (residue.measures['omega'] <= -90):
          percentile = expectations['transpro'].valueAt(cablam_point)
        else:
          percentile = expectations['cispro'].valueAt(cablam_point)
        #percentile = expectations['pro'].valueAt(cablam_point)
      else:
        percentile = expectations['general'].valueAt(cablam_point)

      if percentile < cutoff:
        outliers[resid] = cablam_validation(residue=residue,outlier_level=percentile)
  return outliers

def find_ca_outliers(resdata,expectations,cutoff=0.005):
  outliers = {}
  reskeys = resdata.keys()
  reskeys.sort()
  for resid in reskeys:
    residue = resdata[resid]
    if 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures and 'CA_a' in residue.measures:
      cablam_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CA_a']]
      resname = residue.alts[residue.firstalt('CA')]['resname']
      if resname.upper() == 'GLY':
        percentile = expectations['gly'].valueAt(cablam_point)
      elif resname.upper() == 'PRO':
        if ('omega' not in residue.measures) or (residue.measures['omega'] >= 90) or (residue.measures['omega'] <= -90):
          percentile = expectations['transpro'].valueAt(cablam_point)
        else:
          percentile = expectations['cispro'].valueAt(cablam_point)
        #percentile = expectations['pro'].valueAt(cablam_point)
      else:
        percentile = expectations['general'].valueAt(cablam_point)

      if percentile < cutoff:
        outliers[resid] = cablam_validation(residue=residue,outlier_level=percentile)
  return outliers

def find_all_residue_stats(resdata,peptide_expectations,ca_geom_expectations):
  residue_stats = {}
  reskeys = resdata.keys()
  reskeys.sort()
  for resid in reskeys:
    residue = resdata[resid]
    if 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures and 'CO_d_in' in residue.measures and 'CA_a' in residue.measures:
      peptide_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CO_d_in']]
      ca_geom_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CA_a']]
      resname = residue.alts[residue.firstalt('CA')]['resname']
      if resname.upper() == 'GLY':
        peptide_percentile = peptide_expectations['gly'].valueAt(peptide_point)
        ca_geom_percentile = ca_geom_expectations['gly'].valueAt(ca_geom_point)
      elif resname.upper() == 'PRO':
        if ('omega' not in residue.measures) or (residue.measures['omega'] >= 90) or (residue.measures['omega'] <= -90):
          peptide_percentile = peptide_expectations['transpro'].valueAt(peptide_point)
          ca_geom_percentile = ca_geom_expectations['transpro'].valueAt(ca_geom_point)
        else:
          peptide_percentile = peptide_expectations['cispro'].valueAt(peptide_point)
          ca_geom_percentile = ca_geom_expectations['cispro'].valueAt(ca_geom_point)
      else:
        peptide_percentile = peptide_expectations['general'].valueAt(peptide_point)
        ca_geom_percentile = ca_geom_expectations['general'].valueAt(ca_geom_point)

      residue_stats[resid] = cablam_validation(residue=residue,outlier_level=peptide_percentile,ca_geom_outlier_level=ca_geom_percentile)
  return residue_stats
#-------------------------------------------------------------------------------
#}}}

#{{{ helix_or_sheet (wrapper) and helix_or_sheet_res
#-------------------------------------------------------------------------------
#Labels individual residues as alpha or beta
def helix_or_sheet(residue_stats, motifs):
  alpha_cutoff = 0.001
  beta_cutoff = 0.0001
  threeten_cutoff = 0.001
  for res_stat in residue_stats.values():
    residue = res_stat.residue
    contours = helix_or_sheet_res(residue, motifs)
    if contours:
      res_stat.loose_alpha   = contours['loose_alpha']
      res_stat.regular_alpha = contours['regular_alpha']
      res_stat.loose_beta    = contours['loose_beta']
      res_stat.regular_beta  = contours['regular_beta']
      res_stat.loose_threeten= contours['loose_threeten']
      res_stat.regular_threeten = contours['regular_threeten']

#this single-residue level of the check is intentionally independent from the
#  cablam_validation class in case getting contour levels for a single residue
#  is desirable
def helix_or_sheet_res(residue, motifs):
  contours = {'residue':residue}
  if 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures:
    cablam_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out']]
    for motifname, contour in motifs.items():
      contours[motifname] = contour.valueAt(cablam_point)
    return contours
  else:
    return None
#-------------------------------------------------------------------------------
#}}}

#{{{ assign_sec_struc
def assign_sec_struc(residues):
  alpha_cutoff = 0.001
  beta_cutoff = 0.0001
  threeten_cutoff = 0.001
  ca_geom_cutoff = 0.005
  for resid in residues:
    res = residues[resid]
    res.suggestion = '                 '
    if res.ca_geom_outlier_level < ca_geom_cutoff:
      continue
    if res.loose_alpha >= alpha_cutoff:
      try:
        if residues[res.residue.prevres.resid].loose_alpha >= alpha_cutoff and residues[res.residue.nextres.resid].loose_alpha >= alpha_cutoff:
          res.suggestion = ' try alpha helix '
      except KeyError:
        pass
    if res.regular_beta >= beta_cutoff:
      try:
        if residues[res.residue.prevres.resid].regular_beta >= beta_cutoff and residues[res.residue.nextres.resid].regular_beta >= beta_cutoff:
          res.suggestion = ' try beta sheet  '
      except KeyError:
        pass
    if res.loose_threeten >= threeten_cutoff and res.loose_threeten > res.loose_alpha:
      try:
        if residues[res.residue.prevres.resid].loose_threeten >= threeten_cutoff or residues[res.residue.nextres.resid].loose_threeten >= threeten_cutoff:
          res.suggestion = ' try three-ten   '
      except KeyError:
        pass
#}}}

#{{{ find_partial_sec_struc function
#-------------------------------------------------------------------------------
def find_partial_sec_struc(resdata,ca_outliers={}):

  #these numbers are low-end contour cutoffs for these motifs
  #The numbers are best-guesses from looking at marked-up structures, esp 2o01
  #  and ribosome
  #Values may change during development.
  loose_alpha_cutoff = 0.001 #0.1%
  #reg_beta_cutoff    = 0.001 #0.1%
  reg_beta_cutoff = 0.0001 #this cutoff is a bit generous, beta needs
  #  forthcoming though-space relationships to work properly
  loose_threeten_cutoff = 0.001
  reg_threeten_cutoff = 0.001

  expectations = fetch_peptide_expectations()
  motifs = fetch_motifs()

  allres = find_peptide_outliers(resdata, expectations, cutoff=1.0)
  helix_or_sheet(allres, motifs)

  for resid in resdata:
    resdata[resid].motif_guess = motif_guess()

  #Some alternative logic from earlier versions is included for reference and
  #  development
  for resid in resdata:
    if resid not in allres:
      #couldn't calculate all relevant data for the residue and cannot make
      #further assessment
      continue
    residue = resdata[resid]
    if resid in ca_outliers.keys():
      #if the residue is a ca outlier, it recieves no sec struc assignment
      #(residue.motif guess must be created for each residue, however)
      continue

    #loose alpha
    if allres[resid].loose_alpha >= loose_alpha_cutoff:
      try:
        if allres[residue.nextres.resid].loose_alpha >= loose_alpha_cutoff and allres[residue.prevres.resid].loose_alpha >= loose_alpha_cutoff:
          residue.motif_guess.loose_alpha = True
      except KeyError:
        pass

    #loose beta
    if allres[resid].regular_beta >= reg_beta_cutoff:
      try:
        if allres[residue.nextres.resid].regular_beta >= reg_beta_cutoff and allres[residue.prevres.resid].regular_beta >= reg_beta_cutoff:
          residue.motif_guess.regular_beta = True
          residue.prevres.motif_guess.regular_beta = True
          residue.nextres.motif_guess.regular_beta = True
      except KeyError:
        pass

    #loose 3-10
    if (allres[resid].loose_threeten >= loose_threeten_cutoff) and (allres[resid].loose_threeten > allres[resid].loose_alpha):
      try:
        if allres[residue.nextres.resid].loose_threeten >= loose_threeten_cutoff or allres[residue.prevres.resid].loose_threeten >= loose_threeten_cutoff:
          residue.motif_guess.loose_threeten = True
      except KeyError:
        pass
#-------------------------------------------------------------------------------
#}}}

#{{{ find_whole_sec_struc function
#-------------------------------------------------------------------------------
def find_whole_sec_struc(resdata):
  reskeys = resdata.keys()
  reskeys.sort()
  motifs = []
  current_motif = None
  #search for helix
  for resid in reskeys:
    residue=resdata[resid]
    if residue.motif_guess.loose_alpha:
      if current_motif:
        current_motif.motif_end = residue
      else: #start new motif
        current_motif = motif_chunk()
        motifs.append(current_motif)
        current_motif.motif_type = 'helix'
        current_motif.motif_start = residue
        current_motif.motif_end = residue
      if not residue.nextres:
        current_motif = None
    else:
      if current_motif: #reached end of motif
        current_motif = None
      else:
        pass
  #search for threeten
  for resid in reskeys:
    residue = resdata[resid]
    if residue.motif_guess.loose_threeten:
      if current_motif:
        current_motif.motif_end = residue
      else:
        current_motif = motif_chunk()
        motifs.append(current_motif)
        current_motif.motif_type = 'threeten'
        current_motif.motif_start = residue
        current_motif.motif_end = residue
      if not residue.nextres:
        current_motif = None
    else:
      if current_motif:
        current_motif = None
      else:
        pass
  #search for sheet
  for resid in reskeys:
    residue=resdata[resid]
    if residue.motif_guess.regular_beta:
      if current_motif:
        current_motif.motif_end = residue
      else: #start new motif
        current_motif = motif_chunk()
        motifs.append(current_motif)
        current_motif.motif_type = 'sheet'
        current_motif.motif_start = residue
        current_motif.motif_end = residue
      if not residue.nextres:
        current_motif = None
    else:
      if current_motif: #reached end of motif
        current_motif = None
      else:
        pass
  return motifs
#-------------------------------------------------------------------------------
#}}}

#{{{ check prolines function
#-------------------------------------------------------------------------------
#This is a function in development: cis-proline vs trans-proline seem to have
#  mostly-distinct distributions in CaBLAM space.  This funtion attempts to
#  determine if a proline's CA trace is modeled in agreement with its omega
#  dihedral.  It is clear that errors can be detected.  It is not yet clear how
#  to correct those errors (whether the fault is more likely to be in the CA
#  trace or in the omega dihedral. Use with care.)
def check_prolines(hierarchy,pdbid='pdbid'):
  cis_cutoff   = 0.005
  trans_cutoff = 0.005
  resdata = setup(hierarchy,pdbid)
  #cablam_math.omegacalc(resdata)
  pro_contour = fetch_cis_trans_proline()
  reskeys = resdata.keys()
  reskeys.sort()
  for resid in reskeys:
    residue = resdata[resid]
    if residue.id_with_resname()[0:3].upper() == 'PRO':
      if 'omega' in residue.measures and 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures:
        omega = residue.measures['omega']
        cablam_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out']]
        cislevel = pro_contour['cis'].valueAt(cablam_point)
        translevel = pro_contour['trans'].valueAt(cablam_point)
        if (omega >= -30) and (omega <= 30): #modeled as cis
          if cislevel < cis_cutoff:
            print("bad CIS at ", pdbid, residue.id_with_resname(), "value:%.3f" %cislevel)
            if translevel >= 0.1:#trans_cutoff:
              print("  try TRANS. value:%.3f" %translevel)
            else:
              print("  no suggestion. trans value:%.3f" %translevel)
        elif (omega >=150) or (omega <= -150): #modeled as trans
          if translevel < trans_cutoff:
            print("bad TRANS at ", pdbid, residue.id_with_resname(), "value:%.3f" %translevel)
            if cislevel >= 0.1:#cis_cutoff:
              print("  try CIS. value:%.3f" %cislevel)
            else:
              print("  no suggestion. cis value:%.3f" %cislevel)
        else: #modeled as twisted
          print("TWISTED peptide at ", pdbid, residue.id_with_resname())
          print("  TRANS score: %.3f" %translevel, " CIS score: %.3f" %cislevel)
#}}}

#{{{ print_helix_sheet_records function
#-------------------------------------------------------------------------------
#This function accepts a hierarchy object and returns HELIX and SHEET-style
#  records for major secondary structure elements found in that model by CaBLAM.
#  cablam_multicrit_kin() also returns these records, but as part of a pbd file
def print_helix_sheet_records(hierarchy, ca_cutoff=0.005, pdbid='pdbid',writeto=sys.stdout):
  resdata=setup(hierarchy,pdbid)
  ca_expectations = fetch_ca_expectations()
  motif_contours = fetch_motifs()

  ca_outliers = find_ca_outliers(resdata,ca_expectations,cutoff=ca_cutoff)

  find_partial_sec_struc(resdata,ca_outliers=ca_outliers)
  motifs = find_whole_sec_struc(resdata)
  for motif in motifs:
    motif.print_record(writeto=writeto)
#-------------------------------------------------------------------------------
#}}}

#{{{ print_pdb_with_new_helix_sheet_records function
#-------------------------------------------------------------------------------
def print_pdb_with_new_helix_sheet_records(hierarchy, ca_cutoff=0.005, pdbid='pdbid', writeto=sys.stdout):
  #This function prints CaBLAM-based HELIX and SHEET records as part of a pdb
  #  file.
  print_helix_sheet_records(hierarchy, ca_cutoff, pdbid, writeto)
  writeto.write(hierarchy.as_pdb_string())
#-------------------------------------------------------------------------------
#}}}

#{{{ print_cablam_markup_kin function
#-------------------------------------------------------------------------------
def print_cablam_markup_kin(hierarchy, peptide_cutoff=0.05, peptide_bad_cutoff=0.01, ca_cutoff=0.005, pdbid='pdbid', writeto=sys.stdout):
  #This function prints the three main kinemage markup vectorlists
  resdata=setup(hierarchy,pdbid)
  peptide_expectations = fetch_peptide_expectations()
  ca_expectations = fetch_ca_expectations()
  motif_contours = fetch_motifs()

  peptide_outliers=find_peptide_outliers(resdata,peptide_expectations,cutoff=peptide_cutoff)
  helix_or_sheet(peptide_outliers,motif_contours)
  give_kin(peptide_outliers, peptide_cutoff, color='purple', writeto=writeto)
  give_kin(peptide_outliers, peptide_bad_cutoff, color='magenta', writeto=writeto)

  ca_outliers = find_ca_outliers(resdata,ca_expectations,cutoff=ca_cutoff)
  helix_or_sheet(ca_outliers,motif_contours)
  give_ca_kin(ca_outliers, ca_cutoff, writeto=writeto)
#-------------------------------------------------------------------------------
#}}}

#{{{ cablam_multicrit_kin function
#-------------------------------------------------------------------------------
#This function accepts a hierarchy object and returns a comprehensive validation
#  kinemage, in an open phenix.king window.  Likely to return .kin-format text
#  once I get ahold of the ribbon code
def cablam_multicrit_kin(hierarchy, peptide_cutoff=0.05, peptide_bad_cutoff=0.01, ca_cutoff=0.005, pdbid='pdbid', writeto=sys.stdout):
  outfile_kin = open(pdbid+'_cablam_multi.kin','w')
  outfile_pdb = open(pdbid+'_cablam_multi.pdb','w')
  resdata=setup(hierarchy,pdbid)
  peptide_expectations = fetch_peptide_expectations()
  ca_expectations = fetch_ca_expectations()
  motif_contours = fetch_motifs()

  peptide_outliers=find_peptide_outliers(resdata,peptide_expectations,cutoff=peptide_cutoff)
  helix_or_sheet(peptide_outliers,motif_contours)
  give_kin(peptide_outliers, peptide_cutoff, color='purple', writeto=outfile_kin)
  give_kin(peptide_outliers, peptide_bad_cutoff, color='magenta', writeto=outfile_kin)

  #peptide_bad_outliers=find_peptide_outliers(resdata,peptide_expectations,cutoff=peptide_bad_cutoff)
  #helix_or_sheet(peptide_bad_outliers,motif_contours)
  #give_kin(peptide_bad_outliers, peptide_bad_cutoff, color='magenta', writeto=outfile_kin)

  ca_outliers = find_ca_outliers(resdata,ca_expectations,cutoff=ca_cutoff)
  helix_or_sheet(ca_outliers,motif_contours) #This call for dev purposes, need to know what's flagged and interesting
  give_ca_kin(ca_outliers, ca_cutoff, writeto=outfile_kin)

  find_partial_sec_struc(resdata,ca_outliers=ca_outliers)
  motifs = find_whole_sec_struc(resdata)
  for motif in motifs:
    motif.print_record(writeto=outfile_pdb)
  outfile_pdb.write(hierarchy.as_pdb_string())
  outfile_kin.close()
  outfile_pdb.close()
  king_command = 'phenix.king -m ' + pdbid+'_cablam_multi.pdb ' + pdbid+'_cablam_multi.kin'
  #'king -m' merges the files that follow
  easy_run.fully_buffered(king_command)
#-------------------------------------------------------------------------------
#}}}

#{{{ output functions
#-------------------------------------------------------------------------------
def give_kin(outliers, outlier_cutoff, color='purple', writeto=sys.stdout):
  #The kinemage markup is purple dihedrals following outlier CO angles
  #This angle was chosen because it is indicative of the most likely source of
  #  problems
  #Purple was chosen because it is easily distinguished from (green) Rama
  #  markup, and because the cablam measures share some philosophical
  #  similarity with perp distance in RNA markup
  #writeto.write('\n@kinemage')
  writeto.write('\n@subgroup {cablam out '+str(outlier_cutoff)+'} dominant\n')
  writeto.write('@vectorlist {cablam outliers} color= '+color+' width= 4 master={cablam out '+str(outlier_cutoff)+'}')
  #writeto.write('@vectorlist')
  reskeys = outliers.keys()
  reskeys.sort()
  for resid in reskeys:
    outlier = outliers[resid]
    if outlier.outlier_level >= outlier_cutoff:
      continue
    #Shouldn't have to do checking here, since only residues with calculable values should get to this point
    residue = outlier.residue
    prevres = residue.prevres
    nextres = residue.nextres
    CA_1, O_1 = prevres.getatomxyz('CA'),prevres.getatomxyz('O')
    CA_2, O_2 = residue.getatomxyz('CA'),residue.getatomxyz('O')
    CA_3      = nextres.getatomxyz('CA')

    pseudoC_1 = cablam_math.perptersect(CA_1,CA_2,O_1)
    pseudoC_2 = cablam_math.perptersect(CA_2,CA_3,O_2)

    midpoint = [ (pseudoC_1[0]+pseudoC_2[0])/2.0 , (pseudoC_1[1]+pseudoC_2[1])/2.0 , (pseudoC_1[2]+pseudoC_2[2])/2.0 ]

    stats = ' '.join(['%.5f' %outlier.outlier_level, 'alpha'+'%.5f' %outlier.loose_alpha, 'beta'+'%.5f' %outlier.regular_beta])

    writeto.write('\n{'+stats+'} P '+ str(O_1[0]) +' '+ str(O_1[1]) +' '+ str(O_1[2]))
    writeto.write('\n{'+stats+'} '+ str(pseudoC_1[0]) +' '+ str(pseudoC_1[1]) +' '+ str(pseudoC_1[2]))
    writeto.write('\n{'+stats+'} '+ str(midpoint[0]) +' '+ str(midpoint[1]) +' '+ str(midpoint[2]))
    writeto.write('\n{'+stats+'} '+ str(pseudoC_2[0]) +' '+ str(pseudoC_2[1]) +' '+ str(pseudoC_2[2]))
    writeto.write('\n{'+stats+'} '+ str(O_2[0]) +' '+ str(O_2[1]) +' '+ str(O_2[2]))
  writeto.write('\n')

def give_ca_kin(outliers, outlier_cutoff, writeto=sys.stdout):
  #Kinemage markup as red CA virtual angles, since that is the measure unique to
  #this validation.
  writeto.write('\n@subgroup {ca out '+str(outlier_cutoff)+'} dominant\n')
  writeto.write('@vectorlist {ca outliers} color= red width= 4 master={ca geom outliers}')
  reskeys = outliers.keys()
  reskeys.sort()
  for resid in reskeys:
    outlier = outliers[resid]
    residue = outlier.residue
    prevres = residue.prevres
    nextres = residue.nextres
    CA_1 = prevres.getatomxyz('CA')
    CA_2 = residue.getatomxyz('CA')
    CA_3 = nextres.getatomxyz('CA')

    #Note: 'stats' should be simplified after I'm certain that the cutoffs are in the right places
    stats = ' '.join(['%.5f' %outlier.outlier_level, 'alpha'+'%.5f' %outlier.loose_alpha, 'beta'+'%.5f' %outlier.regular_beta])

    #CA_2[0]+(CA_2[0]-CA_1[0])*0.9

    writeto.write('\n{'+stats+'} P '+str(CA_2[0]-(CA_2[0]-CA_1[0])*0.9)+' '+str(CA_2[1]-(CA_2[1]-CA_1[1])*0.9)+' '+str(CA_2[2]-(CA_2[2]-CA_1[2])*0.9))
    #writeto.write('\n{'+stats+'} P '+str(CA_1[0])+' '+str(CA_1[1])+' '+str(CA_1[2]))
    writeto.write('\n{'+stats+'} '+str(CA_2[0])+' '+str(CA_2[1])+' '+str(CA_2[2]))
    writeto.write('\n{'+stats+'} '+str(CA_2[0]-(CA_2[0]-CA_3[0])*0.9)+' '+str(CA_2[1]-(CA_2[1]-CA_3[1])*0.9)+' '+str(CA_2[2]-(CA_2[2]-CA_3[2])*0.9))
    #writeto.write('\n{'+stats+'} '+str(CA_3[0])+' '+str(CA_3[1])+' '+str(CA_3[2]))
  writeto.write('\n')


def give_points(outliers, outlier_cutoff=0.05, writeto=sys.stdout):
  #This prints a dotlist in kinemage format
  #Each dot is one outlier residue in 3D space for visual comparison to expected
  #  behavior contours
  writeto.write('\n@kinemage\n')
  writeto.write('@group {cablam outliers} dominant\n')
  writeto.write('@dotlist {cablam outliers}')
  reskeys = outliers.keys()
  reskeys.sort()
  for resid in reskeys:
    outlier = outliers[resid]
    if outlier.outlier_level >= outlier_cutoff:
      continue
    residue = outlier.residue
    if 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures and 'CO_d_in' in residue.measures:
      cablam_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CO_d_in']]
      writeto.write('{'+residue.id_with_resname()+'} '+'%.3f' %cablam_point[0]+' '+'%.3f' %cablam_point[1]+' '+'%.3f' %cablam_point[2]+'\n')
  writeto.write('\n')

def give_text(outliers, writeto=sys.stdout):
  #This prints a comma-separated line of data for each outlier residue
  #Intended for easy machine readability
  #writeto.write('\nresidue,contour_level,loose_alpha,regular_alpha,loose_beta,regular_beta,threeten')
  writeto.write('residue : outlier_type : contour_level : ca_contour_level : sec struc recommendation : alpha score : beta score : three-ten score')
  reskeys = outliers.keys()
  reskeys.sort()
  for resid in reskeys:
    outlier = outliers[resid]
    if outlier.ca_geom_outlier_level < 0.005:
      outlier_type = ' CA Geom Outlier    '
    elif outlier.outlier_level < 0.01:
      outlier_type = ' CaBLAM Outlier     '
    elif outlier.outlier_level < 0.05:
      outlier_type = ' CaBLAM Disfavored  '
    else:
      outlier_type =   '                    '
    #outlist = [outlier.residue.mp_id(), '%.5f' %outlier.outlier_level, '%.5f' %outlier.loose_alpha, '%.5f' %outlier.regular_alpha, '%.5f' %outlier.loose_beta, '%.5f' %outlier.regular_beta, '%.5f' %outlier.loose_threeten]
    #try:
    outlist = [outlier.residue.mp_id(), outlier_type, '%.5f' %outlier.outlier_level, '%.5f' %outlier.ca_geom_outlier_level, outlier.suggestion, '%.5f' %outlier.loose_alpha, '%.5f' %outlier.regular_beta, '%.5f' %outlier.loose_threeten]
    #except TypeError:
    #  print '\n',outlier.residue.mp_id(), outlier_type, outlier.outlier_level, outlier.ca_geom_outlier_level, outlier.loose_alpha, outlier.regular_beta, outlier.loose_threeten
    #  sys.exit()
    writeto.write('\n'+':'.join(outlist))
  writeto.write('\n')
#-------------------------------------------------------------------------------
#}}}

#{{{ oneline function
#-------------------------------------------------------------------------------
#Outputs cablam validation information in format similar to molprobity oneline.
#  Contents of this output subject to change: percents of cablam-relevant
#  outliers and annotation.
def oneline_header(writeto=sys.stdout):
  #Write column labels for oneline output
  writeto.write('pdbid:residues:cablam_disfavored_percent:cablam_outlier_percent:ca_geom_outlier_percent\n')

def oneline(hierarchy, peptide_cutoff=0.05, peptide_bad_cutoff=0.01, ca_cutoff=0.005, pdbid='pdbid', writeto=sys.stdout):
  resdata=setup(hierarchy,pdbid)
  peptide_expectations = fetch_peptide_expectations()
  ca_expectations = fetch_ca_expectations()

  residue_count = 0
  peptide_outlier_count = 0
  peptide_bad_outlier_count = 0
  ca_outlier_count = 0
  reskeys = resdata.keys()
  for resid in reskeys:
    residue = resdata[resid]
    if 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures and 'CO_d_in' in residue.measures and 'CA_a' in residue.measures:
      #This check ensures that only calculable, protein residues are considered
      residue_count += 1
      resname = residue.alts[residue.firstalt('CA')]['resname']
      peptide_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CO_d_in']]
      ca_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CA_a']]
      #there are separate contours for gly, for pro, and for other residues
      if resname.upper() == 'GLY':
        peptide_percentile = peptide_expectations['gly'].valueAt(peptide_point)
        ca_percentile = ca_expectations['gly'].valueAt(ca_point)
      elif resname.upper() == 'PRO':
        if ('omega' not in residue.measures) or (residue.measures['omega'] >= 90) or (residue.measures['omega'] <= -90):
          peptide_percentile = peptide_expectations['transpro'].valueAt(peptide_point)
          ca_percentile = ca_expectations['transpro'].valueAt(ca_point)
        else:
          peptide_percentile = peptide_expectations['cispro'].valueAt(peptide_point)
          ca_percentile = ca_expectations['cispro'].valueAt(ca_point)
      #elif resname.upper() == 'PRO':
      #  peptide_percentile = peptide_expectations['pro'].valueAt(peptide_point)
      #  ca_percentile = ca_expectations['pro'].valueAt(ca_point)
      else:
        peptide_percentile = peptide_expectations['general'].valueAt(peptide_point)
        ca_percentile = ca_expectations['general'].valueAt(ca_point)

      if peptide_percentile < peptide_cutoff:
        peptide_outlier_count += 1
      if peptide_percentile < peptide_bad_cutoff:
        peptide_bad_outlier_count += 1
      if ca_percentile < ca_cutoff:
        ca_outlier_count += 1

  try:
    peptide_outlier_percent = peptide_outlier_count/residue_count*100
  except ZeroDivisionError:
    peptide_outlier_percent = 0
  try:
    peptide_bad_outlier_percent = peptide_bad_outlier_count/residue_count*100
  except ZeroDivisionError:
    peptide_bad_outlier_percent = 0
  try:
    ca_outlier_percent = ca_outlier_count/residue_count*100
  except ZeroDivisionError:
    ca_outlier_percent = 0
  #sys.stderr.write(pdbid+'\n') #Doesn't need bc pdbid printed in output
  writeto.write(pdbid.lower()+':'+str(residue_count)+':'+'%.1f'%peptide_outlier_percent+':'+'%.1f'%peptide_bad_outlier_percent+':'+'%.2f'%ca_outlier_percent+'\n')
#-------------------------------------------------------------------------------
#}}}

#{{{ analyze_pdb
#-------------------------------------------------------------------------------
#Returns a list of cablam_validation objects, one for each residue identified as
#  a potential outliers.  This list is the object of choice for internal access
#  to cablam output.
#The default outlier_cutoff = 0.05 catches most outliers of interest, although
#  it may also flag loops and other low-population motifs as potential outliers.
#  Setting outlier_cutoff=1.0 will return a validation object for every residue.
def analyze_pdb(hierarchy, outlier_cutoff=0.05, pdbid='pdbid'):
  resdata=setup(hierarchy,pdbid)
  peptide_expectations = fetch_peptide_expectations()
  ca_geom_expectations = fetch_ca_expectations()
  motifs = fetch_motifs()

  residue_stats = find_all_residue_stats(resdata, peptide_expectations, ca_geom_expectations)
  helix_or_sheet(residue_stats, motifs)
  assign_sec_struc(residue_stats)

  return residue_stats

#Alternative call to find ca-trace outliers
def analyze_pdb_ca(hierarchy, outlier_cutoff=0.005, pdbid='pdbid'):
  resdata=setup(hierarchy,pdbid)
  ca_expectations = fetch_ca_expectations()
  motifs = fetch_motifs()

  outliers = find_ca_outliers(resdata, ca_expectations, cutoff=outlier_cutoff)
  helix_or_sheet(outliers, motifs)

  return outliers

def analyze_pdb_gui(hierarchy, pdbid='pdbid'):
  outlier_cutoff = 0.05
  ca_outlier_cutoff = 0.005
  resdata=setup(hierarchy,pdbid)
  peptide_expectations = fetch_peptide_expectations()
  ca_expectations = fetch_ca_expectations()
  motifs = fetch_motifs()

  outliers = find_peptide_outliers(resdata, peptide_expectations, cutoff=outlier_cutoff)
  ca_outliers = find_ca_outliers(resdata, ca_expectations, cutoff=ca_outlier_cutoff)
  helix_or_sheet(outliers,motifs)
  find_partial_sec_struc(resdata,ca_outliers=ca_outliers)
  for outlier in outliers.values():
    if outlier.residue.motif_guess.regular_beta:
      outlier.suggestion = 'try Beta Strand'
    elif outlier.residue.motif_guess.loose_threeten:
      outlier.suggestion = 'try ThreeTen'
    elif outlier.residue.motif_guess.loose_alpha:
      outlier.suggestion = 'try Alpha Helix'
    else:
      pass #suggestion = ''
  return outliers
#-------------------------------------------------------------------------------
#}}}

#{{{ run
#-------------------------------------------------------------------------------
#Run is intended largely for commandline access.  Its default output is comma-
#  separated text to stdout.
def run(args):
  #{{{ phil parsing
  #-----------------------------------------------------------------------------
  interpreter = libtbx.phil.command_line.argument_interpreter(master_phil=master_phil)
  sources = []
  for arg in args:
    if os.path.isfile(arg): #Handles loose filenames
      input_file = file_reader.any_file(arg)
      if (input_file.file_type == "pdb"):
        sources.append(interpreter.process(arg="pdb_infile=\"%s\"" % arg))
      elif (input_file.file_type == "phil"):
        sources.append(input_file.file_object)
    elif os.path.isdir(arg):
      sources.append(interpreter.process(arg="pdb_infile=\"%s\"" % arg))
    else: #Handles arguments with xxx=yyy formatting
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  params = work_params.cablam_validate
  #-----------------------------------------------------------------------------
  #}}} end phil parsing

  if params.help:
    usage()
    interpretation()
    sys.exit()

  if not params.pdb_infile:
    sys.stdout.write(
      '\nMissing input data, please provide .pdb file\n')
    usage()
    sys.exit()
  if os.path.isdir(params.pdb_infile):
    fileset = []
    dirpath = params.pdb_infile
    for filename in os.listdir(params.pdb_infile):
      fileset.append(os.path.join(dirpath,filename))
  elif os.path.isfile(params.pdb_infile):
    fileset = [params.pdb_infile]
  if params.output=='oneline':
    oneline_header()
  for pdb_infile in fileset:
    if not os.path.isfile(pdb_infile): continue
    pdb_in = file_reader.any_file(pdb_infile)
    if pdb_in.file_type != "pdb":
      sys.stderr.write(pdb_infile +" not id'd as readable file\n")
      continue
    pdbid = os.path.splitext(os.path.basename(pdb_infile))[0]
    pdb_io = pdb.input(pdb_infile)
    hierarchy = pdb_io.construct_hierarchy()

    if params.output=='oneline':
      oneline(hierarchy, pdbid=pdbid)
    elif params.output=='ca_kin':
      ca_outliers = analyze_pdb_ca(
        hierarchy, outlier_cutoff=params.outlier_cutoff, pdbid=pdbid)
      give_ca_kin(ca_outliers,params.outlier_cutoff)
    elif params.output=='full_kin':
      cablam_multicrit_kin(hierarchy,peptide_cutoff=params.outlier_cutoff,pdbid=pdbid)
    elif params.output=='markup_no_ribbons':
      print_cablam_markup_kin(hierarchy, pdbid=pdbid)
    elif params.check_prolines:
      check_prolines(hierarchy,pdbid=pdbid)
    elif params.output=='records':
      print_helix_sheet_records(hierarchy,ca_cutoff=0.005,pdbid='pdbid',writeto=sys.stdout)
    elif params.output=='records_plus_pdb':
      print_pdb_with_new_helix_sheet_records(hierarchy,ca_cutoff=0.005,pdbid='pdbid',writeto=sys.stdout)
    else:
      outliers = analyze_pdb(
        hierarchy, outlier_cutoff=params.outlier_cutoff, pdbid=pdbid)

      if params.output=='kin':
        give_kin(outliers,params.outlier_cutoff)
      if params.output=='points':
        give_points(outliers,params.outlier_cutoff)
      if params.output=='text':
        give_text(outliers)
#-------------------------------------------------------------------------------
#}}}

#class beta_strand(object):
#  def __init__(self):
#    self.members = []
#    self.mates = []
#
#def get_strands(motifs):
#  strands = []
#  for motif in motifs:
#    if motif.motif_type == 'sheet':
#      strand = beta_strand()
#      start = motif.motif_start
#      end = motif.motif_end
#      curres = start
#      strand.members.append(start)
#      while curres is not end:
#        curres = curres.nextres
#        strand.members.append(curres)
#      strands.append(strand)
#  return strands
#
#def pair_beta(strands):
#  #CA-CA interstrand distance for antiparallel alternates between ~4 and ~6
#  #for parallel, distances are consistantly just under 5
#  dist_cutoff = 6.5
#  index1 = 0
#  while index1 < len(strands):
#    strand1 = strands[index1]
#    index2 = index1 + 1
#    while index2 < len(strands):
#      strand2 = strands[index2]
#      for residue1 in strand1.members:
#        res1ca = residue1.getatomxyz('CA')
#        for residue2 in strand2.members:
#          res2ca = residue2.getatomxyz('CA')
#          resdist = cablam_math.veclen(cablam_math.vectorize(res1ca,res2ca))
#          if resdist <= dist_cutoff:
#            strand1.mates.append(strand2)
#            strand2.mates.append(strand1)
#      index2 += 1
#    index1 += 1
#
#
#
##The following is a function in progress to connect beta strands into sheets.
##Do not use without further development
#def stitch_beta(motifs): #pass in an iterable of motif_chunk class instances
#  strands = []
#  for motif in motifs:
#    if motif.motif_type == 'sheet':
#      strands.append(motif)
#    else: pass
#  for strand1 in strands:
#    for strand2 in strands:
#      strand1_dir = cablam_math.vectorize(strand1.motif_start.getatomxyz('CA'),strand1.motif_end.getatomxyz('CA'))
#      strand2_dir = cablam_math.vectorize(strand2.motif_start.getatomxyz('CA'),strand2.motif_end.getatomxyz('CA'))
#      strand_dot_product = cablam_math.dot(strand1_dir,strand2_dir)
#      if strand_dot_product < 0: #antiparallel
#        curres = strand1.motif_start
#        pass
#      elif strand_dot_product >0: #parallel
#        pass
#      pass
#      #check if parallel or antiparallel with dot product
#      #for residue in strand 1:
#      #for residue in strand 2:
#        #if close enough:
#        #walk in direction
#
#
#  return sheets
#
##Find potential mates based on distance
##For each mate, see if there's a good parrallel/ antiparallel local direction (3 residues?)
##Based on par/anti-par, find best alignment of residues minimizing res-res distance
##Check pleat for this alignment
##
##Define units as peptides, not residues
#class beta_strand(object):
#  def __init__(self):
#    self.residues
#    self.peptides = []
#    self.mates = []
#
#class peptide(object):
#  def check_for_mate(self, other_peptide, dist_cutoff = 6.5):
#    dist11 = find_distance(self.res1_CA,other_peptide.res1_CA)
#    dist22 = find_distance(self.res2_CA,other_peptide.res2_CA)
#    dist12 = find_distance(self.res1_CA,other_peptide.res2_CA)
#    dist21 = find_distance(self.res2_CA,other_peptide.res1_CA)
#    if dist11 <= dist_cutoff and dist22 <= dist_cutoff:
#      self.mates.append(peptide_mate(other_peptide))
#    elif dist12 <= dist_cutoff and dist21 <= dist_cutoff:
#      self.mates.append(peptide_mate(other_peptide))
#    else:
#      pass
#
#  def __init__(self, residue):
#    self.dist_cutoff = 6.5
#    self.res1 = residue
#    self.res1_CA = residue.getatomxyz("CA")
#    self.res2 = residue.nextres
#    self.res2_CA = res2.getatomxyz("CA")
#    self.mates = []
#    self.nextpep = None
#    self.prevpep = None
#
#class peptide_mate(object):
#  def __init__(self, peptide):
#    self.peptide = peptide
#
#def find_distance(point1,point2):
#  return cablam_math.veclen(cablam_math.vectorize(point1,point2))
#
#def get_strands(motifs):
#  strands = []
#  for motif in motifs:
#    if motif.motif_type == 'sheet':
#      strand = beta_strand()
#      start = motif.motif_start
#      end = motif.motif_end
#      curres = start
#      strand.members.append(start)
#      prevpep = None
#      while curres is not end:
#        peptide = peptide(curres)
#        if prevpep:
#          peptide.prevpep = prevpep
#          prevpep.nextpep = peptide
#        strand.peptides.append(peptide)
#        curres = curres.nextres
#        prevpep = peptide
#        strand.residues.append(curres)
#      strands.append(strand)
#  return strands
