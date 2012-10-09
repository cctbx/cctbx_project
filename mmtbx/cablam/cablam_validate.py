from __future__ import division
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
import libtbx.load_env

#{{{ phil
#-------------------------------------------------------------------------------
master_phil = libtbx.phil.parse("""
cablam_validate {
  pdb_infile = None
    .type = path
    .help = '''input PDB file'''
  give_kin = False
    .type = bool
    .help = '''print markup kinemage as output to automatically-named file'''
  give_text = False
    .type = bool
    .help = '''print text output to sys.stdout'''
  give_points = False
    .type = bool
    .help = '''print cablam-space points to sys.stdout, in kinemage dotlist format'''
  outlier_cutoff = 0.05
    .type = float
    .help = '''sets the contour level for detecting outliers, e.g. outlier_cutoff=0.01 for accepting the top 99%, defaults to 0.05'''
  help = False
    .type = bool
    .help = '''help and data interpretation messages'''
}
""", process_includes=True)
#-------------------------------------------------------------------------------
#}}}

#{{{ cablam_validation class
#-------------------------------------------------------------------------------
#This object holds information on one residue for purposes of cablam_validate
#It's mostly just a package for data passing and access
class cablam_validation():
  def __init__(self, residue=None,outlier_level=None):
    self.residue = residue #cablam_res.linked_residue object
    self.outlier_level = outlier_level #percentile level in "expectations" contours
    self.loose_alpha  = None #percentile level in motif contour
    self.regular_alpha= None #percentile level in motif contour
    self.loose_beta   = None #percentile level in motif contour
    self.regular_beta = None #percentile level in motif contour
    if self.residue:
      self.rg = self.residue.rg #rg object from a source hierarchy
    else:
      self.rg = None
#-------------------------------------------------------------------------------
#}}}

#{{{ usage
#-------------------------------------------------------------------------------
def usage():
  sys.stderr.write("""
phenix.cablam_validate file.pdb [options ...]

Options:

  pdb_infile=filename      input PDB file
  outlier_cutoff=0.05      sets the contour level for detecting outliers,
                             e.g. outlier_cutoff=0.01 for accepting the top 99%,
                             defaults to 0.05
  give_kin=False           prints markup kinemage to screen
  give_text=False          prints machine-readable columnated data to screen
                           (default output)
  give_points=False        prints cablam-space points to screen,
                             in kinemage dotlist format
  help=False               prints this usage text, plus notes on data
                             interpretation to screen

Example:

phenix.cablam_validate file.pdb give_text=True
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
  resdata=cablam_res.construct_linked_residues(hierarchy,targetatoms=['CA','C','N','O'],pdbid=pdbid)
  cablam_res.prunerestype(resdata, 'HOH')
  cablam_math.cablam_measures(resdata)
  return resdata
#-------------------------------------------------------------------------------
#}}}

#{{{ fetch_expectations
#-------------------------------------------------------------------------------
#This function finds, unpickles, and returns N-Dim Tables of expected residue
#  behavior for use in determining outliers
#The return object is a dict keyed by residue type: 'general','gly','pro'
def fetch_expectations():
  categories = ['general','gly','pro']
  unpickled = {}
  for category in categories:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=("chem_data/cablam_data/cablam.8000.expected."+category+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for category "+category+" in chem_data.\nExiting.\n")
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
  motifs = ['loose_beta','regular_beta','loose_alpha','regular_alpha']
  unpickled = {}
  for motif in motifs:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=("chem_data/cablam_data/cablam.8000.motif."+motif+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for motif "+motif+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[motif] = ndt
  return unpickled
#-------------------------------------------------------------------------------
#}}}

#{{{ find_outliers
#
#These are used by the helix_or_sheet function
def find_outliers(resdata,expectations,cutoff=0.05):
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
        percentile = expectations['pro'].valueAt(cablam_point)
      else:
        percentile = expectations['general'].valueAt(cablam_point)

      if percentile < cutoff:
        #sys.stderr.write(str(percentile)+'\n') #debug
        outliers[resid] = cablam_validation(residue=residue,outlier_level=percentile)
  return outliers
#}}}

#{{{ helix_or_sheet (wrapper) and helix_or_sheet_res
#targets is a list of reskeys to be looked at, targets=resdata.keys to check all
def helix_or_sheet(outliers, motifs):
  for outlier in outliers.values():
    residue = outlier.residue
    contours = helix_or_sheet_res(residue, motifs)
    if contours:
      outlier.loose_alpha   = contours['loose_alpha']
      outlier.regular_alpha = contours['regular_alpha']
      outlier.loose_beta    = contours['loose_beta']
      outlier.regular_beta  = contours['regular_beta']

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
#}}}

#{{{ output functions
#-------------------------------------------------------------------------------
def give_kin(outliers, outlier_cutoff, writeto=sys.stdout):
  #The kinemage markup is purple dihedrals following outlier CO angles
  #This angle was chosen because it is indicative of the most likely source of
  #  problems
  #Purple was chosen because it is easily distinguished from (green) Rama
  #  markup, and because the cablam measures share some philosophical
  #  similarity with perp distance in RNA markup
  #writeto.write('\n@kinemage')
  writeto.write('\n@group {cablam out '+str(outlier_cutoff)+'} dominant\n')
  writeto.write('@vectorlist {cablam outliers} color= purple width= 4 \n')
  reskeys = outliers.keys()
  reskeys.sort()
  for resid in reskeys:
    outlier = outliers[resid]
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

    stats = ' '.join(['%.5f' %outlier.outlier_level, 'a'+'%.5f' %outlier.loose_alpha, 'rega'+'%.5f' %outlier.regular_alpha, 'b'+'%.5f' %outlier.loose_beta, 'regb'+'%.5f' %outlier.regular_beta])

    writeto.write('\n{'+stats+'} P '+ str(O_1[0]) +' '+ str(O_1[1]) +' '+ str(O_1[2]))
    writeto.write('\n{'+stats+'} '+ str(pseudoC_1[0]) +' '+ str(pseudoC_1[1]) +' '+ str(pseudoC_1[2]))
    writeto.write('\n{'+stats+'} '+ str(midpoint[0]) +' '+ str(midpoint[1]) +' '+ str(midpoint[2]))
    writeto.write('\n{'+stats+'} '+ str(pseudoC_2[0]) +' '+ str(pseudoC_2[1]) +' '+ str(pseudoC_2[2]))
    writeto.write('\n{'+stats+'} '+ str(O_2[0]) +' '+ str(O_2[1]) +' '+ str(O_2[2]))
  writeto.write('\n')

def give_points(outliers, writeto=sys.stdout):
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
    residue = outlier.residue
    if 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures and 'CO_d_in' in residue.measures:
      cablam_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CO_d_in']]
      writeto.write('{'+residue.id_with_resname()+'} '+'%.3f' %cablam_point[0]+' '+'%.3f' %cablam_point[1]+' '+'%.3f' %cablam_point[2]+'\n')
  writeto.write('\n')

def give_text(outliers, writeto=sys.stdout):
  #This prints a comma-separated line of data for each outlier residue
  #Intended for easy machine readability
  writeto.write('\nresidue,contour_level,loose_alpha,regular_alpha,loose_beta,regular_beta')
  reskeys = outliers.keys()
  reskeys.sort()
  for resid in reskeys:
    outlier = outliers[resid]
    outlist = [outlier.residue.id_with_resname(), '%.5f' %outlier.outlier_level, '%.5f' %outlier.loose_alpha, '%.5f' %outlier.regular_alpha, '%.5f' %outlier.loose_beta, '%.5f' %outlier.regular_beta]
    writeto.write('\n'+','.join(outlist))
  writeto.write('\n')
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
  expectations = fetch_expectations()
  motifs = fetch_motifs()

  outliers = find_outliers(resdata, expectations, cutoff=outlier_cutoff)
  helix_or_sheet(outliers, motifs)

  return outliers
#-------------------------------------------------------------------------------
#}}}

#{{{ run
#-------------------------------------------------------------------------------
#Run is intended largely for commandline access.  Its default output is comma-
#  sepparated text to stdout.
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
    else: #Handles arguments with xxx=yyy formatting
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  params = work_params.cablam_validate
  #if not work_params.cablam_validate.pdb_infile:
  #  usage()
  #  sys.exit()
  params = work_params.cablam_validate
  #-----------------------------------------------------------------------------
  #}}} end phil parsing

  if params.help:
    usage()
    interpretation()
    sys.exit()

  pdbid = 'pdbid'
  if params.pdb_infile:
    pdbid = os.path.basename(params.pdb_infile)
    pdb_io = pdb.input(params.pdb_infile)
    hierarchy = pdb_io.construct_hierarchy()
  else:
    sys.stdout.write(
      '\nMissing input data, please provide .pdb file\n')
    usage()
    sys.exit()

  outliers = analyze_pdb(
    hierarchy, outlier_cutoff=params.outlier_cutoff, pdbid=pdbid)

  if not (params.give_kin or params.give_points):
    #Set default output as text
    params.give_text = True

  if params.give_kin:
    give_kin(outliers,params.outlier_cutoff)
  if params.give_points:
    give_points(outliers)
  if params.give_text:
    give_text(outliers)
#-------------------------------------------------------------------------------
#}}}

#{{{ "__main__"
#-------------------------------------------------------------------------------
if __name__ == "__main__":
  run(sys.argv[1:])
#-------------------------------------------------------------------------------
#}}}
