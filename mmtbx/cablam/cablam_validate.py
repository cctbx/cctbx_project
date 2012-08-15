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
  give_hierarchy = False
    .type = bool
    .help = '''return a list of dictionaries, each dict has a hierarchy rg and validation on it'''
  outlier_cutoff = 0.05
    .type = float
    .help = '''sets the contour level for detecting outliers, e.g. outlier_cutoff=0.01 for accepting the top 99%, defaults to 0.05'''
}
""", process_includes=True)

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
  give_points=False        prints cablam-space points to screen,
                             in kinemage dotlist format
  give_hierarchy=False     returns a list of validations objects
                             each object contains validation data for one
                             residue, plus the hierarchy 'rg' object for that
                             residue
                           Intended for developer use, this will only work if
                             cablam_validate is given a hierarchy object to work
                             on (not possible from commandline)
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
  cablam_math.COpseudodihedrals(resdata)
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
  outliers = []
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
        outliers.append(cablam_validation(residue=residue,outlier_level=percentile))
  return outliers
#}}}

#{{{ helix_or_sheet (wrapper) and helix_or_sheet_res
#targets is a list of reskeys to be looked at, targets=resdata.keys to check all
def helix_or_sheet(outliers, motifs):
  for outlier in outliers:
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
def give_hierarchy(outliers):
  #Yes, this is a bit trivial at the moment. It exists for symmetry and in case
  #  I need custom treatment for this in the future
  return outliers

def give_kin(outliers, writeto=sys.stdout):
  #The kinemage markup is purple dihedrals following outlier CO angles
  #This angle was chosen because it is indicative of the most likely source of
  #  problems
  #Purple was chosen because it is easily distinguished from (green) Rama
  #  markup, and because the cablam measures share some philosophical
  #  similarity with perp distance in RNA markup
  writeto.write('\n@kinemage\n')
  writeto.write('@group {cablam outliers} dominant\n')
  writeto.write('@vectorlist {cablam outliers} color= purple width= 4 \n')
  for outlier in outliers:
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

def give_points(outliers, writeto=sys.stdout):
  #This prints a dotlist in kinemage format
  #Each dot is one outlier residue in 3D space for visual comparison to expected
  #  behavior contours
  writeto.write('\n@kinemage\n')
  writeto.write('@group {cablam outliers} dominant\n')
  writeto.write('@dotlist {cablam outliers}')
  for outlier in outliers:
    residue = outlier.residue
    if 'CA_d_in' in residue.measures and 'CA_d_out' in residue.measures and 'CO_d_in' in residue.measures:
      cablam_point = [residue.measures['CA_d_in'],residue.measures['CA_d_out'],residue.measures['CO_d_in']]
      writeto.write('{'+residue.id_with_resname(sep=' ')+'} '+'%.3f' %cablam_point[0]+' '+'%.3f' %cablam_point[1]+' '+'%.3f' %cablam_point[2]+'\n')

def give_text(outliers, writeto=sys.stdout):
  #This prints a comma-separated line of data for each outlier residue
  #Intended for easy machine readability
  writeto.write('\nresidue,contour_level,loose_alpha,regular_alpha,loose_beta,regular_beta')
  for outlier in outliers:
    outlist = [outlier.residue.id_with_resname(), '%.5f' %outlier.outlier_level, '%.5f' %outlier.loose_alpha, '%.5f' %outlier.regular_alpha, '%.5f' %outlier.loose_beta, '%.5f' %outlier.regular_beta]
    writeto.write('\n'+','.join(outlist))
#-------------------------------------------------------------------------------
#}}}

#{{{ run
#-------------------------------------------------------------------------------
def run(args, hierarchy=None):
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
  if not work_params.cablam_validate.pdb_infile:
    usage()
    sys.exit()
    #sys.stderr.write('\nNo input file specified or could not find input file\n')
  params = work_params.cablam_validate
  #-----------------------------------------------------------------------------
  #}}} end phil parsing

  pdbid = 'pdbid'
  if not hierarchy:
    if params.pdb_infile:
      pdbid = os.path.basename(params.pdb_infile)
      pdb_io = pdb.input(params.pdb_infile)
      hierarchy = pdb_io.construct_hierarchy()
    else:
      sys.stdout.write(
        '\nMissing input data, please provide .pdb file or hierarchy\n')
      sys.exit()

  resdata=setup(hierarchy,pdbid)
  expectations = fetch_expectations()
  motifs = fetch_motifs()

  outliers = find_outliers(resdata, expectations, cutoff=params.outlier_cutoff)
  helix_or_sheet(outliers, motifs)

  if params.give_kin:
    give_kin(outliers)
  if params.give_points:
    give_points(outliers)
  if params.give_text:
    give_text(outliers)
  if params.give_hierarchy:
    return give_hierarchy(outliers)
#-------------------------------------------------------------------------------
#}}}

#{{{ "__main__"
#-------------------------------------------------------------------------------
if __name__ == "__main__":
  run(sys.argv[1:])
#-------------------------------------------------------------------------------
#}}}
