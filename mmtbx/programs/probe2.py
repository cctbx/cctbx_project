"""Run molprobity. version 2"""
##################################################################################
# Copyright(c) 2021-2023, Richardson Lab at Duke
# Licensed under the Apache 2 license
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, division, print_function
import sys
import math
from datetime import datetime
from pathlib import Path
from libtbx.program_template import ProgramTemplate
from libtbx import group_args, phil
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry
import mmtbx
import mmtbx_probe_ext as probeExt
from mmtbx.probe import Helpers
from iotbx import pdb
from iotbx.pdb import common_residue_names_get_class

version = "4.11.0"

master_phil_str = '''
profile = False
  .type = bool
  .short_caption = Profile the run
  .help = Profile the performance of the entire run

source_selection = "occupancy > 0.33"
  .type = atom_selection
  .short_caption = Source selection
  .help = Source selection description

target_selection = None
  .type = atom_selection
  .short_caption = Target selection
  .help = Target selection description ('=' means same as source)

use_neutron_distances = False
  .type = bool
  .short_caption = Use neutron distances
  .help = Use neutron distances (-nuclear in probe)

approach = *self both once surface count_atoms
  .type = choice
  .short_caption = What to count
  .help = self (src -> src) both (src <=> targ) once (src -> targ) surface (VdW surface) count_atoms (count atoms)

excluded_bond_chain_length = 4
  .type = int
  .short_caption = Excluded bond chain length
  .help = Exclude chain of atoms bonded to source for this many hops (-4H, -3, -2 , -1 in probe).  When set to 4, an atom chain longer than 3 is only excluded when either the first or the last atom in the chain is a Hydrogen.

minimum_water_hydrogen_occupancy = 0.25
  .type = float
  .short_caption = Minimum water hydrogen occupancy
  .help = Minimum occupancy for polar hydrogens (0.66 in original Reduce)

maximum_water_hydrogen_b = 80.0
  .type = float
  .short_caption = Maximum water hydrogen B factor
  .help = Minimum b-factor for polar hydrogens (40.0 in original Reduce)

include_mainchain_mainchain = True
  .type = bool
  .short_caption = Include mainchain-mainchain interactions
  .help = Include mainchain -> mainchain interactions (-mc in probe)

include_water_water = False
  .type = bool
  .short_caption = Include water-water interactions
  .help = Include water-to-water interactions (-wat2wat in probe)

keep_unselected_atoms = True
  .type = bool
  .short_caption = Unselected atoms block dots
  .help = Include atoms that are not selected in the collision neighbor lists (-keep, -drop, -scsurface, -exposed, -asurface, -access in probe)

atom_radius_scale = 1.0
  .type = float
  .short_caption = Scale applied to atom radii
  .help = Atom radius = (r*atom_radius_scale)+atom_radius_offset (-scalevds, -vswscale in probe)

atom_radius_offset = 0.0
  .type = float
  .short_caption = Offset applied to atom radii
  .help = Atom radius = (r*atom_radius_scale)+atom_radius_offset (-addvdw in probe)

minimum_occupancy = 0.02
  .type = float
  .short_caption = Minimum occupancy
  .help = Minimum occupancy for a source atom (-minoccupancy in probe)

overlap_scale_factor = 0.5
  .type = float
  .short_caption = Overlap scale factor
  .help = Fraction of overlap assigned to each atom (-spike in probe)

ignore_lack_of_explicit_hydrogens = False
  .type = bool
  .short_caption = Ignore lack of explicit hydrogens
  .help = For an explicit-hydrogen model, ignore lack of hydrogens (probe behaved this way)

output
  .style = menu_item auto_align
{
  write_files = True
    .type = bool
    .short_caption = Write the output files
    .help = Write the output files(s) when this is True (default). Set to False when harnessing the program.

  file_name = None
    .type = str
    .short_caption = Output file name
    .help = Output file name

  dump_file_name = None
    .type = str
    .short_caption = Dump file name
    .help = Dump file name for regression testing atom characteristics (-DUMPATOMS in probe)

  format = *kinemage raw oneline json
    .type = choice
    .short_caption = Output format
    .help = Type of output to write (-oneline -unformated -kinemage in probe)

  contact_summary = False
    .type = bool
    .short_caption = Summarize contacts
    .help = Report summary of contacts (-oneline, -summary in probe)

  condensed = False
    .type = bool
    .short_caption = Condensed output
    .help = Condensed output format (-condense, -kinemage in probe)

  count_dots = False
    .type = bool
    .short_caption = Count dots, don't list
    .help = Count dots rather than listing all contacts (-countdots in probe)

  record_added_hydrogens = False
    .type = bool
    .short_caption = Record Phantom Hydrogens
    .help = Output hydrogen-bond contacts (-dumph2o in probe)

  report_hydrogen_bonds = True
    .type = bool
    .short_caption = Report hydrogen bonds
    .help = Report hydrogen bonds (-nohbout in probe)

  report_clashes = True
    .type = bool
    .short_caption = Report clashes
    .help = Report clashes (-noclashout in probe)

  report_vdws = True
    .type = bool
    .short_caption = report VdW contacts
    .help = Report van der Waals contects (-novdwout in probe)

  separate_worse_clashes = False
    .type = bool
    .short_caption = Separately report worse clashes
    .help = Separately report worse clashes (-sepworse in probe)

  group_name = ""
    .type = str
    .short_caption = Group name to use
    .help = Specify the group name (-name in probe)

  add_group_name_master_line = False
    .type = bool
    .short_caption = Add group master line
    .help = Add a master=name line on lists (-dotmaster in probe)

  add_group_line = True
    .type = bool
    .short_caption = Add group line
    .help = Add a group line on kinemage output (-nogroup in probe)

  add_kinemage_keyword = False
    .type = bool
    .short_caption = Add kinemage keyword
    .help = Add kinemage 1 to beginning of kin file (-kinemage in probe)

  add_lens_keyword = False
    .type = bool
    .short_caption = Add lens keyword
    .help = Add lens keywoard to kin file (-lens, -nolens in probe)

  color_by_na_base = False
    .type = bool
    .short_caption = Color by nucleic acid base
    .help = Color by nucleic acid base (-basecolor, -colorbase in probe)

  color_by_gap = True
    .type = bool
    .short_caption = Color by gap
    .help = Assign a color to reported gaps (-atomcolor, -gapcolor, -basecolor in probe)

  group_label = ""
    .type = str
    .short_caption = Surface-dots group label
    .help = Label for the surface-dots group (-name, -scsurface, -exposed, -asurface, -access in probe)

  bin_gaps = False
    .type = bool
    .short_caption = Bin the gaps
    .help = Bin the gaps (-gapbins in probe)

  merge_contacts = True
    .type = bool
    .short_caption = Combine wide and close contacts
    .help = Combine wide and close contacts (True in probe)

  only_report_bad_clashes = False
    .type = bool
    .help = Only report bad clashes (-onlybadout in probe)

  atoms_are_masters = False
    .type = bool
    .short_caption = Atoms are masters
    .help = Atoms are listed as masters (-element in probe)

  default_point_color = "gray"
    .type = str
    .short_caption = Default color for points
    .help = Default color for output points (-outcolor in probe)

  compute_scores = True
    .type = bool
    .short_caption = Compute scores rather than counting
    .help = Compute scores rather than just counting dots (-spike, -nospike in probe)

   altid_as_pointmaster = True
    .type = bool
    .short_caption = Add alternate IDs as point masters for atoms that are in alternate conformations
    .help = Add alternate IDs as point masters for atoms that are in alternate conformations
}
''' + Helpers.probe_phil_parameters

program_citations = phil.parse('''
citation {
  authors = Word, et. al.
  journal = J. Mol. Biol.
  volume = 285
  pages = 1711-1733
  year = 1999
  external = True
}
''')

################################################################################
# List of all of the keys for atom classes, including all elements and all
# nucleic acid types.  These are in the order that the original Probe reported
# them.  Based on atomprops.h:INIT_ATOM_TABLE from original probe.
_allAtomClasses = ['ignore',
                      'H','C','N','O','P','S','As','Se','F','Cl','Br','I',
                      'Li','Na','Al','K','Mg','Ca','Mn','Fe','Co','Ni','Cu','Zn',
                      'Rb','Sr','Mo','Ag','Cd','In','Cs','Ba','Au','Hg','Tl','Pb',
                      'V','Cr','Te','Sm','Gd','Yb','W','Pt','U',
                      'He','Be','B','Ne','Se','Ar','Sc','Ti','Ga','Ge','Kr','Y','Zr',
                      'Sn','Sb','Xe','La','Ce','Fr','Ra','Th',
                      'Nb','Tc','Ru','Rh','Pd','Pr','Nd','Pm','Eu','Tb','Dy','Ho','Er',
                      'Tm','Lu','Hf','Ta','Re','Os','Ir','Bi','Po','At','Rn','Ac','Pa',
                      'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
                      'a','c','t/u','g','other na','nonbase']

################################################################################
# Dictionary of dictionaries of lists structure holding lists of DotInfo class objects,
# indexed by atom class and then by interaction type.  Fill in empty lists for all of
# the possible classes and types.
_interactionTypes = [
    probeExt.InteractionType.WideContact,
    probeExt.InteractionType.CloseContact,
    probeExt.InteractionType.WeakHydrogenBond,
    probeExt.InteractionType.SmallOverlap,
    probeExt.InteractionType.Bump,
    probeExt.InteractionType.BadBump,
    probeExt.InteractionType.StandardHydrogenBond
  ]

# ------------------------------------------------------------------------------

def _color_for_gap(gap, interactionType):
  '''
    Report the color associated with a gap (and interaction type).
    :param gap: Size of the gap in Angstroms.
    :param interactionType: InteractionType of the dot.
    :return: Kinemage name of the color associated with the class.
  '''

  if interactionType == probeExt.InteractionType.StandardHydrogenBond:
    return "greentint "
  elif gap > 0.35:
    return "blue "
  elif gap > 0.25:
    return "sky "
  elif gap > 0.15:
    return "sea "
  elif gap > 0.0:
    return "green "
  elif gap > -0.1:
    return "yellowtint "
  elif gap > -0.2:
    return "yellow "
  elif gap > -0.3:
    return "orange "
  elif gap > -0.4:
    return "red "
  else:
    return "hotpink "

# ------------------------------------------------------------------------------

def _color_for_atom_class(c):
  '''
    Report the color associated with an atom class.
    Based on atomprops.h:INIT_ATOM_TABLE from original probe.
    :param c: Class of the atom.
    :return: Kinemage name of the color associated with the class.
  '''

  # Make sure the atom class is one that we know about
  if not c in _allAtomClasses:
    return 'magenta'

  # Check to see if this atom belongs to one of the special colors.
  if c in ['C','Ag','other na']:
    return 'white'
  elif c in ['N','He','t/u']:
    return 'sky'
  elif c in ['O']:
    return 'red'
  elif c in ['P','Ne','a']:
    return 'pink'
  elif c in ['S','c']:
    return 'yellow'
  elif c in ['Se','F','Cl']:
    return 'green'
  elif c in ['Br','I']:
    return 'brown'
  elif c in ['Co']:
    return 'blue'
  elif c in ['Cu','Ar']:
    return 'orange'
  elif c in ['Au']:
    return 'gold'
  elif c in ['Kr']:
    return 'greentint'
  elif c in ['Xe']:
    return 'magenta'
  elif c in ['Rn']:
    return 'pinktint'
  elif c in ['g']:
    return 'sea'

  # Most atom types, the default.
  return 'grey'

# ------------------------------------------------------------------------------

def _condense(dotInfoList, condense):
  '''
    Condensing the list of dots for use in raw or JSON output, sorting and removing
    duplicates.
    :param dotInfoList: List of DotInfo structures to sort and perhaps condense.
    :param condense: Boolean telling whether to condense the output, removing duplicates.
    :return: Condensed dotlist.
  '''

  ret = []

  # Handle all of the dots associated with each source atom as a group.
  # This will be from curAtomIndex to curAtomEndIndex.
  curAtomIndex = 0
  while curAtomIndex < len(dotInfoList):

    # Find the last dot in the current atom, which may be at the end of the list.
    curAtomEndIndex = len(dotInfoList) - 1
    for curAtomEndIndex in range(curAtomIndex+1, len(dotInfoList)):
      if dotInfoList[curAtomIndex].src != dotInfoList[curAtomEndIndex].src:
        curAtomEndIndex -= 1
        break

    # Sort the dots for the same source atom based on characteristics of their target atom.
    # We include the XYZ position in the sort so that we get the same order and grouping each
    # time even though the phantom H? atoms are otherwise identical.
    thisAtom = sorted(dotInfoList[curAtomIndex:curAtomEndIndex+1])

    # Remove duplicates (same target atom) if we've been asked to.
    # We do this by scanning through and accumulating counts as long as the target
    # atom is the same and by appending a new entry when the target atom is different.
    # The result is a single entry for each target atom with a count of the number of
    # dots that were associated with it in the resulting entry.
    if condense and len(thisAtom) > 0:
      thisAtom[0].dotCount = 1
      condensed = [ thisAtom[0] ]
      for i in range(1,len(thisAtom)):
        if thisAtom[i-1].target.memory_id() == thisAtom[i].target.memory_id():
          condensed[-1].dotCount += 1
        else:
          thisAtom[i].dotCount = 1
          condensed.append(thisAtom[i])
      thisAtom = condensed

    # Append the sorted and potentially condensed list to the return list
    ret.extend(thisAtom)

    # Handle the chunk of dots on the next atom
    curAtomIndex = curAtomEndIndex + 1

  return ret

# ------------------------------------------------------------------------------

def _totalInteractionCount(chainCounts):
  '''
    Find the total count of interactions of any type for the specified chain-pair type.
    :param chainCounts: One of the structures that hold the counts of interaction
    types for a given pair of chain types: _MCMCCount, _SCSCCount, _MCSCCount, _otherCount,
    or _sumCount.
    :return: Sum of results across all interaction types.
  '''
  ret = 0
  for v in chainCounts.values():
    ret += v
  return ret

# ------------------------------------------------------------------------------

class DotInfo:
  # Dot class storing information about an individual dot.
  def __init__(self, src, target, loc, spike, overlapType, gap, ptmaster, angle):
    self.src = src                  # Source atom for the interaction
    self.target = target            # Target atom for the interactions
    self.loc = loc                  # Location of the dot start
    self.spike = spike              # Location of the dot end
    self.overlapType = overlapType  # Type of overlap the interaction represents
    self.gap = gap                  # Gap between the atoms
    self.ptmaster = ptmaster        # Main/side chain interaction type
    self.angle = angle              # Angle associated with the bump
    self.dotCount = 1               # Used by _condense and raw/JSON output to count dots on the same source + target

  def _makeName(self, atom):
      # Make the name for an atom, which includes its chain and residue information
      # along with other atom data, and also includes its location to distinguish
      # among H? atoms that are otherwise identical.
      return "{}{:4.4s}{}{} {}{:1s} {:.3f} {:.3f} {:.3f}".format(
        atom.parent().parent().parent().id, # chain
        str(atom.parent().parent().resseq_as_int()), # residue number
        atom.parent().parent().icode, # insertion code
        atom.parent().resname, # residue name
        atom.name, # atom name
        atom.parent().altloc, # alternate location
        atom.xyz[0], atom.xyz[1], atom.xyz[2])

  def __lt__(self, other):
      # Sort dots based on characteristics of their source and target atoms, then their gap
      # (smallest/most negative gap to largest).
      # We include the XYZ position in the sort so that we get the same order and grouping each
      # time even though the phantom H? atoms are otherwise identical.
      # There may be no target atoms specified (may be Python None value), which will
      # be treated as the empty string.

      selfName = self._makeName(self.src)
      if self.target is not None:
        selfName += self._makeName(self.target)
      otherName = other._makeName(other.src)
      if other.target is not None:
        otherName += other._makeName(other.target)

      return (selfName < otherName) or (
        (selfName == otherName) and (self.gap < other.gap)
      )

# ------------------------------------------------------------------------------

def Test():
  '''
    Run tests on the functions that are not part of the program class.
    Throw an assertion error if there is a problem with one of them.
  '''

  #=====================================================================================
  # Test the _condense() method.

  atoms = [ # Different atoms for different indices
    pdb.hierarchy.atom(), pdb.hierarchy.atom(), pdb.hierarchy.atom(), pdb.hierarchy.atom()
  ]
  # Name the atoms distinctly so that they will sort in order.
  for i,a in enumerate(atoms):
    a.name = str(i)
  ag1 = pdb.hierarchy.atom_group()
  for a in atoms:
    ag1.append_atom(a)
  rg1 = pdb.hierarchy.residue_group()
  rg1.append_atom_group(ag1)
  rg1.resseq = 1
  c1 = pdb.hierarchy.chain()
  c1.append_residue_group(rg1)

  sourceTarget = [  # Index of source atom, target atom pairs to add into the dots list
    (1,1), (1,2), (1,1), (1,2),
    (2,1),
    (3,1), (3,1), (3,1), (3,2), (3,2), (3,2)
  ]
  dots = [  # Construct a test dots list based on the sourceTarget tuples.
    DotInfo(atoms[src],atoms[trg],(0,0,0), (0,0,0), probeExt.OverlapType.Ignore, 0.0, ' ', 0.0)
      for (src,trg) in sourceTarget
  ]

  # Test when only sorting
  inorder = _condense(dots, False)
  assert len(inorder) == len(dots), "probe2:Test(): Unexpected length from _condense when not condensing"
  assert inorder[0].target == inorder[1].target, "probe2:Test(): Unexpected sorted value from _condense when not condensing"
  assert inorder[1].target != inorder[2].target, "probe2:Test(): Unexpected sorted value from _condense when not condensing"

  # Test when also condensing
  inorder = _condense(dots, True)
  assert len(inorder) == 5, "probe2:Test(): Unexpected length from _condense when condensing"
  assert inorder[0].target != inorder[1].target, "probe2:Test(): Unexpected sorted value from _condense when condensing"
  assert inorder[-1].dotCount == 3, "probe2:Test(): Unexpected dot count value from _condense when condensing"

  #=====================================================================================
  # Test the _totalInteractionCount() method.  We make stand-in dictionaries using a stand-in
  # list.
  interactionTypes = [0, 1, 2, 3, 4, 5, 6]
  MCMCCount = {}
  for t in interactionTypes:
    MCMCCount[t] = 1

  assert _totalInteractionCount(MCMCCount) == len(interactionTypes), "probe2:Test(): _totalInteractionCount(MCMCCount) failed"

  #=====================================================================================
  # Test the _color_for_gap() method.
  table = [ [0.3, "sky "], [0.1, "green "], [-0.5, "hotpink "]]
  hydro = "greentint "
  for t in table:
    assert _color_for_gap(t[0], probeExt.InteractionType.CloseContact) == t[1], "probe2:Test(): _color_for_gap("+str(t[0])+") failed to return "+t[1]
    assert _color_for_gap(t[0], probeExt.InteractionType.StandardHydrogenBond) == hydro, "probe2:Test(): _color_for_gap() for a hydrogen bond failed to return "+hydro

  #=====================================================================================
  # Test the _color_for_atom_class() method.
  table = [ ["Bob", "magenta"], ["Ag", "white"], ["Cu", "orange"], ["Rn","pinktint"] ]
  for t in table:
    assert _color_for_atom_class(t[0]) == t[1], "probe2:Test(): _color_for_atom_class("+str(t[0])+") failed to return "+t[1]

  print('Success!')

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
probe2 version {}

This program replaces the original "probe" program from the Richarson lab
at Duke University and was developed by them as part of a supplemental award.

It computes the MolProbity Probe score for a file, or a subset of the file,
producing summaries or lists of all contacts, in Kinemage/raw/JSON format, depending
on the Phil parameters.

By default, it compares all atoms in all alternates that meet an occupancy
criterion against themselves and produces a Kinemage-format file showing all of
the dot interactions.  See below for the Phil parameter equivalents to some
original probe command-line arguments.  (Note that the original probe selected
only the a alternate by default, but version 4 of probe2 selects all alternates by default
because it also adds point masters for all alternates.)

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed

Output:
  Kinemage or text file describing the score and other information,
  depending on the parameters.

  If neither output.file_name nor output.filename is specified, it will write
  to a file with the same name as the input model file name but with the
  extension replaced with with '.kin', '.txt', or '.json' depending on the
  parameters (.kin when output.format == kinemage and output.count_dots == False).

  In addition to writing files, this is derived from the Program Template object
  and the run() method returns a dictionary whose key values are atom classes
  (atom names, NA bases, other na and nonbase depending on how the program
  was run).  Each value is a dictionary of dot interaction types (wide contact,
  close contact, weak hydrogen bonds, small overlap, bump, bad bump, hydrogen
  bond) where not all types will be filled in based on the way the program was
  run.  The value for each interaction type entry is an array of DotInfo objects
  that describe all dots of that type found by the run.

Note:
  Some approaches require the target_selection parameter.  Setting the
  target_selection to "=" will re-use the source for the target.  In all
  other cases, the string passed in will be used as a CCTBX selection on
  the model to select a subset of its atoms.

  The original Probe program had two ways to specify whether HET atoms were included
  and whether water atoms were include, in the selection description and as separate
  command-line arguments.  The command-line arguments are not present in Probe2, they
  must be specified as part of the selection criteria.  Also Probe2 does not break out
  aromatic Carbons as Car in a separate category when counting dots, they are treated
  as C for reporting purposes.

  The most simple dotkin:
    mmtbx.probe2 approach=self source_selection="all" output.file_name=out.kin input.pdb

  The probe2 command line to test a ligand named TMP against everything else:
    mmtbx.probe2 approach=both source_selection="resname TMP" target_selection="not resname TMP" PDBfilename

  Equivalent PHIL arguments for original Probe command-line options:
    -defaults:
      source_selection="(altid a or altid '' or altid ' ') and occupancy > 0.33"
      approach=self
      excluded_bond_chain_length=4
      include_mainchain_mainchain=True
    -kinemage:
      output.add_kinemage_keyword=True
      output.count_dots=False
      output.format=kinemage
      output.condensed=False
    -scsurface:
      approach=surface
      source_selection="not water"
      keep_unselected_atoms=False
      probe.radius=1.4
      group_name="SCS"
    -exposed:
      approach=surface
      source_selection="(altid a or altid '' or altid ' ') and occupancy > 0.33"
      keep_unselected_atoms=False
      probe.radius=1.4
      group_name="SCS"
    -asurface:
      approach=surface
      source_selection="not water"
      keep_unselected_atoms=False
      probe.radius=0.0
      group_name="AS"
    -access:
      approach=surface
      source_selection="not water"
      keep_unselected_atoms=False
      atom_radius_offset=1.4
      probe.radius=0.0
      group_name="AS"
    -scan0:
      source_selection="(altid a or altid '' or altid ' ') and bfactor < 40 occupancy > 0.33"
      approach=self
      excluded_bond_chain_length=4
      include_mainchain_mainchain=True
    -scan1:
      approach=once
      excluded_bond_chain_length=4
      source_selection="(altid a or altid '' or altid ' ') and bfactor < 40 and occupancy > 0.33"
      target_selection="((altid a or altid '' or altid ' ') and bfactor < 40 and occupancy > 0.65) or (not water and occupancy > 0.33)"
'''.format(version)
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']
  citations = program_citations
  epilog = '''
  For additional information and help, see http://kinemage.biochem.duke.edu/software/probe
  and http://molprobity.biochem.duke.edu
  '''

# ------------------------------------------------------------------------------

  def _scaled_atom_radius(self, a):
    '''
      Find the scaled and offset radius for the specified atom.  This will be called on each
      atom after their extra information has been loaded to determine the scaled and offset
      value to use for the remainder of the program.
      :param a: Atom whose radius is to be scaled
      :return: Scaled and offset radius of the atom.
    '''
    rad = self._extraAtomInfo.getMappingFor(a).vdwRadius
    if rad <= 0:
      alt = a.parent().altloc
      if alt == "":
        alt = " "
      resName = a.parent().resname.strip().upper()
      resID = str(a.parent().parent().resseq_as_int())
      chainID = a.parent().parent().parent().id
      myFullName = "chain "+str(chainID)+" "+resName+" "+resID+" "+a.name+" "+alt
      raise Sorry("Invalid radius for atom look-up: "+myFullName+"; rad = "+str(rad))
    return self.params.atom_radius_offset + (rad * self.params.atom_radius_scale)


  def _describe_atom_for_debug(self, a):
      resName = a.parent().resname.strip().upper()
      resID = str(a.parent().parent().resseq_as_int())
      chainID = a.parent().parent().parent().id
      iCode = a.parent().parent().icode
      alt = a.parent().altloc
      return "{:>2s}{:>4s}{}{} {}{:1s}".format(chainID, resID, iCode, resName, a.name, alt)

# ------------------------------------------------------------------------------

  def _atom_class_for(self, a):
    '''
      Assign the atom class for a specified atom.
      :param a: Atom whose class is to be specified
      :return: If our parameters have been set to color and sort by NA base,
      then it returns the appropriate base name.  Otherwise, it returns the
      element of the atom with any second letter in the element name lower-case
      to match the values in the _allAtomClasses list.
    '''
    if not self.params.output.color_by_na_base:
      val = a.element
      if len(val) > 1:
        val = val[0] + val[1:].lower()
      return val
    else:
      resName = a.parent().resname
      cl = common_residue_names_get_class(name = resName)
      if cl == "common_rna_dna" or cl == "modified_rna_dna":
        cleanName = resName.upper().strip()
        if cleanName in ['U','URA','UTP','UDP','UMP','UR',
                         'T','THY','TTP','TDP','TMP','5MU','DT','TR']:
          return 't/u'
        elif cleanName in ['A','ADE','ATP','ADP','AMP','1MA','RIA','T6A','DA','AR']:
          return 'a'
        elif cleanName in ['C','CYT','CTP','CDP','CMP','5MC','OMC','DC','CR']:
          return 'c'
        elif cleanName in ['G','GUA','GTP','GDP','GMP','GSP','1MG','2MG','M2G','7MG','OMG','DG','GR']:
          return 'g'
        return 'other na'
      else:
        return "nonbase"

# ------------------------------------------------------------------------------

  def _save_dot(self, src, target, atomClass, loc, spike, overlapType, gap, ptmaster, angle):
    '''
      Generate and store a DotInfo entry with the specified parameters.  It will be stored
      into the self._results data structure.
      :param src: Source atom for the dot.
      :param target: Target atom for the dot, if any.
      :param atomClass: Atom class of this dot, indicates where to store.
      :param loc: Location of the dot start.
      :param spike: Location of the dot end.
      :param overlapType: Type of overlap for the dot.
      :param gap: Gap spacing for the dot.
      :param ptmaster: ptmaster entry for the dot.
      :param angle: angle for the dot.
      :return: As a side effect, this will add a new entry into one of the lists in the
      self._results data structure.
    '''
    self._results[atomClass][self._dotScorer.interaction_type(
        overlapType,gap, self.params.output.separate_worse_clashes)].append(
      DotInfo(src, target, loc, spike, overlapType, gap, ptmaster, angle)
    )

# ------------------------------------------------------------------------------

  def _generate_interaction_dots(self, sourceAtoms, targetSet, spatialQuery, phantomsQuery, bondedNeighborLists):
    '''
      Find all interaction dots for the specified atoms.
      This does not include locations where the probe is inside a bonded neighbor.
      :param sourceAtoms: Sorted list of atoms that can be the source of an interaction.
      :param targetSet: Set of atoms that are targets; others will block dots but not interact
      (can be the same as sourceAtoms for some approaches).
      :param spatialQuery: Spatial-query structure that can be used to look up nearby atoms
      (all atoms, whether or not they are targets).
      :param phantomsQuery: Spatial-query structure that can be used to look up nearby Phantom
      Hydrogens (whether or not they are in the source or target atoms).
      :param bondedNeighborLists: List of bonded neighbors for atoms in sourceAtoms.
      :return: Side effect: Add dots to the self._results data structure by
      atomclass and dot type.
    '''

    # Store constants used frequently
    probeRadius = self.params.probe.probe_radius
    include_mainchain_mainchain = self.params.include_mainchain_mainchain
    minimum_occupancy = self.params.minimum_occupancy
    include_water_water = self.params.include_water_water
    excluded_bond_chain_length = self.params.excluded_bond_chain_length
    maxRadius = 2*self._maximumVDWRadius + 2 * self.params.probe.probe_radius

    for src in sourceAtoms:
      # Find out what class of dot we should place for this atom.
      atomClass = self._atomClasses[src]

      # Generate no dots for ignored atoms.
      if atomClass == 'ignore':
        continue

      # Generate no dots for atoms with too-low occupancy
      if src.occ < minimum_occupancy:
        continue

      # Find atoms that are close enough that they might touch.
      nearby = spatialQuery.neighbors(src.xyz, 0.00001, maxRadius)

      # Find our characteristics
      srcMainChain = self._inMainChain[src]
      srcSideChain = self._inSideChain[src]
      srcHet = self._inHet[src]
      srcInWater = self._inWater[src]
      srcExtra = self._extraAtomInfo.getMappingFor(src)
      srcModel = src.parent().parent().parent().parent().id

      # Select those that are actually within the contact distance based on their
      # particular radius (this query includes only target atoms, so we don't need to check separately for that).
      # Also verify that the potential target atoms meet our criteria based on parameters.
      # Keep a list of nearby Phantom Hydrogens in case we need to exclude them.
      atomSet = set()
      for n in nearby:
        nMainChain = self._inMainChain[n]
        nHet = self._inHet[n]
        nInWater = self._inWater[n]
        nExtra = self._extraAtomInfo.getMappingFor(n)
        nModel = n.parent().parent().parent().parent().id

        d = (Helpers.rvec3(n.xyz) - Helpers.rvec3(src.xyz)).length()
        if (d <= nExtra.vdwRadius + srcExtra.vdwRadius + 2*probeRadius):

          # if both atoms are in the same non-HET chain and on the main chain, then skip
          # if we're not allowing mainchain-mainchain interactions.
          # The atoms must be on the same chain to be skipped.
          if not include_mainchain_mainchain and (
                (srcMainChain and nMainChain) and not (srcHet or nHet) and
                (src.parent().parent().parent().id == n.parent().parent().parent().id) # Same chain
              ):
            continue
          # Skip atoms that are marked to be ignored
          if self._atomClasses[n] == 'ignore':
            continue
          # Skip atoms with too low occupancy
          elif n.occ < minimum_occupancy:
            continue
          # Check for water-water interactions
          elif srcInWater and nInWater:
            # Skip water-water interactions unless they are allowed
            if (not include_water_water):
              continue
            # Ensure that we're not from the same water; we don't interact with ourself
            elif src.parent() == n.parent():
              continue
          # Skip atoms that are in non-compatible alternate conformations
          elif not Helpers.compatibleConformations(src, n):
            continue
          # Skip atoms that are in different models.
          elif srcModel != nModel:
            continue
          atomSet.add(n)

      # Check the atoms for interactions
      if len(atomSet) > 0:
        # Find the atoms that are bonded to the source atom within the specified hop
        # count.  Limit the length of the chain to 3 if neither the source nor the final
        # atom is a Hydrogen.
        excluded = Helpers.getAtomsWithinNBonds(src, bondedNeighborLists, self._extraAtomInfo, probeRadius,
          excluded_bond_chain_length, 3)

        # For Phantom Hydrogens, add any non-Acceptor atom in the atom list into the
        # excluded list and also add nearby Phantom Hydrogens into the excluded list.
        # @todo Consider whether we'd rather handle this by making bonds between the
        # Phantoms and their water Oxygens (both directions), which will shield their
        # contacts from one another and (1) avoid removing sections of hydrogen bond patterns
        # that fall inside atoms that are covalently bonded to acceptors, and (2) remove
        # the inner collision of the water Oxygen with the acceptor that also makes
        # a Hydrogen bond with the Phantom Hydrogen.
        if srcExtra.isDummyHydrogen:
          nearbyPhantomHydrogens = set(phantomsQuery.neighbors(src.xyz, 0.00001, maxRadius))
          newExclusions = set()
          for a in atomSet:
            if not self._extraAtomInfo.getMappingFor(a).isAcceptor:
              newExclusions.add(a)
          excluded = list(set(excluded).union(newExclusions).union(nearbyPhantomHydrogens))

        # Remove all of the excluded atoms from the interaction set so we don't
        # put spurious dots on them.
        for e in excluded:
          atomSet.discard(e)

        # Check each dot to see if it interacts with non-bonded nearby target atoms.
        srcDots = self._dots[src]
        scale = self.params.overlap_scale_factor
        for dotvect in srcDots:

          # Find out if there is an interaction
          res = self._dotScorer.check_dot(src, dotvect, probeRadius, list(atomSet), excluded, scale)

          # Classify the interaction and store appropriate results unless we should
          # ignore the result because there was not valid overlap.
          overlapType = res.overlapType

          # If the overlap type is NoOverlap, check dot to make sure it is not annular.
          # This excludes dots that are further from the contact than dots could be at
          # the ideal just-touched contact.
          if overlapType == probeExt.OverlapType.NoOverlap and res.annular:
            continue

          # Handle any dots that should not be ignored.
          if overlapType != probeExt.OverlapType.Ignore:

            # If the cause of the dot is not in the target set, we ignore the dot.
            if not res.cause in targetSet:
              continue

            # If the overlap type is not a Hydrogen bond, then check the occupancy of atoms where
            # at least one of the pair is on the "" or " " alternate conformation to make sure the
            # sum of their occupancies is greater than 1.
            if overlapType != probeExt.OverlapType.HydrogenBond:
              if (src.parent().altloc in ['',' ']) or (res.cause.parent().altloc in ['',' ']):
                if src.occ + res.cause.occ <= 1:
                  continue

            # See whether this dot is allowed based on our parameters.
            spo = self.params.output
            show = False
            interactionType = self._dotScorer.interaction_type(overlapType,res.gap, self.params.output.separate_worse_clashes)
            if interactionType == probeExt.InteractionType.Invalid:
              print('Warning: Invalid interaction type encountered (internal error)', file=self.logger)
              continue

            # Main branch if we're reporting other than bad clashes
            if (not spo.only_report_bad_clashes):
              # We are reporting other than bad clashes, see if our type is being reported
              if spo.report_hydrogen_bonds and (overlapType == probeExt.OverlapType.HydrogenBond):
                show = True
              elif spo.report_clashes and (overlapType == probeExt.OverlapType.Clash):
                show = True
              elif spo.report_vdws and (overlapType == probeExt.OverlapType.NoOverlap):
                show = True
            else:
              # We are only reporting bad clashes.  See if we're reporting clashes and this is
              # a bad one.
              if (spo.report_clashes and interactionType in [
                    probeExt.InteractionType.Bump, probeExt.InteractionType.BadBump]):
                show = True

            # If we're not showing this one, skip to the next
            if not show:
              continue

            # Determine the ptmaster (main/side chain interaction type) and keep track of
            # counts for each type.
            causeMainChain = self._inMainChain[res.cause]
            causeSideChain = self._inSideChain[res.cause]
            causeHet = self._inHet[res.cause]
            ptmaster = ' '
            if srcMainChain and causeMainChain:
              if (not srcHet) and (not causeHet): # This may be a redundant check
                ptmaster = 'M'
                self._MCMCCount[interactionType] += 1
            elif srcSideChain and causeSideChain:
              if (not srcHet) and (not causeHet): # This may be a redundant check
                ptmaster = 'S'
                self._SCSCCount[interactionType] += 1
            elif ( (srcMainChain and causeSideChain) or (srcSideChain and causeMainChain) ):
              if (not srcHet) and (not causeHet): # This may be a redundant check
                ptmaster = 'P'
                self._MCSCCount[interactionType] += 1
            else:
              ptmaster = 'O'
              self._otherCount[interactionType] += 1

            # Find the locations of the dot and spike by scaling the dot vector by the atom radius and
            # the (negative because it is magnitude) overlap.
            loc = Helpers.rvec3(src.xyz) + Helpers.rvec3(dotvect)
            spikeloc = ( Helpers.rvec3(src.xyz) + Helpers.rvec3(dotvect).normalize() *
                         (self._extraAtomInfo.getMappingFor(src).vdwRadius - res.overlap) )

            # Save the dot
            self._save_dot(src, res.cause, atomClass, loc, spikeloc, overlapType, res.gap, ptmaster, 0)

# ------------------------------------------------------------------------------

  def _generate_surface_dots_for(self, src, nearby):
    '''
      Find all surface dots for the specified atom.
      This does not include locations where the probe is interacting with
      a nearby atom, so it is a subset of the skin dots (for which only the
      dots themselves are outside of the nearby atoms).
      :param src: Atom whose surface dots are to be found.
      :param nearby: Atoms that are nearby to src and might block surface dots.
      :return: Side effect: Add dots on the surface of the atom to the
              self._results data structure by atomclass and dot type.
    '''

    # Generate no dots for ignored atoms.
    if self._atomClasses[src] == 'ignore':
      return

    # Check all of the dots for the atom and see if they should be
    # added to the list.
    srcInWater = self._inWater[src]
    r = self._extraAtomInfo.getMappingFor(src).vdwRadius
    pr = self.params.probe.probe_radius
    srcDots = self._dots[src]
    for dotvect in srcDots:
      # Dot on the surface of the atom, at its radius; both dotloc and spikeloc from original code.
      # This is where the probe touches the surface.
      dotloc = Helpers.rvec3(src.xyz) + Helpers.rvec3(dotvect)
      # Dot that is one probe radius past the surface of the atom, exploring for contact with nearby
      # atoms.  This is the location of the center of the probe.
      exploc = Helpers.rvec3(src.xyz) + Helpers.rvec3(dotvect).normalize() * (r + pr)

      # If the exploring dot is within a probe radius + vdW radius of a nearby atom,
      # we don't add a dot.
      okay = True
      for b in nearby:
        bInWater = self._inWater[b]
        # If we should ignore the nearby atom, we don't check it.
        if self._atomClasses[b] == 'ignore':
          continue
        # If we're ignoring water-water interactions and both src and
        # nearby are in a water, we should ignore this as well (unless
        # both are hydrogens from the same water, in which case we
        # continue on to check.)
        elif ((not self.params.include_water_water) and srcInWater and bInWater
              and src.parent() != b.parent() ):
          continue

        # The nearby atom is one that we should check interaction with, see if
        # we're in range.  If so, mark this dot as not okay because it is inside a
        # nearby atom.
        if ( (Helpers.rvec3(b.xyz) - exploc).length() <=
            pr + self._extraAtomInfo.getMappingFor(b).vdwRadius ):
          okay = False

      # If this dot is okay, add it to the internal data structure based on its
      # atom class and overlap type.
      if okay:
        self._save_dot(src, None, self._atomClasses[src], dotloc, dotloc,
                       probeExt.OverlapType.NoOverlap, 0.0, ' ', 0.0)

# ------------------------------------------------------------------------------

  def _count_skin_dots_for(self, src, bonded):
    '''
      Count all skin dots for the specified atom.
      :param src: Atom whose surface dots are to be found.
      :param bonded: Atoms that are bonded to src by one or more hops.
      :return: Side effect: Add dots on the surface of the atom to the
              self._results data structure by atomclass and dot type.
    '''

    # No dots yet...
    ret = 0

    # Generate no dots for ignored atoms or for phantom hydrogens
    if self._atomClasses[src] == 'ignore' or self._extraAtomInfo.getMappingFor(src).isDummyHydrogen:
      return 0

    # If we should ignore the bonded element, we don't check it.
    # Remove any ignored atoms from the list of bonded atoms to pull this check out of
    # the inner loop.
    srcDots = self._dots[src]
    realBonded = []
    for b in bonded:
      if self._atomClasses[b] != 'ignore':
        realBonded.append(b)

    # Check all of the dots for the atom and see if they should be
    # added to the list.
    return self._dotScorer.count_surface_dots(src, srcDots, realBonded)

# ------------------------------------------------------------------------------

  def _count_skin_dots(self, atoms, bondedNeighborLists):
    '''
      Count all skin dots for the atoms passed in.
      :param atoms: Atoms to check.
      :param bondedNeighborLists: Neighbor list including these atoms.
      This is used to normalize output scores.
      :return: Number of skin dots on any of the atoms in the source selection.
    '''

    ret = 0

    # Store parameters that are used in the inner loop
    excluded_bond_chain_length = self.params.excluded_bond_chain_length

    for src in atoms:
      # Find the atoms that are bonded to the source atom within the specified hop
      # count.  Limit the length of the chain to 3 if neither the source nor the final
      # atom is a Hydrogen.
      # We check only out to a probe radius of 0 (atoms actually overlapping)
      neighbors = Helpers.getAtomsWithinNBonds(src, bondedNeighborLists, self._extraAtomInfo, 0.0,
        excluded_bond_chain_length, 3)

      # Count the skin dots for this atom.
      ret += self._count_skin_dots_for(src, neighbors)

    # Return the total count
    return ret

# ------------------------------------------------------------------------------

  def _writeRawOutput(self, groupName, masterName, writeJSON = False):
    '''
      Describe raw summary counts for data of various kinds.
      :param groupName: Name to give to the group.
      :param masterName: Name for the beginning of each line.
      :param writeJSON: If True, write the output in JSON format.  Otherwise, write one raw entry per line.
      :return: String to be added to the output.
    '''

    if writeJSON:
      ret = '{ "flat_results" : ['
    else:
      ret = ''

    # Provide a short name for each interaction type
    mast = {}
    for t in _interactionTypes:
      mast[t] = probeExt.DotScorer.interaction_type_short_name(t)

    # Store values that we will need often
    density = self.params.probe.density
    gap_weight = self.params.probe.gap_weight
    bump_weight = self.params.probe.bump_weight
    hydrogen_bond_weight = self.params.probe.hydrogen_bond_weight

    # Go through all atom types and contact types and report the contacts.
    first_line = True
    for atomClass in _allAtomClasses:
      for interactionType in _interactionTypes:

        # Condensed report all of the dots of this type.
        condensed = _condense(self._results[atomClass][interactionType], self.params.output.condensed)
        for node in condensed:

          if writeJSON:
            if first_line:
              first_line = False
            else:
              ret += ','
            ret += '\n'
            ret += '{{"master": "{}", "group": "{}", "type": "{}"'.format(masterName, groupName, mast[interactionType])
          else:
            ret += "{}:{}:{}:".format(masterName, groupName, mast[interactionType])

          # Describe the source atom
          a = node.src
          resName = a.parent().resname.strip().upper()
          resID = str(a.parent().parent().resseq_as_int())
          chainID = a.parent().parent().parent().id
          iCode = a.parent().parent().icode
          if iCode == "":
            iCode = " "
          alt = a.parent().altloc
          if writeJSON:
            ret += ', "src": {{"chainID": "{}", "resID": {}, "iCode": "{}", "resName": "{}", "atomName": "{}", "alt": "{}"}}'.format(
              chainID, resID, iCode.strip(), resName, a.name, alt.strip())
          else:
            ret += "{:>2s}{:>4s}{}{} {}{:1s}:".format(chainID, resID, iCode, resName.strip(), a.name, alt)

          # Describe the target atom, if it exists
          t = node.target
          if t is None:
            if not writeJSON:
              ret += ":::::::"
          else:
            resName = t.parent().resname.strip().upper()
            resID = str(t.parent().parent().resseq_as_int())
            chainID = t.parent().parent().parent().id
            iCode = t.parent().parent().icode
            if iCode == "":
              iCode = " "
            alt = t.parent().altloc
            if writeJSON:
              ret += ', "target": {{"chainID": "{}", "resID": {}, "iCode": "{}", "resName": "{}", "atomName": "{}", "alt": "{}"}}'.format(
                chainID, resID, iCode.strip(), resName, t.name, alt.strip())
            else:
              ret += "{:>2s}{:>4s}{}{} {:<3s}{:1s}:".format(chainID, resID, iCode, resName.strip(), t.name, alt)

            r1 = self._extraAtomInfo.getMappingFor(a).vdwRadius
            r2 = self._extraAtomInfo.getMappingFor(t).vdwRadius
            sl = (node.loc-node.spike).length()
            gap = (Helpers.rvec3(a.xyz)-Helpers.rvec3(t.xyz)).length() - (r1 + r2)
            dtgp = node.gap
            score = 0.0

            if interactionType in [probeExt.InteractionType.WideContact, probeExt.InteractionType.WideContact]:
              scaledGap = dtgp / gap_weight
              score = math.exp(-scaledGap*scaledGap)
            elif interactionType in [
              probeExt.InteractionType.WeakHydrogenBond,  # don't know what to do here, because they can be both wc and cc, so will have to check
              probeExt.InteractionType.SmallOverlap,      # small overlap, doing nothing, as before
              probeExt.InteractionType.Bump,
              probeExt.InteractionType.BadBump]:          # worse overlap, same as bad overlap
                score = score = - bump_weight * sl
            else: # Hydrogen bond
              score = hydrogen_bond_weight * sl

            if self.params.output.condensed:
              if writeJSON:
                ret += ', "dotCount": {}'.format(node.dotCount)
              else:
                ret += "{}:".format(node.dotCount)

            if writeJSON:
              if self.params.output.condensed:
                # Don't write elements that are not needed for condensed output
                ret += ', "gap": {:.3f}'.format(gap)
              else:
                ret += ', "gap": {:.3f}, "dotGap": {:.3f}, "spike": [{:.3f},{:.3f},{:.3f}], "spikeLen": {:.3f}, "scoreOverDensity": {:.4f}'.format(
                  gap, dtgp, node.spike[0], node.spike[1], node.spike[2], sl, score/density)
            else:
              ret += "{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{:.4f}".format(gap, dtgp,
                node.spike[0], node.spike[1], node.spike[2], sl, score/density)

          try:
            tName = self._atomClasses[t]
            tBVal = "{:.2f}".format(t.b)
          except Exception:
            tName = ""
            tBVal = ""

          if writeJSON:
            if self.params.output.condensed:
              # Don't write elements that are not needed for condensed output
              ret += ', "srcClass": "{}", "targetClass": "{}"'.format(self._atomClasses[a], tName)
            else:
              ret += ', "srcClass": "{}", "targetClass": "{}", "loc": [{:.3f},{:.3f},{:.3f}]'.format(
                self._atomClasses[a], tName, node.loc[0], node.loc[1], node.loc[2])
          else:
            ret += ":{}:{}:{:.3f}:{:.3f}:{:.3f}".format(self._atomClasses[a], tName,
              node.loc[0], node.loc[1], node.loc[2])

          if writeJSON:
            ret += ', "srcBFactor": {:.2f}, "targetBFactor": {}}}'.format(a.b, tBVal)
          else:
            ret += ":{:.2f}:{}\n".format(a.b, tBVal)

    if writeJSON:
      ret += '\n] }\n'
    return ret

# ------------------------------------------------------------------------------

  def _writeOutput(self, groupName, masterName):
    '''
      Describe contacts for data of various kinds.
      :param groupName: Name to give to the group.
      :param masterName: Name for the master command.
      :return: String to be added to the output.
    '''

    ret = ''

    ptm = ' '
    color = ''
    mast = {}
    for t in _interactionTypes:
      # Probe uses spaces in these names for this function but underscores for others, so we replace
      # underscores with spaces here.
      mast[t] = probeExt.DotScorer.interaction_type_name(t).replace("_"," ")
    extraMaster = ''
    pointid = ''
    ptmast = ''
    gapNames = ['z','y','x','w','v','u','t','g','r','q','f','F','Q','R','G','T','U','V','W','X','Y','Z']
    # std gapbins scope at least -.5 to +.5, wider if probeRad > 0.25 standard
    gaplimit = int(math.floor(((2*(max(self.params.probe.probe_radius,0.25))+0.5)/0.05)+2))
    gapcounts = [0] * gaplimit
    maxgapcounts = 0
    strcName = ''

    # Rename contacts as needed
    if self.params.output.merge_contacts:
      mast[probeExt.InteractionType.WideContact] = mast[probeExt.InteractionType.CloseContact] = 'vdw contact'
    if self.params.approach == 'surface':
      mast[probeExt.InteractionType.CloseContact] = 'surface'

    if self.params.output.add_group_name_master_line:
      extraMaster = ' master={{{}}}'.format(masterName)

    ret += "@subgroup dominant {{{}}}\n".format(groupName)

    if self.params.approach == 'surface':
      ret += "@master {{{}}}\n".format(mast[1])
    else:
      if self.params.output.report_vdws and not self.params.output.only_report_bad_clashes:
        ret += "@master {{{}}}\n".format(mast[probeExt.InteractionType.WideContact])
        if not self.params.output.merge_contacts:
          ret += "@master {{{}}}\n".format(mast[probeExt.InteractionType.CloseContact])
      if self.params.output.report_clashes or self.params.output.only_report_bad_clashes:
        if not self.params.output.only_report_bad_clashes:
          ret += "@master {{{}}}\n".format(mast[probeExt.InteractionType.SmallOverlap])
        ret += "@master {{{}}}\n".format(mast[probeExt.InteractionType.Bump])
        if self.params.output.separate_worse_clashes:
          ret += "@master {{{}}}\n".format(mast[probeExt.InteractionType.BadBump])
      if self.params.output.report_hydrogen_bonds and not self.params.output.only_report_bad_clashes:
        ret += "@master {{{}}}\n".format(mast[probeExt.InteractionType.StandardHydrogenBond])
        if self.params.probe.allow_weak_hydrogen_bonds:
          ret += "@master {{{}}}\n".format(mast[probeExt.InteractionType.WeakHydrogenBond])

    # Report count legend if any counts are nonzero.
    if _totalInteractionCount(self._MCMCCount) > 0:
      ret += "@pointmaster 'M' {McMc contacts}\n"
    if _totalInteractionCount(self._SCSCCount) > 0:
      ret += "@pointmaster 'S' {ScSc contacts}\n"
    if _totalInteractionCount(self._MCSCCount) > 0:
      ret += "@pointmaster 'P' {McSc contacts}\n"
    if _totalInteractionCount(self._otherCount) > 0:
      ret += "@pointmaster 'O' {Hets contacts}\n"

    # Report binned gap legend if we're binning gaps
    if self.params.output.bin_gaps:
      for i in range(gaplimit):
        ret += "@pointmaster '{}' {{gap {:3.2f}}}\n".format(gapNames[i],((i-11.0)/20.0)+0.05)

    # Go through all atom types and contact types and report the contacts.
    for atomClass in _allAtomClasses:
      for interactionType in _interactionTypes:
        # When we switch point types, we need to repeat the point ID.
        lastPointID = ''

        # Write list headers for types that have entries.  Do not write one for weak Hydrogen
        # bonds unless we're separating them out.
        if (len(self._results[atomClass][interactionType]) > 0 and
              (self.params.probe.allow_weak_hydrogen_bonds or
                  interactionType != probeExt.InteractionType.WeakHydrogenBond
              )
            ):
          # The formatting of the header depends on the type
          # of dot it is and whether atoms are masters.  There is a basic line for each, with addition of
          # a lens string for some cases.  Some entries are dots and others are vectors.
          lensDots = ""
          listType = '@dotlist'
          if interactionType in [probeExt.InteractionType.WideContact, probeExt.InteractionType.CloseContact]:
            if self.params.output.add_lens_keyword:
              lensDots = " lens"
          elif interactionType in [probeExt.InteractionType.SmallOverlap, probeExt.InteractionType.Bump,
              probeExt.InteractionType.BadBump]:
            listType = '@vectorlist'
          elif interactionType == probeExt.InteractionType.StandardHydrogenBond:
            # Nothing special
            pass

          # Write the header based on the settings above and whether atoms are masters.
          if self.params.output.atoms_are_masters:
            ret += "{} {{x}} color={} master={{{} dots}} master={{{}}}{}{}\n".format(
                    listType,
                    _color_for_atom_class(atomClass), atomClass, mast[interactionType], extraMaster,
                    lensDots
                   )
          else:
            ret += "{} {{x}} color={} master={{{}}}{}{}\n".format(
                      listType,
                      _color_for_atom_class(atomClass), mast[interactionType], extraMaster,
                      lensDots
                    )

        # Report all of the dots of this type.
        for node in self._results[atomClass][interactionType]:
          a = node.src
          t = node.target
          if self.params.output.bin_gaps:
            # Include trailing space for a gapbin character (second point master)
            ptmast = " '{} ' ".format(node.ptmaster)
          elif node.ptmaster == " ":
            # Blank means no point master
            ptmast = ""
          else:
            ptmast = " '{}' ".format(node.ptmaster)

          pointid = "{}{:1s}{:>3s} {:>3d} {:1s}{}".format(a.name, a.parent().altloc, a.parent().resname,
            a.parent().parent().resseq_as_int(), a.parent().parent().icode,
            a.parent().parent().parent().id)
          if pointid != lastPointID:
            lastPointID = pointid
            ret += '{{{}}}'.format(pointid)
          else:
            ret += '{"}'

          if self.params.output.color_by_gap:
            if t is not None:
              color = _color_for_gap(node.gap, interactionType)
              ret += "{}".format(color)
            else:
              ret += "{} ".format(self.params.output.default_point_color)

          # Handle gap binning if we're doing it
          if self.params.output.bin_gaps:
            Lgotgapbin = False    # until identify which gapbin
            for k in range(gaplimit):
              # pt master intervals of 0.05 from -0.5 to +0.5
              if node.gap < ((k-11.0)/20.0)+0.05:
                # Replace the fourth character of ptmast with the appropriate gap name
                ptmast = ptmast[:3]+gapNames[k]+ptmast[4:]
                gapcounts[k] += 1
                maxgapcounts = max(gapcounts[k], maxgapcounts)
                if k < gaplimit:
                  Lgotgapbin = True
                  break
            if not Lgotgapbin:
              # assign this node, aka dot, to overflow gapbin
              ptmast = ptmast[:3]+gapNames[-1]+ptmast[4:]
              gapcounts[-1] += 1

          # Add alternate conformation masters to the pointmaster unless we've been told not to
          if self.params.output.altid_as_pointmaster:
            alt = ''
            if (not a.parent().altloc in [' ', '']):
              alt = a.parent().altloc
            if (not t is None) and (not t.parent().altloc in [' ', '']):
              alt += t.parent().altloc
            if alt != '':
              # Find the index of the second single quote in the ptmaster string

              idx = ptmast.rfind("'")
              # Insert the alternate conformation into the ptmaster string at this index,
              # using the lower-case version of the alternate conformation character.
              ptmast = ptmast[:idx] + alt.lower() + ptmast[idx:]

          if interactionType in [probeExt.InteractionType.SmallOverlap, probeExt.InteractionType.Bump,
              probeExt.InteractionType.BadBump]:
            ret += 'P {}{:.3f},{:.3f},{:.3f} {{"}}{} {}{:.3f},{:.3f},{:.3f}\n'.format(
                      ptmast, node.loc[0], node.loc[1], node.loc[2],
                      color,
                      ptmast, node.spike[0], node.spike[1], node.spike[2]
                    )
          else: # Contact or H bond
            ret += "{}{:.3f},{:.3f},{:.3f}\n".format(
                      ptmast, node.loc[0], node.loc[1], node.loc[2]
                   )

    # Print the gap bins if we have computed them.
    if self.params.output.bin_gaps:
      ret += "@text\n"
      for k in range(gaplimit):
        ret += "{{{:5.2f}, {:8d} }}\n".format(
                (((k-11.0)/20.0)+0.05),gapcounts[k]
                )

      # kinemage 2
      ret += "@kinemage 2\n"
      ret += "@group {{gapbins}} dominant\n"
      ret += "@vectorlist {gapbins}\n"
      for k in range(gaplimit-1):
        ret += "{{{:5.2f}, {:8d} }} {:5.2f}, {:8f}, 0.00\n".format(
                 (((k-11.0)/20.0)+0.05), gapcounts[k],
                 (((k-11.0)/20.0)+0.05)*maxgapcounts, gapcounts[k]
               )
      ret += "@labellist {gapbins}\n"
      ret += "{0} 0.0, -1.0, 0.0\n"

      # LXHvector output in probe had to do with -oneDotEach, so we don't include it here.

    return ret

# ------------------------------------------------------------------------------

  def _doEnumeration(self, reportSubScores, isSurface, numSkinDots):
    '''
      Compute summary counts for data of various kinds.  Called by _rawEnumerate() and
      _enumerate() to do the shared work.
      :param reportSubScores: Provide reports on different contact subscores.
      :param isSurface: Are these all surface dots?
      :param numSkinDots: The number of dots on atom skins. This is used to normalize output scores.
      :return: Tuple of values: (string_to_output, tgs, ths, thslen, tbs, tbslen, tsas,
              tGscore, tHscore, tBscore, tscore)
    '''
    retString = ''

    # Store values that we will need often
    approach = self.params.approach
    density = self.params.probe.density
    gap_weight = self.params.probe.gap_weight
    bump_weight = self.params.probe.bump_weight
    hydrogen_bond_weight = self.params.probe.hydrogen_bond_weight

    # Compute the counts
    tgs = ths = thslen = tbs = tbslen = tsas = 0
    tGscore = tHscore = tBscore = tscore = 0
    for c in _allAtomClasses:
      for t in _interactionTypes:
        res = self._results[c][t]
        if len(res) > 0:

          # gs stores all of the values unless reportSubScores is True
          gs = hs = hslen = bs = bslen = score = psas = 0

          # Print a line describing the atom class and interaction type.
          label = "external_dots "
          if not isSurface:
            label = probeExt.DotScorer.interaction_type_name(t)
          retString += "{:>3s} {:14s} ".format(c, label)

          for node in self._results[c][t]:
            if reportSubScores:
              if t in [probeExt.InteractionType.WideContact, probeExt.InteractionType.CloseContact,
                  probeExt.InteractionType.WeakHydrogenBond]:
                gs += 1
                dtgp = node.gap
                scaledGap = dtgp/gap_weight
                scoreValue = math.exp(-scaledGap*scaledGap)
                score   += scoreValue
                tGscore += scoreValue
              elif t in [probeExt.InteractionType.SmallOverlap, probeExt.InteractionType.Bump,
                  probeExt.InteractionType.BadBump]:
                bs += 1
                slen = 0.5*abs(node.gap);
                bslen += slen
                scoreValue = - bump_weight * slen
                score   += scoreValue
                tBscore += scoreValue
              else: # Hydrogen bond
                hs += 1
                slen = 0.5*abs(node.gap)
                hslen += slen
                scoreValue = hydrogen_bond_weight * slen
                score   += scoreValue
                tHscore += scoreValue
            else:
              gs += 1

            if approach == 'surface':
              p_radius = self.params.probe.probe_radius
              a_radius = self._extraAtomInfo.getMappingFor(node.src).vdwRadius
              psas += (a_radius + p_radius)*(a_radius + p_radius)/(a_radius * a_radius)

          # Finish reporting by atom class and interaction type
          if reportSubScores:
            if t in [probeExt.InteractionType.WideContact, probeExt.InteractionType.CloseContact,
                probeExt.InteractionType.WeakHydrogenBond]:
              retString += "{:7d} {:5.1f}% {:9.1f} {:9.2f}\n".format(gs, 100.0*gs/numSkinDots, score/density,
                                1000.0*score/numSkinDots)
            elif t in [probeExt.InteractionType.SmallOverlap, probeExt.InteractionType.Bump,
                probeExt.InteractionType.BadBump]:
              retString += "{:7d} {:5.1f}% {:9.1f} {:9.2f}\n".format(bs, 100.0*bs/numSkinDots, score/density,
                                1000.0*score/numSkinDots)
            else: # Hydrogen bond
              retString += "{:7d} {:5.1f}% {:9.1f} {:9.2f}\n".format(hs, 100.0*hs/numSkinDots, score/density,
                                1000.0*score/numSkinDots)
          else:
            retString += "{:7d} {:5.1f}%\n".format(gs, 100.0*gs/numSkinDots)

          # Done computing for this category, calculate totals
          tgs += gs
          ths += hs
          thslen += hslen
          tbs += bs
          tbslen += bslen
          tscore += score
          if approach == 'surface':
            tsas += psas  # tally the solvent accessible surface

    return (retString, tgs, ths, thslen, tbs, tbslen, tsas, tGscore, tHscore, tBscore, tscore)

# ------------------------------------------------------------------------------

  def _rawEnumerate(self, groupName, numberSelected, reportSubScores, isSurface, numSkinDots, masterName):
    '''
      Describe summary counts for data of various kinds.
      :param groupName: Name to give to the group.
      :param numberSelected: Number of atoms in the selection.
      :param reportSubScores: Provide reports on different contact subscores.
      :param isSurface: Are these all surface dots?
      :param numSkinDots: The number of dots on atom skins. This is used to normalize output scores.
      :param masterName: Name for the beginning of each line.
      :return: String to be added to the output.
    '''
    # The C code has a rawName parameter, but it was only nonempty for autobondrot/movingDoCommand
    # The C code has a scoreBias parameter, but it was only nonzero for autobondrot/movingDoCommand

    ret = ""

    # If we have an empty selection, report zero.
    if numberSelected <= 0 or numSkinDots <= 0:
      ret += "{:9.3f}".format(0.0)

    else:
      # Compute and report the score.  Discard anything from the return string in the count
      # routine -- we don't want to print it.
      (retString, tgs, ths, thslen, tbs, tbslen, tsas, tGscore, tHscore, tBscore, tscore
        ) = self._doEnumeration(reportSubScores, isSurface, numSkinDots)

      # Output one line of information.
      if isSurface:
        ret += "{:9.3f}".format( (tgs+tbs+ths)/self.params.probe.density )
      elif reportSubScores:
        ret += "{:9.3f}".format( tscore/self.params.probe.density )
      else:
        ret += "{:9.3f}".format( tgs )

    # Report the same information at the end of the line whether or not we counted the scores.
    if len(groupName) > 0 or len(masterName) > 0:
      ret += "#"
    if len(masterName) > 0:
      ret += " {}".format(masterName)
    if len(groupName) > 0:
      ret += " {}".format(groupName)
    ret += "\n"
    return ret

# ------------------------------------------------------------------------------

  def _count_summary(self, modeName, completed = True):
    '''
      Describe summary counts for chain-vs.-chain counts.
      :param modeName: Description of the mode of operation to report.
      :param completed: This is the last iteration, so print the accumulated values.
      :return: String to be added to the output.
    '''

    ret = ''

    # Keep a running total of values for each chain-vs.-chain list.
    # The first time we're run, fill the values with 0.
    # Clear the global counts once they have been added to the running total so we can run a new count.
    if not hasattr(self,'_MCMCTotal'):
      self._MCMCTotal = {}
      self._SCSCTotal = {}
      self._MCSCTotal = {}
      self._otherTotal = {}
      for t in _interactionTypes:
        self._MCMCTotal[t] = 0
        self._SCSCTotal[t] = 0
        self._MCSCTotal[t] = 0
        self._otherTotal[t] = 0
    for t in _interactionTypes:
      self._MCMCTotal[t] += self._MCMCCount[t]
      self._MCMCCount[t] = 0
      self._SCSCTotal[t] += self._SCSCCount[t]
      self._SCSCCount[t] = 0
      self._MCSCTotal[t] += self._MCSCCount[t]
      self._MCSCCount[t] = 0
      self._otherTotal[t] += self._otherCount[t]
      self._otherCount[t] = 0

    # Compute the sum of all subtypes per interaction type.
    sumTotal = {}
    for t in _interactionTypes:
      sumTotal[t] = self._MCMCTotal[t] + self._SCSCTotal[t] + self._MCSCTotal[t] + self._otherTotal[t]

    # If we're at the last pass, fill in our return string.
    if completed:
      if self.params.output.format == 'oneline':
        # Report the file name that was read along with its summary data on one line
        ret += ": {} ".format(self.data_manager.get_model_names()[0])
        for c in [self._MCMCTotal, self._SCSCTotal, self._MCSCTotal, self._otherTotal]:
          for t in _interactionTypes:
            ret += ":{:9d} ".format(c[t])
        ret += ":\n"
      else:
        ret += "@text\n"
        ret += "probe: {}\n".format(modeName)
        ret += "{}\n".format(self.data_manager.get_model_names()[0])
        ret += ":CONTACT:   WIDE   :  CLOSE   :  weak H-bonds  : SMALL   :   BAD    :  WORSE  :  H-BOND  :\n"
        for (c,name) in [(self._MCMCTotal, "MCMC"), (self._SCSCTotal, "SCSC"),
                         (self._MCSCTotal, "MCSC"), (self._otherTotal, "OTHER"),
                         (sumTotal, "SUM")]:
          ret += ":{:7s}".format(name)
          for t in _interactionTypes:
            ret += ":{:9d} ".format(c[t])
          ret += ":\n"

    return ret

# ------------------------------------------------------------------------------

  def _enumerate(self, groupName, numberSelected, reportSubScores, isSurface, numSkinDots):
    '''
      Describe summary counts for data of various kinds.
      :param groupName: Name to give to the group.
      :param numberSelected: Number of atoms in the selection.
      :param reportSubScores: Provide reports on different contact subscores.
      :param isSurface: Are these all surface dots?
      :param numSkinDots: The number of dots on atom skins. This is used to normalize output scores.
      :return: String to be added to the output.
    '''

    ret = ''

    # Store values that we will need often
    density = self.params.probe.density

    ret += "        \nsubgroup: {}\n".format(groupName)
    ret += "atoms selected: {}\npotential dots: {}\npotential area: {:.1f} A^2\n".format(
      numberSelected, numSkinDots, numSkinDots/density)
    if numberSelected <=0 or numSkinDots <= 0:
      ret += "empty selection\n"
      return

    if reportSubScores:
      ret += "  type                 #      %       score score/A^2 x 1000\n"
    else:
      ret += "  type                 #      %\n"

    # Compute the counts
    (retString, tgs, ths, thslen, tbs, tbslen, tsas, tGscore, tHscore, tBscore, tscore
      ) = self._doEnumeration(reportSubScores, isSurface, numSkinDots)
    ret += retString

    # Report the counts
    if reportSubScores:
      ret += "\n     tot contact:  {:7d} {:5.1f}% {:9.1f} {:9.2f}\n".format(
        tgs, 100.0*tgs/numSkinDots, tGscore/density, 1000.0*tGscore/numSkinDots
      )
      ret += "     tot overlap:  {:7d} {:5.1f}% {:9.1f} {:9.2f}\n".format(
        tbs, 100.0*tbs/numSkinDots, tBscore/density, 1000.0*tBscore/numSkinDots
      )
      ret += "     tot  H-bond:  {:7d} {:5.1f}% {:9.1f} {:9.2f}\n".format(
        ths, 100.0*ths/numSkinDots, tHscore/density, 1000.0*tHscore/numSkinDots
      )

      ret += "\n       grand tot:  {:7d} {:5.1f}% {:9.1f} {:9.2f}\n".format(
        (tgs+tbs+ths), 100.0*(tgs+tbs+ths)/numSkinDots, tscore/density, 1000.0*tscore/numSkinDots
      )
      ret += "\ncontact surface area: {:.1f} A^2\n".format((tgs+tbs+ths)/density)
    else:
      ret += "             tot:  {:7d} {:5.1f}%\n\n".format(tgs, 100.0*tgs/numSkinDots)
      ret += "   contact surface area: {:.1f} A^2\n".format(tgs/density)
      if self.params.approach == 'surface':
        ret += "accessible surface area: {:.1f} A^2\n\n".format(tsas/density)

    return ret


# ------------------------------------------------------------------------------

  def _describe_selection_and_parameters(self, groupLabel, selectionName):
    '''
      Describe the selection type and other parameters for a run.  Called by various run types.
      :param groupLabel: Name to give to the group.
      :param selectionName: Name of the selection mode: 'self', 'once'.
      :return: String to be added to the output.
    '''

    ret = ''
    ret += "selection: {}\nname: {}\n".format(selectionName, groupLabel)
    ret += "density: {:.1f} dots per A^2\nprobeRad: {:.3f} A\nVDWrad: (r * {:.3f}) + {:.3f} A\n".format(
      self.params.probe.density, self.params.probe.probe_radius, self.params.atom_radius_scale,
      self.params.atom_radius_offset)
    ret += "score weights: gapWt={:0g}, bumpWt={:0g}, HBWt={:0g}\n".format(
      self.params.probe.gap_weight, self.params.probe.bump_weight, self.params.probe.hydrogen_bond_weight)
    return ret

# ------------------------------------------------------------------------------

  def _report_single_interaction(self, groupLabel, selectionName, comparisonString, intersectionName,
      numModels, modelIndex, bondedNeighborLists):
    '''
      Print information about a single interaction, either self interaction or once interaction.
      :param groupLabel: Name to give to the group.
      :param selectionName: Name of the selection mode: 'self', 'once'.
      :param comparisonString: String decribing the comparison: '1->1', '1->2'.
      :param intersectionName: Name of the intersection being done: 'SelfIntersect', 'IntersectOnce'.
      :param numModels: Number of models we are running over.
      :param modelIndex: Current model we are running.
      :param bondedNeighborLists: List of bonded neighbors for atoms in sourceAtoms.
      :return: String to be added to the output.
    '''

    ret = ''
    # Count the dots if we've been asked to do so.
    if self.params.output.count_dots:
      numSkinDots = self._count_skin_dots(self._source_atoms_sorted, bondedNeighborLists)
      if self.params.output.format != 'raw':
        ret += self._describe_run("program:","command:")
        ret += self._describe_selection_and_parameters(groupLabel, selectionName)

      nsel = len(self._source_atoms_sorted)
      if self.params.output.format == 'raw':
        ret += self._rawEnumerate("", nsel, self.params.output.compute_scores, False, numSkinDots, groupLabel)
      else:
        ret += self._enumerate("{} dots".format(selectionName), nsel, self.params.output.compute_scores, False, numSkinDots)

    else: # Not counting the dots

      # Check for various output format types.
      # We're not implementing O format or XV format, but we still allow raw and oneline
      if self.params.output.format == 'raw':
        ret += self._writeRawOutput(comparisonString,groupLabel, False)

      elif self.params.output.format == 'json':
        ret += self._writeRawOutput(comparisonString,groupLabel, True)

      elif self.params.output.format == 'oneline':
        ret += self._count_summary(intersectionName)

      elif self.params.output.format == 'kinemage': # Kinemage format
        ret += self._describe_run("@caption"," command:")
        if self.params.output.contact_summary:
          ret += self._count_summary(intersectionName)

        if self.params.output.add_group_line:
          if numModels > 1:
            # doing one of multiple models of an ensemble
            ret += "@group dominant {{{} M{}}} animate\n".format(groupLabel,modelIndex+1)
          else:
            ret += "@group dominant {{{}}}\n".format(groupLabel)

        ret += self._writeOutput("{} dots".format(selectionName), groupLabel)
        # Put the water after so they end up in the right group
        ret += self._phantomHydrogenOutput

      else:
        raise ValueError("Unrecognized output format: "+self.params.output.format+" (internal error)")

    return ret

# ------------------------------------------------------------------------------

  def _clear_results(self):
    # Initialize the results to empty.
    self._results = {}
    for c in _allAtomClasses:
      interactionTypeDicts = {}
      for i in _interactionTypes:
        interactionTypeDicts[i] = []
      self._results[c] = interactionTypeDicts

# ------------------------------------------------------------------------------

  def _describe_run(self, header1, header2):
    '''
      Describe the command-line and other Phil options used for this run so that
      it could be reproduced.
      :param header1: Header for the first output line (the version/time line).
      :param header2: Header for the second output line (the command and its arguments).
      :return: String to be added to the output.
    '''

    global version
    ret = '{} probe2 v.{}, run {}\n'.format(header1, version, datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    ret += header2
    for a in sys.argv:
      ret += ' {}'.format(a)
    ret += '\n'

    return ret

# ------------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)
    if self.params.output.filename is None:
      # If the output file name is not specified, use the same root as the
      # input file and replace the suffix with .kin for Kinemage output or
      # .txt for raw or .json for JSON.
      suffix = '.kin'
      if self.params.output.format == 'json':
        suffix = '.json'
      elif self.params.output.format != 'kinemage' or self.params.output.count_dots:
        suffix = '.txt'
      inName = self.data_manager.get_default_model_name()
      p = Path(inName)
      self.params.output.filename = str(p.with_suffix(suffix))
      print('Setting output.filename Phil parameter to',self.params.output.filename)
    if self.params.source_selection is None:
      raise Sorry("Must specify a source parameter for approach "+self.params.approach)
    if self.params.approach in ['once','both'] and self.params.target_selection is None:
      raise Sorry("Must specify a target parameter for approach "+self.params.approach)
    aScale = self.params.atom_radius_scale
    if aScale < 0.0001 or aScale > 1000:
      raise Sorry("Invalid atom_radius_scale value: {:0g}".format(aScale))
    ao = self.params.atom_radius_offset
    if ao < -10 or ao > 1000:
      raise Sorry("Invalid atom_radius_offset value: {:0g}".format(ao))

    # Ensure consistency among parameters
    if self.params.probe.contact_cutoff < self.params.probe.probe_radius:
      self.params.probe.contact_cutoff = self.params.probe.probe_radius

    # Turn on profiling if we've been asked to in the Phil parameters
    if self.params.profile:
      import cProfile
      self._pr = cProfile.Profile()
      self._pr.enable()

# ------------------------------------------------------------------------------

  def overrideModel(self, model):
    '''This is a hack to let another program harness probe2 without having to write a
    new model file for it to read. After initializing probe2, but before calling
    run(), call this function to override the model that it should use.
    '''
    self.model = model

# ------------------------------------------------------------------------------

  def run(self):
    # String that will be output to the specified file.
    outString = ''

    if (self.params.output.add_kinemage_keyword and not self.params.output.count_dots
        and self.params.output.format == 'kinemage'):
      outString += '@kinemage 1\n'

    make_sub_header('Interpret Model', out=self.logger)

    # Allow the model to be overridden using the overrideModel() method.
    if not hasattr(self, 'model'):
      # Get our model.
      self.model = self.data_manager.get_model()

    ################################################################################
    # Get the bonding information we'll need to exclude our bonded neighbors.
    allAtoms = self.model.get_atoms()
    make_sub_header('Compute neighbor lists', out=self.logger)

    self.model.set_stop_for_unknowns(False)
    p = mmtbx.model.manager.get_default_pdb_interpretation_params()
    p.pdb_interpretation.use_neutron_distances = self.params.use_neutron_distances
    p.pdb_interpretation.allow_polymer_cross_special_position=True
    p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
    p.pdb_interpretation.proceed_with_excessive_length_bonds=True
    p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check=True
    # We need to turn this on because without it the interpretation is
    # renaming atoms to be more correct.  Unfortunately, this causes the
    # dot names to no longer match the input file.
    p.pdb_interpretation.flip_symmetric_amino_acids=False
    try:
      self.model.process(make_restraints=True, pdb_interpretation_params=p) # make restraints
      geometry = self.model.get_restraints_manager().geometry
      sites_cart = self.model.get_sites_cart() # cartesian coordinates
      bondProxies, asu = \
          geometry.get_all_bond_proxies(sites_cart = sites_cart)
    except Exception as e:
      try:
        # Fix up bogus unit cell when it occurs by checking crystal symmetry.
        self.model.add_crystal_symmetry_if_necessary()

        # Retry with the adjusted model
        self.model.process(make_restraints=True, pdb_interpretation_params=p) # make restraints
        geometry = self.model.get_restraints_manager().geometry
        sites_cart = self.model.get_sites_cart() # cartesian coordinates
        bondProxies, asu = \
            geometry.get_all_bond_proxies(sites_cart = sites_cart)

      except Exception as e:
        raise Sorry("Could not get bonding information for input file: " + str(e))

    ################################################################################
    # Get the bonding information we'll need to exclude our bonded neighbors.
    self._allBondedNeighborLists = Helpers.getBondedNeighborLists(allAtoms, bondProxies)

    ################################################################################
    # Get the extra atom information needed to score all of the atoms in the model.
    make_sub_header('Compute extra atom information', out=self.logger)
    ret = Helpers.getExtraAtomInfo(model = self.model,
      bondedNeighborLists = self._allBondedNeighborLists,
      useNeutronDistances = self.params.use_neutron_distances,
      probePhil = self.params.probe)
    self._extraAtomInfo = ret.extraAtomInfo
    if len(ret.warnings) > 0:
      print('Warnings returned by getExtraAtomInfo():\n'+ret.warnings, file=self.logger)

    # Scale and offset the radius values for all atoms based on our command-line arguments.
    for a in allAtoms:
      ei = self._extraAtomInfo.getMappingFor(a)
      ei.vdwRadius = self._scaled_atom_radius(a)
      self._extraAtomInfo.setMappingFor(a, ei)

    ################################################################################
    # Find the maximum VDW radius of any of our atoms, used to limit searches for nearby
    # atoms.
    self._maximumVDWRadius = 1
    for a in allAtoms:
      self._maximumVDWRadius = max(self._maximumVDWRadius, self._extraAtomInfo.getMappingFor(a).vdwRadius)

    ################################################################################
    # Get the extra atom information needed to sort all of the atoms in the model
    # into proper classes for reporting.  These classes may be atom names, when we're
    # sorting by atoms and it can be nucleic acid base names when we're sorting by that.
    # Comes from newAtom() and dotType() functions in probe.c.
    # Rather than a table indexed by type, we directly write the result.
    # Handle all atoms, not only selected atoms.
    self._atomClasses = {}
    for a in allAtoms:
      if not a.element_is_hydrogen():
        # All elements except hydrogen use their own names.
        self._atomClasses[a] = self._atom_class_for(a)
      else:
        # For hydrogen, assign based on what it is bonded to.
        if len(self._allBondedNeighborLists[a]) < 1:
          raise Sorry("Found Hydrogen with no neigbors: " + self._describe_atom_for_debug(a))
        else:
          self._atomClasses[a] = self._atom_class_for(self._allBondedNeighborLists[a][0])

    ################################################################################
    # Get the other characteristics we need to know about each atom to do our work.
    make_sub_header('Getting extra atom characteristics', out=self.logger)
    self._inWater = {}
    self._inHet = {}
    self._inMainChain = {}
    self._inSideChain = {}
    hetatm_sel = self.model.selection("hetatm")
    mainchain_sel = self.model.selection("backbone")  # Will NOT include Hydrogen atoms on the main chain
    sidechain_sel = self.model.selection("sidechain") # Will include Hydrogen atoms on the side chain
    for a in allAtoms:
      self._inWater[a] = common_residue_names_get_class(name=a.parent().resname) == "common_water"
      self._inHet[a] = hetatm_sel[a.i_seq]
      if not a.element_is_hydrogen():
        self._inMainChain[a] = mainchain_sel[a.i_seq]
      else:
        # Check our bonded neighbor to see if it is on the mainchain if we are a Hydrogen
        if len(self._allBondedNeighborLists[a]) < 1:
          raise Sorry("Found Hydrogen with no neigbors: " + self._describe_atom_for_debug(a))
        else:
          self._inMainChain[a] = mainchain_sel[self._allBondedNeighborLists[a][0].i_seq]
      self._inSideChain[a] = sidechain_sel[a.i_seq]

    ################################################################################
    # Ensure that the model we've been passed has at least one Hydrogen bonded to a Carbon
    # and at least one polar Hydrogen (bonded to N, O, or S).  Otherwise, raise a Sorry.
    if not self.params.ignore_lack_of_explicit_hydrogens:
      foundCBonded = False
      foundPolar = False
      for a in allAtoms:
        if Helpers.isPolarHydrogen(a, self._allBondedNeighborLists):
          foundPolar = True
        elif a.element_is_hydrogen():
          if len(self._allBondedNeighborLists[a]) < 1:
            raise Sorry("Found Hydrogen with no neigbors: " + self._describe_atom_for_debug(a))
          else:
            neighbor = self._allBondedNeighborLists[a][0]
            if neighbor.element == 'C':
              foundCBonded = True
      if not (foundCBonded and foundPolar):
        raise Sorry("Did not find both polar and non-polar Hydrogens in model.  For proper operation, "+
                    "Probe requires explicit Hydrogens.  Run Reduce2 or another placement "+
                    "program on the model before running Probe, or else add the Phil "+
                    "parameter ignore_lack_of_explicit_hydrogens=True.")

    ################################################################################
    # Get the source selection (and target selection if there is one).  These will be
    # lists of atoms that are in each selection, a subset of the atoms in the model.
    # If there is no model_id in the selection criteria, these may include atoms from
    # multiple models in the hierarchy.
    make_sub_header('Getting atom selections', out=self.logger)
    source_sel = self.model.selection(self.params.source_selection)
    allSourceAtoms = set()
    for a in allAtoms:
      if source_sel[a.i_seq]:
        allSourceAtoms.add(a)

    allTargetAtoms = set()
    if self.params.target_selection is not None:
      # If the target selection is "=", that means that it should be the same as the source selection.
      if self.params.target_selection == "=":
        allTargetAtoms = allSourceAtoms
      else:
        target_sel = self.model.selection(self.params.target_selection)
        for a in allAtoms:
          if target_sel[a.i_seq]:
            allTargetAtoms.add(a)

    ################################################################################
    # We usually have the selection pick a model, but in the case of SELFINTERSECT with one
    # input file and no model specified in the source and target patterns, we loop over all
    # models in the file.
    # We get lists of all atoms present in each hierarchy model that we're running.
    # This is a list of one when only one is selected and it is all of the available ones
    # when no particular one is selected.
    make_sub_header('Getting atom lists', out=self.logger)
    atomLists = [ self.model.get_atoms() ]
    if (self.params.approach == 'self' and
        (self.params.source_selection is None or 'model' not in self.params.source_selection) and
        (self.params.target_selection is None or 'model' not in self.params.target_selection)):
      # Handle the multiple-model case by looping modelID over all models.
      numModels = self.model.get_hierarchy().models_size()
      atomLists = []
      for i in range(numModels):
        atomLists.append( self.model.get_hierarchy().models()[i].atoms() )

    make_sub_header('Processing atom lists', out=self.logger)
    for modelIndex, atoms in enumerate(atomLists):

      ################################################################################
      # Get the subset of the source selection and target selection for this hierarchy
      # model.

      # Make an acceleration structure for determining whether an atom is in the ones we
      # are considering in the current list. This is a set rather than a list to make it
      # rapid to determine membership.
      atomsInThisModel = set(atoms)
      source_atoms = set()
      for a in allSourceAtoms:
        if a in atomsInThisModel:
          source_atoms.add(a)

      target_atoms = set()
      for a in allTargetAtoms:
        if a in atomsInThisModel:
          target_atoms.add(a)

      ###########################
      # Helper utility function to sort atoms consistently from run to run so that we get
      # the same ordering on coarse angles.
      def atomID(a):
        # Return the ID of the atom, which includes its chain, residue name,
        # residue number, atom name, and alternate separated by spaces. This
        # is used to sort the atoms. This must work in the case where we have
        # test atoms that are not completely fleshed out.
        try:
          return ( a.parent().parent().parent().id + a.parent().resname +
            str(a.parent().parent().resseq_as_int()) + a.name + a.parent().altloc )
        except Exception:
          return ""
      #
      ###########################

      ################################################################################
      # Find a list of all of the selected atoms with no duplicates
      # Get the bonded neighbor lists for the atoms that are in this selection.
      # We have to do this so that when keep_unselected_atoms is set to False we don't
      # follow bonds to neighbor atoms that should not exist.
      # Sort the atoms by an ID that is consistent from run to run so that they end up
      # in our data structures in the same order for each run.

      make_sub_header('Sorting atoms', out=self.logger)
      all_selected_atoms = sorted(source_atoms.union(target_atoms), key=lambda x:atomID(x))
      make_sub_header('Getting bonded-neighbor lists', out=self.logger)
      bondedNeighborLists = Helpers.getBondedNeighborLists(all_selected_atoms, bondProxies)

      ################################################################################
      # Build a spatial-query structure that tells which atoms are nearby.
      # Include all atoms in the structure, not just the ones that have been selected,
      # unless we've been asked not to keep them.
      make_sub_header('Make spatial-query accelerator', out=self.logger)
      if self.params.keep_unselected_atoms:
        self._spatialQuery = Helpers.createSpatialQuery(atoms, self.params.probe)
        # Replace the bonded-neighbor list with all bonded neighbors, even ones that
        # are not selected, so that they will block dots that overlap with bonded atoms.
        bondedNeighborLists = self._allBondedNeighborLists
        selectedAtomsIncludingKept = atoms
      else:
        self._spatialQuery = Helpers.createSpatialQuery(all_selected_atoms, self.params.probe)
        selectedAtomsIncludingKept = all_selected_atoms

      ################################################################################
      # Add Phantom hydrogens to waters and mark
      # the water oxygens as not being donors in atoms that are in the source or target selection.
      # Also clear the donor status of all N, O, S atoms because we have explicit hydrogen donors.
      self._phantomHydrogenOutput = ""
      phantomHydrogens = []
      make_sub_header('Adjusting for explicit hydrogens', out=self.logger)
      if self.params.output.record_added_hydrogens:
        self._phantomHydrogenOutput += "@master {water H?}\n"
        self._phantomHydrogenOutput += '@vectorlist {water H?} color= gray master={water H?}\n'

      # @todo Look up the radius of a water Hydrogen.  This may require constructing a model with
      # a single water in it and asking about the hydrogen radius.  This could also become a
      # Phil parameter.  Also look up the OH bond distance rather than hard-coding it here.
      phantomHydrogenRadius = 1.05
      placedHydrogenDistance = 0.84
      if self.params.use_neutron_distances:
        phantomHydrogenRadius = 1.0
        placedHydrogenDistance = 0.98

      adjustedHydrogenRadius = self.params.atom_radius_offset + (phantomHydrogenRadius * self.params.atom_radius_scale)

      # Check all selected atoms to see if we need to add Phantom Hydrogens to them.
      # Don't add Phantom Hydrogens to atoms that are not selected, even if they are kept.
      maxISeq = Helpers.getMaxISeq(self.model)
      for a in all_selected_atoms:

        # Ignore Hydrogens whose parameters are out of bounds.
        if a.element_is_hydrogen():
          # In the original code, this looks at H atoms with parent N,O,S atoms
          # and marks them as donors.  This is handled for us below in the call
          # to Helpers.fixupExplicitDonors().

          # If we are in a water, make sure our occupancy and temperature (b) factor are acceptable.
          # If they are not, set the class for the atom to 'ignore'.
          # This handles the case where there were explicit Hydrogens on waters and so
          # we won't add Phantom Hydrogens.
          if self._inWater[a] and (a.occ < self.params.minimum_water_hydrogen_occupancy or
              a.b > self.params.maximum_water_hydrogen_b):
            self._atomClasses[a] = 'ignore'

        # If we are the Oxygen in a water, then add phantom hydrogens pointing towards nearby acceptors
        elif self._inWater[a] and a.element == 'O':
          # We're an acceptor and not a donor.
          # @todo Original Probe code only cleared the donor status if it found a bonded
          # Hydrogen in the same conformation whose occupancy was > 0.1.  Here, we're turning
          # it off regardless of the occupancy.
          ei = self._extraAtomInfo.getMappingFor(a)
          ei.isDonor = False
          ei.isAcceptor = True
          self._extraAtomInfo.setMappingFor(a, ei)

          # If we don't yet have Hydrogens attached, add phantom hydrogen(s)
          if len(bondedNeighborLists[a]) == 0:
            # NOTE: The Phantoms have i_seq numbers that are sequential and that are higher than
            # all other atoms in the model.  Each one has a unique i_seq.
            newPhantoms = Helpers.getPhantomHydrogensFor(maxISeq, a, self._spatialQuery, self._extraAtomInfo,
                            0.0, True, adjustedHydrogenRadius, placedHydrogenDistance)
            maxISeq += len(newPhantoms)
            for p in newPhantoms:

              # Put in our list of Phantom Hydrogens
              phantomHydrogens.append(p)

              # Add the atom to the general spatial-query data structure
              self._spatialQuery.add(p)

              # Set the extra atom information for this atom
              ei = probeExt.ExtraAtomInfo(adjustedHydrogenRadius, False, True, True)
              self._extraAtomInfo.setMappingFor(p, ei)

              # Set the atomClass and other data based on the parent Oxygen.
              self._atomClasses[p] = self._atom_class_for(a)
              self._inWater[p] = self._inWater[a]
              self._inMainChain[p] = self._inMainChain[a]
              self._inSideChain[p] = self._inSideChain[a]
              self._inHet[p] = self._inHet[a]

              # Mark the Phantom Hydrogens as being bonded to their Oxygen so that
              # dots on a Phantom Hydrogen within its Oxygen will be excluded.
              bondedNeighborLists[p] = [a]

              # It was thought that in the future, we may add these bonds, but that will cause the
              # Phantom Hydrogens to mask their water Oxygens from close contacts or
              # clashes with the acceptors, which is a change in behavior from the
              # original Probe and would have the undesirable effect of a potential
              # Hydrogen hiding a true collision.
              # Not marking these as bonded requires special-case handling
              # of Phantom Hydrogen interactions in the dot-scoring code.
              # This means that we have a one-way bond, which is unusual but suits our
              # purposes.
              # Not done: bondedNeighborLists[a].append(p)

              # Add the new atom to any selections that the old atom was in.
              if a in source_atoms:
                source_atoms.add(p)
              if a in target_atoms:
                target_atoms.add(p)

              # Report on the creation if we've been asked to
              if self.params.output.record_added_hydrogens:

                resName = a.parent().resname.strip().upper()
                resID = str(a.parent().parent().resseq_as_int())
                chainID = a.parent().parent().parent().id
                iCode = a.parent().parent().icode
                alt = a.parent().altloc
                self._phantomHydrogenOutput += '{{{:4.4s}{:1s}{:>3s}{:>2s}{:>4s}{:1s}}}P {:8.3f}{:8.3f}{:8.3f}\n'.format(
                  a.name, alt, resName, chainID, resID, iCode,
                  a.xyz[0], a.xyz[1], a.xyz[2])

                resName = p.parent().resname.strip().upper()
                resID = str(p.parent().parent().resseq_as_int())
                chainID = p.parent().parent().parent().id
                iCode = p.parent().parent().icode
                alt = p.parent().altloc
                self._phantomHydrogenOutput += '{{{:4.4s}{:1s}{:>3s}{:>2s}{:>4s}{:1s}}}L {:8.3f}{:8.3f}{:8.3f}\n'.format(
                  p.name, alt, resName, chainID, resID, iCode,
                  p.xyz[0], p.xyz[1], p.xyz[2])

      # Fix up the donor status for all of the atoms now that we've added the final explicit
      # Phantom Hydrogens.
      Helpers.fixupExplicitDonors(selectedAtomsIncludingKept, bondedNeighborLists, self._extraAtomInfo)

      ################################################################################
      # Add ionic bonds to the bonded-neighbor list so that we won't count interactions
      # between two atoms that are both bonded to the same ion (such as Nitrogens on
      # Histidine rings around Cu or Zn).  Do this after we've added the Phantom Hydrogens
      # so that we don't see ionic bonds in the Phantom-Hydrogen addition code checks.
      Helpers.addIonicBonds(bondedNeighborLists, selectedAtomsIncludingKept, self._spatialQuery, self._extraAtomInfo)

      # Make a query structure to return the Phantom Hydrogens (if there are any)
      self._phantomHydrogensSpatialQuery = Helpers.createSpatialQuery(phantomHydrogens, self.params.probe)

      ################################################################################
      # Re-fill all_selected_atoms
      all_selected_atoms = sorted(source_atoms.union(target_atoms), key=lambda x:atomID(x))

      ################################################################################
      # Get the dot sets we will need for each atom.  This is the set of offsets from the
      # atom center where dots should be placed.  We use a cache to reduce the calculation
      # time by returning the same answer for atoms that have the same radius.
      # This must be done after we've added all Phantom Hydrogens and adjusted all of
      # the ExtraAtomInfo.
      dotCache = Helpers.createDotSphereCache(self.params.probe)
      self._dots = {}
      for a in all_selected_atoms:
        self._dots[a] = dotCache.get_sphere(self._extraAtomInfo.getMappingFor(a).vdwRadius).dots()

      ################################################################################
      # Construct a DotScorer object.  This must be done after we've added all Phantom
      # Hydrogens and adjusted all of the ExtraAtomInfo.
      make_sub_header('Make dot scorer', out=self.logger)
      self._dotScorer = Helpers.createDotScorer(self._extraAtomInfo, self.params.probe)

      ################################################################################
      # Sums of interaction types of dots based on whether their source and/or target
      # were mainchain, sidechain, both, or neither.  There is another place to store
      # the sum of multiple passes.
      # Each contains an entry for each InteractionType and for the total.
      self._clear_results();
      self._MCMCCount = {}
      self._SCSCCount = {}
      self._MCSCCount = {}
      self._otherCount = {}
      for t in _interactionTypes:
        self._MCMCCount[t] = 0
        self._SCSCCount[t] = 0
        self._MCSCCount[t] = 0
        self._otherCount[t] = 0

      ################################################################################
      # Generate sorted lists of the selected atoms, so that we run them in the same order
      # they appear in the model file.  This will group phantom hydrogens with the oxygens
      # they are associated with because they share the same sequence ID.
      # We add the location to the sorting criteria because the phantom hydrogens have the
      # same sequence ID as their parent O and as each other.
      self._source_atoms_sorted = sorted(source_atoms, key=lambda atom: "{} {:.3f} {:.3f} {:.3f}".format(
        atom.i_seq, atom.xyz[0], atom.xyz[1], atom.xyz[2]))
      self._target_atoms_sorted = sorted(target_atoms, key=lambda atom:  "{} {:.3f} {:.3f} {:.3f}".format(
        atom.i_seq, atom.xyz[0], atom.xyz[1], atom.xyz[2]))

      ################################################################################
      # Find our group label
      if self.params.output.format in ['raw','json']:
        groupLabel = ""
      else:
        groupLabel = "dots"
      if len(self.params.output.group_label) > 0:
        groupLabel = self.params.output.group_label

      ################################################################################
      # Do the calculations; which one depends on the approach and other phil parameters.
      # Append the information to the string that will be written to file.

      if self.params.approach == 'count_atoms':
        make_sub_header('Counting atoms', out=self.logger)
        # Report the number of atoms in the source selection
        outString += 'atoms selected: '+str(len(self._source_atoms_sorted))+'\n'

      elif self.params.approach == 'surface':
        make_sub_header('Find surface dots', out=self.logger)

        # Store constants used frequently
        minimum_occupancy = self.params.minimum_occupancy
        include_water_water = self.params.include_water_water

        # Produce dots on the surfaces of the selected atoms.
        maxRadius = 2*self._maximumVDWRadius + 2 * self.params.probe.probe_radius
        for src in self._source_atoms_sorted:
          srcInWater = self._inWater[src]
          srcModel = src.parent().parent().parent().parent().id

          # Find nearby atoms that might come into contact.  This greatly speeds up the
          # search for touching atoms.
          maxRadius = (self._extraAtomInfo.getMappingFor(src).vdwRadius + self._maximumVDWRadius +
            2 * self.params.probe.probe_radius)
          nearby = self._spatialQuery.neighbors(src.xyz, 0.00001, maxRadius)

          # Select those that are actually within the contact distance based on their
          # particular radius.  Only accept atoms that are in compatible conformations.
          atomList = []
          for n in nearby:
            nInWater = self._inWater[n]
            nModel = n.parent().parent().parent().parent().id

            # Skip atoms that are marked to be ignored
            if self._atomClasses[n] == 'ignore':
              continue
            # Skip water-water interactions unless they are between atoms in the same residue
            elif (not include_water_water) and srcInWater and nInWater and (src.parent() != n.parent()):
              continue
            # Skip atoms that are in non-compatible alternate conformations
            elif not Helpers.compatibleConformations(src, n):
              continue
            # Skip atoms that are in different models.
            elif srcModel != nModel:
              continue
            d = (Helpers.rvec3(n.xyz) - Helpers.rvec3(src.xyz)).length()
            if (d <= self._extraAtomInfo.getMappingFor(n).vdwRadius +
                self._extraAtomInfo.getMappingFor(src).vdwRadius + 2*self.params.probe.probe_radius):
              atomList.append(n)

          # Find out what class of dot we should place for this atom.
          atomClass = self._atomClasses[src]

          # Generate all of the dots for this atom.
          self._generate_surface_dots_for(src, atomList)

        # Count the dots if we've been asked to do so.
        if self.params.output.count_dots:
          numSkinDots = self._count_skin_dots(self._source_atoms_sorted, bondedNeighborLists)
          if self.params.output.format != 'raw':
            outString += self._describe_selection_and_parameters(groupLabel, "external")

          nsel = len(self._source_atoms_sorted)
          if self.params.output.format == 'raw':
            outString += self._rawEnumerate("", nsel, False, True, numSkinDots, groupLabel)
          else:
            outString += self._describe_run("program:","command:")
            outString += self._enumerate("extern dots", nsel, False, True, numSkinDots)

        # Otherwise, produce the dots as output
        else:
          # Check for various output format types other than Kinemage.
          # We're not implementing O format or XV format, but we still allow raw and oneline
          if self.params.output.format == 'raw':
            outString += self._writeRawOutput("1->none",groupLabel)

          elif self.params.output.format == 'oneline':
            # Do nothing for this mode when computing the surface
            pass

          elif self.params.output.format == 'kinemage': # Kinemage format
            outString += self._describe_run("@caption"," command:")
            masterName = "dots"
            if len(self.params.output.group_name) > 0:
              masterName = self.params.output.group_name

            if self.params.output.add_group_line:
              outString += "@group dominant {{{}}}\n".format(masterName)

            outString += self._writeOutput("extern dots", masterName)

          else:
            raise ValueError("Unrecognized output format: "+self.params.output.format+" (internal error)")

      elif self.params.approach == 'self':
        make_sub_header('Find self-intersection dots', out=self.logger)

        # Generate dots for the source atom set against itself.
        self._generate_interaction_dots(self._source_atoms_sorted, source_atoms,
          self._spatialQuery, self._phantomHydrogensSpatialQuery, bondedNeighborLists)

        # Generate our report
        outString += self._report_single_interaction(groupLabel, "self", "1->1", "SelfIntersect",
            len(atomLists), modelIndex, bondedNeighborLists)

      elif self.params.approach == 'once':
        make_sub_header('Find single-direction intersection dots', out=self.logger)

        # Generate dots for the source atom set against the target atom set.
        self._generate_interaction_dots(self._source_atoms_sorted, target_atoms,
          self._spatialQuery, self._phantomHydrogensSpatialQuery, bondedNeighborLists)

        # Generate our report
        outString += self._report_single_interaction(groupLabel, "once", "1->2", "IntersectOnce",
            len(atomLists), modelIndex, bondedNeighborLists)

      elif self.params.approach == 'both':
        make_sub_header('Find both-directions intersection dots', out=self.logger)

        # @todo The code below here is similar to -once but is repeated twice and has different string values.
        # It is also somewhat re-ordered in terms of where the selection is printed.  This keeps us from
        # re-using _report_single_interaction() directly without generalizing it.

        # Preliminary information before running both intersections.
        if self.params.output.count_dots:
          if self.params.output.format != 'raw':
            outString += self._describe_run("program:","command:")
            outString += self._describe_selection_and_parameters(groupLabel, "once")
        else: # Not counting the dots
          if self.params.output.format == 'raw':
            pass
          elif self.params.output.format == 'kinemage':
            outString += self._describe_run("@caption"," command:")
            if self.params.output.add_group_line:
              outString += "@group {{{}}}\n".format(groupLabel)

        # =================== First direction ========================

        # Generate dots for the source atom set against the target atom set.
        self._generate_interaction_dots(self._source_atoms_sorted, target_atoms,
          self._spatialQuery, self._phantomHydrogensSpatialQuery, bondedNeighborLists)

        # Count the dots if we've been asked to do so.
        if self.params.output.count_dots:
          numSkinDots = self._count_skin_dots(self._source_atoms_sorted, bondedNeighborLists)
          nsel = len(self._source_atoms_sorted)
          if self.params.output.format == 'raw':
            outString += self._rawEnumerate("1->2", nsel, self.params.output.compute_scores, False, numSkinDots, groupLabel)
          else:
            outString += self._enumerate("1->2", nsel, self.params.output.compute_scores, False, numSkinDots)

        else: # Not counting the dots

          # Check for various output format types.
          # We're not implementing O format or XV format, but we still allow raw and oneline
          if self.params.output.format == 'raw':
            outString += self._writeRawOutput("1->2",groupLabel)

          elif self.params.output.format == 'oneline':
            # Acculumlate but do not report results
            outString += self._count_summary("IntersectBothWays 1->2", False)

          elif self.params.output.format == 'kinemage': # Kinemage format
            outString += self._writeOutput("1->2", groupLabel)
            if self.params.output.contact_summary:
              # Acculumlate but do not report results
              outString += self._count_summary("IntersectBothWays 1->2", False)

        # =================== Second direction ========================

        # Clear the results before running interactions the other direction.
        self._clear_results();

        # Generate dots for the target atom set against the source atom set.
        self._generate_interaction_dots(self._target_atoms_sorted, source_atoms,
          self._spatialQuery, self._phantomHydrogensSpatialQuery, bondedNeighborLists)

        # Count the dots if we've been asked to do so.
        if self.params.output.count_dots:
          numSkinDots = self._count_skin_dots(self._target_atoms_sorted, bondedNeighborLists)
          nsel = len(self._target_atoms_sorted)
          if self.params.output.format == 'raw':
            outString += self._rawEnumerate("2->1", nsel, self.params.output.compute_scores, False, numSkinDots, groupLabel)
          else:
            outString += self._enumerate("2->1", nsel, self.params.output.compute_scores, False, numSkinDots)

        else: # Not counting the dots

          # Check for various output format types.
          # We're not implementing O format or XV format, but we still allow raw and oneline
          if self.params.output.format == 'raw':
            outString += self._writeRawOutput("2->1",groupLabel)

          elif self.params.output.format == 'oneline':
            # Accumulate and report results
            outString += self._count_summary("IntersectBothWays 2->1", True)

          elif self.params.output.format == 'kinemage': # Kinemage format
            outString += self._writeOutput("2->1", groupLabel)
            if self.params.output.contact_summary:
              # Accumulate and report results
              outString += self._count_summary("IntersectBothWays 2->1", True)

          else:
            raise ValueError("Unrecognized output format: "+self.params.output.format+" (internal error)")

    # Write the output to the specified file.
    if self.params.output.write_files:
      self.data_manager._write_text("Text", outString, self.params.output.filename)

    # If we have a dump file specified, write the atom information into it.
    # We write it at the end because the extra atom info may have been adjusted
    # during the code that handles hydrogen adjustements.
    if self.params.output.write_files and self.params.output.dump_file_name is not None:
      atomDump = Helpers.writeAtomInfoToString(allAtoms, self._extraAtomInfo)
      with open(self.params.output.dump_file_name,"w") as df:
        df.write(atomDump)

    # Report profiling info if we've been asked to in the Phil parameters
    if self.params.profile:
      print('Profile results:')
      import pstats
      profile_params = {'sort_by': 'time', 'num_entries': 20}
      self._pr.disable()
      ps = pstats.Stats(self._pr).sort_stats(profile_params['sort_by'])
      ps.print_stats(profile_params['num_entries'])

    # Return the results object that has all of the dots and the output string.
    return self._results, outString

# ------------------------------------------------------------------------------

  #def get_results(self):
  #  return group_args(model = self.model)

