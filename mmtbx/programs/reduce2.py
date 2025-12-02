"""Run reduce (add hydrogens). version 2"""
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
# See the License for the specific language governing permissionsand
# limitations under the License.

from __future__ import absolute_import, division, print_function
import sys
import os
import time
from datetime import datetime
from libtbx.program_template import ProgramTemplate
from libtbx import group_args, phil
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry, null_out
import mmtbx
from mmtbx.probe import Helpers
from iotbx import pdb, cli_parser
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.reduce import Optimizers
from scitbx.array_family import flex
from iotbx.pdb import common_residue_names_get_class, amino_acid_codes, nucleic_acid_codes
from mmtbx.programs import probe2
import copy
from iotbx.data_manager import DataManager
import csv

version = "2.14.0"

master_phil_str = '''
approach = *add remove
  .type = choice
  .short_caption = Add or remove Hydrogens
  .help = Determines whether Reduce will add (and optimize) or remove Hydrogens from the model
keep_existing_H = False
  .type = bool
  .short_caption = Do not remove Hydrogens in the original model
  .help = Keep existing H atoms in the model
n_terminal_charge = *residue_one first_in_chain no_charge
  .type = choice(multi=False)
  .short_caption = N terminal charge approach
  .help = Mode for placing H3 at terminal nitrogen.
use_neutron_distances = False
  .type = bool
  .short_caption = Use neutron distances
  .help = Use neutron distances (-nuclear in reduce)
preference_magnitude = 1.0
  .type = float
  .short_caption = Rotational-preference magnitude
  .help = Multiplier on the rotational-preference energy for rotatable Movers (-penalty in reduce)
alt_id = None
  .type = str
  .short_caption = Alternate to optimize
  .help = Alternate to optimize.  The default is to optimize all of them.
model_id = None
  .type = int
  .short_caption = Model ID to optimize
  .help = Model ID to optimize.  The default is to optimize all of them.  If one is selected, the others are removed from the output file.
add_flip_movers = False
  .type = bool
  .short_caption = Add flip movers
  .help = Insert flip movers (-flip, -build, -noflip, -demandflipallnhqs in reduce)
non_flip_preference = 0.5
  .type = float
  .short_caption = Preference to not flip
  .help = For flip movers, only do the flip if the score in the flipped orientation is this much better.
skip_bond_fix_up = False
  .type = bool
  .short_caption = Skip fixup step for Movers that are flips
  .help = For debugging purposes, it can be useful to only do flips with no bond fix-up to compare scores. This fixup is done to make the bond angles within expected ranges for structures that are going to be deposited in the PDB. If further refinement is to be done, then it may also be useful to set this to True.
set_flip_states = None
  .type = str
  .short_caption = Comma-separated list of flip Mover states to set
  .help = String with comma-separated entries. Each entry has the form without the single quotes "'1 . A HIS 11H Flipped AnglesAdjusted'". These are space-separated values. The first word is the model number, starting with 1. The second is the lower-case alternate, or '.' for all alternates -- also use this when there are no alternates in the file. The third is the chain ID. The fourth is the residue name. The fifth is the residue id, which may include an insertion code as its last character. The sixth is either Flipped or Unflipped. If it is Flipped, then another word is added -- AnglesAdjusted or AnglesNotAdjusted, specifying whether to do the three-point dock to adjust the bond angles after the flip. An example with several entries is: again, no quotes are included: "'1 a A HIS 11H Unflipped,1 b A ASN 15 Flipped AnglesNotAdjusted,1 . B GLN 27 Flipped AnglesAdjusted'". Any Flip Movers that would be placed at the specified location are instead locked in the specified configuration.
profile = False
  .type = bool
  .short_caption = Profile the entire run
  .help = Profile the performance of the entire run
comparison_file = None
  .type = str
  .short_caption = Compare the Mover scores from this run with those in comparison_file
  .help = Points to a comparison_file that is the result of running Hydrogenate or Reduce or Reduce2 or some other hydrogen-placement program. The Probe2 scores for the Movers found in the current run are compared against the scores for comparison_file and stored in a file with the same name as output.file_name with _comparison.csv appended. If None, no comparison is done.
verbosity = 2
  .type = int
  .short_caption = Level of detail in description file
  .help = Level of detail in description file.
bonded_neighbor_depth = 4
  .type = int
  .short_caption = How many neighbors to consider bonded (>3 accepted only if hydrogen)
  .help = When looking for interactions between atoms, this specifies how many hops we should take when looking for atoms that are excluded because we share a common chain of bonds. Lengths more than 3 only happen when one end is a hydrogen. This value should not be changed from the default except for regression tests against earlier versions.
stop_on_any_missing_hydrogen = False
  .type = bool
  .short_caption = Emit a Sorry and stop when any hydrogen in the model has insufficient restraints.
  .help = Emit a Sorry and stop when any hydrogen in the model has insufficient restraints. It will always emit when there are one or more residues without restraints.
ignore_missing_restraints = False
  .type = bool
  .short_caption = Don't stop if restraints for a residue are missing.
  .help = Don't stop if restraints for a residue are missing.
output
  .style = menu_item auto_align
{
  write_files = True
    .type = bool
    .short_caption = Write the output files
    .help = Write the output files(s) when this is True (default). Set to False when harnessing the program.
  description_file_name = None
    .type = str
    .short_caption = Description output file name
    .help = This file holds a description of the operations that Reduce2 performed when placing and optimize hydrogens. Its creation is required because it can contain important warnings that must be attended to by the person running the program. It also includes a REPORT section that lists the final orientation along with initial and final score for each Mover that can be parsed to determine the results. The REPORT information and some additional information used to be included in the resulting PDB file produced by Reduce.
  flipkin_directory = None
    .type = str
    .short_caption = Where to place the Flipkin Kinemages
    .help = Where to place the Flipkin Kinemages. If None, no Flipkin files are made.
  clique_outline_file_name = None
    .type = str
    .short_caption = Where to save a Kinemage showing all positions for atoms in each Mover in each clique
    .help = Where to save a Kinemage showing all positions for atoms in each Mover for each clique. This enable exploration of the reasons for the cliques. Each clique has the atoms expanded by the probe radius and as each is turned off, the Movers inside are revealed. Clicking on each Mover shows its description. If None, no such Kinemage is made.
  print_atom_info = False
    .type = bool
    .short_caption = Print extra atom info
    .help = Print extra atom info
}
''' + Helpers.probe_phil_parameters.replace("bump_weight = 10.0", "bump_weight = 100.0").replace("hydrogen_bond_weight = 4.0", "hydrogen_bond_weight = 40.0")
# @todo We replace the default weights to avoid the issue described in
# https://github.com/cctbx/cctbx_project/issues/1072 until we can figure out the appropriate
# fix.

program_citations = phil.parse('''
citation {
  authors = Word, et. al.
  journal = J. Mol. Biol.
  volume = 285
  pages = 1735-1747
  year = 1999
  external = True
}
''')

# ------------------------------------------------------------------------------

class _MoverLocation(object):
  # Holds information needed to identify a Mover within a model file.
  def __init__(self, moverType, modelId, altId, chain, resName, resIdWithICode):
    self.moverType = moverType  # String indicating Mover type: AmideFlip or HisFlip
    self.modelId = modelId      # Integer indicating the modelId that the entry corresponds to
    self.altId = altId          # String indicating the altId that the entry corresponds to
    self.chain = chain          # String Chain that the residue is in
    self.resName = resName      # String Name of the residue
    try:                        # String holding the integer ID of the residue and any insertion code
      self.resId = int(resIdWithICode)
      self.iCode = ''
    except Exception:
      self.resId = int(resIdWithICode[:-1])
      self.iCode = resIdWithICode[-1]

  def __str__(self):
    return "{} {} '{}' {} {} {}{}".format(self.moverType, self.modelId, self.altId,
      self.chain, self.resName, self.resId, self.iCode)
  def __repr__(self):
      return "_MoverLocation({})".format(str(self))


def _FindMoversInOutputString(s, moverTypes = ['SingleHydrogenRotator',
      'NH3Rotator', 'AromaticMethylRotator', 'AmideFlip', 'HisFlip']):
  '''Return a list of _MoverLocation items that include all of them found in the
  output string from an Optimizer.
  :param s: String returned from the getInfo() method on an optimizer.
  :param moverTypes: List of names for the movertypes to find.
  :return: List of _MoverLocation objects indicating all flip Movers listed inside
  any BEGIN...END REPORT block in the string, including whether they were marked as
  flipped.
  '''
  ret = []
  modelId = None
  altId = None
  inBlock = False
  for line in s.splitlines():
    words = line.split()
    if inBlock:
      if words[0:2] == ['END','REPORT']:
        inBlock = False
      elif words[0] in moverTypes:
        ret.append(_MoverLocation(words[0], modelId, altId, words[3], words[4], int(words[5])))
    else:
      if words[0:2] == ['BEGIN','REPORT:']:
        modelId = int(words[3])
        # Remove the single-quote and colon characters from the AltId, possibly leaving it empty
        trim = words[5].replace("'", "")
        trim = trim.replace(":", "")
        altId = trim
        inBlock = True
  return ret


def _FindFlipsInOutputString(s, moverType):
  '''Return a list of Optimizers.FlipMoverState items that include all of them found in the
  output string from an Optimizer.
  :param s: String returned from the getInfo() method on an optimizer.
  :param moverType: Either AmideFlip or HisFlip, selects the flip type to report.
  :return: List of Optimizers.FlipMoverState objects indicating all flip Movers listed inside
  any BEGIN...END REPORT block in the string, including whether they were marked as
  flipped.
  '''
  ret = []
  modelId = None
  altId = None
  inBlock = False
  for line in s.splitlines():
    words = line.split()
    if inBlock:
      if words[0:2] == ['END','REPORT']:
        inBlock = False
      elif words[0] == moverType:
        ret.append(Optimizers.FlipMoverState(moverType, modelId, altId, words[3], words[4],
                                   words[5], words[14] == 'Flipped', words[15] == 'AnglesAdjusted') )
    else:
      if words[0:2] == ['BEGIN','REPORT:']:
        modelId = int(words[3])
        # Remove the single-quote and colon characters from the AltId, possibly leaving it empty
        trim = words[5].replace("'", "")
        trim = trim.replace(":", "")
        altId = trim
        inBlock = True
  return ret

def _IsMover(rg, movers):
  '''Is a residue group in the list of residues that have Movers?
  :param rg: residue group
  :param movers: List of _MoverLocation or Optimizers.FlipMoverState objects
  :return: True if the residue is in the list, False if not
  '''
  for m in movers:
    chain = rg.parent()
    modelId = chain.parent().id
    # We must offset the index of the model by 1 to get to the 1-based model ID
    if ( (modelId == m.modelId + 1 or modelId == '') and
         (chain.id == m.chain) and
         (rg.resseq_as_int() == m.resId) and
         (rg.icode.strip() == m.iCode.strip())
       ):
      return True
  return False

def _AltFromFlipOutput(fo):
  '''Return a string that describes the alternate from a record in _FindFlipsInOutputString()
  output.  Reports ' ' for an empty alternate. Returns lower-case representation.
  '''
  if fo.altId in ["", ""]:
    return ' '
  return fo.altId.lower()


def _AddPosition(a, tag, group, partner=None):
  '''Return a string that describes the point at or line to the specified atom or just
  a sphere location.
  This is used when building Kinemages.  Reports the alternate only if it is not empty.
  :param a: Atom to describe.
  :param tag: 'P' for point, 'L' for line, '' for sphere location.
  :param group: The dominant group name the point or line is part of.
  :param partner: Lines connect two atoms, and the alternate of either of them
  causes the line to be in that alternate.  This provides a way to tell about
  a second atom whose alternate should be checked as well.
  '''
  if len(tag) > 0:
    tagString = ' {}'.format(tag)
  else:
    tagString = ''
  if a.parent().altloc in ['', ' ']:
    altLoc = ''
    altTag = ''
  else:
    altLoc = a.parent().altloc.lower()
    altTag = " '{}'".format(a.parent().altloc.lower())
  if (partner is not None) and not (partner.parent().altloc in ['', ' ']):
    altLoc = partner.parent().altloc.lower()
    altTag = " '{}'".format(partner.parent().altloc.lower())
  return '{{{:4s}{:1s}{} {} {:3d} B{:.2f} {}}}{}{} {:.3f}, {:.3f}, {:.3f}'.format(
    a.name.lower(),                       # Atom name
    altLoc,                               # Alternate, if any
    a.parent().resname.strip().lower(),   # Residue name
    a.parent().parent().parent().id,      # chain
    a.parent().parent().resseq_as_int(),  # Residue number
    a.b,                                  # B factor
    group,                                # Dominant group name
    tagString,                            # Tag (P or L)
    altTag,                               # Tag for alternate, if any
    a.xyz[0],                             # location
    a.xyz[1],
    a.xyz[2]
  )


def _DescribeMainchainLink(a0s, a1s, group):
  '''Return a string that describes the links between two atom sets, each of
  which may contain multiple alternates. Only link ones that have compatible
  alternates.
  :param a0s: Alternates of first atom.
  :param a1s: Alternates of second atom.
  :param group: The dominant group name the point or line is part of.
  @todo When we have alternates that end at a common atom, that atom's
  alternative is used for the line, making it appear in both alts.
  '''
  ret = ''
  if (a0s is None) or (a1s is None):
    return ret
  for a0 in a0s:
    for a1 in a1s:
      if Helpers.compatibleConformations(a0, a1):
        ret += _AddPosition(a0, 'P', group) + ' ' + _AddPosition(a1, 'L', group, a0) + '\n'
  return ret


_amino_acid_resnames = sorted(amino_acid_codes.one_letter_given_three_letter.keys())
def _IsStandardResidue(resname):
  return resname.strip().upper() in _amino_acid_resnames


_nucleic_acid_resnames = set(nucleic_acid_codes.rna_one_letter_code_dict.keys()).union(
  set(nucleic_acid_codes.dna_one_letter_code_dict.keys()))
def _IsNucleicAcidResidue(resname):
  return resname.strip().upper() in _nucleic_acid_resnames


def _MainChainAtomsWithHydrogen(resname):
  # Find the main chain atoms with hydrogen bonds
  if _IsStandardResidue(resname):
    return ['N', 'CA', 'C', 'O']
  elif _IsNucleicAcidResidue(resname):
    return ["OP2", "OP3", "C5'", "C4'", "C3'", "C2'", "C1'", "O3'", "O2'"]
  return []


def _DescribeMainchainResidue(r, group, prevCs):
  '''Return a string that describes the mainchain for a specified residue.
  Add the point for the first mainchain atom in the previous residue
  (none for the first) and lines to the N, Ca, C, and O.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param prevCs: Atom(s) that is the mainchain last connection in all conformations
  of the previous residue. For protein, this is C and for nucleic acid, this is O3'.
  '''

  ret = ''

  ############################################################
  # Protein main chain
  if _IsStandardResidue(r.atom_groups()[0].resname):

    # Find the lists of atoms that we might need.
    Ns = [a for a in r.atoms() if a.name.strip().upper() == 'N']
    CAs = [a for a in r.atoms() if a.name.strip().upper() == 'CA']
    Cs = [a for a in r.atoms() if a.name.strip().upper() == 'C']
    Os = [a for a in r.atoms() if a.name.strip().upper() == 'O']
    OXTs = [a for a in r.atoms() if a.name.strip().upper() == 'OXT']

    # Make all of the connections.
    ret += _DescribeMainchainLink(prevCs, Ns, group)
    ret += _DescribeMainchainLink(Ns, CAs, group)
    ret += _DescribeMainchainLink(CAs, Cs, group)
    ret += _DescribeMainchainLink(Cs, Os, group)
    ret += _DescribeMainchainLink(Cs, OXTs, group)

  ############################################################
  # Nucleic acid main chain
  if _IsNucleicAcidResidue(r.atom_groups()[0].resname):

    # Find the lists of atoms that we might need.
    Ps = [a for a in r.atoms() if a.name.strip().upper() == "P"]
    OP3s = [a for a in r.atoms() if a.name.strip().upper() == "OP3"]
    OP2s = [a for a in r.atoms() if a.name.strip().upper() == "OP2"]
    OP1s = [a for a in r.atoms() if a.name.strip().upper() == "OP1"]
    O5Ps = [a for a in r.atoms() if a.name.strip().upper() == "O5'"]
    O4Ps = [a for a in r.atoms() if a.name.strip().upper() == "O4'"]
    O3Ps = [a for a in r.atoms() if a.name.strip().upper() == "O3'"]
    O2Ps = [a for a in r.atoms() if a.name.strip().upper() == "O2'"]
    C5Ps = [a for a in r.atoms() if a.name.strip().upper() == "C5'"]
    C4Ps = [a for a in r.atoms() if a.name.strip().upper() == "C4'"]
    C3Ps = [a for a in r.atoms() if a.name.strip().upper() == "C3'"]
    C2Ps = [a for a in r.atoms() if a.name.strip().upper() == "C2'"]
    C1Ps = [a for a in r.atoms() if a.name.strip().upper() == "C1'"]

    # Make all of the connections.
    ret += _DescribeMainchainLink(prevCs, Ps, group)
    ret += _DescribeMainchainLink(Ps, OP1s, group)
    ret += _DescribeMainchainLink(Ps, OP2s, group)
    ret += _DescribeMainchainLink(Ps, OP3s, group)
    ret += _DescribeMainchainLink(Ps, O5Ps, group)
    ret += _DescribeMainchainLink(O5Ps, C5Ps, group)
    ret += _DescribeMainchainLink(C5Ps, C4Ps, group)
    ret += _DescribeMainchainLink(C4Ps, C3Ps, group)
    ret += _DescribeMainchainLink(C3Ps, C2Ps, group)
    ret += _DescribeMainchainLink(C2Ps, O2Ps, group)

    ret += _DescribeMainchainLink(C4Ps, O4Ps, group)
    ret += _DescribeMainchainLink(O4Ps, C1Ps, group)
    ret += _DescribeMainchainLink(C1Ps, C2Ps, group)

    ret += _DescribeMainchainLink(C3Ps, O3Ps, group)

  return ret


def _DescribeMainchainResidueHydrogens(r, group, bondedNeighborLists):
  '''Return a string that describes the mainchain hydrogens for a specified residue.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  '''
  ret = ''

  # Find all of the Hydrogens in the residue
  Hs = [a for a in r.atoms() if a.element_is_hydrogen()]
  for h in Hs:
    try:
      n = bondedNeighborLists[h][0]
      # If the hydrogen is bonded to a mainchain atom, add it
      if n.name.strip().upper() in _MainChainAtomsWithHydrogen(r.atom_groups()[0].resname):
        ret += _AddPosition(n, 'P', group) + ' ' + _AddPosition(h, 'L', group, n) + '\n'
    except Exception:
      pass

  return ret


def _DescribeSidechainResidue(r, group, bondedNeighborLists):
  '''Return a string that describes the sidechain non-hydrogen portions for a specified residue.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  '''
  ret = ''

  described = []   # A list of sets of two atoms that we have already described bonds between
  queued = []

  ############################################################
  # Protein main chain
  if _IsStandardResidue(r.atom_groups()[0].resname):

    # Start with the CA atom and mark as handled all links that go back to the main chain
    # or to Hydrogens.  Add the CA to the list of atoms queued to be handled.
    try:
      # Get all of the neighbors of CA that are not N or C.  Queue them for testing.
      # Do this for all CAs found because there may be multiple atom groups (alts)
      for aCA in [a for a in r.atoms() if a.name.strip().upper() == 'CA']:
        queued.append(aCA)
        known = [a for a in bondedNeighborLists[aCA]
                 if a.element_is_hydrogen() or a.name.strip().upper() in ['N','C']]
        for a in known:
          described.append({aCA, a})
    except Exception:
      pass

  elif _IsNucleicAcidResidue(r.atom_groups()[0].resname):

    # Start with the C1' atom and mark as handled all links that go back to the main chain
    # or to Hydrogens.  Add the C1' to the list of atoms queued to be handled.
    try:
      # Get all of the neighbors of CA that are not C2' or O4'.  Queue them for testing.
      # Do this for all Cs found because there may be multiple atom groups (alts)
      for aC in [a for a in r.atoms() if a.name.strip().upper() == "C1'"]:
        queued.append(aC)
        known = [a for a in bondedNeighborLists[aC]
                 if a.element_is_hydrogen() or a.name.strip().upper() in ["C2'","O4'"]]
        for a in known:
          described.append({aC, a})
    except Exception:
      pass

  # Cycle through the list of queued atoms until we run out of them.
  # For each, look for a non-hydrogen neighbor that we've not yet described a bond
  # with.  If none are found, remove this entry from the list and cycle again.
  # If one is found, add a point and then a line for the first neighbor found, adding
  # that neighbor and all others to the queued list and continuing to chase to (line to
  # the first, adding it and all others to the queued list) until we find no more links
  # to atoms in the same residue.
  while len(queued) > 0:
    last = queued[0]
    links = [a for a in bondedNeighborLists[last]
              if (not {last, a} in described) and (not a.element_is_hydrogen())
                and (last.parent().parent() == a.parent().parent())
            ]
    if len(links) == 0:
      # First entry on the list yielded no useful neighbors; remove it and check the next
      queued = queued[1:]
      continue
    ret += _AddPosition(last, 'P', group) + ' '
    while len(links) != 0:
      # Put all but the first link into the list to be checked later.
      for a in links[1:]:
        queued.append(a)
      # Add the description for our first one and keep chasing this path
      curr = links[0]
      described.append({last,curr})
      ret += _AddPosition(curr, 'L', group, last) + '\n'
      links = [a for a in bondedNeighborLists[curr]
                if (not {curr, a} in described) and (not a.element_is_hydrogen())
                and (curr.parent().parent() == a.parent().parent())
              ]
      last = curr
  return ret


def _DescribeSidechainResidueHydrogens(r, group, bondedNeighborLists):
  '''Return a string that describes the sidechain hydrogens for a specified residue.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  '''
  ret = ''

  # Find all of the Hydrogens in the residue
  Hs = [a for a in r.atoms() if a.element_is_hydrogen()]
  for h in Hs:
    try:
      n = bondedNeighborLists[h][0]
      # If the hydrogen is bonded to a mainchain atom, add it
      if not n.name.strip().upper() in _MainChainAtomsWithHydrogen(r.atom_groups()[0].resname):
        ret += _AddPosition(n, 'P', group) + ' ' + _AddPosition(h, 'L', group, n) + '\n'
    except Exception:
      pass

  return ret


def _DescribeHet(r, group, bondedNeighborLists):
  '''Return a string that describes the bonds in a Hetatm structure.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  '''
  ret = ''

  # Start with the first atom and no described links.
  described = []   # A list of sets of two atoms that we have already described bonds between
  queued = [r.atoms()[0]]

  # Cycle through the list of queued atoms until we run out of them.
  # For each, look for a non-hydrogen neighbor that we've not yet described a bond
  # with.  If none are found, remove this entry from the list and cycle again.
  # If one is found, add a point and then a line for the first neighbor found, adding
  # that neighbor and all others to the queued list and continuing to chase to (line to
  # the first, adding it and all others to the queued list) until we find no more links
  # to atoms in the same residue.
  while len(queued) > 0:
    last = queued[0]
    links = [a for a in bondedNeighborLists[last]
              if (not {last, a} in described) and (not a.element_is_hydrogen())
                and (last.parent() == a.parent())
            ]
    if len(links) == 0:
      # First entry on the list yielded no useful neightbors; remove it and check the next
      queued = queued[1:]
      continue
    ret += _AddPosition(last, 'P', group) + ' '
    while len(links) != 0:
      # Put all but the first link into the list to be checked later.
      for a in links[1:]:
        queued.append(a)
      # Add the description for our first one and keep chasing this path
      curr = links[0]
      described.append({last,curr})
      ret += _AddPosition(curr, 'L', group, last) + '\n'
      links = [a for a in bondedNeighborLists[curr]
                if (not {curr, a} in described) and (not a.element_is_hydrogen())
                and (curr.parent() == a.parent())
              ]
      last = curr
  return ret


def _DescribeHetHydrogens(r, group, bondedNeighborLists):
  '''Return a string that describes the hydrogens for a specified residue.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  '''
  ret = ''

  # Find all of the Hydrogens in the residue
  Hs = [a for a in r.atoms() if a.element_is_hydrogen()]
  for h in Hs:
    try:
      n = bondedNeighborLists[h][0]
      ret += _AddPosition(n, 'P', group) + ' ' + _AddPosition(h, 'L', group, n) + '\n'
    except Exception:
      pass

  return ret


def _AddFlipkinBase(states, views, fileName, fileBaseName, model, alts, bondedNeighborLists,
    moverList, inSideChain, inWater, inHet):
  '''Return a string that forms the basis for a Flipkin file without the optional positions
  for the specified movers.  This includes the views that will be used to look at them.
  :param states: Return value from _FindFlipsInOutputString() indicating behavior of
  each Mover.
  :param views: List of viewpoints, one per entry in states.
  :param fileName: Name of the optimized output file associated with this Flipkin.
  :param fileBaseName: The base name of the file, without path or extension.
  :param model: The model we're optimizing.
  :param alts: A list of alternates, empty if there are none. Sorted in increasing order.
  :param bondedNeighborLists: List of neighboring atoms bonded to each atom.
  :param moverList: List of Movers, which will be used to exclude residues.
  :param inSideChain: Dictionary looked up by atom telling whether it is in a side chain.
  :param inWater: Dictionary looked up by atom telling whether it is in water.
  :param inHet: Dictionary looked up by atom telling whether it is a hetatm.
  '''
  ret = '@kinemage 1\n'
  ret += '@caption\nfrom file: {}\n'.format(fileName)

  # Compute the views for each Mover as the center of mass of all of the moving atoms and
  # record them, indicating which are flipped in Reduce.
  ret += ' views marked with * are for groups flipped by reduce\n'
  for i, s in enumerate(states):
    # We only have views for the states in the first alternate tried, so we end up with
    # more states than views.  They are repeats, so once we have done all views we are
    # done.
    if i >= len(views):
      break;

    # See whether the state is flipped in Reduce and add a star if so
    star = ' '
    if s.flipped:
      star = '*'

    # Find out the type of the residue, used to determine the type of flip.
    type = '?'
    if s.resName[-3:] == 'ASN':
      type = 'N'
    elif s.resName[-3:] == 'GLN':
      type = 'Q'
    elif s.resName[-3:] == 'HIS':
      type = 'H'

    if i > 0:
      indexString = str(i+1)
    else:
      indexString = ''
    # @todo The original Flipkin generation sometimes reported alternate conformations. We currently always
    # report the average view over all conformations.
    # ret += '@{}viewid {{{}{}{} {} {}}}\n'.format(indexString, star, type, s.resId, _AltFromFlipOutput(s), s.chain)
    ret += '@{}viewid {{{}{}{} {} {}}}\n'.format(indexString, star, type, s.resId, ' ', s.chain)
    ret += '@{}span 12\n'.format(indexString)
    ret += '@{}zslab 100\n'.format(indexString)
    ret += '@{}center{:9.3f}{:9.3f}{:9.3f}\n'.format(indexString, views[i][0], views[i][1], views[i][2])

  # Add the master descriptions
  ret += '@master {mainchain}\n'
  ret += '@master {sidechain}\n'
  ret += "@master {H's}\n"
  ret += '@master {hets}\n'
  ret += '@master {water}\n'

  # Add width and group description, which is the base name of the file without
  # its extension
  ret += '@onewidth\n'
  ret += '@group {{{}}} dominant\n'.format(fileBaseName)

  # Add the masters for the alternates if there are any. The sorted list always starts
  # with '' and then adds the others in increasing order.
  defaultAltSet = False
  for a in alts:
    # The first one is turned on and the others are turned off by default
    if defaultAltSet:
      state = 'off'
    else:
      defaultAltSet = True
      state = 'on'
    ret += "@pointmaster '{}' {{{}}} {}\n".format(a.lower(), 'alt'+a.lower(), state)

  # Add the mainchain (no hydrogens) for all residues (even Movers of the type
  # we're looking at right now). Record the location(s) for the last mainchain
  # atom in the previous residue (None for the first). Handle multiple alternates.
  ret += '@subgroup {{mc {}}} dominant\n'.format(fileBaseName)
  ret += '@vectorlist {mc} color= white  master= {mainchain}\n'
  for c in model.chains():
    prevCs = None
    for rg in c.residue_groups():
      ret += _DescribeMainchainResidue(rg, fileBaseName, prevCs)
      try:
        prevCs = [a for a in rg.atoms() if a.name.strip().upper() == 'C']
        # If not protein, check nucleic acid
        if len(prevCs) == 0:
          prevCs = [a for a in rg.atoms() if a.name.strip().upper() == "O3'"]
      except Exception:
        pass

  # Add the Hydrogens on the mainchain
  ret += "@vectorlist {mc H} color= gray  nobutton master= {mainchain} master= {H's}\n"
  for c in model.chains():
    for rg in c.residue_groups():
      if not inHet[rg.atoms()[0]] and not inWater[rg.atoms()[0]]:
        ret += _DescribeMainchainResidueHydrogens(rg, fileBaseName, bondedNeighborLists)

  # Add the sidechain non-hydrogen atoms for residues that do not have Movers
  ret += '@subgroup {{sc {}}} dominant\n'.format(fileBaseName)
  ret += '@vectorlist {sc} color= cyan  master= {sidechain}\n'
  for c in model.chains():
    for rg in c.residue_groups():
      if not inHet[rg.atoms()[0]] and not inWater[rg.atoms()[0]]:
        if not _IsMover(rg, moverList):
          ret += _DescribeSidechainResidue(rg, fileBaseName, bondedNeighborLists)

  # Add the Hydrogens on the sidechains for residues that do not have Movers
  ret += "@vectorlist {sc H} color= gray  nobutton master= {sidechain} master= {H's}\n"
  for c in model.chains():
    for rg in c.residue_groups():
      if not inHet[rg.atoms()[0]] and not inWater[rg.atoms()[0]]:
        if not _IsMover(rg, moverList):
          ret += _DescribeSidechainResidueHydrogens(rg, fileBaseName, bondedNeighborLists)

  # Describe links between atoms in a sidechain and another residue where neither of the
  # involved residues include Movers.  Don't repeat bonds that have already been
  # described.
  ret += '@vectorlist {SS} color= yellow  master= {sidechain}\n'
  described = []
  for a in model.get_atoms():
    for n in bondedNeighborLists[a]:
      if (a.parent().parent() != n.parent().parent() and inSideChain[a]
          and not _IsMover(a.parent().parent(), moverList)
          and not _IsMover(n.parent().parent(), moverList)
          ):
        if {a,n} not in described:
          ret += _AddPosition(a, 'P', fileBaseName) + ' ' + _AddPosition(n, 'L', fileBaseName, a) + '\n'
          described.append({a,n})

  # Add spheres for ions (was single-atom Het groups in original Flipkins?)
  ret += '@subgroup {het groups} dominant\n'
  ret += '@spherelist {het M} color= gray  radius= 0.5  nubutton master= {hets}\n'
  for a in model.get_atoms():
    if a.element_is_ion():
      ret += _AddPosition(a, '', fileBaseName) + '\n'

  # Add bonded structures for het groups that are not Movers
  ret += '@vectorlist {het} color= orange  master= {hets}\n'
  for c in model.chains():
    for rg in c.residue_groups():
      if inHet[rg.atoms()[0]] and not inWater[rg.atoms()[0]] and not _IsMover(rg, moverList):
         ret += _DescribeHet(rg, fileBaseName, bondedNeighborLists)
  ret += "@vectorlist {ht H} color= gray  master= {hets} master= {H's}\n"
  for c in model.chains():
    for rg in c.residue_groups():
      if inHet[rg.atoms()[0]] and not inWater[rg.atoms()[0]] and not _IsMover(rg, moverList):
         ret += _DescribeHetHydrogens(rg, fileBaseName, bondedNeighborLists)

  # Add waters
  ret += '@subgroup {waters} dominant\n'
  ret += '@balllist {water O} color= pink  radius= 0.15  master= {water}\n'
  for a in model.get_atoms():
    if inWater[a]:
      ret += _AddPosition(a, 'P', fileBaseName) + '\n'

  return ret

def _AddFlipkinMovers(states, fileBaseName, name, color, model, alts, bondedNeighborLists,
    moverList, inSideChain, inWater, inHet):
  '''Return a string that describes the Movers and atoms that are bonded to them.
  :param states: Return value from _FindFlipsInOutputString() indicating behavior of
  each Mover.
  :param fileBaseName: The base name of the file, without path or extension.
  :param name: Name for the master group.
  :param color: Color to use for the residues in states.
  :param model: The model we're optimizing.
  :param alts: A list of alternates, empty if there are none. Sorted in increasing order.
  :param bondedNeighborLists: List of neighboring atoms bonded to each atom.
  :param moverList: List of Movers, which will be used to exclude residues.
  :param inSideChain: Dictionary looked up by atom telling whether it is in a side chain.
  :param inWater: Dictionary looked up by atom telling whether it is in water.
  :param inHet: Dictionary looked up by atom telling whether it is a hetatm.
  '''
  ret = '@group {'+name+'} animate\n'
  ret += '@subgroup {sidechain} nobutton dominant\n'

  # Add balls on the nitrogens and oxygens in the Movers in the states so that
  # we can tell their orientations.  Nitrogens are sky colored and Oxygens are
  # red.
  ret += '@balllist {sc N} color= sky radius= 0.1000  nobutton master= {sidechain}\n'
  for c in model.chains():
    for rg in c.residue_groups():
      if _IsMover(rg, states):
        for a in rg.atoms():
          if a.element == 'N' and inSideChain[a]:
            ret += _AddPosition(a, 'P', fileBaseName) + '\n'
  ret += '@balllist {sc N} color= red radius= 0.1000  nobutton master= {sidechain}\n'
  for c in model.chains():
    for rg in c.residue_groups():
      if _IsMover(rg, states):
        for a in rg.atoms():
          if a.element == 'O' and inSideChain[a]:
            ret += _AddPosition(a, 'P', fileBaseName) + '\n'

  # Add the sidechain non-hydrogen atoms for the Movers that are in the states list.
  ret += '@vectorlist {sc} color= '+color+'  master= {sidechain}\n'
  for c in model.chains():
    for rg in c.residue_groups():
      if _IsMover(rg, states):
        ret += _DescribeSidechainResidue(rg, fileBaseName, bondedNeighborLists)

  # Add the Hydrogens on the sidechains for residues that are in the states list
  ret += "@vectorlist {sc H} color= gray  nobutton master= {sidechain} master= {H's}\n"
  for c in model.chains():
    for rg in c.residue_groups():
      if _IsMover(rg, states):
        ret += _DescribeSidechainResidueHydrogens(rg, fileBaseName, bondedNeighborLists)

  # Add the sidechain non-hydrogen atoms for the Movers that are not in the states list.
  ret += '@vectorlist {sc} color= cyan  master= {sidechain}\n'
  for c in model.chains():
    for rg in c.residue_groups():
      if _IsMover(rg, moverList) and not _IsMover(rg, states):
        ret += _DescribeSidechainResidue(rg, fileBaseName, bondedNeighborLists)

  # Add the Hydrogens on the sidechains for the Movers that are not in the states list
  ret += "@vectorlist {sc H} color= gray  nobutton master= {sidechain} master= {H's}\n"
  for c in model.chains():
    for rg in c.residue_groups():
      if _IsMover(rg, moverList) and not _IsMover(rg, states):
        ret += _DescribeSidechainResidueHydrogens(rg, fileBaseName, bondedNeighborLists)

  # Describe links between atoms in a sidechain and another residue where one of the
  # involved residues include Movers.  Don't repeat bonds that have already been
  # described.
  ret += '@vectorlist {SS} color= yellow  master= {sidechain}\n'
  described = []
  for a in model.get_atoms():
    for n in bondedNeighborLists[a]:
      if (a.parent().parent() != n.parent().parent() and inSideChain[a]
          and (_IsMover(a.parent().parent(), moverList)
          or _IsMover(n.parent().parent(), moverList))
          ):
        if {a,n} not in described:
          ret += _AddPosition(a, 'P', fileBaseName) + ' ' + _AddPosition(n, 'L', fileBaseName, a) + '\n'
          described.append({a,n})

  # Add bonded structures for het groups that are Movers
  ret += '@vectorlist {het} color= orange  master= {hets}\n'
  for c in model.chains():
    for rg in c.residue_groups():
      if inHet[rg.atoms()[0]] and not inWater[rg.atoms()[0]] and _IsMover(rg, moverList):
         ret += _DescribeHet(rg, fileBaseName, bondedNeighborLists)
  ret += "@vectorlist {ht H} color= gray  master= {hets} master= {H's}\n"
  for c in model.chains():
    for rg in c.residue_groups():
      if inHet[rg.atoms()[0]] and not inWater[rg.atoms()[0]] and _IsMover(rg, moverList):
         ret += _DescribeHetHydrogens(rg, fileBaseName, bondedNeighborLists)

  return ret

def _RemoveModelsExceptIndex(model_manager, model_index):
    hierarchy = model_manager.get_hierarchy()
    models = hierarchy.models()
    if model_index < len(models):
        selected_model = models[model_index]
        for model in models:
            if model != selected_model:
                hierarchy.remove_model(model=model)
    return model_manager

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
reduce2 version {}
Add Hydrogens to a model and optimize their placement by adjusting movable groups and
flippable groups of atoms.

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed
Output:
  PDB or mmCIF file with added hydrogens.  If output.filename is specified, then the
  type of file to write will be determined by its suffix (.pdb or .cif).
  If output.filename is not specified, the output file will be
  written into the current working directory with the same base name and type as the
  original file and with FH or H added to the base name (FH when flips are requested);
  1xs0.pdb would be written to ./1xsoH.pdb and 1xso.cif to ./1xsoH.cif by default.

NOTES:
  If multiple alternates are present in the file and a specific one is not specified on the
  command line, they will all be processed in reverse order, such that the lowest-named
  one (perhaps A) will be processed last.  The hydrogen addition and arrangements for
  residues that are not part of any alternate will be left in the configuration that is best
  for the final alternate tested.  This may leave other alternates in sub-optimal configurations.
  When a single alternate is selected using alt_id= on the command line, everything is
  optimized a single time and for that configuration.

  Note that the program also takes probe Phil arguments; run with --show_defaults to see
  all Phil arguments.

  As of 7/14/2025, reduce2 is switching to new default parameters for the relative weighting of
  external contacts (which remains the same) and both hydrogen bonds and collisions (whic are
  increasing by 10x). This is to work as expected with a radius larger than 0 (it is being switched
  to 0.25, which matches the probe2 default). The original Reduce had switched to a radius of 0.0
  and was not considering external contacts; this fixes that without causing hydrogen bonds to be
  broken. See https://github.com/cctbx/cctbx_project/issues/1072 for details. The defaults for
  probe2 are not being changed, so its default behavior will be different from reduce2 until the
  issue can be fully resolved and both set to the same defaults. The probe radius and relative
  weights can be set in both probe2 and reduce2 using the probe2 Phil parameters if different
  behavior is desired.

  Equivalent PHIL arguments for original Reduce command-line options:
    -quiet: No equivalent; metadata is never written to the model file, it is always
            written to the description file, and progress information is always written
            to standard output.
    -trim: approach=remove
    -build: approach=add add_flip_movers=True
    -flip: approach=add add_flip_movers=True
    -allalt: This is the default.
    -penalty200: preference_magnitude=200
    -nobuild9999: approach=add preference_magnitude=9999
    -noflip: approach=add add_flip_movers=True preference_magnitude=9999
    -onlya: alt_id=A
    -nuclear: use_neutron_distances=True
    -demandflipallnhqs: add_flip_movers=True
    -rad0.25: probe.probe_radius=0.25
'''.format(version)
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']
  citations = program_citations
  epilog = '''
  For additional information and help, see http://kinemage.biochem.duke.edu/software/reduce
  and http://molprobity.biochem.duke.edu
  '''

# ------------------------------------------------------------------------------

  def _GetViews(self, movers):
    '''Produce a list of views that will center the sidechains of the specified Movers.
    It ignores the alternate and averages over all of them.  It only makes one entry
    if the same Mover is in more than one alternate.
    :param movers: Movers returned by _FindFlipsInOutputString().
    :return: List of 3-tuples of centers of the sidechains.
    '''
    views = []
    selStrings = []
    for mover in movers:
      # Fill in information needed to construct the view.
      # As of 2/5/2023, the CCTBX selection returns no atoms on a file when the model
      # clause is used unless there is a MODEL statement in the file.  The get_number_of_models()
      # function returns 1 if there are 0 or 1 MODEL statements, so we check to see if there
      # are 2 or more (indicating the need to select) before adding the clause.
      # The model ID that the selection is looking for is 1-based, so we must add 1 to the
      # model index.
      if self.model.get_number_of_models() >= 2:
        modelClause = 'model {} and '.format(mover.modelId + 1)
      else:
        modelClause = ''
      x = 0.0
      y = 0.0
      z = 0.0
      selString = modelClause + "chain {} and resseq {} and sidechain".format(
            mover.chain, mover.resId)
      if selString in selStrings:
        break
      selStrings.append(selString)
      sel = self.model.selection(selString)
      count = 0;
      for a in self.model.get_hierarchy().atoms():
        if sel[a.i_seq]:
          x += a.xyz[0]
          y += a.xyz[1]
          z += a.xyz[2]
          count += 1
      if count > 0:
        x /= count
        y /= count
        z /= count
      views.append( (x, y, z) )
    return views
# ------------------------------------------------------------------------------

  # Create a parser for Probe2 PHIL parameters, overriding specific defaults.
  # Use it to parse the parameters we need to change and return the parser so we
  # can extract values from it.
  def _MakeProbePhilParser(self, movers_to_check, extraArgs = []):

    # Determine the source and target selections based on the Movers
    # that we are checking being tested against everything else.
    source_selection = 'sidechain and ('
    for i, m in enumerate(movers_to_check):
      if i > 0:
        term = 'or ('
      else:
        term = ' ('
      term += 'chain {} and resid {}) '.format(m.chain, m.resId)
      source_selection += term
    source_selection += ')'
    target_selection = 'all'

    parser = cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
    args = [
      "source_selection='{}'".format(source_selection),
      "target_selection='{}'".format(target_selection),
      "use_neutron_distances={}".format(self.params.use_neutron_distances),
      "approach=both",
      "excluded_bond_chain_length={}".format(self._bondedNeighborDepth),
      "minimum_water_hydrogen_occupancy=0.66",
      "maximum_water_hydrogen_b=40.0",
      "minimum_occupancy=0.01",
      "output.write_files=False",
      "ignore_lack_of_explicit_hydrogens=True",
      "output.add_group_line=False"
    ]
    args.extend(extraArgs)

    # Add the probe parameters from the command line to the parser
    for arg in sys.argv:
      if arg.startswith("probe."):
        # Remove the leading "probe." and add it to the parser
        arg = arg[6:]
        args.append(arg)

    parser.parse_args(args)
    return parser

# ------------------------------------------------------------------------------

  def _AddHydrogens(self):
    reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
      model = self.model,
      use_neutron_distances=self.params.use_neutron_distances,
      n_terminal_charge=self.params.n_terminal_charge,
      exclude_water=True,
      stop_for_unknowns=self.params.stop_on_any_missing_hydrogen,
      keep_existing_H=self.params.keep_existing_H
    )
    reduce_add_h_obj.run()
    reduce_add_h_obj.show(self.logger)
    missed_residues = set(reduce_add_h_obj.no_H_placed_mlq)
    if not self.params.ignore_missing_restraints:
      if len(missed_residues) > 0:
        bad = ""
        for res in missed_residues:
          bad += " " + res
        raise Sorry("Restraints were not found for the following residues:"+bad)
    insufficient_restraints = list(reduce_add_h_obj.site_labels_no_para)
    if self.params.stop_on_any_missing_hydrogen and len(insufficient_restraints) > 0:
      bad = insufficient_restraints[0]
      for res in insufficient_restraints[1:]:
        bad += "," + res
      raise Sorry("Insufficient restraints were found for the following atoms:"+bad)

    self.model = reduce_add_h_obj.get_model()

    if not self.model.has_hd():
      raise Sorry("It was not possible to place any H atoms. Is this a single atom model?")

    # We must re-interpret the model when _type_energies is not defined, which happens
    # for example when running on 8zst.cif.
    if not hasattr(self.model, '_type_energies'):
      self._ReinterpretModel(make_restraints=True)

# ------------------------------------------------------------------------------

  def _ReinterpretModel(self, make_restraints=True):
    # Reinterpret the model using the same approach that hydrogen-placement does.
    # :param make_restraints: Should we compute restraints during the interpretation?
    self.model.get_hierarchy().sort_atoms_in_place()
    self.model.get_hierarchy().atoms().reset_serial()

    p = reduce_hydrogen.get_reduce_pdb_interpretation_params(
      self.params.use_neutron_distances)
    p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check=True
    # We need to turn this on because without it 1zz0.txt kept flipping the ring
    # in A TYR 214 every time we re-interpreted. The original interpretation done
    # by Hydrogen placement will have flipped them, so we don't need to do it again.
    p.pdb_interpretation.flip_symmetric_amino_acids=False
    #p.pdb_interpretation.sort_atoms=True
    self.model.set_stop_for_unknowns(self.params.stop_on_any_missing_hydrogen)
    self.model.process(make_restraints=make_restraints, pdb_interpretation_params=p)

# ------------------------------------------------------------------------------

  def _GetAtomCharacteristics(self, bondedNeighborLists):
    '''Determine additional characteristics needed to determine probe dots.
    :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
    structure that the atom from the first parameter interacts with that lists all of the
    bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
    :return: Tuple of maps to Booleans for whether the atom is in water, in Hetatm
    group, in the main chain, and in a side chain.
    '''
    inWater = {}
    inHet = {}
    inMainChain = {}
    inSideChain = {}
    hetatm_sel = self.model.selection("hetatm")
    mainchain_sel = self.model.selection("backbone")  # Will NOT include Hydrogen atoms on the main chain
    sidechain_sel = self.model.selection("sidechain") # Will include Hydrogen atoms on the side chain
    for a in self.model.get_atoms():
      inWater[a] = common_residue_names_get_class(name=a.parent().resname) == "common_water"
      inHet[a] = hetatm_sel[a.i_seq]
      if not a.element_is_hydrogen():
        inMainChain[a] = mainchain_sel[a.i_seq]
      else:
        # Check our bonded neighbor to see if it is on the mainchain if we are a Hydrogen
        if len(bondedNeighborLists[a]) < 1:
          raise Sorry("Found Hydrogen with no neigbors.")
        else:
          inMainChain[a] = mainchain_sel[bondedNeighborLists[a][0].i_seq]
      inSideChain[a] = sidechain_sel[a.i_seq]
    return (inWater, inHet, inMainChain, inSideChain)

# ------------------------------------------------------------------------------

  def _DescribeLockdown(self, flipMover, invertFlip, fixedUp):
    '''Produce a flip-state string describing how to lock down the specified flip mover.
    :param flipMover: The Optimizers.FlipMoverState object being locked down.
    :param invertFlip: Do we want to invert the flip state relative to the one described
    in the mover state?
    :param fixedUp: Should the Fixup be applied after flipping?
    :return: String description suitable for use in the --flip_states command-line option.
    For example: 1 . A HIS 11H Flipped AnglesAdjusted
    '''
    if fixedUp:
      adjustedString = 'AnglesAdjusted'
    else:
      adjustedString = 'AnglesNotAdjusted'
    flipped = flipMover.flipped
    if invertFlip:
      flipped = not flipped
    if flipped:
      flipStateString = 'Flipped'
    else:
      flipStateString = 'Unflipped'
    # Look up the residue that is described and find Nitrogens in it whose names are
    # longer than one character. These will be part of the flipping part of the residue.
    # See if any of these atoms residue groups have the same altid as the flipMover; if
    # so, use that one.  If not, use ''.  If we use '', replace it with '.' so that we
    # can properly parse the string.
    # @todo Atom selection here is will not necessarily work for future flip Movers.
    # @todo How to handle cases where the different alternates resulted in different placement
    # for the flip Movers?  Presumably, we want to do them in backwards order and replace all of
    # the alternates with '.' so that they always get handled and the last one is the one that
    # is set.
    altId = ''
    flipAlt = flipMover.altId.strip()
    resAlt = ''
    for chain in self.model.chains():
      if flipMover.chain.strip() == chain.id.strip():
        for rg in chain.residue_groups():
          if int(flipMover.resId) == rg.resseq_as_int():
            for a in rg.atoms():
              if (a.element == 'N') and (len(a.name.strip()) > 1):
                if a.parent().altloc.lower() == flipAlt.lower():
                  resAlt = flipAlt
    if flipAlt == resAlt:
      altId = flipAlt
    if altId.strip() == '':
      altId = '.'
    return '{} {} {} {} {}{} {} {}'.format(flipMover.modelId+1, altId.lower(), flipMover.chain,
      flipMover.resName, flipMover.resId, flipMover.iCode, flipStateString, adjustedString)

  # ------------------------------------------------------------------------------

  def validate(self):
    # Set the default output file name if one has not been given.
    if self.params.output.filename is None:
      inName = self.data_manager.get_default_model_name()
      if inName is None: raise Sorry('model not found')
      suffix = os.path.splitext(os.path.basename(inName))[1]
      if self.params.add_flip_movers:
        pad = 'FH'
      else:
        pad = 'H'
      base = os.path.splitext(os.path.basename(inName))[0] + pad
      self.params.output.filename = base + suffix
      print('Writing model output to', self.params.output.filename, file=self.logger)

    if os.environ.get('PHENIX_OVERWRITE_ALL', False):
      self.data_manager.set_overwrite(True)
    if not self.params.output.overwrite:
      if os.path.exists(self.params.output.filename):
        print('\n\tOutput filename exists. Use overwrite=True to continue.')

    self.data_manager.has_models(raise_sorry=True)
    if self.params.output.description_file_name is None:
      self.params.output.description_file_name=self.params.output.filename.replace('.pdb',
                                                                                   '.txt')
      self.params.output.description_file_name=self.params.output.description_file_name.replace('.cif',
                                                                                   '.txt')

    # Check the model ID to make sure they didn't set it to 0
    if self.params.model_id == 0:
      raise Sorry("Model ID must be >=1 if specified (None means all models)")

    # Turn on profiling if we've been asked to in the Phil parameters
    if self.params.profile:
      import cProfile
      self._pr = cProfile.Profile()
      self._pr.enable()

  # ------------------------------------------------------------------------------

  def run(self):

    # Set our bonded-neighbor depth
    self._bondedNeighborDepth = self.params.bonded_neighbor_depth

    # String describing the run that will be output to the specified file.
    outString = 'reduce2 v.{}, run {}\n'.format(version, datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for a in sys.argv:
      outString += ' {}'.format(a)
    outString += '\n'

    make_sub_header('Loading Model', out=self.logger)

    # Get our model.
    self.model = self.data_manager.get_model()

    self.model = self.model.select(~self.model.selection('element X'))

    # Use model function to set crystal symmetry if necessary 2025-03-19 TT
    self.model.add_crystal_symmetry_if_necessary()
    if self.data_manager.has_restraints():
      self.model.set_stop_for_unknowns(self.params.stop_on_any_missing_hydrogen)
      self.model.process(make_restraints=False)

    # If we've been asked to only to a single model index from the file, strip the model down to
    # only that index.
    if self.params.model_id is not None:
      make_sub_header('Selecting Model ID ' + str(self.params.model_id), out=self.logger)

      # Select only the current submodel from the hierarchy
      submodel = self.model.deep_copy()
      _RemoveModelsExceptIndex(submodel, self.params.model_id)

      # Construct a hierarchy for the current submodel
      r = pdb.hierarchy.root()
      mdc = submodel.get_hierarchy().models()[0].detached_copy()
      r.append_model(mdc)

      # Make yet another model for the new hierarchy
      subset_model_manager = mmtbx.model.manager(
        model_input       = None,
        pdb_hierarchy     = r,
        stop_for_unknowns = False,
        crystal_symmetry  = submodel.crystal_symmetry(),
        restraint_objects = None,
        log               = None)

      self.model = subset_model_manager

    # Stores the initial coordinates for all of the atoms and the rest of the information
    # about the original model for use by Kinemages.
    initialModel = self.model.deep_copy()

    if self.params.approach == 'add':
      # Add Hydrogens to the model
      make_sub_header('Adding Hydrogens', out=self.logger)
      startAdd = time.time()
      self._AddHydrogens()
      doneAdd = time.time()

      # NOTE: We always optimize all models (leave modelIndex alone) because we've removed all
      # but the desired model ID structure from the model.
      make_sub_header('Optimizing', out=self.logger)
      startOpt = time.time()
      opt = Optimizers.Optimizer(self.params.probe, self.params.add_flip_movers,
        self.model, altID=self.params.alt_id,
        preferenceMagnitude=self.params.preference_magnitude,
        bondedNeighborDepth = self._bondedNeighborDepth,
        nonFlipPreference=self.params.non_flip_preference,
        skipBondFixup=self.params.skip_bond_fix_up,
        flipStates = self.params.set_flip_states,
        verbosity=self.params.verbosity,
        cliqueOutlineFileName=self.params.output.clique_outline_file_name,
        fillAtomDump = self.params.output.print_atom_info)
      doneOpt = time.time()
      warnings = opt.getWarnings()
      # Find the lines from warnings that start with "Warning: " and
      # concatenate them into a single string.
      warnings = '\n'.join([line for line in warnings.split('\n') if line.startswith('Warning:')])
      if len(warnings) > 0:
        print('\nWarnings during optimization:\n'+warnings, file=self.logger)
      outString += opt.getInfo()
      outString += 'Time to Add Hydrogen = {:.3f} sec'.format(doneAdd-startAdd)+'\n'
      outString += 'Time to Optimize = {:.3f} sec'.format(doneOpt-startOpt)+'\n'
      if self.params.output.print_atom_info:
        print('Atom information used during calculations:', file=self.logger)
        print(opt.getAtomDump(), file=self.logger)

    else: # Removing Hydrogens from the model rather than adding them.
      make_sub_header('Removing Hydrogens', out=self.logger)
      sel = self.model.selection("element H")
      for a in self.model.get_atoms():
        if sel[a.i_seq]:
          a.parent().remove_atom(a)

    # Re-process the model because we have removed some atoms that were previously
    # bonded.  Don't make restraints during the reprocessing.
    # We had to do this to keep from crashing on a call to pair_proxies when generating
    # mmCIF files, so we always do it for safety.
    self._ReinterpretModel(False)

    make_sub_header('Writing output', out=self.logger)

    # Skip writing the main output file and description output file if output.write_files is False.
    # This enables a program to harness Reduce2 and not have to deal with handling output files.
    if self.params.output.write_files:
      # Write the description output to the specified file.
      self.data_manager._write_text("description", outString,
        self.params.output.description_file_name)

      # Determine whether to write a PDB or CIF file and write the appropriate text output.
      suffix = os.path.splitext(self.params.output.filename)[1]
      if suffix.lower() == ".pdb":
        txt = self.model.model_as_pdb()
      else:
        txt = self.model.model_as_mmcif()
      self.data_manager._write_text("model", txt, self.params.output.filename)

      print('Wrote', self.params.output.filename,'and',
        self.params.output.description_file_name, file = self.logger)

    # If we've been asked to do a comparison with another program's output, do it.
    if self.params.comparison_file is not None:
      make_sub_header('Comparing with other model', out=self.logger)

      # Construct the file names we'll be using.
      compareOutName = self.params.output.filename + "_comparison.csv"

      # Find the list of all Movers in the model, which will be used to select
      # each in turn for evaluation.
      moverLocations = _FindMoversInOutputString(outString)

      # Make the first line in our table be the header line
      table = [ ["Mover",
                 "Score from Reduce2",
                 "Score from {}".format(self.params.comparison_file),
                 "Difference",
                 "Star if other is higher"] ]

      # Run Probe2 on each of the Movers for our current model and for
      # comparison_file and determine the summary score of the Mover against
      # the rest of the model in each case. Accumulate these into a table.
      for i, m in enumerate(moverLocations):
        print ('Computing scores for Mover',i,'of',len(moverLocations))

        #============================================================
        # Find the score for our output.

        # Make the Probe2 Phil parameters, then overwrite the ones that were
        # filled in with values that we want for our summaries.
        source = [ m ]
        extraArgs = [
          "approach=once",
          "output.format=raw",
          "output.contact_summary=True",
          "output.condensed=True",
          "output.count_dots=True"
          ]
        probeParser = self._MakeProbePhilParser(source, extraArgs)

        # Run Probe2
        p2 = probe2.Program(self.data_manager, probeParser.working_phil.extract(),
                            master_phil=probeParser.master_phil, logger=self.logger)
        p2.overrideModel(self.model)
        dots, output = p2.run()

        # Parse the output to find the scores for each of the two directions
        # and sum them up. There is one per line and the number ends with #.
        lines = output.strip().split('\n')
        values = []
        for line in lines:
          stripped_line = line.strip()
          if stripped_line:
            number = float(stripped_line.split('#')[0].strip())
            values.append(number)
        myScore = sum(values)

        print('My values for Mover', str(m), 'are', values, 'sum is', myScore, file=self.logger)

        #============================================================
        # Find the score for the comparison file.

        # Read the file using a new datamanager
        dm = DataManager()
        dm.process_model_file(self.params.comparison_file)
        #otherModel = dm.get_model(self.params.comparison_file)

        # Make the Probe2 Phil parameters, then overwrite the ones that were
        # filled in with values that we want for our summaries.
        extraArgs = [
          "approach=once",
          "output.format=raw",
          "output.contact_summary=True",
          "output.condensed=True",
          "output.count_dots=True"
          ]
        probeParser = self._MakeProbePhilParser(source, extraArgs)

        # Run Probe2
        p2 = probe2.Program(self.data_manager, probeParser.working_phil.extract(),
                            master_phil=probeParser.master_phil, logger=self.logger)
        p2.overrideModel(dm.get_model(self.params.comparison_file))
        dots, output = p2.run()

        # Parse the output to find the scores for each of the two directions
        # and sum them up. There is one per line and the number ends with #.
        lines = output.strip().split('\n')
        values = []
        for line in lines:
          stripped_line = line.strip()
          if stripped_line:
            number = float(stripped_line.split('#')[0].strip())
            values.append(number)
        otherScore = sum(values)

        print('Other values for Mover', str(m), 'are', values, 'sum is', otherScore, file=self.logger)

        #============================================================
        # Add the line to the table, indicating if the other is better.
        mark = ''
        if otherScore > myScore:
          mark = '*'
        table.append( [str(m), myScore, otherScore, myScore - otherScore, mark] )

      # Write the table to our output CSV file
      with open(compareOutName, 'w', newline='') as csvfile:
          writer = csv.writer(csvfile)
          writer.writerows(table)

    # If we've been asked to make Flipkins, then make each of them.
    if self.params.add_flip_movers and self.params.output.flipkin_directory is not None:
      make_sub_header('Producing flipkins', out=self.logger)

      # Find the base name of the two output files we will produce.
      inName = self.data_manager.get_default_model_name()
      suffix = os.path.splitext(os.path.basename(inName))[1]
      pad = 'FH'
      base = os.path.splitext(os.path.basename(inName))[0] + pad
      flipkinBase = self.params.output.flipkin_directory + "/" + base

      # Find the list of all Movers in the model, which will be used to segment
      # it into parts for the Flipkin.
      moverLocations = _FindMoversInOutputString(outString)

      # Find the list of all alternates in the model, ignoring empty ones.
      # Sort them in increasing alphabetical order.
      alts = Optimizers.AlternatesInModel(self.model)
      alts.discard('')
      alts.discard(' ')
      alts = sorted(list(alts))

      # ===========================================================================

      # Make list of Amides to lock in one flip orientation and then the other,
      # keeping track of which state they are in when Reduce was choosing.
      # We need a different list for the Amide Movers and the Histidine Movers
      # because we generate two different Flipkin files, one for each.
      amides = _FindFlipsInOutputString(outString, 'AmideFlip')

      if len(amides) > 0:
        make_sub_header('Generating Amide Flipkin', out=self.logger)

        # Find the viewpoint locations for each Mover we're going to
        # look at.
        views = self._GetViews(amides)

        # Interpret the model to fill in things we need for determining the neighbor list.
        self._ReinterpretModel()

        carts = flex.vec3_double()
        for a in self.model.get_atoms():
          carts.append(a.xyz)
        bondProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
        bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies)

        # Get the other characteristics we need to know about each atom to do our work.
        inWater, inHet, inMainChain, inSideChain = self._GetAtomCharacteristics(bondedNeighborLists)

        # Write the base information in the Flipkin, not including the moving atoms in
        # the Movers that will be placed, or atoms bonded to the moving atoms.
        flipkinText = _AddFlipkinBase(amides, views, self.params.output.filename, base, self.model,
          alts, bondedNeighborLists, moverLocations, inSideChain, inWater, inHet)

        # Make two configurations, the one that Reduce picked and the one
        # that it did not.
        configurations = ['reduce', 'flipNQ']
        colors = ['sea', 'pink']

        # Run the optimization without fixup on the original orientation of all flippers and
        # again on the flipped orientation for each and combine the info from both of them
        # into the same Flipkin. Handle the re-initialization of atom coordinates before
        # each optimization, remembering the addition and deletion of atoms during placement
        # and optimization.
        for i, c in enumerate(configurations):
          # Restore the model to the state it had before we started adjusting it.
          self.model = initialModel.deep_copy()

          # Rerun hydrogen placement.
          self._AddHydrogens()

          # Run optimization, locking the specified Amides into each configuration.
          # Don't do fixup on the ones that are locked down.  Make sure that we can
          # avoid adding a comma on the first flip state that is added.
          flipStates = self.params.set_flip_states
          if flipStates is None:
            flipStates = ''
          if flipStates.strip() != '':
            flipStates += ','
          if i == 0:
            # For the first configuration, set all Amides to the orientation that Reduce decided
            # on, adding these to the list of locked-down flips.
            for ai, amide in enumerate(amides):
              flipStates += self._DescribeLockdown(amide, invertFlip=False, fixedUp=False)
              if ai < len(amides) - 1:
                flipStates += ','
          else:
            # For the second configuration, set all Amides to the orientation that Reduce did not
            # decide on, adding these to the list of locked-down flips.
            for ai, amide in enumerate(amides):
              flipStates += self._DescribeLockdown(amide, invertFlip=True, fixedUp=False)
              if ai < len(amides) - 1:
                flipStates += ','

          # Optimize the model and then reinterpret it so that we can get all of the information we
          # need for the resulting set of atoms (which may be fewer after Hydrogen removal).
          # NOTE: We always optimize all models (leave modelIndex alone) because we've removed all
          # but the desired model ID structure from the model.
          opt = Optimizers.Optimizer(self.params.probe, self.params.add_flip_movers,
            self.model, altID=self.params.alt_id,
            preferenceMagnitude=self.params.preference_magnitude,
            nonFlipPreference=self.params.non_flip_preference,
            skipBondFixup=self.params.skip_bond_fix_up,
            flipStates = flipStates,
            verbosity=3)
          print('Results of optimization:', file=self.logger)
          print(opt.getInfo(), file=self.logger)
          self._ReinterpretModel()

          # Get the other characteristics we need to know about each atom to do our work.
          # We must do this again here because the atoms change after Hydrogen addition.
          # We also need to re-generate the bonded-neighbor lists for the same reason.
          carts = flex.vec3_double()
          for a in self.model.get_atoms():
            carts.append(a.xyz)
          bondProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
          bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies)
          inWater, inHet, inMainChain, inSideChain = self._GetAtomCharacteristics(bondedNeighborLists)

          # Write the updates to the Flipkin for this configuration, showing the
          # atoms for the amide in the Reduce configuration (i=0) or the other
          # configuration (i=1).
          flipkinText += _AddFlipkinMovers(amides, base, c, colors[i], self.model, alts,
            bondedNeighborLists, moverLocations, inSideChain, inWater, inHet)

          # Compute the dots in both 1->2 and 2->1 directions using probe2.
          # @todo Consider adding methods to probe2 that let us inject the various
          # computed quantities -- neighbor lists and extra atom info and such --
          # before calling the run() method so it does not have to recompute them.

          # Modify the parameters that are passed to include the ones for
          # the harnessed program, including the source and target atom selections.
          probeParser = self._MakeProbePhilParser(amides)

          # Run the program and append its Kinemage output to ours, deleting
          # the temporary file that it produced.
          p2 = probe2.Program(self.data_manager, probeParser.working_phil.extract(),
                              master_phil=probeParser.master_phil, logger=self.logger)

          p2.overrideModel(self.model)
          dots, kinString = p2.run()
          flipkinText += kinString

        # Write the accumulated Flipkin string to the output file.
        with open(flipkinBase+"-flipnq.kin", "w") as f:
          f.write(flipkinText)

      # ===========================================================================
      hists = _FindFlipsInOutputString(outString, 'HisFlip')

      if len(hists) > 0:
        make_sub_header('Generating Histidine Flipkin', out=self.logger)

        # Find the viewpoint locations for each Mover we're going to
        # look at.
        views = self._GetViews(hists)

        # Interpret the model to fill in things we need for determining the neighbor list.
        self._ReinterpretModel()

        carts = flex.vec3_double()
        for a in self.model.get_atoms():
          carts.append(a.xyz)
        bondProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
        bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies)

        # Get the other characteristics we need to know about each atom to do our work.
        inWater, inHet, inMainChain, inSideChain = self._GetAtomCharacteristics(bondedNeighborLists)

        # Write the base information in the Flipkin, not including the moving atoms in
        # the Movers that will be placed, or atoms bonded to the moving atoms.
        flipkinText = _AddFlipkinBase(hists, views, self.params.output.filename, base, self.model,
          alts, bondedNeighborLists, moverLocations, inSideChain, inWater, inHet)

        # Make two configurations, the one that Reduce picked and the one
        # that it did not.
        configurations = ['reduce', 'flipH']
        colors = ['sea', 'pink']

        # Run the optimization without fixup on the original orientation of all flippers and
        # again on the flipped orientation for each and combine the info from both of them
        # into the same Flipkin. Handle the re-initialization of atom coordinates before
        # each optimization, remembering the addition and deletion of atoms during placement
        # and optimization.
        for i, c in enumerate(configurations):
          # Restore the model to the state it had before we started adjusting it.
          self.model = initialModel.deep_copy()

          # Rerun hydrogen placement.
          self._AddHydrogens()

          # Run optimization, locking the specified Histidines into each configuration.
          # Don't do fixup on the ones that are locked down.  Make sure that we can
          # avoid adding a comma on the first flip state that is added.
          flipStates = self.params.set_flip_states
          if flipStates is None:
            flipStates = ''
          if flipStates.strip() != '':
            flipStates += ','
          if i == 0:
            # For the first configuration, set all Histidines to the orientation that Reduce decided
            # on, adding these to the list of locked-down flips.
            for ai, hist in enumerate(hists):
              flipStates += self._DescribeLockdown(hist, invertFlip=False, fixedUp=False)
              if ai < len(hists) - 1:
                flipStates += ','
          else:
            # For the second configuration, set all HIstidines to the orientation that Reduce did not
            # decide on, adding these to the list of locked-down flips.
            for hi, hist in enumerate(hists):
              flipStates += self._DescribeLockdown(hist, invertFlip=True, fixedUp=False)
              if hi < len(hists) - 1:
                flipStates += ','

          # Optimize the model and then reinterpret it so that we can get all of the information we
          # need for the resulting set of atoms (which may be fewer after Hydrogen removal).
          # NOTE: We always optimize all models (leave modelIndex alone) because we've removed all
          # but the desired model ID structure from the model.
          opt = Optimizers.Optimizer(self.params.probe, self.params.add_flip_movers,
            self.model, altID=self.params.alt_id,
            preferenceMagnitude=self.params.preference_magnitude,
            nonFlipPreference=self.params.non_flip_preference,
            skipBondFixup=self.params.skip_bond_fix_up,
            flipStates = flipStates,
            verbosity=3)
          print('Results of optimization:', file=self.logger)
          print(opt.getInfo(), file=self.logger)
          self._ReinterpretModel()

          # Get the other characteristics we need to know about each atom to do our work.
          # We must do this again here because the atoms change after Hydrogen addition.
          # We also need to re-generate the bonded-neighbor lists for the same reason.
          carts = flex.vec3_double()
          for a in self.model.get_atoms():
            carts.append(a.xyz)
          bondProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
          bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies)
          inWater, inHet, inMainChain, inSideChain = self._GetAtomCharacteristics(bondedNeighborLists)

          # Write the updates to the Flipkin for this configuration, showing the
          # atoms for the Histidines in the Reduce configuration (i=0) or the other
          # configuration (i=1).
          flipkinText += _AddFlipkinMovers(hists, base, c, colors[i], self.model, alts,
            bondedNeighborLists, moverLocations, inSideChain, inWater, inHet)

          # Compute the dots in both 1->2 and 2->1 directions using probe2.
          # @todo Consider adding methods to probe2 that let us inject the various
          # computed quantities -- neighbor lists and extra atom info and such --
          # before calling the run() method so it does not have to recompute them.

          # Modify the parameters that are passed to include the ones for
          # the harnessed program, including the source and target atom selections.
          probeParser = self._MakeProbePhilParser(hists)

          # Run the program and append its Kinemage output to ours, deleting
          # the temporary file that it produced.
          p2 = probe2.Program(self.data_manager, probeParser.working_phil.extract(),
                              master_phil=probeParser.master_phil, logger=self.logger)
          p2.overrideModel(self.model)
          dots, kinString = p2.run()
          flipkinText += kinString

        # Write the accumulated Flipkin string to the output file.
        with open(flipkinBase+"-fliphis.kin", "w") as f:
          f.write(flipkinText)

    # Report profiling info if we've been asked to in the Phil parameters
    if self.params.profile:
      print('Profile results:', file=self.logger)
      import pstats
      profile_params = {'sort_by': 'time', 'num_entries': 20}
      self._pr.disable()
      ps = pstats.Stats(self._pr).sort_stats(profile_params['sort_by'])
      ps.print_stats(profile_params['num_entries'])

# ------------------------------------------------------------------------------

  def get_results(self):
    return group_args(model = self.model)

# ------------------------------------------------------------------------------

  def Test(self):
    '''
      Run tests on the methods of the class.  Throw an assertion error if there is a problem with
      one of them and return normally if there is not a problem.
    '''

    #=====================================================================================
    # @todo Unit tests for other methods
