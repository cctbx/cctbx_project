"""Run reduce (add hydrogens). version 2"""
##################################################################################
# Copyright(c) 2021-2026, Richardson Lab at Duke
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
from libtbx.utils import Sorry, null_out, peak_memory_usage, human_readable_size
import mmtbx
from mmtbx.probe import Helpers
from iotbx import pdb, cli_parser
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens import hd_exchange
from mmtbx.hydrogens import capping
from mmtbx.reduce import Optimizers
from scitbx.array_family import flex
from iotbx.pdb import common_residue_names_get_class, amino_acid_codes, nucleic_acid_codes
from mmtbx.programs import probe2
import copy
from iotbx.data_manager import DataManager
import csv

version = "3.5.0"

master_phil_str = '''
approach = *add remove optimize place
  .type = choice
  .short_caption = Add (and optimize), remove, optimize, or place Hydrogens
  .help = Determines whether Reduce will add (and optimize), remove, or optimize Hydrogens from the model, or place them without optimizing (place). Note: if using optimize and add_flip_movers=True, all possible hydrogens must be added to histidines (even ones making collision). remove removes all H and D atoms and renames DOD waters to HOH.
keep_existing_H = False
  .type = bool
  .short_caption = Do not remove Hydrogens in the original model
  .help = Keep existing H and D atoms in the model. Kept D need a neutron_option other than none.
n_terminal_charge = *residue_one first_in_chain all_chain_starts no_charge
  .type = choice(multi=False)
  .short_caption = N terminal charge approach
  .help = Mode for placing H3 at terminal nitrogen: residue_one (resseq 1),\
          first_in_chain (first residue of each chain), all_chain_starts\
          (first residue of each chain and each residue after a chain\
          break), no_charge.
metal_bound_amines_as_nh2 = True
  .type = bool
  .short_caption = Place metal-bound primary amines as NH2
  .help = A primary amine bound to a metal (Lys NZ, an N-terminal NH3, a \
          ligand amine, not with alkali metals or Ca, Sr, Ba, Ra) is placed \
          as NH2, with the metal as the third substituent. Metal-bound chain \
          starts get NH2 whatever n_terminal_charge says. Sites that fail the \
          geometry checks keep their H.
exclude_water = True
  .type = bool
  .short_caption = Add H to water if False
  .help = Add H to water if False
use_neutron_distances = False
  .type = bool
  .short_caption = Use neutron distances
  .help = Use neutron distances (-nuclear in reduce). Switched on by any neutron_option other than none and by qmready=True.
neutron_option = *none all_h all_d hd_and_h hd_and_d all_hd
  .type = choice(multi=False)
  .short_caption = H/D model for neutron refinement
  .help = none (default): no deuterium, distances follow use_neutron_distances \
          and qmready as before. all_h: hydrogens only, neutron distances. \
          all_d: all deuterium. hd_and_h: H/D pairs (altlocs) at exchangeable \
          sites, H elsewhere. hd_and_d: H/D pairs at exchangeable sites, D \
          elsewhere. all_hd: H/D pairs at all sites. Every option other than \
          none uses neutron distances. Waters become DOD in all_d, hd_and_h, \
          hd_and_d and all_hd, with all_h they stay HOH. H whose names do not \
          start with H stay H in every option.
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
  .help = MODEL number to process (as in the file's MODEL records / pdbx_PDB_model_num). Default: all models. If one is selected, the others are removed from the output.
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
  .help = Don't stop if restraints for a residue are missing (otherwise every approach except remove stops). With approach=optimize such a residue is left exactly as in the input, also with a neutron_option (its H/D are not converted), and is listed in the log and the description file. With approach=add or place and keep_existing_H=True, its existing H/D are kept unchanged and it is listed in the log and the description file. With keep_existing_H=False it gets no H and is listed in the log.
skip_optimizations_after_seconds = 0.0
  .type = float
  .short_caption = Skip optimizations after this many seconds
  .help = Once the program has run for this many seconds, skip any further coarse optimizations.  Note that for many runs most of the time is spent on PDB interpretation and Mover placement, which will not be cut short. Fine optimizations are completed after coarse optimizations cease.  If this is 0.0, no time limit is imposed.
qmready = False
  .type = bool
  .short_caption = Run with custom preset options suitable for QM
  .help = Add H to water, use longer XH distances, add OXT and terminal NH3 at all chain ends and breaks, \
          and the terminal hydroxyl H of nucleic acids (HO3-prime at a 3-prime end, HO5-prime at a \
          5-prime OH end)
capping = False
  .type = bool
  .short_caption = Cap dangling bonds with H
  .help = Add capping H to every non-H atom whose valence is unsatisfied \
          because a bonded non-H atom is missing (truncated side chains, \
          incomplete ligands): as many H as the missing bonds need, never \
          a non-H atom. A missing atom between two present atoms of a \
          residue (a ring atom, CG of an ARG with CD) stops with a Sorry, \
          as their caps would overlap: rebuild it. So does a P, or an S \
          with three or more substituents, missing an atom, as an H there \
          would change the chemistry: add the atom (OP3 at a 5-prime end). \
          Done as the very last \
          step. The caps have names unknown to the restraint libraries, so \
          the capped model is for QM use and cannot be processed for \
          restraints. Bond orders come from the CCD Kekule structure \
          (DOUB and TRIP), all other bonds count as single. Requires \
          qmready=True.
capping_guess_for_unknowns = False
  .type = bool
  .short_caption = Guess the valence of residues without a dictionary
  .help = With capping=True, residues without a dictionary (no monomer \
          library, CCD or user CIF entry) cannot be checked. True: guess \
          the valence of every atom from its element and formal charge, bonds \
          from distances, and add H accordingly (reported as guessed). \
          False: stop with a Sorry listing those residues.
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
  ret += '@spherelist {het M} color= gray radius= 0.5 nobutton master= {hets}\n'
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

def _ModelNumbers(hierarchy):
  '''MODEL number (model.id, PDB MODEL / mmCIF pdbx_PDB_model_num) of each
  model, in order; None for an id that is not a number. A single model
  without a number is MODEL 1.
  '''
  models = hierarchy.models()
  if len(models) == 1 and models[0].id.strip() == '':
    return [1]
  result = []
  for m in models:
    try:
      result.append(int(m.id))
    except ValueError:
      result.append(None)
  return result

def _PdbFormatProblems(hierarchy):
  '''The reasons found in the hierarchy why fits_in_pdb_format() (hybrid-36)
  is False, as text; empty if none of them is found.
  '''
  chains = sorted(set(c.id for c in hierarchy.chains() if len(c.id) > 2))
  resnames = sorted(set(ag.resname for ag in hierarchy.atom_groups()
    if len(ag.resname) > 3))
  found = []
  if resnames:
    found.append('residue name%s longer than 3 characters: %s' % (
      's' if len(resnames) > 1 else '', ', '.join(resnames)))
  if chains:
    found.append('chain ID%s longer than 2 characters: %s' % (
      's' if len(chains) > 1 else '', ', '.join(chains)))
  return '; '.join(found)

def _SelectModelNumber(hierarchy, number):
  '''The model of the hierarchy with MODEL number number (in _ModelNumbers()).
  '''
  return hierarchy.models()[_ModelNumbers(hierarchy).index(number)]

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
  A .pdb request for a model that cannot be written in PDB format gives mmCIF
  with the same base name (.cif), and the log says so.
  If output.filename is not specified, the output file will be
  written into the current working directory with the same base name and type as the
  original file and with FH or H added to the base name (FH when flips are requested);
  1xs0.pdb would be written to ./1xsoH.pdb and 1xso.cif to ./1xsoH.cif by default.
  With neutron_option other than none, _<option> follows: ./1xsoH_all_d.pdb.

NOTES:
  If multiple alternates are present in the file and a specific one is not specified on the
  command line, they will all be processed in reverse order, such that the lowest-named
  one (perhaps A) will be processed last.  The hydrogen addition and arrangements for
  residues that are not part of any alternate will be left in the configuration that is best
  for the final alternate tested.  This may leave other alternates in sub-optimal configurations.
  When a single alternate is selected using alt_id= on the command line, everything is
  optimized a single time and for that configuration.

  neutron_option writes a model for neutron refinement. Every option other than none
  switches on use_neutron_distances (qmready=True does too). Exchangeable H are H on a
  non-carbon atom. A split site is written as two altlocs, H and D, each at half the
  parent's occupancy. Per option: non-exchangeable H, exchangeable H, water:
    all_h     H      H      HOH
    all_d     D      D      DOD
    hd_and_h  H      split  DOD
    hd_and_d  D      split  DOD
    all_hd    split  split  DOD
  In an existing alternate conformer, exchangeable H are written as D (and in all_hd the
  other H as H), not split, and are listed. H whose names do not start with H stay H.
  Input D is used only when the H/D are kept (approach=optimize, or keep_existing_H=True).
  It must be unambiguous: there is no repair, and anything else stops with a Sorry that
  includes a hint. Kept D input with neutron_option=none stops. Refined occupancies and
  B of input H/D pairs are kept where the output writes the same pair, new pairs get half
  the parent occupancy each. The log and the description file have an H/D section.
  Limitation: after refinement, H and D of a rotatable group (OH, SH, NH3+) can point in
  different directions. If they are more than 0.2 A apart, the input is refused. Use
  approach=add keep_existing_H=False instead. With ignore_missing_restraints=True,
  approach=optimize leaves residues without restraints exactly as in the input (not
  converted) and lists them. With approach=add or place and keep_existing_H=True, their
  existing H/D are kept unchanged and they are listed. With keep_existing_H=False such a
  residue gets no H. Cannot be combined with approach=remove or qmready=True.

  qmready=True prepares a model for QM: H on water, longer XH distances, OXT and NH3+ at
  chain ends and breaks, HO3' and HO5' at nucleic-acid chain ends. Adding capping=True to it also puts H on every non-H atom left with
  an unsatisfied valence by a missing non-H atom (truncated side chains, incomplete ligands),
  as the very last step; the capped model is for QM use only, because the cap names are
  unknown to the restraint libraries.

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
    # placement queries the server of the model it gets: run() starts it
    # fresh and it holds the run's restraint objects, so it is not stale
    if self.params.qmready:
      reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
        model                 = self.model,
        use_neutron_distances = self.params.use_neutron_distances,
        n_terminal_charge     = "all_chain_starts",
        exclude_water         = False,
        stop_for_unknowns     = self.params.stop_on_any_missing_hydrogen,
        keep_existing_H       = True,
        mon_lib_srv           = self.model.get_mon_lib_srv(),
        na_terminal_h         = True,
        metal_bound_amines_as_nh2 = self.params.metal_bound_amines_as_nh2,
      )
    else:
      reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
        model = self.model,
        use_neutron_distances=self.params.use_neutron_distances,
        n_terminal_charge=self.params.n_terminal_charge,
        exclude_water = self.params.exclude_water,
        stop_for_unknowns=self.params.stop_on_any_missing_hydrogen,
        keep_existing_H=self.params.keep_existing_H,
        mon_lib_srv=self.model.get_mon_lib_srv(),
        metal_bound_amines_as_nh2=self.params.metal_bound_amines_as_nh2
      )
    reduce_add_h_obj.supplied_h_hint = self._RebuildHint()
    # the model is boxed for the calculation: contacts across its own symmetry
    reduce_add_h_obj.input_crystal_symmetry = self._output_cs
    reduce_add_h_obj.run()
    reduce_add_h_obj.show(self.logger, verbose=self.params.verbosity > 2)
    # its report of atom names not in the dictionary, for the description
    self._unknown_atoms_note = ''
    if reduce_add_h_obj.residues_unknown_atoms:
      self._unknown_atoms_note = '\n' + reduce_add_h_obj.unknown_atoms_text(
        verbose=self.params.verbosity > 2) + '\n'
    # its report of thiols and N-H at metals, for the description
    self._metal_thiol_note = ''
    for text in (reduce_add_h_obj.metal_thiol_text(),
                 reduce_add_h_obj.metal_nh_text(),
                 reduce_add_h_obj.metal_amine_text()):
      if text:
        self._metal_thiol_note += '\n' + text + '\n'
    missed_residues = set(reduce_add_h_obj.no_H_placed_mlq)
    if not self.params.ignore_missing_restraints:
      if len(missed_residues) > 0:
        raise self._MissingRestraintsSorry(missed_residues)
    insufficient_restraints = list(reduce_add_h_obj.site_labels_no_para)
    if self.params.stop_on_any_missing_hydrogen and len(insufficient_restraints) > 0:
      bad = insufficient_restraints[0]
      for res in insufficient_restraints[1:]:
        bad += "," + res
      raise Sorry("Insufficient restraints were found for the following atoms:"+bad)

    self.model = reduce_add_h_obj.get_model()
    # input pairs keep their values per site when placement renamed their H
    self._hd_records = hd_exchange.rekey_records(self._hd_records,
      reduce_add_h_obj.renamed_existing_h)

    # Placement replaced the D of an input DOD water by H but kept the name
    self._dod_renamed_note = ''
    n_renamed = 0
    for rg in self.model.get_hierarchy().residue_groups():
      water = [ag for ag in rg.atom_groups()
               if common_residue_names_get_class(name=ag.resname) ==
                  'common_water']
      elements = set(a.element.strip().upper() for ag in water
                     for a in ag.atoms())
      if 'H' not in elements or 'D' in elements: continue
      dod = [ag for ag in water if ag.resname == 'DOD']
      if not dod: continue
      for ag in dod:
        ag.resname = 'HOH'
      n_renamed += 1
    if n_renamed > 0:
      self._dod_renamed_note = ('%d water(s) renamed DOD -> HOH (their D were '
        'replaced by H)\n' % n_renamed)
      print(self._dod_renamed_note, file=self.logger)

    if not self.model.has_hd():
      raise Sorry("It was not possible to place any H atoms. Is this a single atom model?")

    # We must re-interpret the model when _type_energies is not defined, which happens
    # for example when running on 8zst.cif.
    if not hasattr(self.model, '_type_energies'):
      self._ReinterpretModel(make_restraints=True)

# ------------------------------------------------------------------------------

  def _MissingRestraintsSorry(self, resnames):
    return Sorry("Restraints were not found for the following residues:" +
      "".join(" " + r.strip() for r in sorted(resnames)))

# ------------------------------------------------------------------------------

  def _PrepareRestraints(self):
    # Once per run, before the flipkin snapshot: one mon_lib_query() per
    # residue name. CCD-only residues get placement's CCD-generated restraint
    # objects, kept for the whole run (main run, input check, placement,
    # flipkin passes; a server returns them only the first time it sees a
    # residue) and dropped with _DropAutoRestraints(). Residues without
    # restraints (as placement's no_H_placed_mlq) stop the run unless
    # ignore_missing_restraints=True; then all their atoms are left as in
    # the input.
    # :return: group_args(keys (atom keys of the residues left unchanged),
    # labels).
    # the objects are appended to the model's list: its own copy, so that
    # the input's list stays as it is
    ro = self.model.get_restraint_objects()
    if ro is not None:
      self.model.set_restraint_objects(list(ro))
    r = reduce_hydrogen.add_ccd_restraint_objects(self.model)
    self._auto_restraint_names = r.auto_restraint_names
    if r.no_restraints and not self.params.ignore_missing_restraints:
      raise self._MissingRestraintsSorry(r.no_restraints)
    hierarchy = self.model.get_hierarchy()
    n_models = len(hierarchy.models())
    missing = set()   # (residue_group memory_id, resname)
    labels = []
    for model_index, m in enumerate(hierarchy.models()):
      for chain in m.chains():
        for rg in chain.residue_groups():
          for resname in sorted(set(ag.resname for ag in rg.atom_groups())
                                & r.no_restraints):
            missing.add((rg.memory_id(), resname))
            label = '%s %s %s%s' % (resname.strip(), chain.id,
              rg.resseq.strip(), rg.icode.strip())
            if n_models > 1: label += ' (model %d)' % (model_index + 1)
            labels.append(label)
    keys = set()
    if missing:
      for a, k in zip(hierarchy.atoms(), hd_exchange.atom_keys(hierarchy)):
        if (a.parent().parent().memory_id(), a.parent().resname) in missing:
          keys.add(k)
    return group_args(keys=keys, labels=labels)

# ------------------------------------------------------------------------------

  def _ElementXText(self, selection):
    # Atoms with element X, removed at loading: one line, then the atoms,
    # capped at 20 unless verbosity > 2. Printed here.
    # :param selection: The element X atoms of self.model.
    # :return: The description-file text.
    atoms = list(self.model.get_hierarchy().atoms().select(selection))
    cap = None if self.params.verbosity > 2 else 20
    lines = ['', 'Removed %d atom(s) with element X: the element is unknown, so '
      'H cannot be placed around them; they are not in the output.' %
      len(atoms)]
    shown = atoms if cap is None else atoms[:cap]
    for a in shown:
      ag = a.parent()
      rg = ag.parent()
      label = '%s %s %s%s %s' % (ag.resname.strip(), rg.parent().id.strip(),
        rg.resseq.strip(), rg.icode.strip(), a.name.strip())
      if ag.altloc.strip():
        label += ' altloc %s' % ag.altloc
      lines.append('  ' + label)
    if len(shown) < len(atoms):
      lines.append('  (+%d more)' % (len(atoms) - len(shown)))
    text = '\n'.join(lines) + '\n'
    print(text, file=self.logger)
    return text

# ------------------------------------------------------------------------------

  def _UnchangedText(self):
    # Residues without restraints (ignore_missing_restraints=True) that the
    # run leaves as in the input: optimize, and placement keeping the input
    # H/D (keep_existing_H, qmready). Printed here.
    # :return: The description-file text ('' if none).
    if not self._unrestrained.labels: return ''
    if self.params.approach != 'optimize' and \
       not (self.params.keep_existing_H or self.params.qmready):
      return ''
    text = ('\nResidues without restraints, left unchanged:\n' +
      '\n'.join(self._unrestrained.labels) + '\n')
    print(text, file=self.logger)
    return text

# ------------------------------------------------------------------------------

  def _DropAutoRestraints(self, model):
    # the CCD-generated restraint objects of _PrepareRestraints(), off a
    # model that leaves the run (its restraints stay), so that output and
    # downstream validation do not see them (92b7fa5beb)
    reduce_hydrogen.remove_auto_restraint_objects(model,
      self._auto_restraint_names)

# ------------------------------------------------------------------------------

  def _UnrestrainedKeys(self):
    # (also before run(): the checks are called directly in tests)
    unrestrained = getattr(self, '_unrestrained', None)
    if unrestrained is None: return set()
    return unrestrained.keys

# ------------------------------------------------------------------------------

  def _RemoveFloatingHD(self):
    # Remove H/D atoms that have no covalent bond partner in the restraints,
    # except in residues without restraints (left unchanged).
    # :return: Tuple of the description-file text listing the removed atoms
    # ('' if none) and their atom keys (hd_exchange.atom_keys()).
    partners = hd_exchange.covalent_h_partners(self.model)
    atoms = self.model.get_hierarchy().atoms()
    skip = self._UnrestrainedKeys()
    keys = hd_exchange.atom_keys(self.model.get_hierarchy()) if skip else None
    sel = flex.bool(atoms.size(), False)
    for i, a in enumerate(atoms):
      if a.element.strip().upper() in ('H', 'D') and i not in partners and \
         not (skip and keys[i] in skip):
        sel[i] = True
    if sel.count(True) == 0:
      return '', []
    if keys is None:
      keys = hd_exchange.atom_keys(self.model.get_hierarchy())
    removed = [keys[i] for i in sel.iselection()]
    # in residues with atom names not in the dictionary: listed under that
    # cause (_UnknownAtomsText())
    unknown = set(r for r, names in getattr(self, '_unknown_atoms', []))
    labels = []
    for i in sel.iselection():
      label = atoms[i].id_str().replace('pdb=', '').replace('"', '')
      if reduce_hydrogen.residue_label(atoms[i]) in unknown:
        self._unknown_atoms_removed.append(label)
      else:
        labels.append(label)
    text = ''
    if labels:
      text = ('\nH/D without a covalent partner, removed before optimization:'
              '\n' + '\n'.join(labels) + '\n')
      print(text, file=self.logger)
    self.model = self.model.select(~sel)
    return text, removed

# ------------------------------------------------------------------------------

  def _UnknownAtomsText(self):
    # optimize: placement's report of atom names not in the dictionary, with
    # the H/D removed in those residues (at the input check and as floating)
    # instead of in the floating lists. Printed here.
    # :return: The description-file text ('' if none); in the log without
    # its leading blank line (a sub-header or text with one precedes it).
    if not self._unknown_atoms: return ''
    resids = set(r for r, names in self._unknown_atoms)
    moved = [l for l, r in zip(self._input_floating,
             self._input_floating_resids) if r in resids]
    kept = [(l, r) for l, r in zip(self._input_floating,
            self._input_floating_resids) if r not in resids]
    self._input_floating = [l for l, r in kept]
    self._input_floating_resids = [r for l, r in kept]
    labels = moved + self._unknown_atoms_removed
    text = '\n' + reduce_hydrogen.atoms_not_in_dictionary_text(
      self._unknown_atoms, len(labels), 'H/D on them were removed',
      labels) + '\n'
    print(text[1:], file=self.logger)
    return text

# ------------------------------------------------------------------------------

  def _RemoveFromSnapshot(self, snapshot, removed):
    # Remove the atoms with the given atom keys from an unprocessed snapshot
    # of the model, keeping its atom order (keys, not indices: the processed
    # model is re-sorted). Each key must occur exactly once.
    # :return: The snapshot without those atoms.
    keys = hd_exchange.atom_keys(snapshot.get_hierarchy())
    index = dict((k, i) for i, k in enumerate(keys))
    sel = flex.bool(len(keys), False)
    for key in removed:
      if key not in index:
        raise Sorry("Internal error: removed atom not in the initial model: "
          + str(key))
      sel[index[key]] = True
    return snapshot.select(~sel)

# ------------------------------------------------------------------------------

  def _CheckAmbiguousParents(self):
    # neutron_option: an H/D with more than one covalent partner cannot be
    # classified; stop before the Optimizer runs, naming the restraints.
    partners = hd_exchange.covalent_h_partners(self.model)
    atoms = self.model.get_hierarchy().atoms()
    skip = self._UnrestrainedKeys()
    keys = hd_exchange.atom_keys(self.model.get_hierarchy()) if skip else None
    bad = []
    for i in sorted(partners):
      if len(partners[i]) < 2: continue
      if skip and keys[i] in skip: continue
      lines = [' ' + atoms[i].id_str()]
      for j in sorted(partners[i]):
        for info in partners[i][j]:
          lines.append('    %s: %s, ideal %.3f A, %s' % (atoms[j].id_str(),
            info.origin, info.distance_ideal, info.sym_op))
      bad.append('\n'.join(lines))
    if bad:
      raise Sorry('neutron_option=%s: H/D with more than one covalent bond '
        'partner (check the restraints):\n%s' % (self.params.neutron_option,
        '\n'.join(bad)))

# ------------------------------------------------------------------------------

  def _CaptureParentTable(self):
    # Parent of every H/D while the restraints match the hierarchy; the
    # output step uses it after the restraints are gone.
    self._h_parents = hd_exchange.h_parent_table(self.model)

# ------------------------------------------------------------------------------

  def _SuppliedMetalHText(self):
    # approach=optimize keeps every input H: the H on thiols and N bound to
    # metals, and the other N-H next to a metal, as placement lists them
    # (nothing is deleted here).
    # :return: The description-file text ('' if none); also in the log.
    cs = reduce_hydrogen.contact_symmetry(self._output_cs)
    thiols = reduce_hydrogen.thiol_metal_sites(self.model, cs)
    nh = reduce_hydrogen.metal_nh_sites(self.model, cs)
    result = ''
    for text in (
        reduce_hydrogen.metal_thiol_text(
          supplied = reduce_hydrogen.supplied_metal_h(thiols, lambda a: True),
          hint = self._RebuildHint()),
        reduce_hydrogen.metal_nh_text(
          supplied = reduce_hydrogen.supplied_metal_h(nh, lambda a: True),
          reported = reduce_hydrogen.reported_metal_h(nh),
          hint = self._RebuildHint())):
      if not text: continue
      print('\n' + text, file=self.logger)
      result += '\n' + text + '\n'
    return result

# ------------------------------------------------------------------------------

  def _RebuildHint(self):
    # The parameters that discard the input H/D for this run: qmready and
    # optimize always keep them (add also optimizes), add/place keep them
    # only on request.
    if self.params.qmready:
      return 'remove the H/D first (approach=remove), then run qmready=True'
    if self.params.approach == 'optimize':
      return 'approach=add keep_existing_H=False'
    return 'keep_existing_H=False'

# ------------------------------------------------------------------------------

  def _MatchInputHD(self, model=None):
    # Candidates and collapse check on an interpreted model (self.model if
    # None).
    # :return: Tuple (candidates, match).
    if model is None: model = self.model
    table = hd_exchange.h_parent_table(model)
    candidates = hd_exchange.find_hd_candidates(model.get_hierarchy(),
      table, model=model)
    return candidates, hd_exchange.collapse_check(candidates)

# ------------------------------------------------------------------------------

  def _NormalizeInput(self):
    # Input with D that the run keeps: every H/D must be a single or a
    # collapsible pair (no repair). Floating H/D are dropped. Pairs are
    # collapsed to H with their values recorded (self._hd_records), single D
    # renamed to H, DOD -> HOH. Residues without restraints stop the run, or
    # with ignore_missing_restraints=True are left as they are. self.model is
    # rebuilt from the H-only hierarchy, with restraints for optimize (the
    # branch uses them).
    # :return: Unprocessed model of the H-only hierarchy, the snapshot for
    # the flipkin passes.
    make_sub_header('Checking input H/D', out=self.logger)
    if self.params.neutron_option == 'none':
      raise Sorry('input contains deuterium; set neutron_option, or rebuild '
        'H: ' + self._RebuildHint())
    self._ReinterpretModel(make_restraints=True)
    skip = self._unrestrained.keys
    def checked_part(full):
      # the model without the residues left unchanged, and the full model's
      # atom keys of its atoms (its own keys can differ: chain occurrence)
      if not skip: return full, None
      keys = hd_exchange.atom_keys(full.get_hierarchy())
      sel = flex.bool([k not in skip for k in keys])
      return full.select(sel), [k for k, x in zip(keys, sel) if x]
    full = self.model
    part, part_keys = checked_part(full)
    candidates, match = self._MatchInputHD(part)
    clean = lambda id_str: id_str.replace('pdb=', '').replace('"', '')
    titles = {
      'pair': 'H/D pairs that cannot be collapsed',
      'different parents': 'Same-name H and D on different parents',
      'same conformer': 'H and D of one site in one conformer',
      'name collision': 'D whose H name already exists'}
    categories = []
    def add(title, line):
      for t, lines in categories:
        if t == title:
          lines.append(line)
          return
      categories.append((title, [line]))
    for x in match.problems:
      ids = [clean(i) for i in (x.h_id, x.d_id) if i is not None]
      reasons = '; '.join(x.reasons)
      add(titles[x.kind], ' / '.join(ids) + ': ' + reasons if ids else reasons)
    for u in match.unresolved:
      if u.action != 'stop': continue
      add('H/D with more than one covalent partner', clean(u.atom_id) + ': ' +
        '; '.join('%s %s, ideal %.3f A, %s' % (clean(s[0]), s[1], s[2], s[3])
                  for s in u.source))
    if categories:
      text = []
      for title, lines in categories:
        text.append('  %s (%d):' % (title, len(lines)))
        text.extend('    ' + l for l in lines[:20])
        if len(lines) > 20:
          text.append('    (+%d more)' % (len(lines) - 20))
      raise Sorry('The input H/D cannot be used with approach=%s (no repair):'
        '\n%s\nRebuild H: %s' % (self.params.approach, '\n'.join(text),
        self._RebuildHint()))
    # floating: dropped (select() keeps the restraints), then matched again
    floating = [u for u in candidates.unresolved if u.status == 'floating']
    self._input_floating = [clean(u.atom_id) for u in floating]
    # residue_label() from the atom key, for _UnknownAtomsText()
    self._input_floating_resids = ['%s %s %s' % (u.atom_key[6].strip(),
      u.atom_key[1].strip(), u.atom_key[3].strip()) for u in floating]
    if floating:
      drop = set(u.atom_key for u in floating)
      if part_keys is not None:
        drop = set(fk for k, fk in zip(
          hd_exchange.atom_keys(part.get_hierarchy()), part_keys) if k in drop)
      keys = hd_exchange.atom_keys(full.get_hierarchy())
      full = full.select(flex.bool([k not in drop for k in keys]))
      part, part_keys = checked_part(full)
      candidates, match = self._MatchInputHD(part)
    self.model = full
    self._input_notes = [group_args(output=clean(n.h_id) + ' / ' +
      clean(n.d_id), text=n.text) for n in match.notes]
    h = self.model.get_hierarchy().deep_copy()
    # atom_groups left unchanged are taken out while the matched part is
    # normalized (its atom order is then the part's), and put back
    detached = []
    if skip:
      key_of = dict((a.memory_id(), k) for a, k in
                    zip(h.atoms(), hd_exchange.atom_keys(h)))
      for rg in h.residue_groups():
        for index, ag in enumerate(list(rg.atom_groups())):
          if key_of[ag.atoms()[0].memory_id()] in skip:
            detached.append((rg, index, ag))
      for rg, index, ag in detached:
        rg.remove_atom_group(ag)
      assert h.atoms_size() == part.get_hierarchy().atoms_size()
    normalized = hd_exchange.normalize_to_h(h, match)
    for rg, index, ag in detached:
      rg.insert_atom_group(min(index, rg.atom_groups_size()), ag)
    if detached:
      h.sort_atoms_in_place()
      h.atoms().reset_serial()
      h.atoms().reset_i_seq()
    self._hd_records = normalized.records
    c = normalized.counts
    print('Input H/D: %d pairs collapsed to H, %d D renamed to H, %d DOD '
      'renamed to HOH, %d floating dropped' % (c.n_pairs, c.n_d_renamed,
      c.n_dod, len(floating)), file=self.logger)
    snapshot = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = h.deep_copy(),
      crystal_symmetry  = self.model.crystal_symmetry(),
      restraint_objects = self.model.get_restraint_objects(),
      stop_for_unknowns = self.params.stop_on_any_missing_hydrogen,
      log               = null_out())
    self._normalized_with_restraints = self.params.approach == 'optimize'
    self.model = self._RebuildModel(h,
      make_restraints=self._normalized_with_restraints)
    return snapshot

# ------------------------------------------------------------------------------

  def _NeutronDescription(self):
    # H/D section for the description file and the log: distances, counts
    # and the listed sites with their reasons, by their written names. Lists
    # are capped at 20 entries unless verbosity > 2.
    r = self._neutron_result
    cap = None if self.params.verbosity > 2 else 20
    lines = ['', 'H/D model (neutron_option=%s):' % self.params.neutron_option,
      '  Neutron distances used (%s)' % self._neutron_distance_reason,
      '  Exchangeable H sites (not water): %d' % r.n_exchangeable,
      '  Split into H/D pairs: %d (new, at 0.5 x parent occupancy: %d)' % (
        r.n_split, r.n_split_new),
      '  Converted to D: %d' % r.n_to_d,
      '  Waters renamed DOD: %d, HOH: %d' % (r.n_water, r.n_water_to_hoh),
      '  Values restored from the input: %d' % r.n_restored]
    def section(title, entries, detail, remark=None):
      if not entries: return
      lines.append('  %s (%d):' % (title, len(entries)))
      if remark: lines.append('    ' + remark)
      shown = entries if cap is None else entries[:cap]
      for e in shown:
        lines.append(('    %s  %s' % (e.output, detail(e))).rstrip())
      if len(shown) < len(entries):
        lines.append('    (+%d more)' % (len(entries) - len(shown)))
    by_reason = lambda entries, reason: [e for e in entries if e.reason == reason]
    section('Conformer sites written as one isotope',
      by_reason(r.conformer_sites, 'conformer'), lambda e: e.reason)
    section('Not split, no free altloc letters',
      by_reason(r.conformer_sites, 'no free altloc letters'), lambda e: e.reason)
    section('Name not pairable by phenix.refine, kept as H',
      r.unclassified_sites, lambda e: e.reason)
    section('Conformer copies merged into one site', r.merged_sites,
      lambda e: '%s (%d copies)' % (e.reason, len(e.copies)))
    section('Split sites on a partially occupied parent',
      r.partial_parent_sites, lambda e: 'parent occupancy %.2f' % e.parent_occ,
      remark='New pairs get half the parent occupancy each, restored pairs '
             'keep their input values; phenix.refine currently refines '
             'H + D = 1.')
    section('Occupancy notes', r.notes, lambda e: e.text)
    section('Input H/D notes (values not changed)', self._input_notes,
      lambda e: e.text)
    section('H/D without a covalent partner, dropped at the input check',
      [group_args(output=label) for label in self._input_floating],
      lambda e: '')
    def record_label(key):
      # site key: residue part, parent altloc, name
      model_index, chain_id, _, resseq, icode, _, resname, altloc, name = key
      label = '%s %s %s %s%s' % (name.strip(), resname, chain_id,
        resseq.strip(), icode.strip())
      if altloc: label += ' (parent altloc %s)' % altloc
      if model_index > 0: label += ' (model %d)' % (model_index + 1)
      return label
    section('Input values not used (site no longer in the model)',
      [group_args(output=record_label(k)) for k in r.unused_records],
      lambda e: '')
    text = '\n'.join(lines) + '\n'
    # the same section in the log, so it is not lost without files
    print(text, file=self.logger)
    return text

# ------------------------------------------------------------------------------

  def _CheckNeutronOutput(self):
    # The written H/D model against the parent table: provenance is a
    # bijection, every pair collapses, exchangeable_hd_selections() finds
    # exactly the split pairs, no residue_group gained a resname, and water
    # names match their H/D.
    out_h = self._output_model.get_hierarchy()
    provenance = self._neutron_result.provenance
    atoms = out_h.atoms()
    keys = hd_exchange.atom_keys(out_h)
    hd = [k for k, a in zip(keys, atoms)
          if a.element.strip().upper() in ('H', 'D')]
    problems = []
    for k, a in zip(keys, atoms):
      if a.element.strip().upper() in ('H', 'D') and k not in provenance:
        problems.append('no provenance: ' + a.id_str())
    for k in sorted(set(provenance) - set(hd)):
      problems.append('provenance names a missing atom: ' + str(k))
    # residues without restraints: in the provenance, not paired (their
    # atoms are taken out; the keys of the rest can change)
    skip = self._UnrestrainedKeys()
    pair_h, pair_provenance = out_h, provenance
    if not problems and skip:
      sel = flex.bool([k not in skip for k in keys])
      pair_h = out_h.select(sel, copy_atoms=True)
      pair_h.atoms().reset_i_seq()
      pair_provenance = dict((k, provenance[fk]) for k, fk in zip(
        hd_exchange.atom_keys(pair_h), [k for k, x in zip(keys, sel) if x])
        if fk in provenance)
    if not problems:
      match = hd_exchange.collapse_check(hd_exchange.find_hd_candidates(
        pair_h, self._h_parents, provenance=pair_provenance))
      for x in match.problems:
        problems.append('%s %s: %s' % (x.h_id, x.d_id, '; '.join(x.reasons)))
      for x in match.unresolved:
        problems.append('%s: %s' % (x.atom_id, x.status))
      pairs = set(frozenset((c.i_h, c.i_d)) for c in match.pairs)
      written = set(frozenset((s[0][0], s[1][0]))
                    for s in pair_h.exchangeable_hd_selections())
      if len(pairs) != self._neutron_result.n_split or pairs != written:
        problems.append('split pairs: %d written, %d paired, %d found by '
          'exchangeable_hd_selections()' % (self._neutron_result.n_split,
          len(pairs), len(written)))
    rgs_in = list(self.model.get_hierarchy().residue_groups())
    rgs_out = list(out_h.residue_groups())
    if len(rgs_in) != len(rgs_out):
      problems.append('%d residue groups in, %d out' % (len(rgs_in),
        len(rgs_out)))
    else:
      for rg_in, rg_out in zip(rgs_in, rgs_out):
        n_in = len(set(ag.resname for ag in rg_in.atom_groups()))
        n_out = len(set(ag.resname for ag in rg_out.atom_groups()))
        if n_in != n_out:
          problems.append('resnames changed in %s' %
            rg_out.atoms()[0].id_str())
    # water names follow the content: with D DOD, with H not DOD
    for rg in rgs_out:
      water = [ag for ag in rg.atom_groups()
               if common_residue_names_get_class(name=ag.resname) ==
                  'common_water']
      elements = set(a.element.strip().upper() for ag in water
                     for a in ag.atoms())
      names = set(ag.resname for ag in water)
      if ('D' in elements and names != set(['DOD'])) or \
         ('H' in elements and 'DOD' in names):
        problems.append('water name %s does not match its H/D: %s' % (
          '/'.join(sorted(names)), rg.atoms()[0].id_str()))
    if problems:
      raise Sorry('Internal error: H/D output check failed:\n  ' +
        '\n  '.join(problems))

# ------------------------------------------------------------------------------

  def _ReinterpretModel(self, make_restraints=True, model=None):
    # Reinterpret the model using the same approach that hydrogen-placement does.
    # :param make_restraints: Should we compute restraints during the interpretation?
    # :param model: The model to interpret; self.model if None.
    if model is None:
      model = self.model
    model.get_hierarchy().sort_atoms_in_place()
    model.get_hierarchy().atoms().reset_serial()

    p = reduce_hydrogen.get_reduce_pdb_interpretation_params(
      self.params.use_neutron_distances)
    p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check=True
    # We need to turn this on because without it 1zz0.txt kept flipping the ring
    # in A TYR 214 every time we re-interpreted. The original interpretation done
    # by Hydrogen placement will have flipped them, so we don't need to do it again.
    p.pdb_interpretation.flip_symmetric_amino_acids=False
    #p.pdb_interpretation.sort_atoms=True
    model.set_stop_for_unknowns(self.params.stop_on_any_missing_hydrogen)
    # no ADP pair table, as in place_hydrogens
    model.process(make_restraints=make_restraints, pdb_interpretation_params=p,
      plain_pairs_radius=None)

# ------------------------------------------------------------------------------

  def _RebuildModel(self, hierarchy, make_restraints):
    # Turn a modified hierarchy back into a model: working crystal symmetry,
    # user restraint objects and stop_for_unknowns of self.model, interpreted
    # with the same settings as _ReinterpretModel().
    # :return: The new model.
    model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = hierarchy,
      crystal_symmetry  = self.model.crystal_symmetry(),
      restraint_objects = self.model.get_restraint_objects(),
      stop_for_unknowns = self.params.stop_on_any_missing_hydrogen,
      log               = null_out())
    self._ReinterpretModel(make_restraints=make_restraints, model=model)
    return model

# ------------------------------------------------------------------------------

  def _FinalizeModel(self, hierarchy):
    # The model that ends the run: a fresh, unprocessed manager built from the
    # finished hierarchy (crystal symmetry, restraint objects and
    # stop_for_unknowns as _RebuildModel()). Replaces a full re-interpretation
    # with make_restraints=False; the steps below are the parts of it the
    # output depends on. No restraints manager, as before: that was the guard
    # against the crash in pair_proxies() when writing mmCIF after atoms were
    # removed, and it leaves no links to write (C4). No monomer mappings:
    # selection() uses the atom-selection cache.
    # :param hierarchy: The finished hierarchy; changed in place.
    # :return: The new model.
    # the sort and serials of _ReinterpretModel()
    hierarchy.sort_atoms_in_place()
    hierarchy.atoms().reset_serial()
    cs = self.model.crystal_symmetry()
    if cs is not None:
      # process() rounds coordinates and ADPs through an X-ray structure
      # (pdb_interpretation.py ~4270), visible in the written bytes. Extracted
      # here, after the sort: adopt_xray_structure() matches by order.
      hierarchy.adopt_xray_structure(
        hierarchy.extract_xray_structure(crystal_symmetry=cs))
    # process()'s _check_for_capitalised_element(): adopt gives "Ca"
    for a in hierarchy.atoms():
      if a.element != a.element.upper():
        a.element = a.element.upper()
    model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = hierarchy,
      crystal_symmetry  = cs,
      restraint_objects = self.model.get_restraint_objects(),
      stop_for_unknowns = self.params.stop_on_any_missing_hydrogen,
      log               = null_out())
    # the writer's symmetry choice and shift back (_figure_out_cs_to_output())
    if self.model.shift_cart() is not None:
      model.set_unit_cell_crystal_symmetry_and_shift_cart(
        unit_cell_crystal_symmetry = self.model.unit_cell_crystal_symmetry(),
        shift_cart = self.model.shift_cart())
    # as select() does: info() (file names) and the input format, which
    # write_model_file(format=Auto) follows
    model.set_info(copy.deepcopy(self.model.info()))
    model._original_model_format = self.model._original_model_format
    return model

# ------------------------------------------------------------------------------

  def _CheckQmreadyClashes(self, model):
    # qmready=True: Sorry if the finished model has a non-bonded contact
    # below 1.0 A. Checked on a copy with restraints (self.model has none);
    # H included, as geometry_statistics() drops H with the X-ray table.
    # :param model: The rebuilt copy of self.model with restraints.
    from cctbx import geometry_restraints
    sites = model.get_sites_cart()
    nb = model.get_restraints_manager().geometry.pair_proxies(
      sites_cart=sites).nonbonded_proxies
    # simple proxies first, then symmetry (asu) proxies
    d = geometry_restraints.nonbonded_deltas(sites_cart=sites,
      sorted_asu_proxies=nb)
    if d.size() == 0 or flex.min(d) >= 1.0:
      return
    k = flex.min_index(d)
    n = nb.simple.size()
    if k < n:
      i, j = nb.simple[k].i_seqs
      where = ''
    else:
      i, j = nb.asu[k-n].i_seq, nb.asu[k-n].j_seq
      where = ' (symmetry mate)'
    atoms = model.get_hierarchy().atoms()
    raise Sorry("Bad clash(es) prevents qmready=True: %s - %s%s %.3f A"
      % (atoms[i].id_str(), atoms[j].id_str(), where, d[k]))

# ------------------------------------------------------------------------------

  def _ReturnedModel(self, model):
    # The returned model gets the crystal symmetry of the written file, also
    # with output.write_files=False: the input cell, the cryo-EM placeholder
    # (unit cell symmetry only) or none. A model on reduce2's box is rebuilt
    # unprocessed, as if read from that file: a processed model without
    # symmetry fails in selection(). Shifted input: the returned model is put
    # back into the caller's frame (shift, working box, unit cell symmetry).
    # :return: The model to return.
    cs = self._output_cs
    if (cs is not None) and (not cs.is_nonsense()):
      model.set_unit_cell_crystal_symmetry(cs)
      new = model
    else:
      new = mmtbx.model.manager(
        model_input       = None,
        pdb_hierarchy     = model.get_hierarchy(),
        crystal_symmetry  = None,
        restraint_objects = model.get_restraint_objects(),
        stop_for_unknowns = self.params.stop_on_any_missing_hydrogen,
        log               = null_out())
      if cs is not None:
        new.set_unit_cell_crystal_symmetry(cs)
    if self._input_shift is not None:
      shift_cart, box_cs, unit_cell_cs = self._input_shift
      new.shift_model_and_set_crystal_symmetry(shift_cart=shift_cart,
        crystal_symmetry=box_cs, unit_cell_crystal_symmetry=unit_cell_cs)
    return new

# ------------------------------------------------------------------------------

  def _CalculateTotalCharge(self, model, extra_sigma_bonds=None,
                            parity_i_seqs=()):
    # Calculate the total charge of the model (same as Q|R qr.finalise with
    # calculate_charge=True) and, if the output file names were not specified,
    # insert it into them: model.pdb -> modelH_charge-3.pdb
    # :param model: The rebuilt copy of self.model with restraints (the
    # clash-check copy); not changed.
    # :param extra_sigma_bonds: {i_seq: n} of model, the capping H to come
    # (capping=True): the charge is that of the capped model.
    # :param parity_i_seqs: i_seqs of model to check for a leftover electron
    # (every atom in extra_sigma_bonds).
    # :return: Warning lines, one per atom in parity_i_seqs left open-shell.
    make_sub_header('Calculating total charge', out=self.logger)
    from mmtbx.ligands import electrons
    # On reduce2's own rebuilt copy (interpreted with the user CIFs and the
    # CCD objects of _PrepareRestraints()), not a second process() inside
    # electrons; values verified unchanged against that old path.
    distribution = electrons.run(
      model                 = model,
      extra_sigma_bonds     = extra_sigma_bonds,
      return_formal_charges = True,
    )
    self.total_charge = distribution.get_total_charge()
    print('Total charge:', self.total_charge, file=self.logger)
    if self._default_output_filename:
      base, suffix = os.path.splitext(self.params.output.filename)
      self.params.output.filename = "%s_charge%s%s"%(
        base, str(self.total_charge), suffix)
      if self._default_description_file_name:
        base = os.path.splitext(self.params.output.description_file_name)[0]
        self.params.output.description_file_name = "%s_charge%s.txt"%(
          base, str(self.total_charge))
    # a capped atom with an electron left over: the capped fragment is
    # open-shell (an odd number of electrons)
    atoms = model.get_hierarchy().atoms()
    warnings = []
    for i_seq in parity_i_seqs:
      if distribution.get(i_seq, 0):
        warnings.append('open-shell cap on %s' % atoms[i_seq].id_str())
    return warnings

# ------------------------------------------------------------------------------

  def _CapDanglingAtoms(self, model, dangling, parity=()):
    # capping=True: add the capping H of find_dangling_atoms() to model, the
    # model the records come from, as the last change to the hierarchy (the
    # caps have names no restraint library knows). Sorry if a cap clashes.
    # The output model is an unprocessed manager on the capped hierarchy and
    # nothing interprets it later.
    # :param model: The rebuilt copy with restraints; its hierarchy is capped
    # in place, its X-ray structure and restraints are stale afterwards.
    # :param dangling: find_dangling_atoms(model).
    # :param parity: Open-shell warning lines of _CalculateTotalCharge().
    # :return: The Capping section of the description file (also printed to
    # the log). Lists are capped at 20 entries unless verbosity > 2.
    make_sub_header('Capping dangling bonds', out=self.logger)
    result = capping.cap_dangling_atoms(model, dangling,
      use_neutron_distances=self.params.use_neutron_distances)
    self._output_model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = model.get_hierarchy(),
      crystal_symmetry  = self.model.crystal_symmetry(),
      restraint_objects = self.model.get_restraint_objects(),
      stop_for_unknowns = False,
      log               = null_out())
    cap = None if self.params.verbosity > 2 else 20
    lines = ['', 'Capping']
    def section(entries, title=None):
      if not entries: return
      indent = '  '
      if title is not None:
        lines.append('  ' + title)
        indent = '    '
      shown = entries if cap is None else entries[:cap]
      for e in shown:
        lines.append(indent + e)
      if len(shown) < len(entries):
        lines.append('%s(+%d more)' % (indent, len(entries) - len(shown)))
    # cap_dangling_atoms() reports in record order
    section([line for record, line in zip(dangling.records, result.lines)
             if record.source != 'guess'])
    section([line for record, line in zip(dangling.records, result.lines)
             if record.source == 'guess'],
      'Capped by element valence (no dictionary):')
    section(['%s: %s' % (label, reason) for label, reason in dangling.skipped],
      'Not capped:')
    section(list(result.warnings) + list(parity))
    lines.append('  Capped %d atoms in %d residues with %d H' % (
      result.n_atoms, result.n_residues, result.n_added))
    text = '\n'.join(lines) + '\n'
    # the same section in the log, so it is not lost without files
    print(text, file=self.logger)
    return text

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
    # residues without restraints, left unchanged: their H have no bonds
    skip = self._UnrestrainedKeys()
    keys = None
    for a in self.model.get_atoms():
      inWater[a] = common_residue_names_get_class(name=a.parent().resname) == "common_water"
      inHet[a] = hetatm_sel[a.i_seq]
      if not a.element_is_hydrogen():
        inMainChain[a] = mainchain_sel[a.i_seq]
      else:
        # Check our bonded neighbor to see if it is on the mainchain if we are a Hydrogen
        if len(bondedNeighborLists[a]) < 1:
          if skip:
            if keys is None:
              keys = hd_exchange.atom_keys(self.model.get_hierarchy())
            if keys[a.i_seq] in skip:
              # classified by its residue (het or water), not main chain
              inMainChain[a] = False
              inSideChain[a] = sidechain_sel[a.i_seq]
              continue
          raise Sorry("Found Hydrogen with no neigbors.")
        else:
          inMainChain[a] = mainchain_sel[bondedNeighborLists[a][0].i_seq]
      inSideChain[a] = sidechain_sel[a.i_seq]
    return (inWater, inHet, inMainChain, inSideChain)

# ------------------------------------------------------------------------------

  def _ProbeModel(self):
    # probe2 in the flipkin passes stops on H without a bonded neighbour and
    # has no D radius: it gets a copy without the H/D of the residues left
    # unchanged (placement dropped them up to this change).
    # :return: self.model, or that copy.
    skip = self._UnrestrainedKeys()
    if not skip: return self.model
    atoms = self.model.get_hierarchy().atoms()
    keys = hd_exchange.atom_keys(self.model.get_hierarchy())
    sel = flex.bool([not (a.element_is_hydrogen() and k in skip)
                     for a, k in zip(atoms, keys)])
    if sel.all_eq(True): return self.model
    return self.model.select(sel)

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
    p = self.params
    # phil turns an explicit neutron_option=none into None
    if p.neutron_option is None:
      p.neutron_option = 'none'
    if p.neutron_option != 'none':
      if p.approach == 'remove':
        raise Sorry('neutron_option=%s cannot be used with approach=remove.'
          % p.neutron_option)
      if p.qmready:
        raise Sorry('neutron_option=%s cannot be used with qmready=True.'
          % p.neutron_option)

    # capping: the last step of a qmready run that adds H
    if p.capping:
      if not p.qmready:
        raise Sorry('capping=True requires qmready=True.')
      if p.approach == 'remove':
        raise Sorry('capping=True cannot be used with approach=remove.')
      # residues without a dictionary reach the capping step, where
      # capping_guess_for_unknowns decides
      if not p.ignore_missing_restraints:
        p.ignore_missing_restraints = True
        print('capping=True: ignore_missing_restraints set to True (residues '
          'without a dictionary are handled by capping)', file=self.logger)
    elif p.capping_guess_for_unknowns:
      raise Sorry('capping_guess_for_unknowns=True requires capping=True.')

    # One distance flag for the whole run: qmready and neutron_option place H
    # at neutron distances, so every later step must use them
    self._neutron_distance_reason = None
    if p.use_neutron_distances:
      self._neutron_distance_reason = 'use_neutron_distances=True'
    else:
      reason = None
      if p.qmready:
        reason = 'qmready=True'
      elif p.neutron_option != 'none':
        reason = 'neutron_option=%s' % p.neutron_option
      if reason is not None:
        p.use_neutron_distances = True
        self._neutron_distance_reason = reason
        print('%s: use_neutron_distances set to True' % reason, file=self.logger)

    # place: no optimizer, so nothing that needs it
    if self.params.approach == 'place':
      p = self.params
      bad = [name for name, value in [
        ('add_flip_movers', p.add_flip_movers),
        ('output.flipkin_directory', p.output.flipkin_directory),
        ('comparison_file', p.comparison_file),
        ('qmready', p.qmready)] if value]
      if bad:
        raise Sorry('%s requires optimization; it cannot be used with '
          'approach=place.' % ', '.join(bad))

    # Set the default output file name if one has not been given.
    # With qmready=True the total charge is appended to the default name later.
    self._default_output_filename = self.params.output.filename is None
    self._default_description_file_name = \
      self.params.output.description_file_name is None
    if self.params.output.filename is None:
      inName = self.data_manager.get_default_model_name()
      if inName is None: raise Sorry('model not found')
      suffix = os.path.splitext(os.path.basename(inName))[1]
      if self.params.add_flip_movers:
        pad = 'FH'
      else:
        pad = 'H'
      if self.params.neutron_option != 'none':
        pad += '_' + self.params.neutron_option
      base = os.path.splitext(os.path.basename(inName))[0] + pad
      self.params.output.filename = base + suffix
      print('Writing model output to', self.params.output.filename, file=self.logger)

    if os.environ.get('PHENIX_OVERWRITE_ALL', False):
      self.data_manager.set_overwrite(True)

    self.data_manager.has_models(raise_sorry=True)
    if self.params.output.description_file_name is None:
      self.params.output.description_file_name=self.params.output.filename.replace('.pdb',
                                                                                   '.txt')
      self.params.output.description_file_name=self.params.output.description_file_name.replace('.cif',
                                                                                   '.txt')

    # Existing output files stop the run here, not after it. Names not known
    # yet keep the Sorry from _write_text() at write time: the qmready default
    # names (_charge<N> appended later) and the mmCIF written for a .pdb
    # request that does not fit PDB format.
    if self.params.output.write_files and not self.data_manager.get_overwrite():
      charge_pending = (self.params.qmready and self.params.approach == 'add'
        and self._default_output_filename)
      names = []
      if not charge_pending:
        names.append(self.params.output.filename)
      if not (charge_pending and self._default_description_file_name):
        names.append(self.params.output.description_file_name)
      existing = [n for n in names if os.path.exists(n)]
      if existing:
        raise Sorry('%s already exist%s; use overwrite=True to replace %s.' % (
          ' and '.join(existing), 's' if len(existing) == 1 else '',
          'it' if len(existing) == 1 else 'them'))

    # model_id is a MODEL number of the input
    if self.params.model_id is not None:
      numbers = _ModelNumbers(self.data_manager.get_model().get_hierarchy())
      if self.params.model_id <= 0 or self.params.model_id not in numbers:
        raise Sorry("model_id=%d matches no model. Available MODEL numbers: %s"
          % (self.params.model_id,
             ", ".join(str(n) for n in numbers if n is not None) or "none"))

    # Turn on profiling if we've been asked to in the Phil parameters
    if self.params.profile:
      import cProfile
      self._pr = cProfile.Profile()
      self._pr.enable()

  # ------------------------------------------------------------------------------

  def run(self):

    # Set our bonded-neighbor depth
    self._bondedNeighborDepth = self.params.bonded_neighbor_depth

    # Input H/D check (neutron_option): records of collapsed pairs, notes and
    # floating H/D dropped there
    self._hd_records = {}
    self._input_notes = []
    self._input_floating = []
    self._input_floating_resids = []
    # optimize: atoms_not_in_dictionary() and the H/D removed on them
    self._unknown_atoms = []
    self._unknown_atoms_removed = []
    self._normalized_with_restraints = False
    # residues without restraints (optimize, D input): None until checked
    self._unrestrained = None
    self._auto_restraint_names = []

    # String describing the run that will be output to the specified file.
    outString = 'reduce2 v.{}, run {}\n'.format(version, datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for a in sys.argv:
      outString += ' {}'.format(a)
    outString += '\n'

    make_sub_header('Loading Model', out=self.logger)

    # Get our model.
    self.model = self.data_manager.get_model()

    sel_x = self.model.selection('element X')
    if sel_x.count(True) > 0:
      outString += self._ElementXText(sel_x)
    self.model = self.model.select(~sel_x)
    # A model shifted in memory (e.g. boxed with a map) is shifted back on
    # this copy: placement builds a new model without the shift, so the run
    # is the file run. _ReturnedModel() puts the shift back. Any non-zero
    # shift, not shifted() (1e-3 A tolerance): set_unit_cell_crystal_symmetry()
    # asserts exactly (0, 0, 0).
    self._input_shift = None
    shift_cart = self.model.shift_cart()
    if shift_cart is not None and tuple(shift_cart) != (0, 0, 0):
      self._input_shift = (shift_cart,
        self.model.crystal_symmetry(), self.model.unit_cell_crystal_symmetry())
      self.model.shift_model_back()
    # CA and/or P trace; waters and ions do not count
    names = set()
    for ag in self.model.get_hierarchy().atom_groups():
      if common_residue_names_get_class(name=ag.resname) in (
          'common_water', 'common_element'): continue
      names.update(n.strip() for n in ag.atoms().extract_name())
      if not names <= {'CA', 'P'}: break
    if names == {'CA'}:
      raise Sorry("The model contains only CA atoms (CA trace); H atoms cannot be placed.")
    if names == {'P'}:
      raise Sorry("The model contains only P atoms (P trace); H atoms cannot be placed.")
    if names == {'CA', 'P'}:
      raise Sorry("The model contains only CA and P atoms (CA/P trace); H atoms cannot be placed.")

    # Use model function to set crystal symmetry if necessary 2025-03-19 TT
    # The box is for the calculation only and is not written out. Cushion 5:
    # H added later overhang the heavy-atom box; 3 left 4.4 A to the nonbonded
    # cutoff of 3.4 A. Before add_terminal_oxt: a hierarchy change without
    # symmetry leaves an xray_structure that boxing then asserts against.
    self._output_cs = reduce_hydrogen.get_output_crystal_symmetry(
      self.data_manager.get_model())
    self.model.add_crystal_symmetry_if_necessary(box_cushion=5)

    # qmready=True
    if self.params.qmready:
      # add OXT
      self.model.add_terminal_oxt()
      # make sure model is atom-complete. With capping=True a listing only:
      # the atoms left with a dangling bond are capped at the end
      if not self.model.is_atom_complete():
        missing = self.model.missing_atoms_by_residue()
        print("Missing atoms found (this may or may not be an issue):", file=self.logger)
        for k, v in missing.items():
          print("  where: '%s' missing:"%k,v, file=self.logger)
        if not self.params.capping:
          raise Sorry(
            "Model misses non-H atoms. Model must be atom-complete for qmready=true.")
        print("capping=True: continuing, dangling bonds are capped with H at the end",
          file=self.logger)

    if self.data_manager.has_restraints():
      self.model.set_stop_for_unknowns(self.params.stop_on_any_missing_hydrogen)
      self.model.process(make_restraints=False)

    # If we've been asked to only to a single model from the file, strip the model down to
    # only that one (checked in validate()); it keeps its MODEL number.
    if self.params.model_id is not None:
      make_sub_header('Selecting Model ID ' + str(self.params.model_id), out=self.logger)

      # Construct a hierarchy for the selected model
      r = pdb.hierarchy.root()
      mdc = _SelectModelNumber(self.model.get_hierarchy(),
        self.params.model_id).detached_copy()
      r.append_model(mdc)

      # Make yet another model for the new hierarchy, with the user's
      # restraint CIFs (a new list: _PrepareRestraints() appends to it)
      subset_model_manager = mmtbx.model.manager(
        model_input       = None,
        pdb_hierarchy     = r,
        stop_for_unknowns = False,
        crystal_symmetry  = self.model.crystal_symmetry(),
        restraint_objects = list(self.model.get_restraint_objects() or []),
        log               = None)

      self.model = subset_model_manager

    # CCD restraint objects attached once, for the whole run (also the
    # flipkin passes, which restart from initialModel); residues without
    # restraints stop the run here, or are left unchanged
    if self.params.approach != 'remove':
      # A model passed in memory (e.g. from an earlier run) gives the same
      # result as read from a file. Its cached monomer library server returns
      # no CCD object for a CCD-only residue it has seen (no setter in
      # mmtbx.model.manager; get_mon_lib_srv() builds a new one), and its
      # riding-H manager, which setup_riding_h_manager() keeps, has the
      # earlier run's X-H lengths; a scattering table set on it (kept by
      # process()) sets the water O-H length when use_neutron_distances is
      # off; refinement flags restrict the riding-H idealization to their
      # sites. A model read from a file has none of these.
      self.model._mon_lib_srv = None
      self.model.unset_riding_h_manager()
      self.model.get_xray_structure(force=True)
      self.model.set_refinement_flags(None)
      self._unrestrained = self._PrepareRestraints()

    # Input with D that the run keeps (never for remove; qmready forces
    # keep_existing_H=True in _AddHydrogens()): checked and normalized to H;
    # the H-only model is the snapshot for the flipkins
    snapshot = None
    if self.params.approach != 'remove' and \
       (self.params.keep_existing_H or self.params.qmready or
        self.params.approach == 'optimize') and \
       any(e.strip().upper() == 'D'
           for e in self.model.get_hierarchy().atoms().extract_element()):
      snapshot = self._NormalizeInput()

    # Stores the initial coordinates for all of the atoms and the rest of the information
    # about the original model for use by Kinemages.
    if snapshot is None:
      initialModel = self.model.deep_copy()
    else:
      initialModel = snapshot

    neutron = self.params.neutron_option != 'none'
    if self.params.approach == 'place':
      make_sub_header('Adding Hydrogens', out=self.logger)
      startAdd = time.time()
      self._AddHydrogens()
      outString += self._dod_renamed_note
      outString += self._unknown_atoms_note
      outString += self._metal_thiol_note
      outString += self._UnchangedText()
      if neutron:
        self._CheckAmbiguousParents()
        self._CaptureParentTable()
      outString += 'Time to Add Hydrogen = {:.3f} sec'.format(
        time.time()-startAdd)+'\n'

    elif self.params.approach == 'add' or self.params.approach == 'optimize':

      if self.params.approach == 'add':
        # Add Hydrogens to the model
        make_sub_header('Adding Hydrogens', out=self.logger)
        startAdd = time.time()
        self._AddHydrogens()
        doneAdd = time.time()
        outString += self._dod_renamed_note
        outString += self._unknown_atoms_note
        outString += self._metal_thiol_note
        outString += self._UnchangedText()
        if neutron:
          self._CheckAmbiguousParents()
      else:
        # We need restraints on the model for optimization, so we interpret it to get them.
        # (normalized input was rebuilt with them)
        if not self._normalized_with_restraints:
          self._ReinterpretModel(make_restraints=True)
        outString += self._UnchangedText()
        # atom names not in the dictionary, as placement reports them
        self._unknown_atoms = reduce_hydrogen.atoms_not_in_dictionary(
          self.model)
        # H/D without a covalent partner are removed, as placement does; same
        # rule as the H parent table. Also keeps them out of probe scoring.
        text, removed = self._RemoveFloatingHD()
        outString += text
        outString += self._UnknownAtomsText()
        # and from the snapshot the flipkin passes restart from
        if removed:
          initialModel = self._RemoveFromSnapshot(initialModel, removed)
        # nothing to optimize: the Optimizer has no riding H manager then
        if not self.model.has_hd():
          raise Sorry('The model has no H or D atoms to optimize: '
            'approach=optimize does not add H. Use approach=add to add and '
            'optimize them.')
        if neutron:
          self._CheckAmbiguousParents()
        outString += self._SuppliedMetalHText()
        # a supplied metal-bound NH2 stays as it is: out of riding H before
        # the Optimizer builds and idealizes its manager
        if self.params.metal_bound_amines_as_nh2:
          kept, partly = reduce_hydrogen.protect_supplied_metal_amines(
            self.model, reduce_hydrogen.contact_symmetry(self._output_cs),
            use_neutron_distances = True if self.params.use_neutron_distances
              else None)
          text = reduce_hydrogen.metal_amine_text(supplied_nh2 = kept,
            supplied_nh2_partly = partly)
          if text:
            print('\n' + text, file=self.logger)
            outString += '\n' + text + '\n'

      # NOTE: We always optimize all models (leave modelIndex alone) because we've removed all
      # but the desired model ID structure from the model.
      make_sub_header('Optimizing', out=self.logger)
      startOpt = time.time()
      stopTime = 0.0
      if self.params.skip_optimizations_after_seconds > 0.0:
        stopTime = startOpt + self.params.skip_optimizations_after_seconds
      opt = Optimizers.Optimizer(self.params.probe, self.params.add_flip_movers,
        self.model, altID=self.params.alt_id,
        preferenceMagnitude=self.params.preference_magnitude,
        bondedNeighborDepth = self._bondedNeighborDepth,
        useNeutronDistances=self.params.use_neutron_distances,
        nonFlipPreference=self.params.non_flip_preference,
        skipBondFixup=self.params.skip_bond_fix_up,
        flipStates = self.params.set_flip_states,
        verbosity=self.params.verbosity,
        cliqueOutlineFileName=self.params.output.clique_outline_file_name,
        fillAtomDump = self.params.output.print_atom_info,
        stopOptimizingAt = stopTime,
        crystalSymmetry = reduce_hydrogen.contact_symmetry(self._output_cs))
      doneOpt = time.time()
      warnings = opt.getWarnings()
      # Find the lines from warnings that start with "Warning: " and
      # concatenate them into a single string.
      lines = [line for line in warnings.split('\n') if line.startswith('Warning:')]
      # one line, by residue, for the atoms without an energy type (7a5m: 39)
      prefix = 'Warning: Could not find atom info for '
      no_info = [l for l in lines if l.startswith(prefix) and ' (perhaps' in l]
      if no_info:
        lines = [l for l in lines if l not in no_info]
        residues = []
        for l in no_info:
          tokens = l[len(prefix):l.index(' (perhaps')].split()[:-1]
          if len(tokens) == 3: tokens = [tokens[1], tokens[0], tokens[2]]
          if ' '.join(tokens) not in residues: residues.append(' '.join(tokens))
        shown = ', '.join(residues[:12])
        if len(residues) > 12: shown += ' (+%d more)' % (len(residues) - 12)
        lines.append('Warning: %d atoms have no energy type (element radius used, '
          'not H-bond acceptors) - %d residues: %s' % (
          len(no_info), len(residues), shown))
      warnings = '\n'.join(lines)
      if len(warnings) > 0:
        print('\nWarnings during optimization:\n'+warnings, file=self.logger)
      outString += opt.getInfo()

      # the restraints still match the hierarchy here: before the deletions
      # and before the final reinterpretation without restraints
      if neutron:
        self._CaptureParentTable()

      # Delete any hydrogens that we've been asked to delete, adding this to the info.
      hToDelete = opt.getHydrogensToDelete()
      if len(hToDelete) > 0:
        outString += ' Deleting hydrogens requested for deletion by optimization:\n'
        for a in hToDelete:
          aName = a.name.strip().upper()
          chainID = a.parent().parent().parent().id
          resName = a.parent().resname.strip().upper()
          resID = str(a.parent().parent().resseq_as_int())
          altLoc = a.parent().altloc
          # Don't print the code if it is a space (blank).
          insertionCode = a.parent().parent().icode.strip()
          resNameAndID = "chain "+str(chainID)+" "+altLoc+resName+" "+resID+insertionCode
          outString += "  Deleting {} {}\n".format(resNameAndID, aName)

          a.parent().remove_atom(a)

      if self.params.approach == 'add':
        outString += 'Time to Add Hydrogen = {:.3f} sec'.format(doneAdd-startAdd)+'\n'
      outString += 'Time to Optimize = {:.3f} sec'.format(doneOpt-startOpt)+'\n'
      if self.params.output.print_atom_info:
        print('Atom information used during calculations:', file=self.logger)
        print(opt.getAtomDump(), file=self.logger)

    else: # Removing Hydrogens from the model rather than adding them.
      # Own path, independent of the isotope logic: any H/D layout, no
      # matcher, keep_existing_H ignored
      make_sub_header('Removing Hydrogens', out=self.logger)
      sel = self.model.selection("element H or element D")
      for a in self.model.get_atoms():
        if sel[a.i_seq]:
          a.parent().remove_atom(a)
      # H/D-only atom_groups (and residue_groups, chains) are now empty;
      # waters without D are HOH
      for m in self.model.get_hierarchy().models():
        for chain in list(m.chains()):
          for rg in list(chain.residue_groups()):
            for ag in list(rg.atom_groups()):
              if ag.atoms_size() == 0:
                rg.remove_atom_group(ag)
              elif ag.resname == 'DOD':
                ag.resname = 'HOH'
            if rg.atom_groups_size() == 0:
              chain.remove_residue_group(rg)
          if chain.residue_groups_size() == 0:
            m.remove_chain(chain)

    # Rebuild the model without restraints because we have removed some atoms
    # that were previously bonded. Like the re-interpretation without restraints
    # this replaces, it keeps us from crashing on a call to pair_proxies when
    # generating mmCIF files, so we always do it for safety.
    self.model = self._FinalizeModel(self.model.get_hierarchy())

    # qmready=True (remove adds no H): one copy with restraints for the clash
    # check (the finished model must be clash-free), the dangling atoms
    # (capping=True) and the total charge, as qr.finalise does with
    # calculate_charge=True, added to the default output file name
    # (approach=add). With capping=True the caps go on that copy as the
    # last change to the hierarchy and it is the output model.
    self.total_charge = None
    parity = []
    if self.params.qmready and self.params.approach != 'remove':
      qm_model = self._RebuildModel(self.model.get_hierarchy().deep_copy(),
        make_restraints=True)
      self._CheckQmreadyClashes(qm_model)
      dangling = None
      if self.params.capping:
        dangling = capping.find_dangling_atoms(qm_model,
          self.params.capping_guess_for_unknowns, log=self.logger)
        if dangling.unknown_residues and \
           not self.params.capping_guess_for_unknowns:
          raise Sorry('capping=True: no dictionary for %s; set '
            'capping_guess_for_unknowns=True to cap them by element valence'
            % ', '.join(dangling.unknown_residues))
      if self.params.approach == 'add':
        extra = None
        if dangling is not None:
          # the caps, and the bonds guessing found by distance
          extra = capping.extra_sigma_bonds(dangling)
        parity = self._CalculateTotalCharge(qm_model, extra,
          sorted(extra) if extra else ())
      if self.params.capping:
        outString += self._CapDanglingAtoms(qm_model, dangling, parity)

    # The model that is written and returned; self.model stays H-only for
    # the comparison and the flipkins.
    # With capping=True it is set by _CapDanglingAtoms(), not interpreted
    # again.
    self._neutron_result = None
    if not self.params.capping:
      if not neutron:
        self._output_model = self.model
      else:
        make_sub_header('Writing H/D model (neutron_option=%s)' %
          self.params.neutron_option, out=self.logger)
        h = self.model.get_hierarchy().deep_copy()
        self._neutron_result = hd_exchange.apply_neutron_option(h,
          self.params.neutron_option, parent_table=self._h_parents,
          records=self._hd_records, unchanged_keys=self._UnrestrainedKeys())
        self._output_model = self._FinalizeModel(h)
        self._DropAutoRestraints(self._output_model)
        self._CheckNeutronOutput()
        outString += self._NeutronDescription()

    make_sub_header('Writing output', out=self.logger)

    # Skip writing the main output file and description output file if output.write_files is False.
    # This enables a program to harness Reduce2 and not have to deal with handling output files.
    if self.params.output.write_files:
      # Write the description output to the specified file.
      self.data_manager._write_text("description", outString,
        self.params.output.description_file_name)

      # Determine whether to write a PDB or CIF file and write the appropriate text output.
      suffix = os.path.splitext(self.params.output.filename)[1]
      output_cs = self._output_cs is not None
      if output_cs:
        self._output_model.set_unit_cell_crystal_symmetry(self._output_cs)
      if suffix.lower() == ".pdb" and not self._output_model.can_be_output_as_pdb():
        # model_as_pdb() would return only header records; write mmCIF instead
        requested = self.params.output.filename
        self.params.output.filename = self.data_manager.write_model_file(
          self._output_model, os.path.splitext(requested)[0] + '.cif',
          format='cif', output_cs=output_cs)
        cause = _PdbFormatProblems(self._output_model.get_hierarchy())
        print('Model cannot be represented in PDB format%s; writing mmCIF to %s '
          'instead of %s.' % (' (%s)' % cause if cause else '',
          self.params.output.filename, requested), file=self.logger)
      else:
        if suffix.lower() == ".pdb":
          txt = self._output_model.model_as_pdb(output_cs=output_cs)
        else:
          txt = self._output_model.model_as_mmcif(output_cs=output_cs)
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

        # Interpret the model to fill in things we need for determining the neighbor list
        # and for the sidechain selection of the views (the finished model is not interpreted).
        self._ReinterpretModel()

        # Find the viewpoint locations for each Mover we're going to
        # look at.
        views = self._GetViews(amides)

        carts = flex.vec3_double()
        for a in self.model.get_atoms():
          carts.append(a.xyz)
        bondProxies, asuProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)
        bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies, asuProxies)

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
            useNeutronDistances=self.params.use_neutron_distances,
            nonFlipPreference=self.params.non_flip_preference,
            skipBondFixup=self.params.skip_bond_fix_up,
            flipStates = flipStates,
            verbosity=3,
            crystalSymmetry = reduce_hydrogen.contact_symmetry(self._output_cs))
          print('Results of optimization:', file=self.logger)
          print(opt.getInfo(), file=self.logger)
          self._ReinterpretModel()

          # Get the other characteristics we need to know about each atom to do our work.
          # We must do this again here because the atoms change after Hydrogen addition.
          # We also need to re-generate the bonded-neighbor lists for the same reason.
          carts = flex.vec3_double()
          for a in self.model.get_atoms():
            carts.append(a.xyz)
          bondProxies, asuProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)
          bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies, asuProxies)
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

          p2.overrideModel(self._ProbeModel())
          dots, kinString = p2.run()
          flipkinText += kinString

        # Write the accumulated Flipkin string to the output file.
        with open(flipkinBase+"-flipnq.kin", "w") as f:
          f.write(flipkinText)

      # ===========================================================================
      hists = _FindFlipsInOutputString(outString, 'HisFlip')

      if len(hists) > 0:
        make_sub_header('Generating Histidine Flipkin', out=self.logger)

        # Interpret the model to fill in things we need for determining the neighbor list
        # and for the sidechain selection of the views (the finished model is not interpreted).
        self._ReinterpretModel()

        # Find the viewpoint locations for each Mover we're going to
        # look at.
        views = self._GetViews(hists)

        carts = flex.vec3_double()
        for a in self.model.get_atoms():
          carts.append(a.xyz)
        bondProxies, asuProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)
        bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies, asuProxies)

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
            useNeutronDistances=self.params.use_neutron_distances,
            nonFlipPreference=self.params.non_flip_preference,
            skipBondFixup=self.params.skip_bond_fix_up,
            flipStates = flipStates,
            verbosity=3,
            crystalSymmetry = reduce_hydrogen.contact_symmetry(self._output_cs))
          print('Results of optimization:', file=self.logger)
          print(opt.getInfo(), file=self.logger)
          self._ReinterpretModel()

          # Get the other characteristics we need to know about each atom to do our work.
          # We must do this again here because the atoms change after Hydrogen addition.
          # We also need to re-generate the bonded-neighbor lists for the same reason.
          carts = flex.vec3_double()
          for a in self.model.get_atoms():
            carts.append(a.xyz)
          bondProxies, asuProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)
          bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies, asuProxies)
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
          p2.overrideModel(self._ProbeModel())
          dots, kinString = p2.run()
          flipkinText += kinString

        # Write the accumulated Flipkin string to the output file.
        with open(flipkinBase+"-fliphis.kin", "w") as f:
          f.write(flipkinText)

    # self.model kept them for the flipkin passes; with neutron_option=none
    # the output model is the self.model of before the flipkin passes
    self._DropAutoRestraints(self.model)
    self._DropAutoRestraints(self._output_model)
    self._output_model = self._ReturnedModel(self._output_model)

    # Report profiling info if we've been asked to in the Phil parameters
    if self.params.profile:
      print('Profile results:', file=self.logger)
      import pstats
      profile_params = {'sort_by': 'time', 'num_entries': 20}
      self._pr.disable()
      ps = pstats.Stats(self._pr).sort_stats(profile_params['sort_by'])
      ps.print_stats(profile_params['num_entries'])

    print('', file=self.logger)
    print('peak memory:', human_readable_size(peak_memory_usage()),
      file=self.logger)

# ------------------------------------------------------------------------------

  def get_results(self):
    '''
    :return: group_args with model (the output model) and total_charge (qmready
    add, otherwise None). The returned model is not guaranteed to be
    interpreted. Callers needing interpretation-dependent operations (e.g.
    backbone or sidechain selections) should check processed() and call
    process(make_restraints=False) if necessary.
    '''
    return group_args(
      model        = getattr(self, '_output_model', self.model),
      total_charge = getattr(self, 'total_charge', None))

# ------------------------------------------------------------------------------

  def Test(self):
    '''
      Run tests on the methods of the class.  Throw an assertion error if there is a problem with
      one of them and return normally if there is not a problem.
    '''

    #=====================================================================================
    # @todo Unit tests for other methods
