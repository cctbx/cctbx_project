##################################################################################
#                Copyright 2021-2023  Richardson Lab at Duke University
#
# Licensed under the Apache License, Version 2.0 (the "License"],
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

##################################################################################
# This module exports functions that are helpful to create data structures
# needed by Probe.

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import
import sys
import math
import traceback
import iotbx.map_model_manager
import iotbx.data_manager
import cctbx.maptbx.box
import mmtbx
import scitbx.matrix
from scitbx.array_family import flex
from mmtbx.probe import AtomTypes
from iotbx import pdb
from libtbx.utils import null_out, Sorry

import boost_adaptbx.boost.python as bp
bp.import_ext("mmtbx_probe_ext")
import mmtbx_probe_ext as probeExt

##################################################################################
# Helper string and functions for using Probe PHIL parameters.  This lets Python
# code that want to use Probe PHIL parameters pass these to the various object
# constructors in the C++ ext library.

# PHIL parameters to be added to the master PHIL string.  They make a subobject
# 'probe' and put all the parameters inside it.
probe_phil_parameters = """
probe
  .style = menu_item auto_align
{
  probe_radius = 0.25
    .type = float
    .short_caption = Probe radius
    .help = Probe radius (half distance between touched atoms) (-radius in probe)

  density = 16.0
    .type = float
    .short_caption = Dot density
    .help = Probe dots per square angstrom on atom surface (-density in probe)

  worse_clash_cutoff = 0.5
    .type = float
    .short_caption = Worse clash cutoff
    .help = Cutoff for worse clashes, a positive (-divworse in probe)

  clash_cutoff = 0.4
    .type = float
    .short_caption = Clash cutoff
    .help = Cutoff for the clashes, a positive number (-divlow in probe)

  contact_cutoff = 0.25
    .type = float
    .short_caption = Contact cutoff
    .help = Cutoff for the contact (-divhigh in probe)

  uncharged_hydrogen_cutoff = 0.6
    .type = float
    .short_caption = Uncharged hydrogen cutoff
    .help = Cutoff for uncharged hydrogen overlap (-hbregular in probe)

  charged_hydrogen_cutoff = 0.8
    .type = float
    .short_caption = Charged hydroen cutoff
    .help = Cutoff for charged hydrogen overlap (-hbcharged in probe)

  bump_weight = 10.0
    .type = float
    .short_caption = Bump weight
    .help = Weight applied to bump score (-bumpweight in probe)

  hydrogen_bond_weight = 4.0
    .type = float
    .short_caption = Hydrogen bond weight
    .help = Weight applied to hydrogen bond score (-hbweight in probe)

  gap_weight = 0.25
    .type = float
    .short_caption = Gap weight
    .help = Weight applied to gap score (-gapweight in probe). Should be no smaller than probe_radius.

  allow_weak_hydrogen_bonds = False
    .type = bool
    .short_caption = Count weak hydrogen bonds
    .help = Separately account for weak hydrogen bonds (-LweakHbonds in probe)

  ignore_ion_interactions = False
    .type = bool
    .short_caption = Ignore ion interactions
    .help = Ignore interactions with ions

  set_polar_hydrogen_radius = True
    .type = bool
    .short_caption = Override polar hydrogen radius
    .help = Override the radius of polar Hydrogens with to ensure it has the expected value. (-usepolarh in probe)
}
"""

def createSpatialQuery(atoms, probePhil):
  """
    Helper function to create a SpatialQuery object, passing it all required Phil parameters
    in addition to the other constructor parameters.
    :param atoms: Flex array of atoms to insert in the structure.  This will determine the
    bounds of the structure.
    :param probePhi: Subobject of PHIL parameters for probe.  Can be obtained using
    self.params.probe from a Program Template program that includes the probe_phil_parameters
    from above in its master PHIL parameters string.
    :returns A mmtbx_probe_ext.SpatialQuery object.
  """
  return probeExt.SpatialQuery(atoms)

def createDotSphereCache(probePhil):
  """
    Helper function to create a DotSphereCache object, passing it all required Phil parameters
    in addition to the other constructor parameters.
    :param probePhi: Subobject of PHIL parameters for probe.  Can be obtained using
    self.params.probe from a Program Template program that includes the probe_phil_parameters
    from above in its master PHIL parameters string.
    :returns A mmtbx_probe_ext.DotSphereCache object.
  """
  return probeExt.DotSphereCache(probePhil.density)

def createDotScorer(extraAtomInfo, probePhil):
  """
    Helper function to create a DotScorer object, passing it all required Phil parameters
    in addition to the other constructor parameters. Also enforces the constraint that the
    gap_weight is at least as large as the probe_radius.
    :param extraAtomInfo: Can be obtained from getExtraAtomInfo().  Holds extra atom information
    needed for scoring atoms.
    :param probePhil: Subobject of PHIL parameters for probe.  Can be obtained using
    self.params.probe from a Program Template program that includes the probe_phil_parameters
    from above in its master PHIL parameters string.
    :returns A mmtbx_probe_ext.DotScorer object.
  """
  gap_weight = max(probePhil.gap_weight, probePhil.probe_radius)
  return probeExt.DotScorer(extraAtomInfo, gap_weight,
        probePhil.bump_weight, probePhil.hydrogen_bond_weight,
        probePhil.uncharged_hydrogen_cutoff, probePhil.charged_hydrogen_cutoff,
        probePhil.clash_cutoff, probePhil.worse_clash_cutoff,
        probePhil.contact_cutoff, probePhil.allow_weak_hydrogen_bonds,
        probePhil.ignore_ion_interactions)

##################################################################################
# Other helper functions.

def getBondedNeighborLists(atoms, bondProxies):
  """
    Helper function to produce a dictionary of lists that contain all bonded
    neighbors for each atom in a set of atoms.
    :param atoms: Flex array of atoms (could be obtained using model.get_atoms() if there
    are no chains with multiple conformations, must be a subset of the atoms including
    all in the base conformation and in a particular conformation otherwise).
    :param bondProxies: Flex array of bond proxies for the atoms.  This could be obtained
    using model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart =
    model.get_sites_cart())[0] if the model has only a single conformation.  Otherwise,
    it should be a flex array of atom positions for the atoms that are in the first parameter.
    It can include atoms that are not in the first parameter, but they will not be added
    to the lists.
    :returns a dictionary with one entry for each atom, indexed by i_seq, that contains a
    list of all of the atoms (within the atoms list) that are bonded to it.
  """
  atomDict = {}
  for a in atoms:
    atomDict[a.i_seq] = a
  bondedNeighbors = {}
  for a in atoms:
    bondedNeighbors[a] = []
  for bp in bondProxies:
    try:
      # These lookups will fail if the atoms are not in the list of atoms passed in.
      first = atomDict[bp.i_seqs[0]]
      second = atomDict[bp.i_seqs[1]]
      bondedNeighbors[first].append(second)
      bondedNeighbors[second].append(first)
    except Exception:
      # When an atom is bonded to an atom is not in our atom list (in a different conformer or not
      # in our selection) we just ignore it.
      pass
  return bondedNeighbors

def addIonicBonds(bondedNeighborLists, atoms, spatialQuery, extraAtomInfo):
  """
    Helper function to add ionic bonds to the list of bonded neighbors.
    @todo: Be sure never to make an ionic bond to a water (fix in old and new probe).
    :param bondedNeighborLists: Object to have the ionic bonds added to in place.
    :param atoms: Flex array of atoms (could be obtained using model.get_atoms() if there
    are no chains with multiple conformations, must be a subset of the atoms including
    :param spatialQuery: mmtbx_probe_ext.SpatialQuery structure to rapidly determine which atoms
    are within a specified distance of a location.
    :param extraAtomInfo: mmtbx_probe_ext.ExtraAtomInfo mapper that provides radius and other
    information about atoms beyond what is in the pdb.hierarchy.  Used here to determine
    which atoms may be acceptors.
  """
  for a in atoms:
    if a.element_is_positive_ion():
      myRad = extraAtomInfo.getMappingFor(a).vdwRadius
      minDist = myRad
      maxDist = 0.25 + myRad + 3  # overestimate so we don't miss any
      neighbors = spatialQuery.neighbors(a.xyz, minDist, maxDist)
      for n in neighbors:
        # Never add bonds with Phantom Hydrogens.
        if extraAtomInfo.getMappingFor(n).isDummyHydrogen:
          continue
        # See if we're within range for an ionic bond between the two atoms.
        dist = (rvec3(a.xyz) - rvec3(n.xyz)).length()
        expected = myRad + extraAtomInfo.getMappingFor(n).vdwRadius
        if dist >= (expected - 0.6) and dist <= (expected + 0.2):
          # We're in range; bond each of us to the other.
          try:
            bondedNeighborLists[a].append(n)
          except Exception:
            bondedNeighborLists[a] = [n]
          try:
            bondedNeighborLists[n].append(a)
          except Exception:
            bondedNeighborLists[n] = [a]

def compatibleConformations(a1, a2):
  '''
    Returns True if the two atoms are in compatible conformations, False if not.
    :param a1: First atom.
    :param a2: Second atom.
    :return: True if either atom is in the empty conformation or if both are in the
    same conformation.
  '''
  alt1 = a1.parent().altloc
  alt2 = a2.parent().altloc
  if alt1 in ['', ' ']:
    return True
  if alt2 in ['', ' ']:
    return True
  return alt1 == alt2

def getAtomsWithinNBonds(atom, bondedNeighborLists, extraAtomInfo, probeRad, N, nonHydrogenN = 1e10):
  """
    Helper function to produce a list of all of the atoms that are bonded to the
    specified atom, or to one of the atoms bonded to the specified atom, recursively
    to a depth of N.  The atom itself will not be included in the list, so an atom that
    has no bonded neighbors will have an empty result.  This can be used to
    produce a list of excluded atoms for dot scoring.  It checks to ensure that all of
    the bonded atoms are from compatible conformations (if the original atom
    is in the empty configuration then this will return atoms from all conformations that
    are in the bonded set).
    For Phantom Hydrogens only ever check to a depth of one (their parent Oxygen atom)
    to avoid spurious bonds found through ions.
    :param atom: The atom to be tested.
    :param bondedNeighborLists: Dictionary of lists that contain all bonded neighbors for
    each atom in a set of atoms.  Should be obtained using getBondedNeighborLists() and
    perhaps augmented by calling addIonicBonds().
    :param extraAtomInfo: mmtbx_probe_ext.ExtraAtomInfo mapper that provides radius and other
    information about atoms beyond what is in the pdb.hierarchy.  Used here to determine
    radii of the potential bonded neighbors.  No atom that is further than twice the probe
    radius from the source atom is part of the chain of neighbors.
    :param probeRad: Probe radius.  No atom that is further than twice the probe
    radius from the source atom is part of the chain of neighbors.
    :param N: Depth of recursion.  N=1 will return the atoms bonded to atom.  N=2 will
    also return those bonded to these neighbors (but not the atom itself).
    :param nonHydrogenN: When neither the original atom nor the bonded atom is a Hydrogen,
    limit the depth to this value (if this value is less than N).
    :returns a list of all atoms that are bonded to atom within a depth of N.  The original
    atom is never on the list.
  """
  aLoc = atom.xyz
  aRad = extraAtomInfo.getMappingFor(atom).vdwRadius
  atomIsHydrogen = atom.element_is_hydrogen()
  if extraAtomInfo.getMappingFor(atom).isDummyHydrogen:
    # Only ever allow a Phantom Hydrogen to be bonded to its parent Oxygen
    N = 1
  # Find all atoms to the specified depth
  atoms = {atom}            # Initialize the set with the atom itself
  for i in range(N):        # Repeat the recursion this many times
    current = list(atoms)   # Make a copy so we're not modifying the list we are traversing
    for a in current:       # Add all neighbors of all atoms in the current level
      for n in bondedNeighborLists[a]:
        # If we find a hydrogen, we no longer use the non-Hydrogen N limit.
        if i < nonHydrogenN or atomIsHydrogen or n.element_is_hydrogen():
          # Ensure that the new atom is in a compatible conformation with the original atom.
          if compatibleConformations(atom, n):
            # Ensure that the atom is not too far away from the original atom.
            nLoc = n.xyz
            nRad = extraAtomInfo.getMappingFor(n).vdwRadius
            d = (rvec3(nLoc) - rvec3(aLoc)).length()
            if d <= aRad + nRad + 2 * probeRad:
              atoms.add(n)

  # Remove the original atom from the result and turn the result into a list.
  atoms.discard(atom)
  return list(atoms)

def isPolarHydrogen(atom, bondedNeighborLists):
  '''
    Returns True if the atom is a polar hydrogen, False if not.
    :param atom: The atom to be tested.
    :param bondedNeighborLists: Dictionary of lists that contain all bonded neighbors for
    each atom in a set of atoms.  Should be obtained using
    mmtbx.probe.Helpers.getBondedNeighborLists() and perhaps augmented by calling addIonicBonds().
    :return: True if the atom is a polar hydrogen.
  '''
  if atom.element_is_hydrogen():
    neighbors = bondedNeighborLists[atom]
    if len(neighbors) == 1 and neighbors[0].element in ['N', 'O', 'S']:
      return True
  return False

def getMaxISeq(model):
  """Return the maximum i_seq value for all atoms in a model.
  """
  maxISeq = 0
  for a in model.get_atoms():
    maxISeq = max(maxISeq, a.i_seq)
  return maxISeq

class getExtraAtomInfoReturn(object):
  """
    Return type from getExtraAtomInfo() call.
      extraAtomInfo: ExtraAtomInfoMap with an entry for every atom in the model suitable for
                     passing to the scoring functions.
      warnings: a string that if not empty lists warnings that the person running the program
                might want to know about.  Suitable for printing or logging.
  """
  def __init__(self, extraAtomInfo, warnings):
    self.extraAtomInfo = extraAtomInfo
    self.warnings = warnings

def getExtraAtomInfo(model, bondedNeighborLists, useNeutronDistances = False, probePhil = None):
  """
    Helper function to provide a mapper for ExtraAtomInfo needed by Probe when scoring
    models.
    :param model: Map Model Manager's Model containing all of the atoms to be described.
    PDB interpretation must have been done on the model, perhaps by calling
    model.process(make_restraints=True), with useNeutronDistances matching
    the parameter to this function.
    :param bondedNeighborLists: Lists of atoms that are bonded to each other.
    Can be obtained by calling getBondedNeighborLists() and
    perhaps augmented by calling addIonicBonds().
    :param useNeutronDistances: Default is to use x-ray distances, but setting this to
    True uses neutron distances instead.  This must be set consistently with the
    PDB interpretation parameter used on the model.
    :param probePhil: None or subobject of PHIL parameters for probe.  Can be obtained using
    self.params.probe from a Program Template program that includes the probe_phil_parameters
    from above in its master PHIL parameters string.  If None, local defaults will be used.
    The following are used:
      set_polar_hydrogen_radius(bool): Default is to override the radius for polar
      Hydrogen atoms with 1.05.  Setting this to false uses the CCTBX value.
    :returns a ExtraAtomInfoMap with an entry for every atom in the model suitable for
    passing to the scoring functions.
  """

  warnings = ""

  # Traverse the hierarchy and look up the extra data to be filled in.
  extras = probeExt.ExtraAtomInfoMap([],[])
  mon_lib_srv = model.get_mon_lib_srv()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  ph = model.get_hierarchy()
  for m in ph.models():
    for chain in m.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():

          for a in ag.atoms():
            extra = probeExt.ExtraAtomInfo()
            try:
              # Get the first non-blank character of the alternate for this atom.
              # If there is none, set the value to the empty string.  If there is, set
              # it to a string with just that character.
              alt = a.parent().altloc.strip()
              extra.altLoc = alt

              extra.charge = probeExt.atom_charge(a)

              # For ions, the Richardsons determined in discussion with
              # Michael Prisant that we want to use the ionic radius rather than the
              # larger radius for all purposes.
              # @todo Once the CCTBX radius determination discussion and upgrade is
              # complete (ongoing as of March 2022), this check might be removed
              # and we'll just use the CCTBX radius.
              extra.isIon = a.element_is_ion()
              if extra.isIon:
                extra.vdwRadius = model.get_specific_ion_radius(a.i_seq)
                warnings += ("Using ionic radius for "+a.name.strip()+": "+str(extra.vdwRadius)+
                              " (rather than "+
                              str(model.get_specific_vdw_radius(a.i_seq, False))+
                              ")\n")
              else:
                extra.vdwRadius = model.get_specific_vdw_radius(a.i_seq, False)

              # Mark aromatic ring N and C atoms as acceptors as a hack to enable the
              # ring itself to behave as an acceptor.
              # @todo Remove this once we have a better way to model the ring itself
              # as an acceptor, perhaps making it a cylinder or a sphere in the center
              # of the ring.
              if a.element in ['C','N']:
                if AtomTypes.IsAromaticAcceptor(ag.resname, a.name):
                  extra.isAcceptor = True
                  warnings += "Marking "+a.name.strip()+" as an aromatic-ring acceptor\n"

              # Mark Nitrogens that do not have an attached Hydrogen as acceptors.
              # We only do this in HET atoms because CCTBX routines seem to be working
              # properly in standard residues.
              # We determine whether it is in a het atom by seeing if the residue name
              # is in the amino-acid mapping structure.
              if not a.parent().resname in iotbx.pdb.amino_acid_codes.one_letter_given_three_letter:
                if a.element == 'N':
                  # See if the atom has no Hydrogens covalently bonded to it.
                  found = False
                  for n in bondedNeighborLists[a]:
                    if n.element_is_hydrogen():
                      found = True
                      break
                  if not found and not extra.isAcceptor:
                    extra.isAcceptor = True
                    warnings += "Marking "+a.parent().resname.strip()+" "+a.name.strip()+" as a non-Hydrogen HET acceptor\n"

              # Mark all Carbonyl's with the Probe radius while the Richardsons and
              # the CCTBX decide how to handle this.
              # @todo After 2021, see if the CCTBX has the same values (1.65 and 1.80)
              # for Carbonyls and remove this if so.  It needs to stay with these values
              # to avoid spurious collisions per experiments run by the Richardsons in
              # September 2021.
              if (a.name.strip().upper() == 'C'
                  or AtomTypes.IsSpecialAminoAcidCarbonyl(a.parent().resname.strip().upper(),
                      a.name.upper()) ):
                expected = 1.65
                if extra.vdwRadius != expected:
                  warnings += ("Overriding radius for "+a.name.strip()+": "+str(expected)+
                                " (was "+str(extra.vdwRadius)+")\n")
                  extra.vdwRadius = expected

              # If we've been asked to ensure polar hydrogen radius, do so here.
              if probePhil.set_polar_hydrogen_radius and isPolarHydrogen(a, bondedNeighborLists):
                if extra.vdwRadius != 1.05:
                  warnings += ("Overriding radius for "+a.name.strip()+": 1.05 (was "+
                                str(extra.vdwRadius)+")\n")
                  extra.vdwRadius = 1.05

              # Find the hydrogen-bond type for this atom.
              # Do this after the above checks because it can fail for unknown HETs and we want to have
              # filled in the rest of the information first.
              hb_type = model.get_specific_h_bond_type(a.i_seq)
              if hb_type in ['A', 'B', 'D', 'N', 'H']: #isinstance(hb_type, str):
                if hb_type == "A" or hb_type == "B":
                  extra.isAcceptor = True
                if hb_type == "D" or hb_type == "B":
                  extra.isDonor = True

                extras.setMappingFor(a, extra)

              # Did not find hydrogen-bond information for this atom.
              else:
                fullName = (chain.id + ' ' + a.parent().resname.strip() + ' ' +
                  str(a.parent().parent().resseq_as_int()) + ' ' + a.name.strip())
                warnings += ("Warning: Unrecognized specific H bond type for "+fullName+", got "+hb_type+
                             ": keeping default value\n")
                extras.setMappingFor(a, extra)

            except Exception as e:
              fullName = (chain.id + ' ' + a.parent().resname.strip() + ' ' +
                str(a.parent().parent().resseq_as_int()) + ' ' + a.name.strip())
              warnings += ("Warning: Could not find atom info for "+fullName+
                " (perhaps interpretation was not run on the model?):"+
                " keeping some default values: "+str(e)+"\n")
              extras.setMappingFor(a, extra)

  return getExtraAtomInfoReturn(extras, warnings)

def writeAtomInfoToString(atoms, extraAtomInfo):
  """
    Write information about atoms, including their radius and donor/acceptor status,
    into a string that matches the form used by the new options added to the original
    Probe and Reduce code to enable comparisons between versions.
    :param atoms: The atoms to be described.
    :param extraAtomInfo: mmtbx_probe_ext.ExtraAtomInfo mapper that provides radius and other
    information about atoms beyond what is in the pdb.hierarchy.  Used here to determine
    which atoms may be acceptors.
    :return: String suitable for writing to a file to enable comparison with the dump
    files from the original Probe or Reduce programs.
  """

  #################################################################################
  # Dump information about all of the atoms in the model into a string.
  ret = ""
  acceptorChoices = ["noAcceptor","isAcceptor"]
  donorChoices = ["noDonor","isDonor"]
  metallicChoices = ["noMetallic","isMetallic"]
  for a in atoms:
    chainID = a.parent().parent().parent().id
    resName = a.parent().resname.upper()
    resID = str(a.parent().parent().resseq_as_int())
    alt = a.parent().altloc
    if alt == " " or alt == "":
      alt = "-"
    ret += (
      " "+str(chainID)+" "+resName+" {:3d} ".format(int(resID))+a.name+" "+alt+
      " {:7.3f}".format(a.xyz[0])+" {:7.3f}".format(a.xyz[1])+" {:7.3f}".format(a.xyz[2])+
      " {:5.2f}".format(extraAtomInfo.getMappingFor(a).vdwRadius)+
      " "+acceptorChoices[extraAtomInfo.getMappingFor(a).isAcceptor]+
      " "+donorChoices[extraAtomInfo.getMappingFor(a).isDonor]+
      " "+metallicChoices[a.element_is_positive_ion()]+
      "\n"
    )
  return ret

def compareAtomInfoFiles(fileName1, fileName2, distanceThreshold):
  """
    Write information about atom differences in two files into a string that describes
    their differences.  This is used to perform a difference on files that contain the
    output of writeAtomInfoToString() run on different versions of Probe or Reduce.
    The lines in the files are sorted so that reordering of the atoms between files will
    not cause this function to report differences.
    :param fileName1: The first file.
    :param fileName2: The second file.
    :param distanceThreshold: If the same atom in each file is further apart than this
    between files, report it as a difference; otherwise, not.
    :return: String suitable for writing to a file to describe the differences.
  """

  ###########################
  # Helper utility functions usable only inside this function
  def atomID(s):
    # Return the ID of the atom, which includes its chain, residue name,
    # residue number, atom name, and alternate separated by spaces
    w = s.split()
    return w[0]+" "+w[1]+" "+w[2]+" "+w[3]+" "+w[4]

  def position(s):
    # Return a tuple that has the location of the atom.
    w = s.split()
    return ( float(w[5]), float(w[6]), float(w[7]) )

  def radius(s):
    # Return the radius for the atom.
    w = s.split()
    return ( w[8] )

  def flag(s,f):
    # Return a flag for the atom, indexed by 0-2.
    w = s.split()
    return ( w[9+f] )

  def distance(a, b):
    # Return the distance between 3-space tuples a and b
    return math.sqrt( (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 )

  ###########################
  # Initialize our return value to empty.
  ret = ""

  # Read the the data from each file.
  # Sort the lines from the file so that we'll get chain, residue, atom name sorting.
  with open(fileName1) as f:
    m1 = f.readlines()
  m1.sort(key=lambda x:atomID(x))

  with open(fileName2) as f:
    m2 = f.readlines()
  m2.sort(key=lambda x:atomID(x))

  # Compare the two and report on any differences.
  # Both lists are sorted, so we'll know which has a missing element compared
  # to the other.  We walk from the beginning to the end of each list.
  m1ID = 0
  m2ID = 0
  while m1ID < len(m1) and m2ID < len(m2):

    # Atom names are unique within residues, so we can use that to
    # match atoms in one file from those in the other and print out
    # missing atoms and differences in positions.
    if atomID(m1[m1ID]) < atomID(m2[m2ID]):
      ret += 'Only in first: '+atomID(m1[m1ID])+'\n'
      m1ID += 1
      continue

    if atomID(m1[m1ID]) > atomID(m2[m2ID]):
      ret += 'Only in second: '+atomID(m2[m2ID])+'\n'
      m2ID += 1
      continue

    # Check the radius to see if they differ
    oldRadius = radius(m1[m1ID])
    newRadius = radius(m2[m2ID])
    if oldRadius != newRadius:
      ret += atomID(m1[m1ID])+' radii differ: old is {}, new is {}'.format(oldRadius,newRadius)+'\n'

    # Check the flag fields to see if any of them differ
    for f in range(3):
      oldFlag = flag(m1[m1ID], f)
      newFlag = flag(m2[m2ID], f)
      if oldFlag != newFlag:
        ret += atomID(m1[m1ID])+' flags differ: old is {}, new is {}'.format(oldFlag,newFlag)+'\n'

    # Check the distance between the atoms.
    dist = distance(position(m1[m1ID]),position(m2[m2ID]))
    if dist > distanceThreshold:
      ret += atomID(m1[m1ID])+' Distance between runs: {:.2f}'.format(dist)+'\n'

    # Go on to next line
    m1ID += 1
    m2ID += 1

  # Finish up either list that is not done
  while m1ID < len(m1):
    ret += 'Only in first: '+atomID(m1[m1ID])+'\n'
    m1ID += 1
  while m2ID < len(m2):
    ret += 'Only in second: '+atomID(m2[m2ID])+'\n'
    m2ID += 1

  return ret

def getPhantomHydrogensFor(largestISeq, atom, spatialQuery, extraAtomInfo, minOccupancy,
      acceptorOnly = False, placedHydrogenRadius = 1.05, placedHydrogenDistance = 0.84):
  """
    Get a list of phantom Hydrogens for the atom specified, which is asserted to be an Oxygen
    atom for a water.
    :param largestISeq: The largest i_seq number in the model so far. Each new Phantom will be
    given a different number, starting with one larger than this. This should be incremented
    by the length of the returned list if this function is called multiple times.
    :param atom: The Oxygen that is to have phantoms added to it.
    :param spatialQuery: mmtbx_probe_ext.SpatialQuery structure to rapidly determine which atoms
    are within a specified distance of a location.
    :param extraAtomInfo: mmtbx_probe_ext.ExtraAtomInfo mapper that provides radius and other
    information about atoms beyond what is in the pdb.hierarchy.  Used here to determine
    which atoms may be acceptors.
    :param minOccupancy: Minimum occupancy for a nearby atom to be considered.
    :param acceptorOnly: Only allow bonds with atoms that are acceptors when this is True.
    This is false by default because Reduce needs to check whether the bonded atom is either
    an acceptor or a possible flipped position of an acceptor, and that is not something that
    can be determined at the time we're placing phantom hydrogens.  In that case, we want to
    include all possible interactions and weed them out during optimization.
    :param placedHydrogenRadius: Radius of the Phantom Hydrogen to be placed.  Default is
    for electron-cloud distances.
    :param placedHydrogenDistance: Maximum distance from the placed Phantom Hydrogen to the
    Water Oxygen. The Phantom Hydrogens are placed at the optimal overlap distance so may be
    closer than this.  Default is for electron-cloud distances.
    :return: List of new atoms that make up the phantom Hydrogens, with only their name and
    element type and xyz positions filled in.  They will have i_seq as specified and they
    should not be inserted into a structure.
  """

  ret = []

  # Get the list of nearby atoms.  The center of the search is the water atom
  # and the search radius is 4 (these values are pulled from the Reduce C++ code).
  # Sort these by i_seq so that we get repeatable results from run to run.
  maxDist = 4.0
  nearby = sorted(spatialQuery.neighbors(atom.xyz, 0.001, maxDist), key=lambda x:x.i_seq)

  # Candidates for nearby atoms.  We use this list to keep track of one ones we
  # have already found so that we can compare against them to only get one for each
  # aromatic ring.
  class Candidate(object):
    def __init__(self, atom, overlap):
      self._atom = atom
      self._overlap = overlap
  candidates = []

  atomModel = atom.parent().parent().parent().parent().id
  for a in nearby:
    aModel = a.parent().parent().parent().parent().id

    # Only check atoms in compatible conformations.
    if not compatibleConformations(atom, a):
      continue

    # Skip atoms that are in different models.
    if atomModel != aModel:
      continue

    # Check to ensure the occupancy of the neighbor is above threshold and that it is
    # close enough to potentially bond to the atom.  We want a minimum overlap of 0.1
    # before we consider possible Hs.
    overlap = ( (rvec3(atom.xyz) - rvec3(a.xyz)).length()  -
                (placedHydrogenRadius + extraAtomInfo.getMappingFor(a).vdwRadius + placedHydrogenDistance) )
    if overlap <= -0.1 and a.occ >= minOccupancy and a.element != "H":
      if not acceptorOnly or extraAtomInfo.getMappingFor(a).isAcceptor:
        # If we have multiple atoms in the same Aromatic ring (part of the same residue)
        # we only point at the closest one.  To ensure this, we check all current candidates
        # and if we find one that is on the same aromatic ring then we either ignore this new
        # atom (if it is further) or replace the existing one (if it is closer).
        skip = False
        if AtomTypes.IsAromaticAcceptor(a.parent().resname.strip().upper(), a.name.strip().upper()):
          for c in candidates:
            # See if we belong to the same atom group and are both ring acceptors.  If so, we need to replace
            # or else squash this atom.
            if (AtomTypes.IsAromaticAcceptor(c._atom.parent().resname.strip().upper(), c._atom.name.strip().upper()) and
                a.parent() == c._atom.parent()):
              if overlap < c._overlap:
                # Replace the further atom with this atom.
                c._atom = a
                c._overlap = overlap
                skip = True
                break
              else:
                # This is further away, so we don't insert it.
                skip = True
                break

        # Add the Candidate
        if not skip:
          candidates.append(Candidate(a, overlap))

  # Generate phantoms pointing toward all of the remaining candidates.
  # Make most of their characteristics copied from the source Oxygen.
  # The element, name, and location are modified.
  for c in candidates:
    h = pdb.hierarchy.atom(atom.parent(), atom)
    h.element = "H"
    h.name = " H?"

    # Place the Phantom Hydrogen pointing from the Oxygen towards the candidate at a distance
    # of the standard H bond length plus an offset that is clamped to the range -H bond length..0 that
    # is the sum of the overlap and the best hydrogen-bonding overlap.  This is an approximation
    # to the situation where the Phantom Hydrogen would rotate around the Oxygen to maintain a proper
    # distance from the acceptor that does not involve trying to select a rotation direction.
    # Because Phantom Hydrogens do not block the Oxygen from collisions with their neighbors, and
    # because Phantom Hydrogens cannot clash with any atom, this will not interfere with clashes.
    # Note that the Phantom Hydrogen also will not block any collision between the Oxygen atom
    # in the Water and a nearby acceptor, so those collisions will still show up.
    # The bond length will never be lengthened, but might be shortened if there is more than
    # the ideal amount of overlap, and it will never be shorted to less than 0 (the location
    # of the Oxygen).
    BEST_HBOND_OVERLAP=0.6
    distance = placedHydrogenDistance + max(
      -placedHydrogenDistance, min(0.0, c._overlap + BEST_HBOND_OVERLAP) )
    try:
      normOffset = (rvec3(c._atom.xyz) - rvec3(atom.xyz)).normalize()
      h.xyz = rvec3(atom.xyz) + distance * normOffset
      largestISeq += 1
      probeExt.set_atom_i_seq(h, largestISeq)
      ret.append(h)
    except Exception:
      # If we have overlapping atoms (normalize() fails), don't add.
      pass

  return ret

def fixupExplicitDonors(atoms, bondedNeighborLists, extraAtomInfo):
  """
    Fix up the donor status for models that have explicit hydrogens.  All Nitrogens, Oxygens
    and Sulphur atoms are stripped of their donor status because they will have explicit Hydrogens
    added or else had some other form of covalent bond added.  All hydrogens that are bonded to
    Nitrogens, Oxygens, or Sulphur atoms are marked as donors.  This does not handle any Phantom
    Hydrogens unless those Hydrogens are marked as bonded to their Water Oxygen.
    :param atoms: The list of atoms to adjust.
    :param bondedNeighborLists: Dictionary of lists that contain all bonded neighbors for
    each atom in a set of atoms.  Should be obtained using
    mmtbx.probe.Helpers.getBondedNeighborLists() and perhaps augmented by calling addIonicBonds().
    :param extraAtomInfo: mmtbx_probe_ext.ExtraAtomInfo mapper that provides radius and other
    information about atoms beyond what is in the pdb.hierarchy.  Used here to determine
    which atoms may be acceptors.  This information is modified in place to adjust the donor
    status of atoms.
    :return: None.  As a side effect, the extraAtomInfo is adjusted.
  """

  for a in atoms:
    # If we are a hydrogen that is bonded to a nitrogen, oxygen, or sulfur then we're a donor
    # and our bonded neighbor is not.
    if a.element_is_hydrogen():
      for n in bondedNeighborLists[a]:
        if n.element in ['N','O','S']:
          # Copy the value, set the new values, then copy the new one back in.
          # We are a donor and may have our radius adjusted
          ei = extraAtomInfo.getMappingFor(a)
          ei.isDonor = True
          extraAtomInfo.setMappingFor(a, ei)

          # Set our neigbor to not be a donor, since we are the donor
          ei = extraAtomInfo.getMappingFor(n)
          ei.isDonor = False
          extraAtomInfo.setMappingFor(n, ei)

    # Otherwise, if we're an N, O, or S then remove our donor status because
    # the hydrogens will be the donors.  Because we're not doing implicit
    # hydrogens (and thus are doing explicit hydrogens), if we have a leftover
    # atom that did not have a hydrogen attached we assume that this is because
    # there is some other bonding and we still need to remove the donor status.
    elif a.element in ['N','O','S']:
      ei = extraAtomInfo.getMappingFor(a)
      ei.isDonor = False
      extraAtomInfo.setMappingFor(a, ei)

def dihedralChoicesForRotatableHydrogens(hydrogens, hParams, potentials):
  """
    Determine which of the hydrogens passed in is the one that should be used to
    define the conventional dihedral for a rotatable hydrogen bond, and which of
    the potential third-link neighbors should be used as the other atom for the
    calculation.
    :param hydrogens: List of atoms in a set of rotatable hydrogens.
    :param hParams: List indexed by sequence ID that stores the riding
    coefficients for hydrogens that have associated dihedral angles.  This can be
    obtained by calling model.setup_riding_h_manager() and then model.get_riding_h_manager().
    :param potentials: List of atoms; potential third-bonded neighbor atoms to the hydrogens.
    It can be obtained by walking the covalent bonds in the structure from any of
    the hydrogens.
    :returns: A tuple whose first entry is the hydrogen to use, whose second entry
    is the atom to use to compute the dihedral.
  """
  conventionalH = None
  conventionalFriend = None
  for h in hydrogens:
    # Find the Hydrogen whose "n" entry is 0 and fill in it and the associated
    # potential atom.
    try:
      item = hParams[h.i_seq]
      if item.n == 0:
        whichPotential = item.a2
        for p in potentials:
          if p.i_seq == whichPotential:
            conventionalH = h
            conventionalFriend = p
    except Exception as e:
      pass

  # Throw a Sorry if we were unable to find an answer.
  if conventionalH is None or conventionalFriend is None:
    raise Sorry("mmtbx.probe.Helpers.dihedralChoicesForRotatableHydrogens(): Could not determine atoms to use")

  return conventionalH, conventionalFriend

##################################################################################
# Helper functions to make things that are compatible with vec3_double so
# that we can do math on them.  We need a left-hand and right-hand one so that
# we can make both versions for multiplication.
def rvec3(xyz) :
  return scitbx.matrix.rec(xyz, (3,1))
def lvec3(xyz) :
  return scitbx.matrix.rec(xyz, (1,3))

def Test(inFileName = None):
  """
  Run tests on all of our functions.  Throw an assertion failure if one fails.
  """

  #========================================================================
  # Model file snippets used in the tests.
  pdb_1xso_his_61 = (
"""
ATOM    442  N   HIS A  61      26.965  32.911   7.593  1.00  7.19           N
ATOM    443  CA  HIS A  61      27.557  32.385   6.403  1.00  7.24           C
ATOM    444  C   HIS A  61      28.929  31.763   6.641  1.00  7.38           C
ATOM    445  O   HIS A  61      29.744  32.217   7.397  1.00  9.97           O
ATOM    446  CB  HIS A  61      27.707  33.547   5.385  1.00  9.38           C
ATOM    447  CG  HIS A  61      26.382  33.956   4.808  1.00  8.78           C
ATOM    448  ND1AHIS A  61      26.168  34.981   3.980  1.00  9.06           N
ATOM    449  CD2 HIS A  61      25.174  33.397   5.004  1.00 11.08           C
ATOM    450  CE1 HIS A  61      24.867  35.060   3.688  1.00 12.84           C
ATOM    451  NE2 HIS A  61      24.251  34.003   4.297  1.00 11.66           N
HETATM 4333 CU    CU A   1      22.291  33.388   3.996  1.00 13.22          Cu
END
"""
    )

  pdb_4fenH_C_26 = (
"""
ATOM    234  P     C B  26      25.712  13.817   1.373  1.00 10.97           P
ATOM    235  OP1   C B  26      26.979  13.127   1.023  1.00 12.57           O
ATOM    236  OP2   C B  26      25.641  14.651   2.599  1.00 12.59           O
ATOM    237  O5'   C B  26      25.264  14.720   0.142  1.00 10.58           O
ATOM    238  C5'   C B  26      25.128  14.166  -1.162  1.00  9.37           C
ATOM    239  C4'   C B  26      24.514  15.183  -2.092  1.00  9.05           C
ATOM    240  O4'   C B  26      23.142  15.447  -1.681  1.00 10.46           O
ATOM    241  C3'   C B  26      25.159  16.555  -2.048  1.00  8.49           C
ATOM    242  O3'   C B  26      26.338  16.599  -2.839  1.00  7.94           O
ATOM    243  C2'   C B  26      24.043  17.427  -2.601  1.00  8.44           C
ATOM    244  O2'   C B  26      23.885  17.287  -3.999  1.00  9.33           O
ATOM    245  C1'   C B  26      22.835  16.819  -1.888  1.00 10.06           C
ATOM    246  N1    C B  26      22.577  17.448  -0.580  1.00  9.44           N
ATOM    247  C2    C B  26      21.940  18.695  -0.552  1.00  9.61           C
ATOM    248  O2    C B  26      21.613  19.223  -1.623  1.00 10.13           O
ATOM    249  N3    C B  26      21.704  19.292   0.638  1.00  9.82           N
ATOM    250  C4    C B  26      22.082  18.695   1.771  1.00 10.10           C
ATOM    251  N4    C B  26      21.841  19.328   2.923  1.00 10.94           N
ATOM    252  C5    C B  26      22.728  17.423   1.773  1.00  9.51           C
ATOM    253  C6    C B  26      22.952  16.840   0.586  1.00 10.02           C
ATOM      0  H5'   C B  26      24.573  13.371  -1.128  1.00  9.37           H   new
ATOM      0 H5''   C B  26      25.996  13.892  -1.498  1.00  9.37           H   new
ATOM      0  H4'   C B  26      24.618  14.792  -2.973  1.00  9.05           H   new
ATOM      0  H3'   C B  26      25.463  16.833  -1.170  1.00  8.49           H   new
ATOM      0  H2'   C B  26      24.192  18.375  -2.460  1.00  8.44           H   new
ATOM      0 HO2'   C B  26      23.422  16.605  -4.162  1.00  9.33           H   new
ATOM      0  H1'   C B  26      22.041  16.954  -2.429  1.00 10.06           H   new
ATOM      0  H41   C B  26      22.073  18.968   3.668  1.00 10.94           H   new
ATOM      0  H42   C B  26      21.454  20.096   2.919  1.00 10.94           H   new
ATOM      0  H5    C B  26      22.984  17.013   2.568  1.00  9.51           H   new
ATOM      0  H6    C B  26      23.369  16.009   0.556  1.00 10.02           H   new
"""
    )

  class philLike:
    def __init__(self):
      self.set_polar_hydrogen_radius = True

  #========================================================================
  # Run unit test on getExtraAtomInfo().  We use a specific PDB snippet
  # for which we know the answer and then we verify that the results are what
  # we expect.

  # Spot check the values on the atoms for standard, neutron distances,
  # and original Probe results.
  standardChecks = [
    # Name, vdwRadius, isAcceptor, isDonor, isDummyHydrogen, isIon, charge, altLoc
    ["CU",  0.72, False, False, False,  True, 0, ''],
    ["N",   1.55, False,  True, False, False, 0, ''],
    ["ND1", 1.55, True,   True, False, False, 0, 'A'],
    ["C",   1.65, False, False, False, False, 0, ''],
    ["CB",  1.7,  False, False, False, False, 0, ''],
    ["O",   1.4,  True,  False, False, False, 0, ''],
    ["CD2", 1.75, False, False, False, False, 0, '']
  ]
  neutronChecks = [
    # Name, vdwRadius, isAcceptor, isDummyHydrogen, isDonor, isIon, charge, altLoc
    ["CU",  0.72, False, False, False,  True, 0, ''],
    ["N",   1.55, False,  True, False, False, 0, ''],
    ["ND1", 1.55, True,   True, False, False, 0, 'A'],
    ["C",   1.65, False, False, False, False, 0, ''],
    ["CB",  1.7,  False, False, False, False, 0, ''],
    ["O",   1.4,  True,  False, False, False, 0, ''],
    ["CD2", 1.75, False, False, False, False, 0, '']
  ]

  # Situations to run the test in and expected results:
  cases = [
    # Use neutron distances, expected results
    [False, standardChecks],
    [True,  neutronChecks]
  ]

  for cs in cases:
    useNeutronDistances = cs[0]
    checks = cs[1]
    runType = "; neutron = "+str(useNeutronDistances)

    dm = iotbx.data_manager.DataManager(['model'])
    dm.process_model_str("1xso_snip.pdb",pdb_1xso_his_61)
    model = dm.get_model()
    p = mmtbx.model.manager.get_default_pdb_interpretation_params()
    p.pdb_interpretation.use_neutron_distances = useNeutronDistances
    model.process(make_restraints=True, pdb_interpretation_params = p)

    # Get the atoms for the first model in the hierarchy.
    atoms = model.get_hierarchy().models()[0].atoms()

    # Get the Cartesian positions of all of the atoms we're considering for this alternate
    # conformation.
    carts = flex.vec3_double()
    for a in atoms:
      carts.append(a.xyz)

    # Get the bond proxies for the atoms in the model and conformation we're using and
    # use them to determine the bonded neighbor lists.
    bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
    bondedNeighborLists = getBondedNeighborLists(atoms, bondProxies)

    # Get the extra atom information for the model using default parameters.
    # Make a PHIL-like structure to hold the parameters.
    extras = getExtraAtomInfo(model,bondedNeighborLists,
      useNeutronDistances=useNeutronDistances,probePhil=philLike()).extraAtomInfo

    # Get the atoms for the first model in the hierarchy.
    atoms = model.get_hierarchy().models()[0].atoms()

    for a in atoms:
      e = extras.getMappingFor(a)
      for c in checks:
        if a.name.strip() == c[0]:
          if hasattr(math, 'isclose'):
            assert math.isclose(e.vdwRadius, c[1]), ("Helpers.Test(): Bad radius for "+a.name+runType+": "
            +str(e.vdwRadius)+" (expected "+str(c[1])+")")
          assert e.isAcceptor == c[2], "Helpers.Test(): Bad Acceptor status for "+a.name+": "+str(e.isAcceptor)+runType
          assert e.isDonor == c[3], "Helpers.Test(): Bad Donor status for "+a.name+": "+str(e.isDonor)+runType
          # Check the ability to set and check Dummy/Phantom Hydrogen status
          assert e.isDummyHydrogen == False, "Helpers.Test(): Bad Dummy Hydrogen status for "+a.name+runType
          e.isDummyHydrogen = True
          assert e.isDummyHydrogen == True, "Helpers.Test(): Can't set DummyHydrogen status for "+a.name+runType
          assert e.isIon == c[5], "Helpers.Test(): Bad Ion status for "+a.name+": "+str(e.isIon)+runType
          assert e.charge == c[6], "Helpers.Test(): Bad charge for "+a.name+": "+str(e.charge)+runType
          assert e.altLoc == c[7], "Helpers.Test(): Bad altLoc for "+a.name+": "+str(e.altLoc)+", wanted "+str(c[7])+runType

  #========================================================================
  # Run unit test on getPhantomHydrogensFor().

  # Generate a Water Oxygen atom at the origin with a set of atoms around it.
  # The surrounding ones will include atoms that are acceptors, atoms that are not
  # and atoms that are on an aromatic ring (all on the same ring).  The distance to
  # each atom will be such that it is within range for a Phantom Hydrogen.  The
  # directions will be towards integer lattice points and we'll round-robin the type.
  # NOTE: The atom names and types and radii and such are not correct, just made up
  # to test our functions.
  try:
    # Prepare our extraAtomInfoMap
    atoms = []
    extras = []

    radius = 1.1  # All atoms have the same radius
    rg = pdb.hierarchy.residue_group()
    water = pdb.hierarchy.atom_group()
    water.resname = 'HOH'
    rg.append_atom_group(water)
    o = pdb.hierarchy.atom()
    o.element = "O"
    o.xyz = [ 0.0, 0.0, 0.0 ]
    atoms.append(o)
    extras.append(probeExt.ExtraAtomInfo(radius, False))
    water.append_atom(o)
    type = 0
    # They are all in a residue that has the specified atom name listed as Aromatic.
    # We change the names and types below to make some not match.
    ag = pdb.hierarchy.atom_group()
    ag.resname = 'ADE'
    for x in range(-1,2,2):
      for y in range(-1,2,2):
        for z in range(-1,2,2):
          dist = math.sqrt(x*x+y*y+z*z)
          d = (radius + 1.0) / dist
          a = pdb.hierarchy.atom()
          a.xyz = [ x*d,y*d,z*d ]
          a.occ = 0.5
          if type % 3 == 0: # Acceptor
            a.element = 'C'
            a.name = 'CA'
            extras.append(probeExt.ExtraAtomInfo(radius, True))
          elif type % 3 == 1: # Not Acceptor
            a.element = 'C'
            a.name = 'CB'
            extras.append(probeExt.ExtraAtomInfo(radius, False))
          else: # Aromatic ring atom, also an acceptor
            a.name = 'N1'
            a.element = 'N'
            extras.append(probeExt.ExtraAtomInfo(radius, True))
          type += 1

          atoms.append(a)
          ag.append_atom(a)

    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    # Prepare our extra-atom information mapper.
    extrasMap = probeExt.ExtraAtomInfoMap(atoms, extras)

    # Prepare our spatial-query structure.
    sq = probeExt.SpatialQuery(atoms)

    # Run Phantom placement with different settings for occupancy and acceptorOnly
    # and make sure the atom counts match what is expected.
    # @todo Test changing the radius and distance for the Hydrogen.
    largestISeq = len(m.atoms())
    for occThresh in [0.4,1.0]:
      for acceptorOnly in [False, True]:
        # Check that we get the expected number of contacts
        if acceptorOnly:
          expected = 1 + 3 # One for the aromatics, three others
        else:
          expected = 7 # Only one of the acceptor Aromatics but all other atoms
        if occThresh > 0.5:
          expected = 0
        ret = getPhantomHydrogensFor(largestISeq, o, sq, extrasMap, occThresh, acceptorOnly, 1.0, 1.0)
        assert len(ret) == expected, "Helpers.Test() Unexpected count during Phantom Hydrogen placement: "+str(len(ret))+" (expected "+str(expected)+")"
        largestISeq += len(ret)

        # The location of the each Hydrogen should point towards one of the non-Oxygen atoms.
        # Here we check that we get as many matching directions as we have atoms.
        if len(ret) > 0:
          numMatch = 0
          for h in ret:
            for a in atoms[1:]:
              hLen = rvec3(h.xyz).length()
              aLen = rvec3(a.xyz).length()
              dot = ( h.xyz[0] * a.xyz[0] +
                      h.xyz[1] * a.xyz[1] +
                      h.xyz[2] * a.xyz[2] )
              if hasattr(math, 'isclose'):
                if math.isclose(hLen*aLen, dot):
                  numMatch += 1
          assert numMatch == len(ret), "Helpers.Test(): Direction of Phantom Hydrogen placement incorrect: "+str(numMatch)

  except Exception as e:
    assert len(str(e)) == 0, "Helpers.Test() Exception during test of Phantom Hydrogen placement: "+str(e)+"\n"+traceback.format_exc()

  #========================================================================
  # Run unit test on getBondedNeighborLists().  We use a specific PDB snippet
  # for which we know the answer and then we verify that the results are what
  # we expect.

  dm = iotbx.data_manager.DataManager(['model'])
  dm.process_model_str("1xso_snip.pdb",pdb_1xso_his_61)
  model = dm.get_model()
  model.process(make_restraints=True) # make restraints

  # Get the atoms for the first model in the hierarchy.
  atoms = model.get_hierarchy().models()[0].atoms()

  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms and use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = getBondedNeighborLists(atoms, bondProxies)

  # Check the counts in the neighbor lists to make sure they match what we expect
  neighborCounts = {"N": 1, "CA": 3, "C": 2, "O": 1, "CB": 2,
                    "CG": 3, "ND1": 2, "CD2": 2, "CE1":2, "NE2": 2, "CU": 0}
  for a in atoms:
    assert len(bondedNeighborLists[a]) == neighborCounts[a.name.strip()], (
        "Helpers.Test(): Neighbor count for "+a.name.strip()+" was "+
        str(len(bondedNeighborLists[a]))+", expected "+str(neighborCounts[a.name.strip()]))

  #=====================================================================================
  # Run unit test on compatibleConformations().
  a1 = pdb.hierarchy.atom()
  ag1 = pdb.hierarchy.atom_group()
  ag1.append_atom(a1)
  a2 = pdb.hierarchy.atom()
  ag2 = pdb.hierarchy.atom_group()
  ag2.append_atom(a2)
  ag1.altloc = ""
  ag2.altloc = "A"
  assert compatibleConformations(a1,a2), "Helpers:Test(): altloc expected True for empty first"
  ag1.altloc = "A"
  ag2.altloc = "A"
  assert compatibleConformations(a1,a2), "Helpers:Test(): altloc expected True for compatible"
  ag1.altloc = "A"
  ag2.altloc = "B"
  assert not compatibleConformations(a1,a2), "Helpers:Test(): altloc expected False for incompatible"
  ag1.altloc = "A"
  ag2.altloc = " "
  assert compatibleConformations(a1,a2), "Helpers:Test(): altloc expected True for blank second"
  ag1.altloc = ""
  ag2.altloc = " "
  assert compatibleConformations(a1,a2),  "Helpers:Test(): altloc expected True for empty first and blank second"

  #========================================================================
  # Run unit test on isPolarHydrogen().  We use a specific PDB snippet
  # for which we know the answer and then we verify that the results are what
  # we expect.

  dm = iotbx.data_manager.DataManager(['model'])
  dm.process_model_str("pdb_4fenH_C_26.pdb",pdb_4fenH_C_26)
  model = dm.get_model()
  model.process(make_restraints=True)

  # Get the Cartesian positions of all of the atoms.
  carts = flex.vec3_double()
  atoms = model.get_hierarchy().models()[0].atoms()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms and
  # use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = getBondedNeighborLists(atoms, bondProxies)

  model = dm.get_model()

  for a in model.get_hierarchy().models()[0].atoms():
    if a.name.strip() in ["H41","H42","HO2'"]:
      if not isPolarHydrogen(a, bondedNeighborLists):
        return "Optimizers.Test(): Polar Hydrogen not identified: " + a.name
    if a.name.strip() in ["H5","H6"]:
      if isPolarHydrogen(a, bondedNeighborLists):
        return "Optimizers.Test(): Polar Hydrogen improperly identified: " + a.name

  #========================================================================
  # Run unit test on getAtomsWithinNBonds().
  # Get the atoms within N bounds for a range for the "N" atom and verify that the
  # counts match what is expected.  Do this for the case where we clamp the non-
  # hydrogen ones to the 3 and when we use the default of very large to count
  # them all.
  # NOTE: This re-uses the bondedNeighborLists test results from above.
  N4 = None
  for a in atoms:
    if a.name.strip().upper() == 'N4':
      N4 = a
  assert N4 is not None, ("Helpers.Test(): Could not find N4 (internal failure)")

  # First test with a huge probe radius where all of the atoms are within range.
  # When clamped, we can't go further than the non-Hydrogen bound except for Hydrogens
  nestedNeighborsForN4 = [ None, 3, 5, 8, 9, 9, 9]
  extraInfo = getExtraAtomInfo(model, bondedNeighborLists,
    useNeutronDistances=False,probePhil=philLike()).extraAtomInfo
  hugeRadius = 1000
  for N in range(1,7):
    count = len(getAtomsWithinNBonds(N4, bondedNeighborLists, extraInfo, hugeRadius, N, 3))
    assert count == nestedNeighborsForN4[N], ("Helpers.Test(): Nested clamped count for "+N4.name.strip()+
        " for N = "+str(N)+" was "+str(count)+", expected "+str(nestedNeighborsForN4[N]))
  # When unclamped, we can traverse all the way to the end in all cases
  nestedNeighborsForN4 = [ None, 3, 5, 8, 11, 12, 15]
  for N in range(1,7):
    count = len(getAtomsWithinNBonds(N4, bondedNeighborLists, extraInfo, hugeRadius, N))
    assert count == nestedNeighborsForN4[N], ("Helpers.Test(): Nested unclamped count for "+N4.name.strip()+
        " for N = "+str(N)+" was "+str(count)+", expected "+str(nestedNeighborsForN4[N]))
  # When we start with a Hydrogen, we should never get clamped off.
  H41 = None
  for a in atoms:
    if a.name.strip().upper() == 'H41':
      H41 = a
  assert H41 is not None, ("Helpers.Test(): Could not find H41 (internal failure)")
  nestedNeighborsForH41 = [ None, 1, 3, 5, 8, 11, 12]
  for N in range(1,7):
    count = len(getAtomsWithinNBonds(H41, bondedNeighborLists, extraInfo, hugeRadius, N, 3))
    assert count == nestedNeighborsForH41[N], ("Helpers.Test(): Nested clamped count for "+H41.name.strip()+
        " for N = "+str(N)+" was "+str(count)+", expected "+str(nestedNeighborsForH41[N]))

  # Then test with a small probe radius which limits how far we can go before not finding
  # any more neighbors.
  nestedNeighborsForN4Small = [ None, 3, 5, 6, 6, 6, 6]
  smallRadius = 0.1
  for N in range(1,7):
    count = len(getAtomsWithinNBonds(N4, bondedNeighborLists, extraInfo, smallRadius, N))
    assert count == nestedNeighborsForN4Small[N], ("Helpers.Test(): Nested small-radius count for "+N4.name.strip()+
        " for N = "+str(N)+" was "+str(count)+", expected "+str(nestedNeighborsForN4Small[N]))

  # Test with Phantom Hydrogens to ensure we don't see more than just the nearest atom.
  # @todo

  #========================================================================
  # Test addIonicBonds().
  # @todo

  #========================================================================
  # Generate an example data model with a small molecule in it or else read
  # from the specified file.
  if inFileName is not None and len(inFileName) > 0:
    # Read a model from a file using the DataManager
    dm = iotbx.data_manager.DataManager()
    dm.process_model_file(inFileName)
    model = dm.get_model(inFileName)
  else:
    # Generate a small-molecule model using the map model manager
    # get an initialized instance of the map_model_manager
    mmm=iotbx.map_model_manager.map_model_manager()
    mmm.generate_map()     #   get a model from a generated small library model and calculate a map for it
    model = mmm.model()    #   get the model

  # Fix up bogus unit cell when it occurs by checking crystal symmetry.
  cs = model.crystal_symmetry()
  if (cs is None) or (cs.unit_cell() is None):
    model = cctbx.maptbx.box.shift_and_box_model(model = model)

  # Run PDB interpretation on the model to fill in the required CCTBX information.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = False
  model.process(make_restraints=True, pdb_interpretation_params = p) # make restraints

  ret = getExtraAtomInfo(model, bondedNeighborLists, False, philLike())

  #========================================================================
  # Run unit tests on the dihedralChoicesForRotatableHydrogens class.  Both
  # with and without being able to find the results.
  pdb_dihedral ='''
CRYST1   15.957   14.695   20.777  90.00  90.00  90.00 P 1
ATOM      1  N   VAL A   4      10.552   7.244  12.226  1.00 11.32           N
ATOM      2  CA  VAL A   4       9.264   6.605  12.502  1.00 11.54           C
ATOM      3  C   VAL A   4       8.169   7.396  11.772  1.00 11.25           C
ATOM      4  O   VAL A   4       7.912   8.570  12.087  1.00 12.48           O
ATOM      5  CB  VAL A   4       8.972   6.612  14.048  1.00 12.54           C
ATOM      6  CG1 VAL A   4       7.640   5.915  14.308  1.00 11.37           C
ATOM      7  CG2 VAL A   4      10.139   5.959  14.827  1.00 12.59           C
ATOM      8  H   VAL A   4      10.592   8.064  12.483  1.00 11.32           H
ATOM      9  HA  VAL A   4       9.282   5.685  12.195  1.00 11.54           H
ATOM     10  HB  VAL A   4       8.903   7.525  14.367  1.00 12.54           H
ATOM     11 HG11 VAL A   4       7.454   5.916  15.260  1.00 11.37           H
ATOM     12 HG12 VAL A   4       6.932   6.385  13.840  1.00 11.37           H
ATOM     13 HG13 VAL A   4       7.686   5.000  13.989  1.00 11.37           H
ATOM     14 HG21 VAL A   4       9.942   5.972  15.777  1.00 12.59           H
ATOM     15 HG22 VAL A   4      10.251   5.041  14.533  1.00 12.59           H
ATOM     16 HG23 VAL A   4      10.957   6.454  14.660  1.00 12.59           H
ATOM     17  N   TYR A   5       7.521   6.791  10.776  1.00  9.67           N
ATOM     18  CA  TYR A   5       6.396   7.390  10.072  1.00  9.24           C
ATOM     19  C   TYR A   5       5.146   6.935  10.810  1.00 10.02           C
ATOM     20  O   TYR A   5       5.000   5.738  11.127  1.00 10.62           O
ATOM     21  CB  TYR A   5       6.289   6.896   8.621  1.00  9.83           C
ATOM     22  CG  TYR A   5       7.382   7.426   7.724  1.00 12.79           C
ATOM     23  CD1 TYR A   5       8.649   6.852   7.769  1.00 13.60           C
ATOM     24  CD2 TYR A   5       7.105   8.491   6.868  1.00 12.13           C
ATOM     25  CE1 TYR A   5       9.651   7.354   6.948  1.00 14.77           C
ATOM     26  CE2 TYR A   5       8.116   8.986   6.053  1.00 14.21           C
ATOM     27  CZ  TYR A   5       9.375   8.421   6.104  1.00 14.05           C
ATOM     28  OH  TYR A   5      10.398   8.932   5.306  1.00 17.48           O
ATOM     29  H   TYR A   5       7.729   6.008  10.488  1.00  9.67           H
ATOM     30  HA  TYR A   5       6.507   8.353  10.049  1.00  9.24           H
ATOM     31  HB2 TYR A   5       6.316   5.926   8.614  1.00  9.83           H
ATOM     32  HB3 TYR A   5       5.428   7.159   8.260  1.00  9.83           H
ATOM     33  HD1 TYR A   5       8.823   6.141   8.342  1.00 13.60           H
ATOM     34  HD2 TYR A   5       6.254   8.866   6.843  1.00 12.13           H
ATOM     35  HE1 TYR A   5      10.501   6.978   6.964  1.00 14.77           H
ATOM     36  HE2 TYR A   5       7.945   9.695   5.476  1.00 14.21           H
ATOM     37  HH  TYR A   5      10.168   9.680   5.000  1.00 17.48           H
TER
'''
  pdb_inp = iotbx.pdb.input(lines=pdb_dihedral.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.set_log(null_out()) # don't print out stuff
  ph = model.get_hierarchy()
  model.process(make_restraints=True) # get restraints_manager
  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()
  h_parameterization = riding_h_manager.h_parameterization
  atoms = model.get_atoms()
  assert atoms[10].name.strip().upper() == 'HG11', "Helpers.Test(): Internal error"
  hydrogens = [ atoms[10], atoms[11], atoms[12] ] # HG1, HG2, HG3
  potentials = [ atoms[1], atoms[6] ]  # CA, CG2

  myH, myFriend = dihedralChoicesForRotatableHydrogens(hydrogens, h_parameterization, potentials)
  assert myH == atoms[12], "Helpers.Test(): Unexpected H from dihedralChoicesForRotatableHydrogens"
  assert myFriend == atoms[1], "Helpers.Test(): Unexpected friend from dihedralChoicesForRotatableHydrogens"

  gotSorry = False
  try:
    myH, myFriend = dihedralChoicesForRotatableHydrogens(hydrogens, None, potentials)
  except Exception:
    gotSorry = True
  assert gotSorry, "Helpers.Test(): Unexpected lack of exception from dihedralChoicesForRotatableHydrogens"

  #========================================================================
  # Run unit tests on rvec3 and lvec3.
  v1 = rvec3([0, 0, 0])
  v2 = rvec3([1, 0, 0])
  if hasattr(math, 'isclose'):
    assert math.isclose((v2-v1).length(), 1), "Helpers.Test(): rvec3 test failed"

  v1 = lvec3([0, 0, 0])
  v2 = lvec3([1, 0, 0])
  if hasattr(math, 'isclose'):
    assert math.isclose((v2-v1).length(), 1), "Helpers.Test(): lvec3 test failed"

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB/CIF file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    fileName = sys.argv[i]

  # This will raise an assertion failure if there is a problem
  Test(fileName)
  print('Success!')
