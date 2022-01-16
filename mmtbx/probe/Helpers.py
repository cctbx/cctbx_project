##################################################################################
#                Copyright 2021  Richardson Lab at Duke University
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
  radius = 0.25
    .type = float
    .help = Probe radius (half distance between touched atoms) (-radius in probe)

  density = 16.0
    .type = float
    .help = Probe dots per square angstrom on atom surface (-density in probe)

  worse_clash_cutoff = 0.5
    .type = float
    .help = Cutoff for worse clashes, a positive (-divworse in probe)

  clash_cutoff = 0.4
    .type = float
    .help = Cutoff for the clashes, a positive number (-divlow in probe)

  contact_cutoff = 0.25
    .type = float
    .help = Cutoff for the contact (-divhigh in probe)

  uncharged_hydrogen_cutoff = 0.6
    .type = float
    .help = Cutoff for uncharged hydrogen overlap (-hbregular in probe)

  charged_hydrogen_cutoff = 0.8
    .type = float
    .help = Cutoff for charged hydrogen overlap (-hbcharged in probe)

  bump_weight = 10.0
    .type = float
    .help = Weight applied to bump score (-bumpweight in probe)

  hydrogen_bond_weight = 4.0
    .type = float
    .help = Weight applied to hydrogen bond score (-hbweight in probe)

  gap_weight = 0.25
    .type = float
    .help = Weight applied to gap score (-gapweight in probe)

  allow_weak_hydrogen_bonds = False
    .type = bool
    .help = Separately account for weak hydrogen bonds (-LweakHbonds in probe)

  implicit_hydrogens = False
    .type = bool
    .help = Use implicit hydrogens, no water proxies (-implicit in probe)

  use_original_probe_tables = False
    .type = bool
    .help = Use the original Probe tables rather than CCTBX tables by default (for regression tests)
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
    in addition to the other constructor parameters.
    :param extraAtomInfo: Can be obtained from getExtraAtomInfo().  Holds extra atom information
    needed for scoring atoms.
    :param probePhi: Subobject of PHIL parameters for probe.  Can be obtained using
    self.params.probe from a Program Template program that includes the probe_phil_parameters
    from above in its master PHIL parameters string.
    :returns A mmtbx_probe_ext.DotScorer object.
  """
  return probeExt.DotScorer(extraAtomInfo, probePhil.gap_weight,
        probePhil.bump_weight, probePhil.hydrogen_bond_weight,
        probePhil.uncharged_hydrogen_cutoff, probePhil.charged_hydrogen_cutoff,
        probePhil.clash_cutoff, probePhil.worse_clash_cutoff,
        probePhil.contact_cutoff, probePhil.allow_weak_hydrogen_bonds)

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
    it should be a flex array of atom positions for the atoms that are in the first argument.
    :returns a dictionary with one entry for each atom that contains a list of all of
    the atoms (within the atoms list) that are bonded to it.
  """
  atomDict = {}
  for a in atoms:
    atomDict[a.i_seq] = a
  bondedNeighbors = {}
  for a in atoms:
    bondedNeighbors[a] = []
  for bp in bondProxies:
    try:
      first = atomDict[bp.i_seqs[0]]
      second = atomDict[bp.i_seqs[1]]
      bondedNeighbors[first].append(second)
      bondedNeighbors[second].append(first)
    except Exception:
      # When an atom is bonded to an atom in a different conformer (not in our atom list)
      # we just ignore it.
      pass
  return bondedNeighbors

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

def getAtomsWithinNBonds(atom, bondedNeighborLists, N, nonHydrogenN = 1e10):
  """
    Helper function to produce a list of all of the atoms that are bonded to the
    specified atoms, or to one of the atoms bonded to the specified atom, recursively
    to a depth of N.  The atom itself will not be included in the list, so an atom that
    has no bonded neighbors will always have an empty result.  This can be used to
    produce a list of excluded atoms for dot scoring.  It checks to ensure that all of
    the bonded atoms are from compatible conformations (if the original atom
    is in the empty configuration then this will return atoms from all conformations that
    are in the bonded set).
    :param atom: The atom to be tested.
    :param bondedNeighborLists: Dictionary of lists that contain all bonded neighbors for
    each atom in a set of atoms.  Should be obtained using
    mmtbx.probe.Helpers.getBondedNeighborLists().
    :param N: Depth of recursion.  N=1 will return the atoms bonded to atom.  N=2 will
    also return those bonded to these neighbors (but not the atom itself).
    :param nonHydrogenN: When neither the original atom nor the bonded atom is a Hydrogen,
    limit the depth to this value (if this value is less than N).
    :returns a list of all atoms that are bonded to atom within a depth of N.  The original
    atom is never on the list.
  """
  atomIsHydrogen = atom.element_is_hydrogen()
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
    mmtbx.probe.Helpers.getBondedNeighborLists().
    :return: True if the atom is a polar hydrogen.
  '''
  if atom.element_is_hydrogen():
    neighbors = bondedNeighborLists[atom]
    if len(neighbors) == 1 and neighbors[0].element in ['N', 'O', 'S']:
      return True
  return False

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
    models.  It first tries to find the information in CCTBX.  If it cannot, it looks
    the information up using the original C-code Probe tables and algorithms.
    :param model: Map Model Manager's Model containing all of the atoms to be described.
    PDB interpretation must have been done on the model, perhaps by calling
    model.process(make_restraints=True), with useNeutronDistances matching
    the parameter to this function.
    :param bondedNeighborLists: Lists of atoms that are bonded to each other.
    Can be obtained by calling getBondedNeighborLists().
    :param useNeutronDistances: Default is to use x-ray distances, but setting this to
    True uses neutron distances instead.  This must be set consistently with the
    PDB interpretation parameter used on the model.
    :param probePhil: None or subobject of PHIL parameters for probe.  Can be obtained using
    self.params.probe from a Program Template program that includes the probe_phil_parameters
    from above in its master PHIL parameters string.  If None, local defaults will be used.
    The following are used:
      implicit_hydrogens (bool): Default is to use distances consistent with
      explicitly-listed Hydrgoens, but setting this to True implicit-Hydrogen distances instead.
      This must be set consistently with the hydrogens in the model.
      use_original_probe_tables (bool): Do not attempt to read the data from CCTBX, use the
      original Probe tables.  This is normally the fall-back when it cannot find the data in
      CCTBX.  The Probe tables do not have accurate data on HET atoms, only standard residues.
      The values in the tables may differ from the current CCTBX values, and the Probe tables
      are not being maintained.
    :returns a ExtraAtomInfoMap with an entry for every atom in the model suitable for
    passing to the scoring functions.
  """

  warnings = ""

  # Pull parameters from PHIL parameters, if they are present.  Otherwise, set to
  # defaults
  useImplicitHydrogenDistances = False
  useProbeTablesByDefault = False
  if probePhil is not None:
    useImplicitHydrogenDistances = probePhil.implicit_hydrogens
    useProbeTablesByDefault = probePhil.use_original_probe_tables

  # Construct the AtomTypes object we're going to use, telling it whether to use neutron distances.
  at = AtomTypes.AtomTypes(useNeutronDistances, useImplicitHydrogenDistances)

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
            if not useProbeTablesByDefault:
              # See if we can find out about its Hydrogen-bonding status from the
              # model.  If so, we fill it and the vdwRadius information from
              # CCTBX.
              try:
                hb_type = model.get_specific_h_bond_type(a.i_seq)
                if isinstance(hb_type, str):
                  if hb_type == "A" or hb_type == "B":
                    extra.isAcceptor = True
                  if hb_type == "D" or hb_type == "B":
                    extra.isDonor = True

                  # For ions, the Richardsons determined in discussion with
                  # Michael Prisant that we want to use the ionic radius rather than the
                  # larger radius for all purposes.
                  # @todo Once the CCTBX radius determination discussion and upgrade is
                  # complete (ongoing as of September 2021), this check might be removed
                  # and we'll just use the CCTBX radius.
                  if a.element_is_ion():
                    warnings += "Using ionic radius for "+a.name.strip()+"\n"
                    extra.vdwRadius = model.get_specific_ion_radius(a.i_seq)
                  else:
                    extra.vdwRadius = model.get_specific_vdw_radius(a.i_seq, useImplicitHydrogenDistances)

                  # Mark aromatic ring N and C atoms as acceptors as a hack to enable the
                  # ring itself to behave as an acceptor.
                  # @todo Remove this once we have a better way to model the ring itself
                  # as an acceptor, perhaps making it a cylinder or a sphere in the center
                  # of the ring.
                  if a.element in ['C','N']:
                    if AtomTypes.IsAromatic(ag.resname, a.name):
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

                  # Mark all Carbonyl's with the Probe radius while the Richarsons and
                  # the CCTBX decide how to handle this.
                  # @todo After 2021, see if the CCTBX has the same values (1.65 and 1.80)
                  # for Carbonyls and remove this if so.  It needs to stay with these values
                  # to avoid spurious collisions per experiments run by the Richardsons in
                  # September 2021.
                  if a.name.strip().upper() == 'C':
                    if useImplicitHydrogenDistances:
                      extra.vdwRadius = 1.80
                    else:
                      extra.vdwRadius = 1.65
                    warnings += "Overriding radius for "+a.name.strip()+": "+str(extra.vdwRadius)+"\n"

                  extras.setMappingFor(a, extra)
                  continue

                # Did not find the information from CCTBX, so look it up using
                # the original Probe approach by dropping through to below
                else:
                  warnings += "Could not find "+a.name.strip()+" in CCTBX, using Probe tables\n"
              except Exception as e:
                # Warn and drop through to below.
                warnings += ("Could not look up "+a.name.strip()+" in CCTBX "+
                  "(perhaps interpretation was not run on the model?), using Probe tables"+
                  ": "+str(e)+"\n")

            # Did not find what we were looking for in CCTBX, so drop through to Probe.
            # Probe always returns the result we want as the VdW radius, even for ions.
            extra, warn = at.FindProbeExtraAtomInfo(a)
            if len(warn) > 0:
              warnings += "  Probe says: "+warn+"\n"

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
  for a in atoms:
    chainID = a.parent().parent().parent().id
    resName = a.parent().resname.upper()
    resID = str(a.parent().parent().resseq_as_int())
    acceptorChoices = ["noAcceptor","isAcceptor"]
    donorChoices = ["noDonor","isDonor"]
    metallicChoices = ["noMetallic","isMetallic"]
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

def getPhantomHydrogensFor(atom, spatialQuery, extraAtomInfo, minOccupancy, acceptorOnly = False,
      placedHydrogenRadius = 1.05):
  """
    Get a list of phantom Hydrogens for the atom specified, which is asserted to be an Oxygen
    atom for a water.
    :param atom: The Oxygen that is to have phantoms added to it.
    :param spatialQuery: mmtbx_probe_ext.SpatialQuery structure to rapidly determine which atoms
    are within a specified distance of a location.
    :param extraAtomInfo: mmtbx_probe_ext.ExtraAtomInfo mapper that provides radius and other
    information about atoms beyond what is in the pdb.hierarchy.  Used here to determine
    which atoms may be acceptors.
    :param minOccupancy: Minimum occupancy for an atom to be considered.
    :param acceptorOnly: Only allow bonds with atoms that are acceptors when this is True.
    This is false by default because Reduce needs to check whether the bonded atom is either
    an acceptor or a possible flipped position of an acceptor, and that is not something that
    can be determined at the time we're placing phantom hydrogens.  In that case, we want to
    include all possible interactions and weed them out during optimization.
    :param placedHydrogenRadius: Maximum radius to use for placed Phantom Hydrogen atoms.
    The Hydrogens are placed at the optimal overlap distance so may be closer than this.
    :return: List of new atoms that make up the phantom Hydrogens, with only their name and
    element type and xyz positions filled in.  They will have i_seq 0 and they should not be
    inserted into a structure.
  """

  ret = []

  # Get the list of nearby atoms.  The center of the search is the water atom
  # and the search radius is 4 (these values are pulled from the Reduce C++ code).
  maxDist = 4.0
  nearby = spatialQuery.neighbors(atom.xyz, 0.001, maxDist)

  # Candidates for nearby atoms.  We use this list to keep track of one ones we
  # have already found so that we can compare against them to only get one for each
  # aromatic ring.
  class Candidate(object):
    def __init__(self, atom, overlap):
      self._atom = atom
      self._overlap = overlap
  candidates = []

  for a in nearby:
    # Only check atoms in compatible conformations.
    if not compatibleConformations(atom, a):
      continue

    # Check to ensure the occupancy of the neighbor is above threshold and that it is
    # close enough to potentially bond to the atom.
    OH_BOND_LENGTH = 1.0
    overlap = ( (rvec3(atom.xyz) - rvec3(a.xyz)).length()  -
                (placedHydrogenRadius + extraAtomInfo.getMappingFor(a).vdwRadius + OH_BOND_LENGTH) )
    if overlap < -0.1 and a.occ > minOccupancy and a.element != "H":
      if not acceptorOnly or extraAtomInfo.getMappingFor(a).isAcceptor:
        # If we have multiple atoms in the same Aromatic ring (part of the same residue)
        # we only point at the closest one.  To ensure this, we check all current candidates
        # and if we find one that is on the same aromatic ring then we either ignore this new
        # atom (if it is further) or replace the existing one (if it is closer).
        skip = False
        if AtomTypes.IsAromatic(a.parent().resname.strip().upper(), a.name.strip().upper()):
          for c in candidates:
            # See if we belong to the same atom group and are both ring acceptors.  If so, we need to replace
            # or else squash this atom.
            if (AtomTypes.IsAromatic(c._atom.parent().resname.strip().upper(), c._atom.name.strip().upper()) and
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
  # Make most of their characteristics (including i_seq) copied from the source Oxygen.
  # The element, name, and location are modified.
  for c in candidates:
    h = pdb.hierarchy.atom(atom.parent(), atom)
    h.element = "H"
    h.name = " H?"

    # Place the hydrogen pointing from the Oxygen towards the candidate at a distance
    # of 1 plus an offset that is clamped to the range -1..0 that is the sum of the overlap
    # and the best hydrogen-bonding overlap.
    BEST_HBOND_OVERLAP=0.6
    distance = 1.0 + max(-1.0, min(0.0, c._overlap + BEST_HBOND_OVERLAP))
    try:
      normOffset = (rvec3(c._atom.xyz) - rvec3(atom.xyz)).normalize()
      h.xyz = rvec3(atom.xyz) + distance * normOffset
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
    mmtbx.probe.Helpers.getBondedNeighborLists().
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

##################################################################################
# Helper functions to make things that are compatible with vec3_double so
# that we can do math on them.  We need a left-hand and right-hand one so that
# we can make both versions for multiplication.
def rvec3 (xyz) :
  return scitbx.matrix.rec(xyz, (3,1))
def lvec3 (xyz) :
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
ATOM    448  ND1 HIS A  61      26.168  34.981   3.980  1.00  9.06           N
ATOM    449  CD2 HIS A  61      25.174  33.397   5.004  1.00 11.08           C
ATOM    450  CE1 HIS A  61      24.867  35.060   3.688  1.00 12.84           C
ATOM    451  NE2 HIS A  61      24.251  34.003   4.297  1.00 11.66           N
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

  #========================================================================
  # Run unit test on getExtraAtomInfo().  We use a specific PDB snippet
  # for which we know the answer and then we verify that the results are what
  # we expect.

  # Spot check the values on the atoms for standard, neutron distances, implicit hydrogen distances,
  # and original Probe results.
  standardChecks = [
    # Name, vdwRadius, isAcceptor, isDonor
    ["N",   1.55, False, True ],
    ["ND1", 1.55, True,  True ],
    ["C",   1.65, False, False],
    ["CB",  1.7,  False, False],
    ["O",   1.4,  True,  False],
    ["CD2", 1.75, False, False]
  ]
  neutronChecks = [
    # Name, vdwRadius, isAcceptor, isDonor
    ["N",   1.55, False, True ],
    ["ND1", 1.55, True,  True ],
    ["C",   1.65, False, False],
    ["CB",  1.7,  False, False],
    ["O",   1.4,  True,  False],
    ["CD2", 1.75, False, False]
  ]
  implicitChecks = [
    # Name, vdwRadius, isAcceptor, isDonor
    ["N",   1.6,  False, True ],
    ["ND1", 1.6,  True,  True ],
    ["C",   1.8,  False, False],
    ["CB",  1.92,  False, False],
    ["O",   1.52,  True,  False],
    ["CD2", 1.74,  False, False]
  ]
  probeChecks = [
    # Name, vdwRadius, isAcceptor, isDonor
    ["N",   1.55, False, False],
    ["ND1", 1.55, True,  False],
    ["C",   1.65, False, False],
    ["CB",  1.7,  False, False],
    ["O",   1.4,  True,  False],
    ["CD2", 1.7,  False, False]
  ]

  # Situations to run the test in and expected results:
  cases = [
    # Use neutron distances, use implicit distances, use probe values, expected results
    [False, False, False, standardChecks],
    [True,  False, False, neutronChecks],
    [False, True,  False, implicitChecks],
    [False, False, True,  probeChecks]
  ]

  for cs in cases:
    useNeutronDistances = cs[0]
    useImplicitHydrogenDistances = cs[1]
    useProbe = cs[2]
    checks = cs[3]
    runType = "; neutron,implicit,probe = "+str(useNeutronDistances)+","+str(useImplicitHydrogenDistances)+","+str(useProbe)

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
    class philLike:
      def __init__(self, useImplicitHydrogenDistances = False, useProbe = False):
        self.implicit_hydrogens = useImplicitHydrogenDistances
        self.use_original_probe_tables = useProbe
    philArgs = philLike(useImplicitHydrogenDistances, useProbe)
    extras = getExtraAtomInfo(model,bondedNeighborLists,
      useNeutronDistances=useNeutronDistances,probePhil=philArgs).extraAtomInfo

    # Get the atoms for the first model in the hierarchy.
    atoms = model.get_hierarchy().models()[0].atoms()

    for a in atoms:
      e = extras.getMappingFor(a)
      for c in checks:
        if a.name.strip() == c[0]:
          if hasattr(math, 'isclose'):
            assert math.isclose(e.vdwRadius, c[1]), ("Helpers.Test(): Bad radius for "+a.name+": "
            +str(e.vdwRadius)+" (expected "+str(c[1])+")"+runType)
          assert e.isAcceptor == c[2], "Helpers.Test(): Bad Acceptor status for "+a.name+": "+str(e.isAcceptor)+runType
          assert e.isDonor == c[3], "Helpers.Test(): Bad Donor status for "+a.name+": "+str(e.isDonor)+runType
          # Check the ability to set and check Dummy/Phantom Hydrogen status
          assert e.isDummyHydrogen == False, "Helpers.Test(): Bad Dummy Hydrogen status for "+a.name+runType
          e.isDummyHydrogen = True
          assert e.isDummyHydrogen == True, "Helpers.Test(): Can't set DummyHydrogen status for "+a.name+runType

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
    ag.resname = 'HIS'
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
            a.name = 'ND1'
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
    for occThresh in [0.4,1.0]:
      for acceptorOnly in [False, True]:
        # Check that we get the expected number of contacts
        if acceptorOnly:
          expected = 1 + 3 # One for the aromatics, three others
        else:
          expected = 7 # Only one of the acceptor Aromatics but all other atoms
        if occThresh > 0.5:
          expected = 0
        ret = getPhantomHydrogensFor(o, sq, extrasMap, occThresh, acceptorOnly)
        assert len(ret) == expected, "Helpers.Test() Unexpected count during Phantom Hydrogen placement: "+str(len(ret))

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

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = getBondedNeighborLists(atoms, bondProxies)

  # Check the counts in the neighbor lists to make sure they match what we expect
  neighborCounts = {"N": 1, "CA": 3, "C": 2, "O": 1, "CB": 2,
                    "CG": 3, "ND1": 2, "CD2": 2, "CE1":2, "NE2": 2}
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

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = getBondedNeighborLists(atoms, bondProxies)

  model = dm.get_model().get_hierarchy().models()[0]

  for a in model.atoms():
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
  # NOTE: This re-uses the bondedNeighborLists test results from above
  N4 = None
  for a in atoms:
    if a.name.strip().upper() == 'N4':
      N4 = a
  assert N4 is not None, ("Helpers.Test(): Could not find N4 (internal failure)")
  # When clamped, we can't go further than the non-Hydrogen bound except for Hydrogens
  nestedNeighborsForN4 = [ None, 3, 5, 8, 9, 9, 9]
  for N in range(1,7):
    count = len(getAtomsWithinNBonds(N4, bondedNeighborLists, N, 3))
    assert count == nestedNeighborsForN4[N], ("Helpers.Test(): Nested clamped count for "+N4.name.strip()+
        " for N = "+str(N)+" was "+str(count)+", expected "+str(nestedNeighborsForN4[N]))
  # When unclamped, we can traverse all the way to the end in all cases
  nestedNeighborsForN4 = [ None, 3, 5, 8, 11, 12, 15]
  for N in range(1,7):
    count = len(getAtomsWithinNBonds(N4, bondedNeighborLists, N))
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
    count = len(getAtomsWithinNBonds(H41, bondedNeighborLists, N, 3))
    assert count == nestedNeighborsForH41[N], ("Helpers.Test(): Nested clamped count for "+H41.name.strip()+
        " for N = "+str(N)+" was "+str(count)+", expected "+str(nestedNeighborsForH41[N]))

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

  ret = getExtraAtomInfo(model, bondedNeighborLists)
  # User code should check for and print any warnings.
  #if len(ret.warnings) > 0:
  #  print('Warnings returned by getExtraAtomInfo():\n'+ret.warnings)

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
  print('OK')
