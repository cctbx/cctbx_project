##################################################################################
#                Copyright 2021-2023 Richardson Lab at Duke University
#
# Licensed under the Apache License, Version 2.0 (the "License");
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

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import

import argparse, re

from boost_adaptbx import graph
from boost_adaptbx.graph import connected_component_algorithm as cca

from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
from iotbx.pdb import common_residue_names_get_class
import mmtbx
from scitbx.array_family import flex

# To enable addition of Hydrogens
# @todo See if we can remove the shift and box once reduce_hydrogen is complete
from cctbx.maptbx.box import shift_and_box_model
from mmtbx.hydrogens import reduce_hydrogen

import boost_adaptbx.boost.python as bp
bp.import_ext("mmtbx_probe_ext")
import mmtbx_probe_ext as probeExt
from mmtbx.probe import Helpers
from mmtbx.reduce import Movers, InteractionGraph

bp.import_ext("mmtbx_reduce_ext")
from mmtbx_reduce_ext import OptimizerC, Optimizers_test

##################################################################################
# This file includes a set of functions and classes that implement placement and optimization of
# Reduce's "Movers".

##################################################################################
# Module variables that affect the way computations are handled

_DoCliqueOptimizationInC = True

##################################################################################
# Helper functions

def _VerboseCheck(verbosity, level, message):
  # Returns "" if the level is less than verbosity level, message if it is at or above
  if verbosity >= level:
    return " "*level + message
  else:
    return ""

_lastTime = None
def _ReportTiming(verbosity, message):
  """Use None message to start the timer without printing.
  """
  import time

  global _lastTime
  if message is None:
    _lastTime = time.time()
    return
  curTime = time.time()
  diff = curTime - _lastTime
  _lastTime = curTime
  return _VerboseCheck(verbosity, 2,"Time to {}: {:0.3f}".format(message,diff)+"\n")

def AlternatesInModel(model):
  """Returns a set of altloc names of all conformers in all chains.
  :param model: pdb.hierarchy.model to search for alternates.
  :return: Set of strings.  The set is will include only the empty string if no
  chains in the model have alternates.  It has an additional entry for every altloc
  found in every chain.
  """
  ret = set(['']) # We always add this so that it will be present even when there are alternates
  for c in model.chains():
    for alt in c.conformers():
      ret.add(alt.altloc)
  return ret

def GetAtomsForConformer(model, conf):
  """Returns a list of atoms in the named conformer.  It also includes atoms from
  the first conformation in chains that do not have this conformation, to provide a
  complete model.  For a model whose chains have only one conformer, it will return
  all of the atoms.
  :param model: pdb.hierarchy.model to search for atoms.
  :param conf: String name of the conformation to find.  Can be "", which finds the
  default conformer; if there is no empty conformation, then it will
  pick the first available conformation for each atom group.
  :return: List of atoms consistent with the specified conformer.
  """
  ret = []
  for ch in model.chains():
    confs = ch.conformers()
    which = 0
    for i in range(1,len(confs)):
      if confs[i].altloc == conf:
        which = i
        break
    ret += confs[which].atoms()
  return ret

def _ResNameAndID(a):
  """Make a string describing the residue and chain for the specfied atom.
  """
  chainID = a.parent().parent().parent().id
  resName = a.parent().resname.strip().upper()
  resID = str(a.parent().parent().resseq_as_int())
  altLoc = a.parent().altloc
  # Don't print the code if it is a space (blank).
  insertionCode = a.parent().parent().icode.strip()
  return "chain "+str(chainID)+" "+altLoc+resName+" "+resID+insertionCode

##################################################################################
# Helper classes

class FlipMoverState(object):
  # Holds information needed to identify a flip Mover within a model file and
  # whether or not it is flipped and has its angles adjusted.
  def __init__(self, moverType, modelId, altId, chain, resName, resIdWithICode, flipped, fixedUp):
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
    self.flipped = flipped      # Boolean Whether the Mover is flipped in this configuration
    self.fixedUp = fixedUp      # Boolean Whether the fixup has been done on the angles

  def __str__(self):
    return "{} {} '{}' {} {} {}{} {} {}".format(self.moverType, self.modelId, self.altId,
      self.chain, self.resName, self.resId, self.iCode, self.flipped, self.fixedUp)
  def __repr__(self):
      return "Optimizers.FlipMoverState({})".format(str(self))

##################################################################################
# Optimizer, which wraps the OptimizerC class to do the optimization:

class Optimizer(object):

  def __init__(self, probePhil, addFlipMovers, model, modelIndex = 0, altID = None,
                bondedNeighborDepth = 4,
                useNeutronDistances = False,
                minOccupancy = 0.02,
                preferenceMagnitude = 1.0,
                nonFlipPreference = 0.5,
                skipBondFixup = False,
                flipStates = '',
                verbosity = 1,
                cliqueOutlineFileName = None,
                fillAtomDump = True
              ):
    """Constructor.  This is the wrapper class for the C++ OptimizerC and
    it implements the machinery that finds and optimizes Movers.
    :param probePhil: Phil parameters to Probe to be passed on when subroutines are called.
    :param addFlipMovers: Do we add flip Movers along with other types?
    :param model: iotbx model (a group of hierarchy models).  Can be obtained using
    iotbx.map_model_manager.map_model_manager.model().  The model must have Hydrogens,
    which can be added using mmtbx.hydrogens.reduce_hydrogen.place_hydrogens().get_model().
    It must have a valid unit cell, which can be helped by calling
    cctbx.maptbx.box.shift_and_box_model().  It must have had PDB interpretation run on it,
    which can be done using model.process(make_restraints=True) with PDB
    interpretation parameters and hydrogen placement matching the value of the
    useNeutronDistances parameter described below.
    :param modelIndex: Identifies which index from the hierarchy is to be selected.
    If this value is None, optimization will be run sequentially on every model in the
    hierarchy. This is 1-based, the first model is 1.
    :param altID: The conformer alternate location specifier to use.  The value "" will
    cause it to run on the first conformer found in each model.  If this is set to None,
    optimization will be run sequentially for every conformer in the model, starting with
    the last and ending with the first.  This will leave the initial conformer's values as the
    final location for atoms that are not inside a conformer or are in the first conformer.
    :param bondedNeighborDepth: How many hops to ignore bonding when doing Probe calculations.
    A depth of 3 will ignore my bonded neighbors and their bonded
    neighbors and their bonded neighbors.
    :param useNeutronDistances: False will use X-ray/electron cloud distances.  If set to
    True, it will use neutron (nuclear) distances.  This must be set consistently with the
    values used to generate the hydrogens and to run PDB interpretation.
    :param minOccupancy: Minimum occupancy for an atom to be considered in the Probe score.
    :param preferenceMagnitude: Multiplier for the preference energies expressed
    by some Movers for particular orientations.
    :param nonFlipPreference: Preference for not flipping Movers that are flips.  This keeps
    them from being flipped unless their score is significantly better than in the original position.
    :param skipBondFixup: Should we do fixup on Movers or just leave them flipped?  We always do Hydrogen
    removal and fixup for Histidines that are not placed as Movers, but for flips that
    are Movers we don't adjust the bond angles.  This is overridden for Flips that are
    specified in flipStates.
    :param flipStates: String with comma-separated entries. Each entry has the form
    (without the single quotes) '1 . A HIS 11H Flipped AnglesAdjusted'. These are space-separated values.
    The first word is the model number, starting with 1. The second is the lower-case alternate, or
    '.' for all alternates (also use this when there are no alternates in the file).
    The third is the chain ID. The fourth is the residue name. The fifth is the residue id,
    which may include an insertion code as its last character. The sixth is either Flipped or Unflipped.
    If it is Flipped, then another word is added -- AnglesAdjusted or AnglesNotAdjusted,
    specifying whether to do the three-point dock to adjust the bond angles after the flip.
    An example with several entries is (again, no quotes are included:
    '1 a A HIS 11H Unflipped,1 b A ASN 15 Flipped AnglesNotAdjusted,1 . B GLN 27 Flipped AnglesAdjusted'. Any
    Flip Movers that would be placed at the specified location are instead locked in the
    specified configuration.
    :param verbosity: Default value of 1 reports standard information.
    Value of 2 reports timing information.
    Setting it to 0 removes all informational messages.
    Setting it above 1 provides additional debugging information.
    :param cliqueOutlineFileName: Name of file to write Kinemage with outlines of Movers to.
    This file holds spheres showing all possible locations for each atom in each Mover in different
    colors as one master. It shows the outlines expanded by the probe radius, with a single color
    for each clique, as another master. These are useful for determining why the cliques are as
    they are.
    :param fillAtomDump: If true, fill in the atomDump string with the atom information.
    This can take a long time to do, so the caller may want to turn it off if they don't need it.
    """

    ################################################################################
    # Store the parameters that will be accessed by other methods
    self._probePhil = probePhil
    self._bondedNeighborDepth = bondedNeighborDepth
    self._probeRadius = probePhil.probe_radius
    self._useNeutronDistances = useNeutronDistances
    self._probeDensity = probePhil.density
    self._minOccupancy = minOccupancy
    self._preferenceMagnitude = preferenceMagnitude
    self._nonFlipPreference = nonFlipPreference
    self._skipBondFixup = skipBondFixup
    self._flipStates = flipStates
    self._verbosity = verbosity
    self._cliqueOutlineFileName = cliqueOutlineFileName

    ################################################################################
    # Initialize internal variables.
    self._infoString = ""
    self._warningString = ""
    self._atomDump = ""
    self._waterOccCutoff = 0.66 # @todo Make this a parameter, -WaterOCCcutoff param in reduce
    self._waterBCutoff = 40.0   # @todo Make this a parameter, -WaterBcutoff param in reduce
    self._numCalculated = 0
    self._numCached = 0

    ################################################################################
    # Get the Cartesian positions of all of the atoms in the entire model and find
    # the bond proxies for all of them.  The proxies will not attempt to bond atoms
    # that are in different model indices.
    _ReportTiming(self._verbosity, None) # Reset timer
    carts = flex.vec3_double()
    for a in model.get_atoms():
      carts.append(a.xyz)
    self._infoString += _ReportTiming(self._verbosity, "get coordinates")
    bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
    self._infoString += _ReportTiming(self._verbosity, "compute bond proxies")

    ################################################################################
    # Get the bonded neighbor lists for all of the atoms in the model, so we don't get
    # failures when we look up an atom from another in Helpers.getExtraAtomInfo().
    # We won't get bonds between atoms in different conformations.
    bondedNeighborLists = Helpers.getBondedNeighborLists(model.get_atoms(), bondProxies)
    self._infoString += _ReportTiming(self._verbosity, "compute bonded neighbor lists")

    ################################################################################
    # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
    # This is done for all atoms in the model.
    # This must operate on the entire model, not just the current model index.
    ret = Helpers.getExtraAtomInfo(
      model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=self._useNeutronDistances, probePhil=self._probePhil)
    self._extraAtomInfo = ret.extraAtomInfo
    self._infoString += ret.warnings
    self._infoString += _ReportTiming(self._verbosity, "get extra atom info")
    self._warningString += ret.warnings

    ################################################################################
    # Run optimization for every desired conformer and every desired model, calling
    # placement and then a derived-class single optimization routine for each.  When
    # the modelIndex or altID is None, that means to run over all available cases.
    # For alternates, if there is a non-empty (not "" or " ") case, then we run backwards
    # from the last to the first but do not run for the empty case; if there are only
    # empty cases, then we run just once.  We run the models in order, all of them when
    # None is specified and the specified one if it is specified.

    try:
      model.setup_riding_h_manager()
    except Exception:
      # If an optimized has already been run on the model, our riding H manager will
      # have already been set up.  When one has already been set up, this causes
      # an exception and it doesn't make a new one.
      pass
    riding_h_manager = model.get_riding_h_manager()
    h_parameterization = riding_h_manager.h_parameterization

    # Find the single-hydrogen rotators.
    # Dorothee provided this faster approach that uses the riding_h_manager.
    rotatableHydrogens = flex.size_t()
    for p in h_parameterization:
      if p is not None:
        if p.htype  == 'alg1b':
          rotatableHydrogens.append(p.ih)

    self._infoString += _ReportTiming(self._verbosity, "select rotatable hydrogens")

    startModelIndex = 0
    stopModelIndex = len(model.get_hierarchy().models())
    if modelIndex is not None:
      # The command-line parameter matches the name of the model in the model file, which
      # starts with 1. The internal indexing starts with 0. So we subtract one.
      startModelIndex = (modelIndex - 1)
      stopModelIndex = startModelIndex + 1
    for mi in range(startModelIndex, stopModelIndex):
      # Get the specified model from the hierarchy.
      myModel = model.get_hierarchy().models()[mi]

      ################################################################################
      # Store the states (position and extra atom info) of all of the atoms in this model
      # so that we can restore it for atoms in a given alternate configuration before optimizing
      # each new alternate. Looked up by i_seq.
      initialAtomPositions = {}
      initialExtraAtomInfos = {}
      for a in myModel.atoms():
        initialAtomPositions[a.i_seq] = a.xyz;
        initialExtraAtomInfos[a.i_seq] = probeExt.ExtraAtomInfo(self._extraAtomInfo.getMappingFor(a))

      # Get the list of alternate conformation names present in all chains for this model.
      # If there is more than one result, remove the empty results and then sort them
      # in reverse order so we finalize all non-alternate ones to match the first.
      alts = AlternatesInModel(myModel)
      if len(alts) > 1:
        alts.discard("")
        alts.discard(" ")
      alts = sorted(list(alts), reverse=True)
      self._infoString += _ReportTiming(self._verbosity, "compute alternates")

      # If there is a specified alternate, use it.
      if altID is not None:
        alts = [altID]

      # Clear the Movers list for each model.  It will be retained from one alternate to the next.
      self._movers = []
      firstAlt = True
      for alt in alts:
        # If we are doing the second or later alternate, place all atoms that are in a compatible alternate
        # back into their initial configuration so we start from the same state we would have if this were
        # the only alternate being tested.  This will ensure that we get compatible outputs when run either
        # way.
        if firstAlt:
          firstAlt = False
        else:
          for a in myModel.atoms():
            if a.parent().altloc in ['', ' ', alt]:
              a.xyz = initialAtomPositions[a.i_seq]
              self._extraAtomInfo.setMappingFor(a, initialExtraAtomInfos[a.i_seq])

        # Tell about the run we are currently doing.
        self._infoString += _VerboseCheck(self._verbosity, 1,"Running Reduce optimization on model index "+str(mi)+
          ", alternate '"+alt+"'\n")
        self._infoString += _VerboseCheck(self._verbosity, 1,"  bondedNeighborDepth = "+str(self._bondedNeighborDepth)+"\n")
        self._infoString += _VerboseCheck(self._verbosity, 1,"  probeRadius = "+str(self._probeRadius)+"\n")
        self._infoString += _VerboseCheck(self._verbosity, 1,"  useNeutronDistances = "+str(self._useNeutronDistances)+"\n")
        self._infoString += _VerboseCheck(self._verbosity, 1,"  probeDensity = "+str(self._probeDensity)+"\n")
        self._infoString += _VerboseCheck(self._verbosity, 1,"  minOccupancy = "+str(self._minOccupancy)+"\n")
        self._infoString += _VerboseCheck(self._verbosity, 1,"  preferenceMagnitude = "+str(self._preferenceMagnitude)+"\n")

        # Get the atoms from the specified conformer in the model (the empty string is the name
        # of the first conformation in the model; if there is no empty conformation, then it will
        # pick the first available conformation for each atom group.
        self._atoms = GetAtomsForConformer(myModel, alt)

        ################################################################################
        # Reset the timer
        _ReportTiming(self._verbosity, None)

        ################################################################################
        # Construct the spatial-query information needed to quickly determine which atoms are nearby
        self._spatialQuery = Helpers.createSpatialQuery(self._atoms, self._probePhil)
        self._infoString += _ReportTiming(self._verbosity, "construct spatial query")

        ################################################################################
        # Find the radius of the largest atom we'll have to deal with.
        maxVDWRad = 1
        for a in self._atoms:
          maxVDWRad = max(maxVDWRad, self._extraAtomInfo.getMappingFor(a).vdwRadius)
        self._maximumVDWRadius = maxVDWRad

        ################################################################################
        # Make a set of atoms that are to be deleted based on analysis.  It is initially
        # empty, keeping all of the added Hydrogens in the model.
        self._deleteMes = set()

        ################################################################################
        # Placement of water phantom Hydrogens, including adding them to our 'atoms' list
        # and the spatial query but not adding them to the hierarchy.  This must be done after
        # the bond proxies are constructed and the Movers have been placed so that these
        # fake atoms do not confuse placement.
        # These atoms will not be part of any atom group, and we'll fill in their extra atom
        # info directly.

        # @todo Look up the radius of a water Hydrogen.  This may require constructing a model with
        # a single water in it and asking about the hydrogen radius.  This could also become a
        # Phil parameter.  Also look up the OH bond distance rather than hard-coding it here.
        phantomHydrogenRadius = 1.05
        placedHydrogenDistance = 0.84
        if useNeutronDistances:
          phantomHydrogenRadius = 1.0
          placedHydrogenDistance = 0.98

        # Find every Oxygen that is part of a water and get the Phantom Hydrogens for it
        # unless it has out-of-bounds parameter values.  If the water Oxygen has out-of-bounds
        # parameters, then we remove it from the atoms to be considered (but not from the
        # model) by removing them from the atom list and from the spatial query structure.
        phantoms = [] # List of tuples, the Phantom and its parent Oxygen
        watersToDelete = []
        maxISeq = Helpers.getMaxISeq(model)
        for a in self._atoms:
          if a.element == 'O' and common_residue_names_get_class(name=a.parent().resname) == "common_water":
            if a.occ >= self._waterOccCutoff and a.b < self._waterBCutoff:

              # We're an acceptor and not a donor.
              ei = self._extraAtomInfo.getMappingFor(a)
              ei.isDonor = False
              ei.isAcceptor = True
              self._extraAtomInfo.setMappingFor(a, ei)

              newPhantoms = Helpers.getPhantomHydrogensFor(maxISeq, a, self._spatialQuery, self._extraAtomInfo, self._minOccupancy,
                              False, phantomHydrogenRadius, placedHydrogenDistance)
              if len(newPhantoms) > 0:
                resNameAndID = _ResNameAndID(a)
                self._infoString += _VerboseCheck(self._verbosity, 3,"Added {} phantom Hydrogens on {}\n".format(len(newPhantoms), resNameAndID))
                for p in newPhantoms:
                  self._infoString += _VerboseCheck(self._verbosity, 5,"Added phantom Hydrogen at "+str(p.xyz)+"\n")
                  phantoms.append( (p,a) )

            else:
              # Occupancy or B factor are out of bounds, so remove this atom from consideration.
              self._infoString += _VerboseCheck(self._verbosity, 3,"Ignoring "+
                a.name.strip()+" "+a.parent().resname.strip()+" "+str(a.parent().parent().resseq_as_int())+
                " "+str(a.parent().parent().parent().id)+
                " with occupancy "+str(a.occ)+" and B factor "+str(a.b)+"\n")
              watersToDelete.append(a)

        if len(watersToDelete) > 0:
          self._infoString += _VerboseCheck(self._verbosity, 1,"Ignored "+str(len(watersToDelete))+" waters due to occupancy or B factor\n")
          for a in watersToDelete:
            self._atoms.remove(a)
            self._spatialQuery.remove(a)

        if len(phantoms) > 0:

          # Add these atoms to the list of atoms we deal with.
          # Add these atoms to the spatial-query structure.
          # Insert ExtraAtomInfo for each of these atoms, marking each as a dummy and as a donor.
          # Add to the bondedNeighborList with their parent Oxygen as bonded one way so that dots on
          # a Phantom Hydrogen within its Oxygen will be excluded.  Do not mark the Oxygen as being
          # bonded to the Phantom Hydrogen to avoid having it mask collisions between the Oxygen and
          # other atoms.
          origCount = len(self._atoms)
          for p in phantoms:
            self._atoms.append(p[0])
            self._spatialQuery.add(p[0])
            eai = probeExt.ExtraAtomInfo(phantomHydrogenRadius, False, True, True)
            self._extraAtomInfo.setMappingFor(p[0], eai)
            bondedNeighborLists[p[0]] = [p[1]]

          self._infoString += _VerboseCheck(self._verbosity, 1,"Added "+str(len(phantoms))+" phantom Hydrogens on waters")
          self._infoString += _VerboseCheck(self._verbosity, 1," (Old total "+str(origCount)+", new total "+str(len(self._atoms))+")\n")
        self._infoString += _ReportTiming(self._verbosity, "place water phantom Hydrogens")

        ################################################################################
        # Fix up the donor status for all of the atoms now that we've added the final explicit
        # Phantom Hydrogens.
        Helpers.fixupExplicitDonors(self._atoms, bondedNeighborLists, self._extraAtomInfo)
        self._infoString += _ReportTiming(self._verbosity, "fixup explicit doners")

        ################################################################################
        # Get the list of Movers using the _PlaceMovers private function.
        # The list of rotatable hydrogens comes from the global model, not just the current
        # model index.  However, we only place on atoms that are also in self._atoms, which
        # only includes those from the current model index.
        # NOTE: We must do this after fixupExplicitDonors() so that the Hydrogens are properly
        # marked as donors.
        deleteAtoms = self._PlaceMovers(self._atoms, rotatableHydrogens, bondedNeighborLists, h_parameterization,
                           addFlipMovers)
        self._infoString += _VerboseCheck(self._verbosity, 1,"Inserted "+str(len(self._movers))+" Movers\n")
        self._infoString += _VerboseCheck(self._verbosity, 1,'Marked '+str(len(deleteAtoms))+' atoms for deletion\n')
        self._infoString += _ReportTiming(self._verbosity, "place movers")

        ################################################################################
        # Add the atoms that were unconditionally marked for deletion during placement
        # to the set of atoms to delete.
        self._deleteMes = self._deleteMes.union(deleteAtoms)

        ################################################################################
        # Initialize the Movers to their starting coarse positions.
        for m in self._movers:
          pr = m.CoarsePositions()
          self._setMoverState(pr, 0)
        self._infoString += _ReportTiming(self._verbosity, "initialize Movers")

        ################################################################################
        # Compute the interaction graph, of which each connected component is a Clique.
        # Get a list of singleton Cliques and a list of other Cliques.  Keep separate lists
        # of the singletons and the groups.
        self._interactionGraph, self._atomMoverLists = InteractionGraph.InteractionGraphAllPairs(self._movers,
          self._extraAtomInfo, probeRadius=self._probeRadius)
        components = cca.connected_components( graph = self._interactionGraph )
        maxLen = 0
        singletonCliques = []   # Each entry is a list of integer indices into models with one entry
        groupCliques = []       # Each entry is a list of integer indices into models with >1 entry
        for c in components:
          if len(c) == 1:
            singletonCliques.append(c)
          else:
            groupCliques.append(c)
          if len(c) > maxLen:
            maxLen = len(c)
        self._infoString += _VerboseCheck(self._verbosity, 1,"Found "+str(len(components))+" Cliques ("+
            str(len(singletonCliques))+" are singletons); largest Clique size = "+
            str(maxLen)+"\n")
        self._infoString += _ReportTiming(self._verbosity, "compute interaction graph")

        # If we've been asked to write an interaction graph file, do so.
        if self._cliqueOutlineFileName:
          with open(self._cliqueOutlineFileName, 'w') as f:
            f.write(self._InteractionKinemage(groupCliques))

        ################################################################################
        # Determine excluded atoms to a specified hop count for each atom that will
        # be moved.  Make a dictionary of lists that includes all atoms in all Movers.

        # Get the set of all atoms that can be returned from all conformations of all Movers.
        moverAtoms = set()
        for m in self._movers:
          for a in m.CoarsePositions().atoms:
            moverAtoms.add(a)
          for a in m.FixUp(0).atoms:
            moverAtoms.add(a)

        # Get the excluded list for each atom in the set, making a dictionary.
        # We go at most 3 hops unless one end of the chain has a hydrogen.
        # Look up excluded atoms by i_seq.
        self._excludeDict = {}
        for a in moverAtoms:
          self._excludeDict[a.i_seq] = mmtbx.probe.Helpers.getAtomsWithinNBonds(a,
            bondedNeighborLists, self._extraAtomInfo, self._probeRadius, self._bondedNeighborDepth, 3)
        self._infoString += _ReportTiming(self._verbosity, "determine excluded atoms")

        ################################################################################
        # Construct dot-sphere cache.
        # This must be done after the phantom Hydrogens have been added so that they will be included.
        dotSphereCache = Helpers.createDotSphereCache(self._probePhil)

        ################################################################################
        # Contruct the DotScorer object we'll use to score the dots.
        self._dotScorer = Helpers.createDotScorer(self._extraAtomInfo, self._probePhil)
        self._infoString += _ReportTiming(self._verbosity, "construct dot scorer")

        ################################################################################
        # Construct C++ optimizer.
        optC = OptimizerC(self._verbosity, self._preferenceMagnitude,
                          self._maximumVDWRadius, self._minOccupancy, self._probeRadius, self._probeDensity,
                          self._atoms,
                          self._excludeDict, self._dotScorer, dotSphereCache, self._atomMoverLists,
                          self._spatialQuery, self._extraAtomInfo, self._deleteMes)
        self._infoString += _ReportTiming(self._verbosity, "construct OptimizerC")

        ################################################################################
        # Compute and record the initial score for each Mover in its info
        optC.Initialize(self._movers)
        for m in self._movers:
          self._moverInfo[m] += " Initial score: {:.2f}".format(optC.GetHighScore(m))

        ################################################################################
        # Call internal methods to optimize the single-element Cliques and then to optimize
        # the multi-element Cliques and then to do independent fine adjustment of all
        # Cliques.  Subclasses should overload the called routines, but the global approach
        # taken here will be the same for all of them.  If we want to change the recipe
        # so that we can do global fine optimization, we'll do that here rather than in the
        # subclasses.

        # Do coarse optimization on the singleton Movers.  Record the selected coarse
        # index.
        for s in singletonCliques:
          mover = self._interactionGraph.vertex_label(s[0])
          (bestScore, infoString) = optC.OptimizeSingleMoverCoarse(mover)
          self._infoString += infoString
          ret = bestScore
          self._infoString += _VerboseCheck(self._verbosity, 1,"Singleton optimized with score {:.2f}\n".format(ret))
        self._infoString += _ReportTiming(self._verbosity, "optimize singletons (coarse)")

        # Do coarse optimization on the multi-Mover Cliques.
        for g in groupCliques:
          movers = [self._interactionGraph.vertex_label(i) for i in g]
          subset = _subsetGraph(self._interactionGraph, movers)

          # Find all of the Movers in the clique so that we know which ones to capture the state for.
          movers = [self._interactionGraph.vertex_label(v) for v in subset.vertices()]

          # Turn the graph edges into a versa array that holds the edges. The first index is the
          # number of edge and the second is 2D, listing the index of each mover.
          vertexList = list(subset.vertices())
          edges = flex.int(flex.grid(len(list(subset.edges())), 2))
          for i,e in enumerate(subset.edges()):
            edges[(i,0)] = vertexList.index(subset.source(e))
            edges[(i,1)] = vertexList.index(subset.target(e))

          (bestScore, infoString) = optC.OptimizeCliqueCoarse(movers, edges)
          self._infoString += infoString
          ret = bestScore

          self._infoString += _VerboseCheck(self._verbosity, 1,"Clique optimized with score {:.2f}\n".format(ret))
        self._infoString += _ReportTiming(self._verbosity, "optimize cliques (coarse)")

        # Do fine optimization on the Movers.  This is done independently for
        # each of them, whether they are part of a multi-Mover Clique or not.
        self._infoString += _VerboseCheck(self._verbosity, 1,"Fine optimization on all Movers\n")
        for m in self._movers:
          (bestScore, infoString) = optC.OptimizeSingleMoverFine(m)
          self._infoString += infoString
        self._infoString += _ReportTiming(self._verbosity, "optimize all Movers (fine)")

        ################################################################################
        # Print the final state and score for all Movers
        def _scoreMoverReportClash(self, m, index):
          coarse = m.CoarsePositions()
          score = self._preferenceMagnitude * coarse.preferenceEnergies[index]
          clash = False
          # There may not be as many atoms moved as there are atoms, so use the proper length
          for i in range(len(coarse.positions[index])):
            atom = coarse.atoms[i]
            maxRadiusWithoutProbe = self._extraAtomInfo.getMappingFor(atom).vdwRadius + self._maximumVDWRadius
            res = self._dotScorer.score_dots(atom, self._minOccupancy, self._spatialQuery,
              maxRadiusWithoutProbe, self._probeRadius, self._excludeDict[atom.i_seq],
              optC.GetDots(atom.i_seq), self._probeDensity, False, False)
            score += res.totalScore()
            if res.hasBadBump:
              clash = True
          return score, clash

        def _printPose(self, m):
          description = m.PoseDescription(optC.GetCoarseLocation(m), optC.GetFineLocation(m),
                                          not self._skipBondFixup)

          # If the Mover is a flip of some kind, then the substring "lipped " will be present
          # in the description.
          # When that happens, we check the final state and the other flip
          # state (which is half of the coarse states away) to see if both have clashes or
          # if they are close in energy. If so, then we annotate the output.
          # We add the same number of words to the output string in all cases to make things
          # easier for a program to parse.
          if "lipped " in description:
            coarse = m.CoarsePositions()
            numPositions = len(coarse.positions)
            final = optC.GetCoarseLocation(m)
            other = (final + numPositions//2) % numPositions
            self._setMoverState(coarse, other)
            otherScore, otherBump = _scoreMoverReportClash(self, m, other)
            self._setMoverState(coarse, final)
            finalScore, finalBump = _scoreMoverReportClash(self, m, final)
            if otherBump and finalBump:
              description += " BothClash"
            elif "Unflipped" in description and (
                (otherScore > finalScore) and (otherScore - finalScore <= self._nonFlipPreference)):
              description += " Uncertain"
            else:
              description += " ."
          else:
            description += " ."

          self._infoString += _VerboseCheck(self._verbosity, 1,"  {} final score: {:.2f} pose {}\n".format(
            self._moverInfo[m], optC.GetHighScore(m), description))

        self._infoString += _VerboseCheck(self._verbosity, 1,"BEGIN REPORT: Model "+str(mi)+" Alt '"+alt+"':\n")
        sortedGroups = sorted(groupCliques, key=len, reverse=True)
        for g in sortedGroups:
          self._infoString += _VerboseCheck(self._verbosity, 1," Set of "+str(len(g))+" Movers:")
          movers = [self._interactionGraph.vertex_label(i) for i in g]
          # Parse the record for each mover and pull out the initial score.  Sum the initial and
          # final scores across all Movers in the group and report this.
          initial = 0.0
          final = 0.0
          for m in movers:
            initial += float(self._moverInfo[m].split()[9])
            final += optC.GetHighScore(m)
          self._infoString += _VerboseCheck(self._verbosity, 1," Totals: initial score {:.2f}, final score {:.2f}\n".format(initial, final))
          for m in movers:
            _printPose(self, m)
        self._infoString += _VerboseCheck(self._verbosity, 1," Singleton Movers:\n")
        for s in singletonCliques:
          m = self._interactionGraph.vertex_label(s[0])
          _printPose(self, m)
        self._infoString += _VerboseCheck(self._verbosity, 1,"END REPORT\n")

        ################################################################################
        # Do FixUp on the final coarse orientations.  Set the positions, extra atom info
        # and deletion status for all atoms that have entries for each.
        # This must be done after we print the scores because the print methods move the
        # coarse state to see how much it changed.
        if not self._skipBondFixup:
          self._infoString += _VerboseCheck(self._verbosity, 1,"FixUp on all Movers\n")
          for m in self._movers:
            loc = optC.GetCoarseLocation(m)
            self._infoString += _VerboseCheck(self._verbosity, 3,"FixUp on {} coarse location {}\n".format(
            self._moverInfo[m],loc))
            self._doFixup(m.FixUp(loc))
          self._infoString += _ReportTiming(self._verbosity, "fix up Movers")

        #################################################################################
        # Record the fraction of atoms that were calculated and the fraction that were cached.
        self._numCalculated += optC.GetNumCalculatedAtoms()
        self._numCached += optC.GetNumCachedAtoms()

      ################################################################################
      # Deletion of atoms (Hydrogens) that were requested by Histidine FixUp()s,
      # both in the initial setup and determined during optimization.  Phantom Hydrogens
      # on waters do not need to be adjusted because they were never added to the
      # structure.
      # We only do this after the last-checked alternate configuration for a given model.
      self._infoString += _VerboseCheck(self._verbosity, 1,"Deleting Hydrogens tagged by Histidine Movers\n")
      for a in self._deleteMes:
        aName = a.name.strip().upper()
        resNameAndID = _ResNameAndID(a)
        self._infoString += _VerboseCheck(self._verbosity, 5,"Deleting {} {}\n".format(resNameAndID, aName))
        a.parent().remove_atom(a)
      self._infoString += _ReportTiming(self._verbosity, "delete Hydrogens")

      #################################################################################
      # Dump information about all of the atoms in the model into a string.
      if fillAtomDump:
        self._atomDump = Helpers.writeAtomInfoToString(myModel.atoms(), self._extraAtomInfo)
        self._infoString += _ReportTiming(self._verbosity, "dump atom info to string")

      #################################################################################
      # Report the fraction of atoms that were calculated and the fraction that were cached.
      if self._numCalculated > 0:
        self._infoString += _VerboseCheck(self._verbosity, 1,
            'Calculated : cached atom scores: {} : {}; fraction calculated {:.2f}\n'.format(
            self._numCalculated, self._numCached,
            self._numCalculated/(self._numCalculated + self._numCached)))

  def getInfo(self):
    """
      Returns information that the user may care about regarding the processing.  The level of
      detail on this can be set by setting the verbosity parameter during Optimizer construction.
      :return: the information so far collected in the string.  Calling this method also clears
      the information, so that later calls will not repeat it.
      If the object was constructed with fillAtomDump = False, then the atom dump will always be
      empty.
    """
    ret = self._infoString
    self._infoString = ""
    return ret

  def getWarnings(self):
    """
      Returns warnings that the user may care about regarding the processing.
      :return: the information so far collected in the string.  Calling this method also clears
      the information, so that later calls will not repeat it.
      If the object was constructed with fillAtomDump = False, then the atom dump will always be
      empty.
    """
    ret = self._warningString
    self._warningString = ""
    return ret

  def getAtomDump(self):
    """
      Returns information about the final status of each atom in each model, including VdW radius
      and various flags.  Useful for debugging and for regression testing.
      :return: the information so far collected in the string.  Calling this method also clears
      the information, so that later calls will not repeat it.
    """
    ret = self._atomDump
    self._atomDump = ""
    return ret

  def _setMoverState(self, positionReturn, index):
    """
      Move the atoms to their new positions, updating the spatial query structure
      by removing the old and adding the new location.
    """
    # there may be fewer moved atoms than there are total atoms
    numMoved = len(positionReturn.positions[index])
    for i in range(numMoved):
      a = positionReturn.atoms[i]

      self._spatialQuery.remove(a)
      # Make a slice here so that we get a copy of the location rather than a reference to it
      a.xyz = positionReturn.positions[index][i][:]
      self._spatialQuery.add(a)

    # Update the extra atom information associated with each atom.
    # Note that there may be fewer entries than atoms, but they all correspond.
    for i, e in enumerate(positionReturn.extraInfos[index]):
      self._extraAtomInfo.setMappingFor(positionReturn.atoms[i], e)

    # Manage the deletion status of each atom, including ensuring
    # consistency with the spatial-query structure.
    # Note that there may be fewer entries than atoms, but they all correspond.
    for i, doDelete in enumerate(positionReturn.deleteMes[index]):
      if doDelete:
        self._spatialQuery.remove(positionReturn.atoms[i])
        self._deleteMes.add(positionReturn.atoms[i])
        self._infoString += _VerboseCheck(self._verbosity, 10,"Deleting atom\n")
      else:
        self._spatialQuery.add(positionReturn.atoms[i])
        self._deleteMes.discard(positionReturn.atoms[i])
        self._infoString += _VerboseCheck(self._verbosity, 10,"Ensuring deletable atom is present\n")

  def _doFixup(self, fixUp):
    """
      Move the atoms to their fixup positions, updating their extra atom info
      and deletion status
    """
    myAtoms = fixUp.atoms
    for i, p in enumerate(fixUp.positions):
      self._infoString += _VerboseCheck(self._verbosity, 5,"Moving atom to {}\n".format(p))
      myAtoms[i].xyz = p
    for i, e in enumerate(fixUp.extraInfos):
      self._infoString += _VerboseCheck(self._verbosity, 5,"Atom info to {}\n".format(e))
      self._extraAtomInfo.setMappingFor(myAtoms[i], e)
    for i, d in enumerate(fixUp.deleteMes):
      # Either ensure that it is deleted or ensure that it is not depending on the
      # value of the deletion result.
      self._infoString += _VerboseCheck(self._verbosity, 5,"Atom deleted is {}\n".format(d))
      if d:
        self._deleteMes.add(myAtoms[i])
      else:
        self._deleteMes.discard(myAtoms[i])

  ##################################################################################
  # Placement

  def _PlaceMovers(self, atoms, rotatableHydrogenIDs, bondedNeighborLists, hParameters,
                    addFlipMovers):
    """Produce a list of Movers for atoms in a pdb.hierarchy.conformer that has added Hydrogens.
    :param atoms: flex array of atoms to search.  This must have all Hydrogens needed by the
    Movers present in the structure already.
    :param rotatableHydrogenIDs: List of sequence IDs that include those of single hydrogens
    that are rotatable. These are hydrogens that are the only one bound to their neighbor
    and they are rotatable.
    :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
    structure that the atom from the first parameter interacts with that lists all of the
    bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
    :param hParameters: List indexed by sequence ID that stores the riding
    coefficients for hydrogens that have associated dihedral angles.  This can be
    obtained by calling model.setup_riding_h_manager() and then model.get_riding_h_manager().
    :param addFlipMovers: Do we add flip Movers along with other types?
    :return: List of atoms that should be unconditionally marked for deletion as a result
    of the analysis done during placement.
    """

    # Parse the flipStates parameter to get a list of FlipMoverState objects.
    fs = _ParseFlipStates(self._flipStates)

    # List of Movers
    self._movers = []

    # List of atoms to delete
    deleteAtoms = []

    # Dictionary mapping from Mover to information about the Mover, so that we can keep track of
    # everything from its placement to its final state in the same string.
    self._moverInfo = {}

    # The radius of a Polar Hydrogen (one that forms Hydrogen bonds)
    # @todo Can we look this up from inside CCTBX?
    polarHydrogenRadius = 1.05

    # For each atom, check to see if it is the main atom from a Mover.  If so, attempt to
    # construct the appropriate Mover.  The construction may fail because not all required
    # atoms are present (especially the Hydrogens).  If the construction succeeds, add the
    # Mover to the list of those to return.  When they fail, add the output to the information
    # string. We check single-hydrogen rotators later.

    for a in atoms:
      # Find the stripped upper-case atom and residue names and the residue ID,
      # and an identifier string
      aName = a.name.strip().upper()
      resName = a.parent().resname.strip().upper()
      resNameAndID = _ResNameAndID(a)

      # See if we should construct a MoverNH3Rotator here.
      # Find any Nitrogen that has four total bonded neighbors, three of which are Hydrogens.
      # @todo Is this a valid way to search for them?
      if a.element == 'N' and len(bondedNeighborLists[a]) == 4:
        numH = 0
        for n in bondedNeighborLists[a]:
          if n.element_is_hydrogen():
            numH += 1
        if numH == 3:
          try:
            self._movers.append(Movers.MoverNH3Rotator(a, bondedNeighborLists, hParameters))
            self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverNH3Rotator "+str(len(self._movers))+" to "+resNameAndID+"\n")
            self._moverInfo[self._movers[-1]] = "NH3Rotator at "+resNameAndID+" "+aName;
          except Exception as e:
            self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverNH3Rotator to "+resNameAndID+": "+str(e)+"\n")

      # See if we should construct a MoverAromaticMethylRotator or MoverTetrahedralMethylRotator here.
      # Find any Carbon that has four total bonded neighbors, three of which are Hydrogens.
      # @todo Is this a valid way to search for them?
      if a.element == 'C' and len(bondedNeighborLists[a]) == 4:
        numH = 0
        neighbor = None
        for n in bondedNeighborLists[a]:
          if n.element_is_hydrogen():
            numH += 1
          else:
            neighbor = n
        if numH == 3:
          # See if the Carbon's other neighbor is attached to two other atoms (3 total).  If so,
          # then insert a MoverAromaticMethylRotator and if not, generate a MoverTetrahedralMethylRotator
          # so that the Hydrogens will be staggered but do not add it to those to be optimized.
          if len(bondedNeighborLists[neighbor]) == 3:
            try:
              self._movers.append(Movers.MoverAromaticMethylRotator(a, bondedNeighborLists, hParameters))
              self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverAromaticMethylRotator "+str(len(self._movers))+" to "+resNameAndID+" "+aName+"\n")
              self._moverInfo[self._movers[-1]] = "AromaticMethylRotator at "+resNameAndID+" "+aName;
            except Exception as e:
              self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverAromaticMethylRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")
          else:
            try:
              ignored = Movers.MoverTetrahedralMethylRotator(a, bondedNeighborLists, hParameters)
              self._infoString += _VerboseCheck(self._verbosity, 1,"Used MoverTetrahedralMethylRotator to stagger "+resNameAndID+" "+aName+"\n")
            except Exception as e:
              self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverTetrahedralMethylRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")

      # See if we should insert a MoverAmideFlip here.
      # @todo Is there a more general way than looking for specific names?
      if addFlipMovers and ((aName == 'XD2' and resName == 'ASX') or (aName == 'XE2' and resName == 'GLX')):
        self._infoString += _VerboseCheck(self._verbosity, 1,"Not attempting to adjust "+resNameAndID+" "+aName+"\n")
      if addFlipMovers and ((aName == 'ND2' and resName == 'ASN') or (aName == 'NE2' and resName == 'GLN')):
        # Find the Oxygen and see if it is within a range of ideal bonding distances to a positive ion.
        # Do this in two steps; find the Carbon bonded to the Nitrogen and then the Oxygen bonded to the
        # Carbon.
        foundIon = False
        oxygen = None
        for b in bondedNeighborLists[a]:
          if b.element.upper() == 'C':
            for b2 in bondedNeighborLists[b]:
              if b2.element.upper() == 'O':
                oxygen = b2
        # If we have a close-enough ion, we skip adding a Mover.
        if oxygen is not None:
          # @todo Check both flips and lock down the appropriate one rather than only checking the first.
          myRad = self._extraAtomInfo.getMappingFor(a).vdwRadius
          minDist = myRad
          maxDist = 0.25 + myRad + self._maximumVDWRadius
          neighbors = self._spatialQuery.neighbors(oxygen.xyz, minDist, maxDist)
          for n in neighbors:
            if n.element_is_positive_ion():
              dist = (Helpers.rvec3(oxygen.xyz) - Helpers.rvec3(n.xyz)).length()
              expected = myRad + self._extraAtomInfo.getMappingFor(n).vdwRadius
              self._infoString += _VerboseCheck(self._verbosity, 5,'Checking AmideFlip against '+n.name.strip()+' at '+str(n.xyz)+' from '+str(oxygen.xyz)+
                ' dist = '+str(dist)+', expected = '+str(expected)+'; N rad = '+str(myRad)+
                ', '+n.name.strip()+' rad = '+str(self._extraAtomInfo.getMappingFor(n).vdwRadius)+'\n')
              # @todo Why are we using -0.65 here and -0.55 for Histidine?
              if dist >= (expected - 0.65) and dist <= (expected + 0.25):
                foundIon = True

        if not foundIon:
          try:
            # Check to see if the state of this Mover has been specified. If so, place it in
            # the requested state and don't insert the Mover. If not, then insert the Mover.
            s = _FindFlipState(a, fs)
            if s is not None:
              self._infoString += _VerboseCheck(self._verbosity, 1,"Setting MoverAmideFlip for "+resNameAndID+": flipped = "+
                str(s.flipped)+", angles adjusted = "+str(s.fixedUp)+"\n")
              flip = Movers.MoverAmideFlip(a, "CA", bondedNeighborLists, self._nonFlipPreference)
              cp = flip.CoarsePositions()
              if s.flipped:
                index = 1
              else:
                index = 0
              self._setMoverState(cp, index)
              if s.fixedUp:
                self._doFixup(flip.FixUp(index))
            else:
              self._movers.append(Movers.MoverAmideFlip(a, "CA", bondedNeighborLists, self._nonFlipPreference))
              self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverAmideFlip "+str(len(self._movers))+" to "+resNameAndID+"\n")
              self._moverInfo[self._movers[-1]] = "AmideFlip at "+resNameAndID+" "+aName
          except Exception as e:
            self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverAmideFlip to "+resNameAndID+": "+str(e)+"\n")

      # See if we should insert a MoverHisFlip here.
      # @todo Is there a more general way than looking for specific names?
      if aName == 'NE2' and resName == 'HIS':
        try:
          # Get a potential Mover and test both of its Nitrogens in the original and flipped
          # locations.  If one or both of them are near enough to be ionically bound to an
          # ion, then we remove the Hydrogen(s) and lock the Histidine at that orientation
          # rather than inserting the Mover into the list of those to be optimized.
          # @todo Consider checking both configurations to see if either one has two bonds.
          hist = Movers.MoverHisFlip(a, bondedNeighborLists, self._extraAtomInfo, self._nonFlipPreference)

          # Find the four positions to check for Nitrogen ionic bonds
          # The two atoms are NE2 (0th atom with its Hydrogen at atom 1) and
          # ND1 (4th atom with its Hydrogen at atom 5).
          cp = hist.CoarsePositions()
          ne2Orig = cp.positions[0][0]
          nd1Orig = cp.positions[0][4]
          ne2Flip = cp.positions[4][0]
          nd1Flip = cp.positions[4][4]

          # See if any are close enough to a positive ion to be ionically bonded.
          # If any are, record whether it is the original or
          # the flipped configuration.  Check the original configuration first.
          # Check out to the furthest distance of any atom's VdW radius.
          myRad = self._extraAtomInfo.getMappingFor(a).vdwRadius
          minDist = myRad
          maxDist = 0.25 + myRad + self._maximumVDWRadius
          bondedConfig = None
          for i,pos in enumerate([ne2Orig, nd1Orig, ne2Flip, nd1Flip]):
            neighbors = self._spatialQuery.neighbors(pos, minDist, maxDist)
            for n in neighbors:
              if n.element_is_positive_ion():
                dist = (Helpers.rvec3(pos) - Helpers.rvec3(n.xyz)).length()
                expected = myRad + self._extraAtomInfo.getMappingFor(n).vdwRadius
                self._infoString += _VerboseCheck(self._verbosity, 5,'Checking Histidine '+str(i)+' against '+n.name.strip()+' at '+str(n.xyz)+' from '+str(pos)+
                  ' dist = '+str(dist)+', expected = '+str(expected)+'; N rad = '+str(myRad)+
                  ', '+n.name.strip()+' rad = '+str(self._extraAtomInfo.getMappingFor(n).vdwRadius)+'\n')
                if dist >= (expected - 0.55) and dist <= (expected + 0.25):
                  # The first two elements in the enumeration come from Histidine configuration 0
                  # and the second two from configuration 4, so we figure out which one we should
                  # be in.
                  bondedConfig = (i // 2) * 4
                  self._infoString += _VerboseCheck(self._verbosity, 5,'Found ionic bond in coarse configuration '+str(bondedConfig)+'\n')
                  break
            if bondedConfig is not None:
              # We want the first configuration that is found, so we don't flip if we don't need to.
              break

          # If one of the bonded configurations has at least one Ionic bond, then check each of
          # the Nitrogens in that configuration, removing its Hydrogen if it is bonded to an ion.
          # Set the histidine to that flip state; it will not be inserted as a Mover.
          if bondedConfig is not None:
            # Set the histidine in that flip state
            fixUp = hist.FixUp(bondedConfig)
            coarsePositions = hist.CoarsePositions().positions[bondedConfig]
            for i,a in enumerate(fixUp.atoms):
              if i < len(fixUp.positions):
                a.xyz = fixUp.positions[i]
              if i < len(fixUp.extraInfos):
                self._extraAtomInfo.setMappingFor(a, fixUp.extraInfos[i])

            # See if we should remove the Hydrogen from each of the two potentially-bonded
            # Nitrogens and make each an acceptor if we do remove its Hydrogen.  The two atoms
            # are NE2 (0th atom with its Hydrogen at atom 1) and ND1 (4th atom with
            # its Hydrogen at atom 5).
            # NOTE that we must check the position of the atom in its coarse configuration rather
            # than its FixUp configuration because that is the location that we use to determine if
            # we have an ionic bond.
            def _modifyIfNeeded(nitro, coarseNitroPos, hydro):
              # Helper function to check and change things for one of the Nitrogens.
              result = ""
              myRad = self._extraAtomInfo.getMappingFor(nitro).vdwRadius
              minDist = myRad
              maxDist = 0.25 + myRad + self._maximumVDWRadius
              neighbors = self._spatialQuery.neighbors(coarseNitroPos, minDist, maxDist)
              for n in neighbors:
                if n.element_is_positive_ion():
                  dist = (Helpers.rvec3(coarseNitroPos) - Helpers.rvec3(n.xyz)).length()
                  expected = myRad + self._extraAtomInfo.getMappingFor(n).vdwRadius
                  if dist >= (expected - 0.55) and dist <= (expected + 0.25):
                    result += _VerboseCheck(self._verbosity, 1,'Not adding Hydrogen to '+resNameAndID+nitro.name+' and marking as an acceptor '+
                      '(ionic bond to '+n.name.strip()+')\n')
                    extra = self._extraAtomInfo.getMappingFor(nitro)
                    extra.isAcceptor = True
                    self._extraAtomInfo.setMappingFor(nitro, extra)
                    deleteAtoms.append(hydro)
                    break
              return result

            self._infoString += _modifyIfNeeded(fixUp.atoms[0], coarsePositions[0], fixUp.atoms[1])
            self._infoString += _modifyIfNeeded(fixUp.atoms[4], coarsePositions[4], fixUp.atoms[5])

            self._infoString += _VerboseCheck(self._verbosity, 1,"Set MoverHisFlip on "+resNameAndID+" to state "+str(bondedConfig)+"\n")
          elif addFlipMovers: # Add a Histidine flip Mover if we're adding flip Movers
            # Check to see if the state of this Mover has been specified. If so, place it in
            # the requested state and insert the Mover as a non-flipping Histidine.
            # If not, then insert the Histidine flip Mover.
            s = _FindFlipState(a, fs)
            if s is not None:
              if s.flipped:
                enabledFlips = 2
              else:
                enabledFlips = 1
              hist = Movers.MoverHisFlip(a, bondedNeighborLists, self._extraAtomInfo, self._nonFlipPreference,
                enabledFlips, s.fixedUp)
              self._movers.append(hist)
              self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverHisPlace "+str(len(self._movers))+" to "+resNameAndID+"\n")
              # Name is HisPlace rather than HisFlip, so we don't try and report uncertain or clash for it.
              self._moverInfo[self._movers[-1]] = "HisPlace at "+resNameAndID+" "+aName;
              # Set the position to the first coarse position, which will cause a flip if required
              coarse = hist.CoarsePositions()
              self._setMoverState(coarse, 0)
            else:
              self._movers.append(hist)
              self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverHisFlip "+str(len(self._movers))+" to "+resNameAndID+"\n")
              self._moverInfo[self._movers[-1]] = "HisFlip at "+resNameAndID+" "+aName;
        except Exception as e:
          self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverHisFlip to "+resNameAndID+": "+str(e)+"\n")

    # Make a dictionary looked up up by i_seq that returns the relevant atom. We use this to
    # look up the single-hydrogen rotators. We only place atoms that are in our current model
    # index and alternate.
    atomDict = {}
    for a in atoms:
      atomDict[a.i_seq] = a

    # Insert all single-hydrogen rotators, which we have a list of. Do this after we've
    # done all of the other Movers so that any lock-down placements will have been done
    # before we check.
    for i in rotatableHydrogenIDs:
      # Find out the associated atom
      try:
        a = atomDict[i]
      except KeyError:
        # This atom is not in our current model index or alternate, so skip it.
        continue

      aName = a.name.strip().upper()
      resName = a.parent().resname.strip().upper()
      resNameAndID = _ResNameAndID(a)
      try:
        # Skip Hydrogens that are not bound to an atom that only has a single other
        # atom bound to it -- we will handle those in other cases.
        neighbor = bondedNeighborLists[a][0]
        if len(bondedNeighborLists[neighbor]) == 2:

          # Get the list of atoms that are bonded in a chain to the hydrogen. They
          # cannot be either hydrogen-bond acceptors nor touches/clashes with this atom.
          bonded = Helpers.getAtomsWithinNBonds(a, bondedNeighborLists, self._extraAtomInfo,
                                                self._probeRadius, self._bondedNeighborDepth, 3)

          # Construct a list of nearby atoms that are potential acceptors and potential touches.
          # The potential acceptors are also potential touches.
          potentialAcceptors = []
          potentialTouches = []

          # Get the list of nearby atoms.  The center of the search is the atom that
          # the Hydrogen is bound to and its radius is 4 (these values are pulled from
          # the Reduce C++ code). We skip ones whose radius is inside the neighbor so
          # that we don't count the hydrogen itself.
          maxDist = 4.0
          nearby = self._spatialQuery.neighbors(neighbor.xyz,
                                                self._extraAtomInfo.getMappingFor(neighbor).vdwRadius,
                                                maxDist)

          # Check each nearby atom to see if its distance from the neighbor is within
          # the sum of the hydrogen-bond length of the neighbor atom, the radius of
          # a polar Hydrogen, and the radius of the nearby atom, indicating potential
          # overlap.
          # O-H & N-H bond len == 1.0, S-H bond len == 1.3
          XHbondlen = 1.0
          if neighbor.element == "S":
            XHbondlen = 1.3
          candidates = []
          for n in nearby:
            # We don't count bonded atoms.
            if not n in bonded:
              # We treat them all as potential touches/clashes.
              potentialTouches.append(n)
              d = (Helpers.rvec3(neighbor.xyz) - Helpers.rvec3(n.xyz)).length()
              if d <= XHbondlen + self._extraAtomInfo.getMappingFor(n).vdwRadius + polarHydrogenRadius:
                candidates.append(n)

          # See if each candidate atom is a potential acceptor or a flip partner from
          # Histidine or NH2 flips (which may be moved into that atom's position during a flip).
          # We check the partner (N's for GLN And ASN, C's for HIS) because if it is the O or
          # N atom, we will already be checking it as an acceptor now.
          # In any case, it is a potential touch.
          for c in candidates:
            cName = c.name.strip().upper()
            resName = c.parent().resname.strip().upper()
            flipPartner = (
              (cName == 'ND2' and resName == 'ASN') or
              (cName == 'NE2' and resName == 'GLN') or
              (cName == 'CE1' and resName == 'HIS') or (cName == 'CD2' and resName == 'HIS') )
            acceptor = self._extraAtomInfo.getMappingFor(c).isAcceptor
            if acceptor or flipPartner:
              potentialAcceptors.append(c)

          self._movers.append(Movers.MoverSingleHydrogenRotator(a, bondedNeighborLists, self._extraAtomInfo,
                                                                hParameters, potentialAcceptors,
                                                                potentialTouches))
          self._infoString += _VerboseCheck(self._verbosity, 1,"Added MoverSingleHydrogenRotator "+str(len(self._movers))+" to "+resNameAndID+" "+aName+
            " with "+str(len(potentialAcceptors))+" potential nearby acceptors\n")
          self._moverInfo[self._movers[-1]] = "SingleHydrogenRotator at "+resNameAndID+" "+aName;
      except Exception as e:
        self._infoString += _VerboseCheck(self._verbosity, 0,"Could not add MoverSingleHydrogenRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")

    return deleteAtoms


  def _DescribeMover(self, m):
    type_str = str(type(m))
    type_str = type_str.split("'")[1]  # Extract the type name between the quotes
    type_str = type_str.split(".")[-1]  # Extract the class name from the module path
    info = _ResNameAndID(m.CoarsePositions().atoms[0])
    return type_str + ' ' + info

  def _DescribeMoverAtom(self, m, a, rad, loc):
    return ( ' {' + self._DescribeMover(m) + " " + a.name.strip().upper() + '}' +
            ' r= '+str(rad)+' {:0.3f} {:0.3f} {:0.3f}\n'.format(loc[0], loc[1], loc[2]))

  def _MoverSpheres(self, m, extraRadius=0.0):
    """Returns a Kinemage string listing all of the positions of the movable atoms
    in the specified mover over all possible locations. The radius of each is as
    specified in the atom plus the specified extraRadius.
    """
    ret = ''
    for i, cp in enumerate(m.CoarsePositions().positions):
      for j, loc in enumerate(cp):
        a = m.CoarsePositions().atoms[j]
        rad = str(extraRadius + self._extraAtomInfo.getMappingFor(a).vdwRadius)
        ret += self._DescribeMoverAtom(m, a, rad, loc)
      for fp in m.FinePositions(i).positions:
        for j, loc in enumerate(fp):
          a = m.CoarsePositions().atoms[j]
          rad = str(extraRadius + self._extraAtomInfo.getMappingFor(a).vdwRadius)
          ret += self._DescribeMoverAtom(m, a, rad, loc)
    return ret


  def _InteractionKinemage(self, components):
      # Used to round-robin among a set of differentiable colors
      COLORS = ['gray', 'pink', 'sea', 'sky', 'cyan', 'magenta', 'yellow']
      curColor = 0

      # Write all of the Movers in all of the cliques, with a different color per Mover
      # using the actual atom radius.
      ret = '@kinemage 1\n'
      ret += '@master {movers}\n'
      for c in components:
        movers = [self._interactionGraph.vertex_label(i) for i in c]
        for m in movers:
          color = COLORS[curColor]
          curColor = (curColor + 1) % len(COLORS)
          # Radius will be overridden per atom
          ret += '@spherelist color= '+color+' nobutton master= {movers}\n'
          ret += self._MoverSpheres(m)

      # Write all of the cliques, with a different color per clique
      # using the atom radius plus the probe radius.
      for i, c in enumerate(components):
        cliqueName = 'clique '+str(i)+' size '+str(len(c))
        ret += '@master {'+cliqueName+'}\n'
        color = COLORS[curColor]
        curColor = (curColor + 1) % len(COLORS)
        movers = [self._interactionGraph.vertex_label(i) for i in c]
        for m in movers:
          # Radius will be overridden per atom
          ret += '@spherelist color= '+color+' nobutton master= {'+cliqueName+'}\n'
          ret += self._MoverSpheres(m, self._probeRadius)

      return ret


##################################################################################
# Private helper functions.

def _ParseFlipStates(fs):
  """Produce a list of FlipMoverState objects by parsing a comma-separated flipStates
  string.
  """
  ret = []
  if fs is not None:
    for state in fs.split(','):
      words = state.split()
      if len(words) == 6 or len(words) == 7:
        modelId = int(words[0]) + 1
        altId = words[1].lower()
        if altId == '.':
          altId = ''
        chain = words[2].upper()
        resName = words[3].upper()
        if resName == 'HIS':
          t = 'HisFlip'
        else:
          t = 'AmideFlip'
        resIdWithICode = words[4]
        flipped = words[5] == 'Flipped'
        fixedUp = False
        if len(words) == 7:
          fixedUp = words[6] == 'AnglesAdjusted'
        ret.append(FlipMoverState(t, modelId, altId, chain, resName, resIdWithICode, flipped, fixedUp))
  return ret


def _FindFlipState(a, flipStates):
  '''Is an atom in the list of residues that have flip states? If so, return the associated state
  :param a: atom
  :param flipStates: List of FlipMoverState objects
  :return: FlipMoverState if the residue is in the list, None if not
  '''
  ag = a.parent()
  rg = ag.parent()
  chain = rg.parent()
  for fs in flipStates:
    modelId = chain.parent().id
    # We must offset the index of the model by 1 to get to the 1-based model ID
    if ( (modelId == fs.modelId + 1 or modelId == '') and
          (chain.id == fs.chain) and
          (rg.resseq_as_int() == fs.resId) and
          (ag.altloc.strip() == '' or fs.altId == '.' or ag.altloc.strip().lower() == fs.altId.lower()) and
          (rg.icode.strip() == fs.iCode.strip())
        ):
      return fs
  return None


def _subsetGraph(g, keepLabels):
  """
  Return a new graph that is a subset of g that includes only
  the vertices whose labels are listed, and edges between these.
  :param g: Graph to be subsetted.
  :param keepLabels: Labels on the vertices that should end up in the resulting
  graph.  Must be a subset of the vertices in the graph, but it can include all
  of the vertices, in which case it will make a full copy of the graph.
  :return: Boost graph that is a subset of the original graph.
  """

  # We keep a dictionary from labels to the vertex in the new graph that points to
  # that label so that we can construct edges in the new graph.
  vertexForLabel = {}

  # Construct a subgraph that consists only of vertices to be kept and
  # edges both of whose ends are on these vertices.
  ret = graph.adjacency_list(
        vertex_type = "list",   # List so that deletions do not invalidate iterators and descriptors
        )
  for v in g.vertices():
    label = g.vertex_label(v)
    if label in keepLabels:
      vertexForLabel[label] = ret.add_vertex(label)
  for e in g.edges():
    sourceLabel = g.vertex_label( g.source(e) )
    targetLabel = g.vertex_label( g.target(e) )
    if sourceLabel in keepLabels and targetLabel in keepLabels:
      ret.add_edge(vertex1 = vertexForLabel[sourceLabel], vertex2 = vertexForLabel[targetLabel])

  return ret


##################################################################################
# Test function and associated data and helpers to verify that all functions behave properly.

# Class to pass default Probe parameters as if they were in a probePhil structure
class _philLike:
  def __init__(self):
    self.probe_radius = 0.25
    self.density = 16.0
    self.worse_clash_cutoff = 0.5
    self.clash_cutoff = 0.4
    self.contact_cutoff = 0.25
    self.uncharged_hydrogen_cutoff = 0.6
    self.charged_hydrogen_cutoff = 0.8
    self.bump_weight = 10.0
    self.hydrogen_bond_weight = 4.0
    self.gap_weight = 0.25
    self.allow_weak_hydrogen_bonds = False
    self.ignore_ion_interactions = False
    self.set_polar_hydrogen_radius = True
    self.excluded_bond_chain_length = 4

def _optimizeFragment(pdb_raw, bondedNeighborDepth = 4):
  """Returns an optimizer constructed based on the raw PDB data snippet passed in.
  :param pdb_raw: A string that includes a snippet of a PDB file. Should have no alternates.
  :param bondedNeighborDepth: How many levels of bonded neighbors to include in the
  neighbor lists for each atom. Some tests were done with a value of 3, which differs
  from the newer default of 4.
  :return: Optimizer initialized based on the snippet.
  """
  dm = DataManager(['model'])
  dm.process_model_str("my_data.pdb", pdb_raw)
  model = dm.get_model()

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check=True
  model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints

  # Optimization will place the movers.
  probePhil = _philLike()
  return Optimizer(probePhil, True, model, bondedNeighborDepth=bondedNeighborDepth)

def Test(inFileName = None, dumpAtoms = False):

  """Test function for all functions provided above.
  :param inFileName: Name of a PDB or CIF file to load (default makes a small molecule)
  :return: Empty string on success, string describing the problem on failure.
  """

  #========================================================================
  # Test the imported C++ functions
  ret = Optimizers_test()
  if len(ret) > 0:
    return "Optimizers.Test(): " + ret

  #========================================================================
  # Test model to use for validating the alternate/conformer-selection functions.

  alternates_test = (
"""
ATOM      1  N   HIS A  1       26.965  32.911   7.593  1.00  7.19           N
ATOM      2  N  AHIS A  2       26.965  32.911   7.593  1.00  7.19           N
ATOM      3  N  BHIS A  2       26.965  32.911   7.593  1.00  7.19           N
ATOM      4  N  CHIS A  2       26.965  32.911   7.593  1.00  7.19           N
END
"""
    )

  #========================================================================
  # Test the AlternatesInModel() function.
  # Test the GetAtomsForConformer() function.
  dm = DataManager(['model'])
  dm.process_model_str("alternates_test.pdb",alternates_test)
  model = dm.get_model()
  model.process(make_restraints=False)
  model = dm.get_model().get_hierarchy().models()[0]
  alts = AlternatesInModel(model)
  if alts != set(['','A','B','C']):
    return "Optimizers.Test(): Incorrect set of alternates for AlternatesInModel() test: " + str(alts)

  for a in alts:
    count = len(GetAtomsForConformer(model, a))
    if count != 2:
      return "Optimizers.Test(): Incorrect atom count for GetAtomsForConformer() test: " + str(count)

  #========================================================================
  # Unit tests for each type of Optimizer.

  ################################################################################
  # Test using a snippet from 7c31 to make sure the single-hydrogen rotator sets
  # the orientation of the Hydrogen to make good contact with a nearby Oxygen. The B
  # alternate has been removed and the A moved to the main alternate. This is not
  # close enough for a good Hydrogen bond, so it does not point straight at the Oxygen.

  pdb_7c31_two_residues = """\
CRYST1   27.854   27.854   99.605  90.00  90.00  90.00 P 43          8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.035901  0.000000  0.000000        0.00000
SCALE2      0.000000  0.035901  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010040        0.00000
ATOM     68  N   SER A   5     -31.155  49.742   0.887  1.00 10.02           N
ATOM     69  CA  SER A   5     -32.274  48.937   0.401  1.00  9.76           C
ATOM     70  C   SER A   5     -33.140  49.851  -0.454  1.00  9.47           C
ATOM     71  O   SER A   5     -33.502  50.939  -0.012  1.00 10.76           O
ATOM     72  CB  SER A   5     -33.086  48.441   1.599  1.00 12.34           C
ATOM     73  OG  SER A   5     -34.118  47.569   1.179  1.00 19.50           O
ATOM    758  N   VAL B   2     -34.430  42.959   3.043  1.00 19.95           N
ATOM    759  CA  VAL B   2     -32.994  42.754   2.904  0.50 18.09           C
ATOM    760  C   VAL B   2     -32.311  44.083   2.629  1.00 17.65           C
ATOM    761  O   VAL B   2     -32.869  44.976   1.975  1.00 18.97           O
ATOM    762  CB  VAL B   2     -32.703  41.738   1.778  0.50 19.70           C
ATOM    763  CG1 VAL B   2     -33.098  42.307   0.419  0.50 19.61           C
ATOM    764  CG2 VAL B   2     -31.224  41.334   1.755  0.50 16.49           C
TER    1447      CYS B  47
END
"""
  opt = _optimizeFragment(pdb_7c31_two_residues)
  movers = opt._movers
  if len(movers) != 3:
    return "Optimizers.Test(): Incorrect number of Movers for single-hydrogen rotator test: " + str(len(movers))

  # See what the pose angle is on the Mover. It should be 62 degrees, and is reported
  # after 'pose Angle '.
  # NOTE: This test was generated when the SingleHydrogenRotator was testing all coarse
  # angles in addition to those towards acceptors. When we switched to the original Reduce
  # behavior of only doing a single non-clashing angle along with the acceptor angles, this
  # test no longer passes. Leaving it in for now in case we want to switch back to the
  # original behavior.
  #angle = int(re.search(r'(?<=pose Angle )[-+]?\d+', opt.getInfo()).group(0))
  #if angle != 62:
  #  return "Optimizers.Test(): Unexpected angle ("+str(angle)+") for single-hydrogen rotator, expected 62"

  ################################################################################
  # Test using a modified snippet from 7c31 to make sure the single-hydrogen rotator sets
  # the orientation of the Hydrogen to make a good hydrogen bond with a nearby Oxygen.
  # The above example has been stripped down and has the Oxygen moved closer.

  pdb_7c31_two_residues_close_O = """\
CRYST1   27.854   27.854   99.605  90.00  90.00  90.00 P 43          8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.035901  0.000000  0.000000        0.00000
SCALE2      0.000000  0.035901  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010040        0.00000
ATOM     68  N   SER A   5     -31.155  49.742   0.887  1.00 10.02           N
ATOM     69  CA  SER A   5     -32.274  48.937   0.401  1.00  9.76           C
ATOM     70  C   SER A   5     -33.140  49.851  -0.454  1.00  9.47           C
ATOM     71  O   SER A   5     -33.502  50.939  -0.012  1.00 10.76           O
ATOM     72  CB  SER A   5     -33.086  48.441   1.599  1.00 12.34           C
ATOM     73  OG  SER A   5     -34.118  47.569   1.179  1.00 19.50           O
ATOM    761  O   VAL B   2     -33.462  45.185   1.977  1.00 18.97           O
TER    1447      CYS B  47
END
"""
  opt = _optimizeFragment(pdb_7c31_two_residues_close_O)
  movers = opt._movers
  if len(movers) != 2:
    return "Optimizers.Test(): Incorrect number of Movers for single-hydrogen rotator H-bond test: " + str(len(movers))

  # See what the pose angle is on the Mover. It is reported after 'pose Angle '.
  expected = 113
  angle = int(re.findall(r'(?<=pose Angle )[-+]?\d+', opt.getInfo())[1])
  if angle != expected:
    return "Optimizers.Test(): Unexpected angle ("+str(angle)+") for single-hydrogen rotator H-bond, expected "+str(expected)+", found "+str(angle)

  ################################################################################
  # Test using a modified snippet from 7c31 to make sure the single-hydrogen rotator sets
  # the orientation of the Hydrogen touching a nearby Carbon.
  # The above example has the Oxygen replaced by a Carbon.

  pdb_7c31_two_residues_close_C = """\
CRYST1   27.854   27.854   99.605  90.00  90.00  90.00 P 43          8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.035901  0.000000  0.000000        0.00000
SCALE2      0.000000  0.035901  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010040        0.00000
ATOM     68  N   SER A   5     -31.155  49.742   0.887  1.00 10.02           N
ATOM     69  CA  SER A   5     -32.274  48.937   0.401  1.00  9.76           C
ATOM     70  C   SER A   5     -33.140  49.851  -0.454  1.00  9.47           C
ATOM     71  O   SER A   5     -33.502  50.939  -0.012  1.00 10.76           O
ATOM     72  CB  SER A   5     -33.086  48.441   1.599  1.00 12.34           C
ATOM     73  OG  SER A   5     -34.118  47.569   1.179  1.00 19.50           O
ATOM    761  C   VAL B   2     -33.462  45.185   1.977  1.00 18.97           C
TER    1447      CYS B  47
END
"""
  opt = _optimizeFragment(pdb_7c31_two_residues_close_C)
  movers = opt._movers
  if len(movers) != 2:
    return "Optimizers.Test(): Incorrect number of Movers for single-hydrogen rotator clash test: " + str(len(movers))

  # See what the pose angle is on the Mover. It is reported after 'pose Angle '.
  expected = -64
  angle = int(re.findall(r'(?<=pose Angle )[-+]?\d+', opt.getInfo())[1])
  if angle != expected:
    return "Optimizers.Test(): Unexpected angle ("+str(angle)+") for single-hydrogen rotator clash, expected "+str(expected)+", found "+str(angle)

  ################################################################################
  # Test using a modified snippet from 1ubq to make sure the NH3 rotator sets
  # the orientation as expected to align with a nearby Oxygen and waters.

  pdb_1ubq_two_residues_close_O = """\
CRYST1   50.840   42.770   28.950  90.00  90.00  90.00 P 21 21 21    4
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.019670  0.000000  0.000000        0.00000
SCALE2      0.000000  0.023381  0.000000        0.00000
SCALE3      0.000000  0.000000  0.034542        0.00000
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N
ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C
ATOM      3  C   MET A   1      26.913  26.639   3.531  1.00  9.62           C
ATOM      4  O   MET A   1      27.886  26.463   4.263  1.00  9.62           O
ATOM      5  CB  MET A   1      25.112  24.880   3.649  1.00 13.77           C
ATOM      6  CG  MET A   1      25.353  24.860   5.134  1.00 16.29           C
ATOM      7  SD  MET A   1      23.930  23.959   5.904  1.00 17.17           S
ATOM      8  CE  MET A   1      24.447  23.984   7.620  1.00 16.11           C
ATOM    127  N   VAL A  17      30.310  25.458   5.384  1.00  8.99           N
ATOM    128  CA  VAL A  17      30.288  24.245   6.193  1.00  8.85           C
ATOM    129  C   VAL A  17      29.279  23.227   5.641  1.00  8.04           C
ATOM    130  O   VAL A  17      28.478  23.522   4.725  1.00  8.99           O
ATOM    131  CB  VAL A  17      29.903  24.590   7.665  1.00  9.78           C
ATOM    132  CG1 VAL A  17      30.862  25.496   8.389  1.00 12.05           C
ATOM    133  CG2 VAL A  17      28.476  25.135   7.705  1.00 10.54           C
TER     603      GLY A  76
HETATM  625  O   HOH A  98      25.928  21.774   2.325  1.00 13.70           O
HETATM  637  O   HOH A 110      28.824  25.094   0.886  0.77 36.99           O
END
"""
  # This test is done with a bonded-neighbor depth of 3 because that is how we were
  # doing the calculations when we set it up.
  opt = _optimizeFragment(pdb_1ubq_two_residues_close_O, 3)
  movers = opt._movers
  if len(movers) != 2:
    return "Optimizers.Test(): Incorrect number of Movers for NH3 rotator test: " + str(len(movers))

  # See what the pose angle is on the Mover. It should be 163 degrees, and is reported
  # after 'pose Angle '.
  angle = int(re.findall(r'(?<=pose Angle )[-+]?\d+', opt.getInfo())[0])
  if angle != 163:
    return "Optimizers.Test(): Unexpected angle ("+str(angle)+") for NH3 rotator, expected 163"

  ################################################################################
  # Test using snippet from 1dfu to ensure that the Histidine optimzization code will
  # properly skip hydrogens that are going to be deleted when scoring.
  pdb_1dfu_his80 = (
"""
CRYST1   75.600   76.600   95.100  90.00  90.00  90.00 C 2 2 21      8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.013257  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013063  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010515        0.00000
ATOM   1443  N   HIS P  80      27.889  23.654 -11.281  1.00 41.74           N
ATOM   1444  CA  HIS P  80      28.738  24.315 -10.299  1.00 45.36           C
ATOM   1445  C   HIS P  80      29.467  25.487 -10.941  1.00 46.80           C
ATOM   1446  O   HIS P  80      29.927  25.397 -12.080  1.00 47.32           O
ATOM   1447  CB  HIS P  80      29.765  23.342  -9.714  1.00 47.26           C
ATOM   1448  CG  HIS P  80      30.457  23.865  -8.492  1.00 48.91           C
ATOM   1449  ND1 HIS P  80      29.947  23.705  -7.221  1.00 49.73           N
ATOM   1450  CD2 HIS P  80      31.593  24.589  -8.351  1.00 49.34           C
ATOM   1451  CE1 HIS P  80      30.739  24.306  -6.351  1.00 49.69           C
ATOM   1452  NE2 HIS P  80      31.745  24.850  -7.011  1.00 50.50           N
ATOM   1453  N   PRO P  81      29.581  26.606 -10.211  1.00 48.50           N
ATOM   1454  CA  PRO P  81      30.255  27.810 -10.698  1.00 49.93           C
ATOM   1455  C   PRO P  81      31.662  27.557 -11.236  1.00 51.19           C
ATOM   1456  O   PRO P  81      32.041  28.115 -12.265  1.00 51.83           O
ATOM   1457  CB  PRO P  81      30.264  28.714  -9.468  1.00 50.10           C
ATOM   1458  CG  PRO P  81      29.003  28.336  -8.774  1.00 49.68           C
ATOM   1459  CD  PRO P  81      29.026  26.831  -8.864  1.00 49.20           C
ATOM   1460  N   TYR P  82      32.432  26.716 -10.549  1.00 51.94           N
ATOM   1461  CA  TYR P  82      33.797  26.440 -10.985  1.00 52.99           C
ATOM   1462  C   TYR P  82      34.257  24.985 -10.876  1.00 52.34           C
ATOM   1463  O   TYR P  82      35.456  24.712 -10.890  1.00 52.96           O
ATOM   1464  CB  TYR P  82      34.769  27.349 -10.225  1.00 54.88           C
ATOM   1465  CG  TYR P  82      34.757  27.162  -8.722  1.00 56.85           C
ATOM   1466  CD1 TYR P  82      35.404  26.078  -8.127  1.00 57.60           C
ATOM   1467  CD2 TYR P  82      34.096  28.070  -7.895  1.00 57.47           C
ATOM   1468  CE1 TYR P  82      35.394  25.903  -6.743  1.00 58.41           C
ATOM   1469  CE2 TYR P  82      34.080  27.904  -6.510  1.00 58.47           C
ATOM   1470  CZ  TYR P  82      34.730  26.819  -5.942  1.00 58.59           C
ATOM   1471  OH  TYR P  82      34.718  26.650  -4.577  1.00 59.33           O
HETATM 1782  O   HOH P 134      28.550  22.420  -4.280  1.00 33.59           O
END
"""
    )
  opt = _optimizeFragment(pdb_1dfu_his80)
  movers = opt._movers
  if len(movers) != 3:
    return "Optimizers.Test(): Incorrect number of Movers for His placement test: " + str(len(movers))
  info = opt.getInfo()
  if not 'HD1NotPlaced' in info:
    return "Optimizers.Test(): Unexpected HD1 placement for His."
  if not 'HE2Placed' in info:
    return "Optimizers.Test(): Missing HE2 placement for His."


  ################################################################################
  # Test using snippet from 1xso to ensure that the Histidine placement code will lock down the
  # Histidine, set its Nitrogen states, and mark its Hydrogens for deletion.  To
  # ensure that hydrogen placement puts them there in the first place, we first
  # move the CU and ZN far from the Histidine before adding Hydrogens, then move
  # them back before building the spatial hierarchy and testing.
  pdb_1xso_his_61_and_ions = (
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
HETATM 2190 CU    CU A   1      22.291  33.388   3.996  1.00 13.22          CU
HETATM 2191 ZN    ZN A 152      27.539  36.010   2.881  1.00  9.34          ZN
END
"""
    )
  dm = DataManager(['model'])
  dm.process_model_str("1xso_snip.pdb",pdb_1xso_his_61_and_ions)
  model = dm.get_model()

  # Make sure we have a valid unit cell.  Do this before we add hydrogens to the model
  # to make sure we have a valid unit cell.
  model = shift_and_box_model(model = model)

  # Find and move the Copper and Zinc atoms far from the Histidine so that
  # Hydrogen placement and bond proxies won't consider them to be bonded.
  # This tests the algorithm behavior for the case of all Hydrogens present
  # in all cases.
  for a in model.get_hierarchy().models()[0].atoms():
    if a.element.upper() == "CU":
      origPositionCU = a.xyz
      a.xyz = (1000, 1000, 1000)
    if a.element.upper() == "ZN":
      origPositionZN = a.xyz
      a.xyz = (1000, 1000, 1000)

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints

  # Get the first model in the hierarchy.
  firstModel = model.get_hierarchy().models()[0]

  # Get the list of alternate conformation names present in all chains for this model.
  alts = AlternatesInModel(firstModel)

  # Get the atoms from the first conformer in the first model (the empty string is the name
  # of the first conformation in the model; if there is no empty conformation, then it will
  # pick the first available conformation for each atom group.
  atoms = GetAtomsForConformer(firstModel, "")

  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = Helpers.getBondedNeighborLists(atoms, bondProxies)

  # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
  probePhil = _philLike()
  ret = Helpers.getExtraAtomInfo(model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=False, probePhil=probePhil)
  extra = ret.extraAtomInfo

  # Also compute the maximum VDW radius among all atoms.
  maxVDWRad = 1
  for a in atoms:
    maxVDWRad = max(maxVDWRad, extra.getMappingFor(a).vdwRadius)

  # Put the Copper and Zinc back in their original positions before we build the
  # spatial-query structure.  This will make them close enough to be bonded to
  # the Nitrogens and should cause Hydrogen removal and marking of the Nitrogens
  # as acceptors.
  for a in model.get_hierarchy().models()[0].atoms():
    if a.element.upper() == "CU":
      a.xyz = origPositionCU
    if a.element.upper() == "ZN":
      a.xyz = origPositionZN

  # Optimization will place the movers, which should be none because the Histidine flip
  # will be constrained by the ionic bonds. There is a rotatable hydrogen placed at the
  # terminus, but no Flip Mover.
  probePhil = _philLike()
  opt = Optimizer(probePhil, True, model)
  movers = opt._movers
  if len(movers) != 1:
    return "Optimizers.Test(): Incorrect number of Movers for 1xso Histidine test: " + str(len(movers))

  # Make sure that the two ring Nitrogens have been marked as acceptors.
  # Make sure that the two hydrogens have been marked for deletion.
  for a in model.get_hierarchy().models()[0].atoms():
    name = a.name.strip()
    if name in ["ND1", "NE2"]:
      if not opt._extraAtomInfo.getMappingFor(a).isAcceptor:
        return 'Optimizers.Test(): '+name+' in 1xso Histidine test was not an acceptor'
    if name in ["HD1", "HE2"]:
      if not a in opt._deleteMes:
        return 'Optimizers.Test(): '+name+' in 1xso Histidine test was not set for deletion'

  #========================================================================
  # Check a clique with multiple elements to be sure that it was properly
  # globally optimized. This is a set of ACT residues that are offset
  # such that they want to line up the same way.  A carbon is placed to
  # force the oriention away from the initial solution.
  pdb_multi_act = (
"""
CRYST1   93.586  127.886  251.681  90.00  90.00  90.00 I 2 2 2
SCALE1      0.010685  0.000000  0.000000        0.00000
SCALE2      0.000000  0.007819  0.000000        0.00000
SCALE3      0.000000  0.000000  0.003973        0.00000
HETATM    1  C   ACT A   1       6.494 -47.273 -37.006  1.00 16.65           C
HETATM    2  O   ACT A   1       5.654 -47.981 -37.645  1.00 15.33           O
HETATM    3  CH3 ACT A   1       6.790 -47.410 -35.548  1.00 17.13           C
HETATM    4  OXT ACT A   1       7.110 -46.425 -37.699  1.00 17.05           O
HETATM  101  C   ACT A 101       8.994 -49.773 -37.606  1.00 16.65           C
HETATM  102  O   ACT A 101       8.154 -50.481 -38.245  1.00 15.33           O
HETATM  103  CH3 ACT A 101       9.290 -49.910 -36.148  1.00 17.13           C
HETATM  104  OXT ACT A 101       9.610 -48.925 -38.299  1.00 17.05           O
HETATM  111  C   ACT A 111      11.494 -52.273 -38.206  1.00 16.65           C
HETATM  112  O   ACT A 111      10.654 -52.981 -38.845  1.00 15.33           O
HETATM  113  CH3 ACT A 111      11.790 -52.410 -36.748  1.00 17.13           C
HETATM  114  OXT ACT A 111      12.110 -51.425 -38.899  1.00 17.05           O
HETATM  121  C   ACT A 121      13.994 -54.773 -38.806  1.00 16.65           C
HETATM  122  O   ACT A 121      13.154 -55.481 -39.445  1.00 15.33           O
HETATM  123  CH3 ACT A 121      14.290 -54.910 -37.348  1.00 17.13           C
HETATM  124  OXT ACT A 121      14.610 -53.925 -39.499  1.00 17.05           O
HETATM  131  C   ACT A 131      16.494 -57.273 -39.406  1.00 16.65           C
HETATM  132  O   ACT A 131      15.654 -57.981 -40.045  1.00 15.33           O
HETATM  133  CH3 ACT A 131      16.790 -57.410 -37.948  1.00 17.13           C
HETATM  134  OXT ACT A 131      17.110 -56.425 -40.099  1.00 17.05           O
HETATM  141  C   ACT A 141      18.994 -59.773 -40.006  1.00 16.65           C
HETATM  142  O   ACT A 141      18.154 -60.481 -40.645  1.00 15.33           O
HETATM  143  CH3 ACT A 141      19.290 -59.910 -38.548  1.00 17.13           C
HETATM  144  OXT ACT A 141      19.610 -58.925 -40.699  1.00 17.05           O
HETATM  200  C   ACT A 200       4.500 -45.500 -35.000  1.00 16.65           C
END
"""
    )
  dm = DataManager(['model'])
  dm.process_model_str("pdb_multi_act.pdb",pdb_multi_act)
  model = dm.get_model()

  # Make sure we have a valid unit cell.  Do this before we add hydrogens to the model
  # to make sure we have a valid unit cell.
  model = shift_and_box_model(model = model)

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints

  # Get the first model in the hierarchy.
  firstModel = model.get_hierarchy().models()[0]

  # Get the list of alternate conformation names present in all chains for this model.
  alts = AlternatesInModel(firstModel)

  # Get the atoms from the first conformer in the first model (the empty string is the name
  # of the first conformation in the model; if there is no empty conformation, then it will
  # pick the first available conformation for each atom group.
  atoms = GetAtomsForConformer(firstModel, "")

  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = Helpers.getBondedNeighborLists(atoms, bondProxies)

  # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
  probePhil = _philLike()
  ret = Helpers.getExtraAtomInfo(model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=False,probePhil=probePhil)
  extra = ret.extraAtomInfo

  # Also compute the maximum VDW radius among all atoms.
  maxVDWRad = 1
  for a in atoms:
    maxVDWRad = max(maxVDWRad, extra.getMappingFor(a).vdwRadius)

  # Optimization will place the movers. Make sure we got as many as we expected.
  # Make sure that the orientation for all of the movers is correct.
  # Test with each type of optimizer, from the base to the more derived, so
  # that we find out about failures on the base classes first.
  probePhil = _philLike()

  global _DoCliqueOptimizationInC

  print('Testing Optimizer')
  opt = Optimizer(probePhil, True, model, modelIndex = 0, altID = None,
                bondedNeighborDepth = 4,
                useNeutronDistances = False,
                minOccupancy = 0.02,
                preferenceMagnitude = 1.0,
                nonFlipPreference = 0.5,
                skipBondFixup = False,
                flipStates = '',
                verbosity = 1)
  movers = opt._movers
  if len(movers) != 6:
    return "Optimizers.Test(): Incorrect number of Movers for Optimizer multi-ACT test: " + str(len(movers))
  res = re.findall(r'pose Angle [-+]?\d+', opt.getInfo())
  for r in res:
    if not 'pose Angle 90' in r:
      return "Optimizers.Test(): Unexpected angle ("+str(r)+") for Optimizer multi-ACT test"

  #========================================================================
  # Check a case where an AmideFlip would be locked down and have its Hydrogen removed.
  # @todo

  #========================================================================
  # Check that the occupancy and B-factor cut-offs for water Oxygens are causing them
  # to be ignored in the calculations by putting an atom in the way and making sure it
  # is ignored.
  pdb_b_factor = (
"""
CRYST1   93.586  127.886  251.681  90.00  90.00  90.00 I 2 2 2
SCALE1      0.010685  0.000000  0.000000        0.00000
SCALE2      0.000000  0.007819  0.000000        0.00000
SCALE3      0.000000  0.000000  0.003973        0.00000
HETATM    1  C   ACT A   1       6.494 -47.273 -37.006  1.00 16.65           C
HETATM    2  O   ACT A   1       5.654 -47.981 -37.645  1.00 15.33           O
HETATM    3  CH3 ACT A   1       6.790 -47.410 -35.548  1.00 17.13           C
HETATM    4  OXT ACT A   1       7.110 -46.425 -37.699  1.00 17.05           O
HETATM  200  O   HOH A 200       8.200 -48.610 -35.148  1.00 50.00           O
HETATM  201  O   HOH A 201       5.800 -46.810 -34.548  1.00 20.00           O
END
"""
    )
  dm = DataManager(['model'])
  dm.process_model_str("pdb_b_factor.pdb",pdb_b_factor)
  model = dm.get_model()

  # Make sure we have a valid unit cell.  Do this before we add hydrogens to the model
  # to make sure we have a valid unit cell.
  model = shift_and_box_model(model = model)

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints

  # Get the first model in the hierarchy.
  firstModel = model.get_hierarchy().models()[0]

  # Get the list of alternate conformation names present in all chains for this model.
  alts = AlternatesInModel(firstModel)

  # Get the atoms from the first conformer in the first model (the empty string is the name
  # of the first conformation in the model; if there is no empty conformation, then it will
  # pick the first available conformation for each atom group.
  atoms = GetAtomsForConformer(firstModel, "")

  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = Helpers.getBondedNeighborLists(atoms, bondProxies)

  # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
  probePhil = _philLike()
  ret = Helpers.getExtraAtomInfo(model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=False,probePhil=probePhil)
  extra = ret.extraAtomInfo

  # Also compute the maximum VDW radius among all atoms.
  maxVDWRad = 1
  for a in atoms:
    maxVDWRad = max(maxVDWRad, extra.getMappingFor(a).vdwRadius)

  # Optimization will place the movers. Make sure we got as many as we expected.
  # Make sure that the orientation for all of the movers is correct.
  probePhil = _philLike()
  opt = Optimizer(probePhil, True, model, modelIndex = 0, altID = None,
                bondedNeighborDepth = 4,
                useNeutronDistances = False,
                minOccupancy = 0.02,
                preferenceMagnitude = 1.0,
                nonFlipPreference = 0.5,
                skipBondFixup = False,
                flipStates = '',
                verbosity = 1)
  movers = opt._movers
  if len(movers) != 1:
    return "Optimizers.Test(): Incorrect number of Movers for B-factor test: " + str(len(movers))
  res = re.findall(r'pose Angle [-+]?\d+', opt.getInfo())
  for r in res:
    if not 'pose Angle 90' in r:
      return "Optimizers.Test(): Unexpected angle ("+str(r)+") for B-factor test"

  #========================================================================
  # Generate an example data model with a small molecule in it or else read
  # from the specified file.
  infoString = ""
  if inFileName is not None and len(inFileName) > 0:
    # Read a model from a file using the DataManager
    print('Reading model from',inFileName)
    dm = DataManager()
    dm.process_model_file(inFileName)
    model = dm.get_model(inFileName)
  else:
    # Generate a small-molecule model using the map model manager
    print('Generating model')
    mmm=map_model_manager()         #   get an initialized instance of the map_model_manager
    mmm.generate_map()              #   get a model from a generated small library model and calculate a map for it
    model = mmm.model()             #   get the model

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Interpret the model after shifting and adding Hydrogens to it so that
  # all of the needed fields are filled in when we use them below.
  # @todo Remove this once place_hydrogens() does all the interpretation we need.
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  model.process(make_restraints=True, pdb_interpretation_params=p) # make restraints

  # Run the optimizer on the model to make sure it doesn't crash.
  opt = Optimizer(probePhil, True, model)

  # Write debugging output if we've been asked to
  if dumpAtoms:
    fn = "deleteme.pdb"
    fn = model.get_hierarchy().write_pdb_or_mmcif_file(
        target_format = 'pdb',
        target_filename = fn)
    print("Wrote model to '%s'" %fn)
    f = open("atomDump.pdb","w")
    f.write(opt.getAtomDump())

  #========================================================================
  # @todo Unit test a multi-model case, a multi-alternate case, and singles of each.

  return ""

##################################################################################
# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  parser = argparse.ArgumentParser(description='Test mmtbx.reduce.Optimizers.')
  parser.add_argument("--dumpAtoms", help="dump the atoms into PDB files to help debug", action="store_true")
  parser.add_argument('inputFile', nargs='?', default="")
  args = parser.parse_args()

  ret = Test(args.inputFile, args.dumpAtoms)
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
