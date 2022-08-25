##################################################################################
#                Copyright 2021  Richardson Lab at Duke University
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

import argparse, time, itertools

from mmtbx.reduce import Movers
from mmtbx.reduce import InteractionGraph
from boost_adaptbx import graph
from boost_adaptbx.graph import connected_component_algorithm as cca

from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
from iotbx import pdb
from iotbx.pdb import common_residue_names_get_class
import mmtbx
from scitbx.array_family import flex

from mmtbx.probe import Helpers
import mmtbx_probe_ext as probeExt

# To enable addition of Hydrogens
# @todo See if we can remove the shift and box once reduce_hydrogen is complete
from cctbx.maptbx.box import shift_and_box_model
from mmtbx.hydrogens import reduce_hydrogen

##################################################################################
# This file includes a set of functions and classes that implement placement and optimization of
# Reduce's "Movers".

##################################################################################
# Module-scoped attributes that can be set to modify behavior.

# Default value of 1 reports standard information.
# Value of 2 reports timing information.
# Setting it to 0 removes all informational messages.
# Setting it above 1 provides additional debugging information.
verbosity = 3

# Probe PHIL parameters that are passed to Probe methods.  These enable modification of
# Probe behavior without having to pass individual parameters through the various Optimizer
# classes.
# @todo Consider making this a constructor parameter rather than a global.
probePhil = None

##################################################################################
# Helper functions

def _VerboseCheck(level, message):
  # Returns "" if the level is less than global verbosity level, message if it is at or above
  if verbosity >= level:
    return " "*level + message
  else:
    return ""

_lastTime = None
def _ReportTiming(message):
  """Use None message to start the timer without printing.
  """
  from libtbx.development.timers import work_clock

  global _lastTime
  if message is None:
    _lastTime = work_clock()
    return
  curTime = work_clock()
  diff = curTime - _lastTime
  _lastTime = curTime
  return _VerboseCheck(2,"Time to {}: {:0.3f}".format(message,diff)+"\n")

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
  # Don't print the code if it is a space (blank).
  insertionCode = a.parent().parent().icode.strip()
  return "chain "+str(chainID)+" "+resName+" "+resID+insertionCode

##################################################################################
# Optimizers:
#     _SingletonOptimizer: This is a base Optimizer class that implements the placement and scoring
# of Movers and optimizes all Movers independently, not taking into account their impacts
# on each other.  The other Optimizers will all be derived from this base class and will
# override the multi-Mover Clique optimization routine.  This base class does produce an
# InteractionGraph telling which Movers might possibly interact during their motion but it
# treats all Movers as if they were singletons rather than Cliques.  It should not be used
# by client code.
#     _BruteForceOptimizer: This runs brute-force optimization on each Clique independently.
# This is much to slow for many molecules that have a large number of interacting Movers
# in their largest clique.  It should not be used by client code, and is here mainly to
# support unit testing by providing a baseline for comparison against faster algorithms.
#     _CliqueOptimizer: This uses a recursive divide-and-conquer approach to break
# each clique into separate non-interacting subgraphs so that it only has to compute complete
# interactions for a small number of interacting Movers at a time.  It should not be used
# by client code.
#     FastOptimizer: This is the optimizer that should be used by client code.  This uses
# a dictionary telling which Movers a given atom may interact with to cache scores for atoms
# in each combination of coarse locations for these Movers, greatly reducing the number of
# atom scores that must be computed and speeding up calculations.

class _SingletonOptimizer(object):

  def __init__(self, addFlipMovers, model, modelIndex, altID,
                bondedNeighborDepth,
                probeRadius,
                useNeutronDistances,
                probeDensity,
                minOccupancy,
                preferenceMagnitude
              ):
    """Constructor for _SingletonOptimizer.  This is the base class for all optimizers and
    it implements the machinery that finds and optimized Movers.
    This class optimizes all Movers independently,
    ignoring their impact on one another.  This means that it will not find the global
    minimum for any multi-Mover Cliques.  Derived classes should override the
    _optimizeCliqueCoarse() and other methods as needed to improve correctness and speed.
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
    hierarchy.
    :param altID: The conformer alternate location specifier to use.  The value "" will
    cause it to run on the first conformer found in each model.  If this is set to None,
    optimization will be run sequentially for every conformer in the model, starting with
    the last and ending with the first.  This will leave the initial conformer's values as the
    final location for atoms that are not inside a conformer or are in the first conformer.
    :param bondedNeighborDepth: How many hops to ignore bonding when doing Probe calculations.
    A depth of 3 will ignore my bonded neighbors and their bonded
    neighbors and their bonded neighbors.
    :param probeRadius: Radius of the probe to be used in Probe calculations (Angstroms).
    :param useNeutronDistances: False will use X-ray/electron cloud distances.  If set to
    True, it will use neutron (nuclear) distances.  This must be set consistently with the
    values used to generate the hydrogens and to run PDB interpretation.
    :param probeDensity: How many dots per sq Angstroms in VDW calculations.
    :param minOccupancy: Minimum occupancy for an atom to be considered in the Probe score.
    :param preferenceMagnitude: Multiplier for the preference energies expressed
    by some Movers for particular orientations.
    """

    ################################################################################
    # Store the parameters that will be accessed by other methods
    self._bondedNeighborDepth = bondedNeighborDepth
    self._probeRadius = probeRadius
    self._useNeutronDistances = useNeutronDistances
    self._probeDensity = probeDensity
    self._minOccupancy = minOccupancy
    self._preferenceMagnitude = preferenceMagnitude

    ################################################################################
    # Initialize internal variables.
    self._infoString = ""
    self._atomDump = ""
    self._waterOccCutoff = 0.66 # @todo Make this a parameter, -WaterOCCcutoff param in reduce
    self._waterBCutoff = 40.0   # @todo Make this a parameter, -WaterBcutoff param in reduce

    ################################################################################
    # Get the Cartesian positions of all of the atoms in the entire model and find
    # the bond proxies for all of them.  The proxies will not attempt to bond atoms
    # that are in different model indices.
    _ReportTiming(None) # Reset timer
    carts = flex.vec3_double()
    for a in model.get_atoms():
      carts.append(a.xyz)
    self._infoString += _ReportTiming("get coordinates")
    bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
    self._infoString += _ReportTiming("compute bond proxies")

    ################################################################################
    # Get the bonded neighbor lists for all of the atoms in the model, so we don't get
    # failures when we look up an atom from another in Helpers.getExtraAtomInfo().
    # We won't get bonds between atoms in different conformations.
    bondedNeighborLists = Helpers.getBondedNeighborLists(model.get_atoms(), bondProxies)
    self._infoString += _ReportTiming("compute bonded neighbor lists")

    ################################################################################
    # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
    # This is done for all atoms in the model.
    global probePhil
    # This must operate on the entire model, not just the current model index.
    ret = Helpers.getExtraAtomInfo(
      model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=self._useNeutronDistances, probePhil=probePhil)
    self._extraAtomInfo = ret.extraAtomInfo
    self._infoString += ret.warnings
    self._infoString += _ReportTiming("get extra atom info")

    ################################################################################
    # Run optimization for every desired conformer and every desired model, calling
    # placement and then a derived-class single optimization routine for each.  When
    # the modelIndex or altID is None, that means to run over all available cases.
    # For alternates, if there is a non-empty (not "" or " ") case, then we run backwards
    # from the last to the first but do not run for the empty case; if there are only
    # empty cases, then we run just once.  We run the models in order, all of them when
    # None is specified and the specified one if it is specified.

    model.setup_riding_h_manager()
    riding_h_manager = model.get_riding_h_manager()
    h_parameterization = riding_h_manager.h_parameterization
    allRotatableHydrogens = model.rotatable_hd_selection(iselection=True)
    startModelIndex = 0
    stopModelIndex = len(model.get_hierarchy().models())
    if modelIndex is not None:
      startModelIndex = modelIndex
      stopModelIndex = modelIndex + 1
    for mi in range(startModelIndex, stopModelIndex):
      # Get the specified model from the hierarchy.
      myModel = model.get_hierarchy().models()[mi]

      # Get the list of alternate conformation names present in all chains for this model.
      # If there is more than one result, remove the empty results and then sort them
      # in reverse order so we finalize all non-alternate ones to match the first.
      alts = AlternatesInModel(myModel)
      if len(alts) > 1:
        alts.discard("")
        alts.discard(" ")
      alts = sorted(list(alts), reverse=True)
      self._infoString += _ReportTiming("compute alternates")

      # If there is a specified alternate, use it.
      if altID is not None:
        alts = [altID]

      # Clear the Movers list for each model.  It will be retained from one alternate to the next.
      self._movers = []
      for alt in alts:
        # If we are doing the second or later alternate, place all Movers that are in a compatible alternate
        # back into their initial configuration so we start from the same state we would have if this were
        # the only alternate being tested.  This will ensure that we get compatible outputs when run either
        # way (we may end up with equivalent but different results, like 120 degree rotations for 3 hydrogens).
        for m in self._movers:
          coarse = m.CoarsePositions()
          if coarse.atoms[0].parent().altloc in ['', ' ', alt]:
            self._setMoverState(coarse, 0)

            # Apply any location and information fixups needed for the initial configuration.
            # This will in all cases put things back into their original configuration.
            self._doFixup(m.FixUp(0))

        # Tell about the run we are currently doing.
        self._infoString += _VerboseCheck(1,"Running Reduce optimization on model index "+str(mi)+
          ", alternate '"+alt+"'\n")
        self._infoString += _VerboseCheck(1,"  bondedNeighborDepth = "+str(self._bondedNeighborDepth)+"\n")
        self._infoString += _VerboseCheck(1,"  probeRadius = "+str(self._probeRadius)+"\n")
        self._infoString += _VerboseCheck(1,"  useNeutronDistances = "+str(self._useNeutronDistances)+"\n")
        self._infoString += _VerboseCheck(1,"  probeDensity = "+str(self._probeDensity)+"\n")
        self._infoString += _VerboseCheck(1,"  minOccupancy = "+str(self._minOccupancy)+"\n")
        self._infoString += _VerboseCheck(1,"  preferenceMagnitude = "+str(self._preferenceMagnitude)+"\n")

        # Get the atoms from the specified conformer in the model (the empty string is the name
        # of the first conformation in the model; if there is no empty conformation, then it will
        # pick the first available conformation for each atom group.
        self._atoms = GetAtomsForConformer(myModel, alt)

        ################################################################################
        # Reset the timer
        _ReportTiming(None)

        ################################################################################
        # Construct the spatial-query information needed to quickly determine which atoms are nearby
        self._spatialQuery = probeExt.SpatialQuery(self._atoms)
        self._infoString += _ReportTiming("construct spatial query")

        ################################################################################
        # Initialize any per-alternate data structures now that we have the atoms and
        # other data structures initialized.
        self._initializeAlternate()

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
        # Get the list of Movers using the _PlaceMovers private function.
        # The list of rotatable hydrogens comes from the global model, not just the current
        # model index.  However, we only place on atoms that are also in self._atoms, which
        # only includes those from the current model index.
        self._infoString += _ReportTiming("select rotatable hydrogens")
        ret = _PlaceMovers(self._atoms, allRotatableHydrogens,
                           bondedNeighborLists, h_parameterization, self._spatialQuery,
                           self._extraAtomInfo, self._maximumVDWRadius, addFlipMovers)
        self._infoString += ret.infoString
        self._movers = ret.moverList
        self._infoString += _VerboseCheck(1,"Inserted "+str(len(self._movers))+" Movers\n")
        self._infoString += _VerboseCheck(1,'Marked '+str(len(ret.deleteAtoms))+' atoms for deletion\n')
        self._infoString += _ReportTiming("place movers")
        self._moverInfo = ret.moverInfo

        ################################################################################
        # Add the atoms that were unconditionally marked for deletion during placement
        # to the set of atoms to delete.
        self._deleteMes = self._deleteMes.union(set(ret.deleteAtoms))
        for m in self._movers:
          pr = m.CoarsePositions()
          self._setMoverState(pr, 0)
        self._infoString += _ReportTiming("initialize Movers")

        ################################################################################
        # Initialize high score for each Mover, filling in a None value for each.
        self._highScores = {}
        for m in self._movers:
          self._highScores[m] = None

        ################################################################################
        # Compute the interaction graph, of which each connected component is a Clique.
        # Get a list of singleton Cliques and a list of other Cliques.  Keep separate lists
        # of the singletons and the groups.
        self._interactionGraph, self._atomMoverSets = InteractionGraph.InteractionGraphAllPairs(self._movers, self._extraAtomInfo,
          probeRadius=probeRadius)
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
        self._infoString += _VerboseCheck(1,"Found "+str(len(components))+" Cliques ("+
            str(len(singletonCliques))+" are singletons); largest Clique size = "+
            str(maxLen)+"\n")
        self._infoString += _ReportTiming("compute interaction graph")

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

        # Get the excluded list for each atom in the set, making a dictionary
        self._excludeDict = {}
        for a in moverAtoms:
          self._excludeDict[a] = mmtbx.probe.Helpers.getAtomsWithinNBonds(a,
            bondedNeighborLists, self._extraAtomInfo, probeRadius, self._bondedNeighborDepth)
        self._infoString += _ReportTiming("determine excluded atoms")

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
        for a in self._atoms:
          if a.element == 'O' and common_residue_names_get_class(name=a.parent().resname) == "common_water":
            if a.occ >= self._waterOccCutoff and a.b < self._waterBCutoff:

              # We're an acceptor and not a donor.
              ei = self._extraAtomInfo.getMappingFor(a)
              ei.isDonor = False
              ei.isAcceptor = True
              self._extraAtomInfo.setMappingFor(a, ei)

              newPhantoms = Helpers.getPhantomHydrogensFor(a, self._spatialQuery, self._extraAtomInfo, self._minOccupancy,
                              False, phantomHydrogenRadius, placedHydrogenDistance)
              if len(newPhantoms) > 0:
                resNameAndID = _ResNameAndID(a)
                self._infoString += _VerboseCheck(3,"Added {} phantom Hydrogens on {}\n".format(len(newPhantoms), resNameAndID))
                for p in newPhantoms:
                  self._infoString += _VerboseCheck(5,"Added phantom Hydrogen at "+str(p.xyz)+"\n")
                  phantoms.append( (p,a) )

            else:
              # Occupancy or B factor are out of bounds, so remove this atom from consideration.
              self._infoString += _VerboseCheck(3,"Ignoring "+
                a.name.strip()+" "+a.parent().resname.strip()+" "+str(a.parent().parent().resseq_as_int())+
                " "+str(a.parent().parent().parent().id)+
                " with occupancy "+str(a.occ)+" and B factor "+str(a.b)+"\n")
              watersToDelete.append(a)

        if len(watersToDelete) > 0:
          self._infoString += _VerboseCheck(1,"Ignored "+str(len(watersToDelete))+" waters due to occupancy or B factor\n")
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

          self._infoString += _VerboseCheck(1,"Added "+str(len(phantoms))+" phantom Hydrogens on waters")
          self._infoString += _VerboseCheck(1," (Old total "+str(origCount)+", new total "+str(len(self._atoms))+")\n")
        self._infoString += _ReportTiming("place water phantom Hydrogens")

        ################################################################################
        # Fix up the donor status for all of the atoms now that we've added the final explicit
        # Phantom Hydrogens.
        Helpers.fixupExplicitDonors(self._atoms, bondedNeighborLists, self._extraAtomInfo)
        atomDump = Helpers.writeAtomInfoToString(self._atoms, self._extraAtomInfo)

        ################################################################################
        # Construct dot-spheres for each atom that we may need to find interactions for.
        # This must be done after the phantom Hydrogens have been added so that they will be included.
        dotSphereCache = probeExt.DotSphereCache(self._probeDensity)
        self._dotSpheres = {}
        for a in self._atoms:
          self._dotSpheres[a] = dotSphereCache.get_sphere(self._extraAtomInfo.getMappingFor(a).vdwRadius)
        self._infoString += _ReportTiming("compute dot spheres")

        ################################################################################
        # Contruct the DotScorer object we'll use to score the dots.
        self._dotScorer = probeExt.DotScorer(self._extraAtomInfo)
        self._infoString += _ReportTiming("construct dot scorer")

        ################################################################################
        # Compute and store the initial score for each Mover in its info
        for m in self._movers:
          coarse = m.CoarsePositions()
          score = coarse.preferenceEnergies[0] * self._preferenceMagnitude
          self._setMoverState(coarse, 0)
          for a in coarse.atoms:
            score += self._scoreAtom(a)
          self._moverInfo[m] += " Initial score: {:.2f}".format(score)

        ################################################################################
        # Call internal methods to optimize the single-element Cliques and then to optimize
        # the multi-element Cliques and then to do indepenedent fine adjustment of all
        # Cliques.  Subclasses should overload the called routines, but the global approach
        # taken here will be the same for all of them.  If we want to change the recipe
        # so that we can do global fine optimization, we'll do that here rather than in the
        # subclasses.

        # Do coarse optimization on the singleton Movers.  Record the selected coarse
        # index.
        self._coarseLocations = {}
        for m in self._movers:
          self._coarseLocations[m] = 0
        for s in singletonCliques:
          mover = self._interactionGraph.vertex_label(s[0])
          ret = self._optimizeSingleMoverCoarse(mover)
          self._infoString += _VerboseCheck(1,"Singleton optimized with score {:.2f}\n".format(ret))
        self._infoString += _ReportTiming("optimize singletons (coarse)")

        # Do coarse optimization on the multi-Mover Cliques.
        for g in groupCliques:
          movers = [self._interactionGraph.vertex_label(i) for i in g]
          subset = _subsetGraph(self._interactionGraph, movers)
          ret = self._optimizeCliqueCoarse(subset)
          self._infoString += _VerboseCheck(1,"Clique optimized with score {:.2f}\n".format(ret))
        self._infoString += _ReportTiming("optimize cliques (coarse)")

        # Do fine optimization on the Movers.  This is done independently for
        # each of them, whether they are part of a multi-Mover Clique or not.
        self._fineLocations = {}
        for m in self._movers:
          self._fineLocations[m] = 0
        self._infoString += _VerboseCheck(1,"Fine optimization on all Movers\n")
        for m in self._movers:
          self._optimizeSingleMoverFine(m)
        self._infoString += _ReportTiming("optimize all Movers (fine)")

        ################################################################################
        # Print the final state and score for all Movers
        def _scoreMoverReportClash(self, m, index):
          coarse = m.CoarsePositions()
          score = coarse.preferenceEnergies[index] * self._preferenceMagnitude
          clash = False
          for atom in coarse.atoms:
            maxRadiusWithoutProbe = self._extraAtomInfo.getMappingFor(atom).vdwRadius + self._maximumVDWRadius
            res = self._dotScorer.score_dots(atom, self._minOccupancy, self._spatialQuery,
              maxRadiusWithoutProbe, self._probeRadius, self._excludeDict[atom], self._dotSpheres[atom].dots(),
              self._probeDensity, False)
            score += res.totalScore()
            if res.hasBadBump:
              clash = True
          return score, clash

        def _printPose(self, m):
          description = m.PoseDescription(self._coarseLocations[m], self._fineLocations[m])

          # If the Mover is a flip of some kind, then the substring "lipped" will be present
          # in the description (Flipped and Unflipped both have this subtring).
          # When that happens, we check the final state and the other flip state
          # state (which is half of the coarse states away) to see if both have clashes or
          # if they are close in energy. If so, then we annotate the output.
          # We add the same number of words to the output string in all cases to make things
          # easier for a program to parse.
          if "lipped" in description:
            coarse = m.CoarsePositions()
            numPositions = len(coarse.positions)
            final = self._coarseLocations[m]
            other = (final + numPositions//2) % numPositions
            self._setMoverState(coarse, other)
            otherScore, otherBump = _scoreMoverReportClash(self, m, other)
            self._setMoverState(coarse, final)
            finalScore, finalBump = _scoreMoverReportClash(self, m, final)
            if otherBump and finalBump:
              description += " BothClash"
            elif abs(finalScore - otherScore) <= 0.5:
              description += " Uncertain"
            else:
              description += " ."
          else:
            description += " ."

          self._infoString += _VerboseCheck(1,"  {} final score: {:.2f} pose {}\n".format(
            self._moverInfo[m], self._highScores[m], description ))

        self._infoString += _VerboseCheck(1,"BEGIN REPORT: Model "+str(mi)+" Alt '"+alt+"':\n")
        sortedGroups = sorted(groupCliques, key=len, reverse=True)
        for g in sortedGroups:
          self._infoString += _VerboseCheck(1," Set of "+str(len(g))+" Movers:")
          movers = [self._interactionGraph.vertex_label(i) for i in g]
          # Parse the record for each mover and pull out the initial score.  Sum the initial and
          # final scores across all Movers in the group and report this.
          initial = 0.0
          final = 0.0
          for m in movers:
            initial += float(self._moverInfo[m].split()[9])
            final += self._highScores[m]
          self._infoString += _VerboseCheck(1," Totals: initial score {:.2f}, final score {:.2f}\n".format(initial, final))
          for m in movers:
            _printPose(self, m)
        self._infoString += _VerboseCheck(1," Singleton Movers:\n")
        for s in singletonCliques:
          m = self._interactionGraph.vertex_label(s[0])
          _printPose(self, m)
        self._infoString += _VerboseCheck(1,"END REPORT\n")

        ################################################################################
        # Do FixUp on the final coarse orientations.  Set the positions, extra atom info
        # and deletion status for all atoms that have entries for each.
        # This must be done after we print the scores because the print methods move the
        # coarse state to see how much it changed.
        self._infoString += _VerboseCheck(1,"FixUp on all Movers\n")
        for m in self._movers:
          loc = self._coarseLocations[m]
          self._infoString += _VerboseCheck(3,"FixUp on {} coarse location {}\n".format(
          self._moverInfo[m],loc))
          self._doFixup(m.FixUp(loc))
        self._infoString += _ReportTiming("fix up Movers")

      ################################################################################
      # Deletion of atoms (Hydrogens) that were requested by Histidine FixUp()s,
      # both in the initial setup and determined during optimization.  Phantom Hydrogens
      # on waters do not need to be adjusted because they were never added to the
      # structure.
      # We only do this after the last-checked alternate configuration for a given model.
      self._infoString += _VerboseCheck(1,"Deleting Hydrogens tagged by Histidine Movers\n")
      for a in self._deleteMes:
        aName = a.name.strip().upper()
        resNameAndID = _ResNameAndID(a)
        self._infoString += _VerboseCheck(5,"Deleting {} {}\n".format(resNameAndID, aName))
        a.parent().remove_atom(a)
      self._infoString += _ReportTiming("delete Hydrogens")

      #################################################################################
      # Dump information about all of the atoms in the model into a string.
      self._atomDump = Helpers.writeAtomInfoToString(myModel.atoms(), self._extraAtomInfo)

  def getInfo(self):
    """
      Returns information that the user may care about regarding the processing.  The level of
      detail on this can be set by setting Optimizers.verbosity before creating or calling methods.
      :return: the information so far collected in the string.  Calling this method also clears
      the information, so that later calls will not repeat it.
    """
    ret = self._infoString
    self._infoString = ""
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

  def _initializeAlternate(self):
    """
    Override this method in derived classes.
    This is a place for derived classes to perform any operations that should be done at
    the start of every alternate.  For example, initializing per-atom caches.
    """
    return

  def _setMoverState(self, positionReturn, index):
    """
      Move the atoms to their new positions, updating the spatial query structure
      by removing the old and adding the new location.
    """
    for i, a in enumerate(positionReturn.atoms):
      self._spatialQuery.remove(a)
      a.xyz = positionReturn.positions[index][i]
      self._spatialQuery.add(a)
    # Update the extra atom information associated with the atom.
    for i, e in enumerate(positionReturn.extraInfos[index]):
      self._extraAtomInfo.setMappingFor(positionReturn.atoms[i], e)
    # Manage the deletion status of each atom, including ensuring
    # consistency with the spatial-query structure.
    for i, doDelete in enumerate(positionReturn.deleteMes[index]):
      if doDelete:
        self._spatialQuery.remove(positionReturn.atoms[i])
        self._deleteMes.add(positionReturn.atoms[i])
        self._infoString += _VerboseCheck(10,"Deleting atom\n")
      else:
        self._spatialQuery.add(positionReturn.atoms[i])
        self._deleteMes.discard(positionReturn.atoms[i])
        self._infoString += _VerboseCheck(10,"Ensuring deletable atom is present\n")

  def _doFixup(self, fixUp):
    """
      Move the atoms to their fixup positions, updating their extra atom info
      and deletion status
    """
    myAtoms = fixUp.atoms
    for i, p in enumerate(fixUp.positions):
      self._infoString += _VerboseCheck(5,"Moving atom to {}\n".format(p))
      myAtoms[i].xyz = p
    for i, e in enumerate(fixUp.extraInfos):
      self._infoString += _VerboseCheck(5,"Atom info to {}\n".format(e))
      self._extraAtomInfo.setMappingFor(myAtoms[i], e)
    for i, d in enumerate(fixUp.deleteMes):
      # Either ensure that it is deleted or ensure that it is not depending on the
      # value of the deletion result.
      self._infoString += _VerboseCheck(5,"Atom deleted is {}\n".format(d))
      if d:
        self._deleteMes.add(myAtoms[i])
      else:
        self._deleteMes.discard(myAtoms[i])

  def _scoreAtom(self, atom):
    maxRadiusWithoutProbe = self._extraAtomInfo.getMappingFor(atom).vdwRadius + self._maximumVDWRadius
    return self._dotScorer.score_dots(atom, self._minOccupancy, self._spatialQuery,
      maxRadiusWithoutProbe, self._probeRadius, self._excludeDict[atom], self._dotSpheres[atom].dots(),
      self._probeDensity, False).totalScore()

  def _optimizeSingleMoverCoarse(self, mover):
    # Find the coarse score for the Mover in all orientations by moving each atom into the
    # specified position and summing the scores over all of them.  Determine the best
    # orientation by selecting the highest scorer.
    # Add the preference energy to the sum for each orientation scaled by our preference
    # magnitude.
    # :return: the score for the Mover in its optimal state.
    # :side_effect: self._setMoverState() is called to put the Mover's atoms into its best state.
    # :side_effect: self._coarseLocations is set to the Mover's best state.
    # :side_effect: Changes the value of self._highScores[mover] to the score at the coarse position
    # selected
    coarse = mover.CoarsePositions()
    scores = coarse.preferenceEnergies[:]
    for i in range(len(scores)):
      scores[i] *= self._preferenceMagnitude
    for i in range(len(coarse.positions)):
      self._setMoverState(coarse, i)
      self._coarseLocations[mover] = i

      for a in coarse.atoms:
        scores[i] += self._scoreAtom(a)
      self._infoString += _VerboseCheck(5,"Single Mover score at orientation {} = {:.2f}\n".format(i,scores[i]))

    # Find the maximum score, keeping track of the best score and its index.
    maxScore = scores[0]
    maxIndex = 0
    for i in range(1,len(coarse.positions)):
      if scores[i] > maxScore:
        maxScore = scores[i]
        maxIndex = i;

    # Put the Mover into its final position (which may be back to its initial position)
    self._infoString += _VerboseCheck(3,"Setting single Mover to coarse orientation {}".format(maxIndex)+
      ", max score = {:.2f} (initial score {:.2f})\n".format(maxScore, scores[0]))
    self._setMoverState(coarse, maxIndex)
    self._coarseLocations[mover] = maxIndex

    # Record and return the best score for this Mover.
    self._highScores[mover] = maxScore
    return maxScore

  def _optimizeCliqueCoarse(self, clique):
    """
    Override this method in derived classes.
    The _SingletonOptimizer class just calls the single-Mover optimimization for each
    of the elements in the Clique and returns the vector of their results.  This should
    be overridden in derived classes to actually check for the joint maximum score over
    all of the Movers simultaneously.
    :param clique: Boost graph whose vertices contain all of the Movers to be optimized
    and whose edges contain a description of all pairwise interections between them.
    :return: the score for the Movers in their optimal state.
    :side_effect: self._setMoverState() is called to put the Movers into the best combined state.
    :side_effect: self._coarseLocations is set to the Mover's best state.
    :side_effect: self._highScores is set to the individual score for each of the Movers.
    """
    self._infoString += _VerboseCheck(3,"Optimizing clique of size {} as singletons\n".format(len(list(clique.vertices()))))
    ret = 0.0
    for v in clique.vertices():
      ret += self._optimizeSingleMoverCoarse(clique.vertex_label(v))
    return ret

  def _optimizeSingleMoverFine(self, mover):
    """
    Find the score for the Mover in all fine orientations by moving each atom into the
    specified position and summing the scores over all of them.  Determine the best
    orientation by selecting the highest scorer.
    Add the preference energy to the sum for each orientation scaled by our preference
    :return: the score for the Mover in its optimal state.
    :side effect: Changes the value of self._highScores[mover] to the score at the fine position
    selected if one is selected.
    """
    maxScore = 0.0
    coarse = mover.CoarsePositions()  # Record in case we need to put it back
    fine = mover.FinePositions(self._coarseLocations[mover])
    if len(fine.positions) > 0:
      scores = fine.preferenceEnergies[:]
      for i in range(len(scores)):
        scores[i] *= self._preferenceMagnitude
      for i in range(len(fine.positions)):
        self._setMoverState(fine, i)

        for a in fine.atoms:
          scores[i] += self._scoreAtom(a)
        self._infoString += _VerboseCheck(5,"Single Mover score at orientation {} = {:.2f}\n".format(i,scores[i]))

      # Find the maximum score, keeping track of the best score and its index.
      maxScore = scores[0]
      maxIndex = 0
      for i in range(1,len(fine.positions)):
        if scores[i] > maxScore:
          maxScore = scores[i]
          maxIndex = i;

      # Put the Mover into its final position (which may be back to its initial position)
      # and update the high score.
      if maxScore > self._highScores[mover]:
        self._infoString += _VerboseCheck(3,"Setting single Mover to fine orientation {}".format(maxIndex)+
          ", max score = {:.2f} (coarse score {:.2f})\n".format(maxScore,self._highScores[mover]))
        self._setMoverState(fine, maxIndex)
        self._fineLocations[mover] = maxIndex

        # Record the best score for this Mover.
        self._highScores[mover] = maxScore
      else:
        # Put us back to the initial coarse location and don't change the high score.
        self._infoString += _VerboseCheck(3,"Leaving single Mover at coarse orientation\n")
        self._setMoverState(coarse, self._coarseLocations[mover])
        self._fineLocations[mover] = 0
    return maxScore

class _BruteForceOptimizer(_SingletonOptimizer):
  def __init__(self, addFlipMovers, model, modelIndex, altID,
                bondedNeighborDepth,
                probeRadius,
                useNeutronDistances,
                probeDensity,
                minOccupancy,
                preferenceMagnitude
              ):
    """Constructor for _BruteForceOptimizer.  This tries all combinations of Mover positions
    within a Clique.  It will be too slow for many files, but it provides a baseline against
    which to compare the results from faster optimizers.
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
    hierarchy.
    :param altID: The conformer alternate location specifier to use.  The value "" will
    cause it to run on the first conformer found in each model.  If this is set to None,
    optimization will be run sequentially for every conformer in the model, starting with
    the last and ending with the first.  This will leave the initial conformer's values as the
    final location for atoms that are not inside a conformer or are in the first conformer.
    :param bondedNeighborDepth: How many hops to ignore bonding when doing Probe calculations.
    A depth of 3 will ignore my bonded neighbors and their bonded
    neighbors and their bonded neighbors.
    :param probeRadius: Radius of the probe to be used in Probe calculations (Angstroms).
    :param useNeutronDistances: False will use X-ray/electron cloud distances.  If set to
    True, it will use neutron (nuclear) distances.  This must be set consistently with the
    values used to generate the hydrogens and to run PDB interpretation.
    :param probeDensity: How many dots per sq Angstroms in VDW calculations.
    :param minOccupancy: Minimum occupancy for an atom to be considered in the Probe score.
    :param preferenceMagnitude: Multiplier for the preference energies expressed
    by some Movers for particular orientations.
    """
    super(_BruteForceOptimizer, self).__init__(addFlipMovers, model, modelIndex = modelIndex, altID = altID,
                bondedNeighborDepth = bondedNeighborDepth,
                probeRadius = probeRadius, useNeutronDistances = useNeutronDistances, probeDensity = probeDensity,
                minOccupancy = minOccupancy, preferenceMagnitude = preferenceMagnitude)

  def _optimizeCliqueCoarse(self, clique):
    """
    The _BruteForceOptimizer class checks for the joint maximum score over
    all of the Movers simultaneously.  It tries all Movers in all possible positions against all
    other Movers in all combinations of positions.
    :param clique: Boost graph whose vertices contain all of the Movers to be optimized
    and whose edges contain a description of all pairwise interections between them.
    :return: the score for the Movers in their optimal state.
    :side_effect: self._setMoverState() is called to put the Movers into the best combined state.
    :side_effect: self._coarseLocations is set to the Mover's best state.
    :side_effect: self._highScores is set to the individual score for each of the Movers.
    """
    self._infoString += _VerboseCheck(3,"Optimizing clique of size {} using brute force\n".format(len(list(clique.vertices()))))

    # Find the value for the current set of states, compare it against the max, and store it if
    # it is the best so far.
    # We will cycle the states[] list through all possible states for each Mover;
    # use the generating function to cycle through all possible states, keeping track of the best.
    bestState = None
    bestScore = -1e100  # Any score will be better than this
    states = []         # Coarse position state return for each Mover
    numStates = []      # Number of positions for each Mover
    movers = [self._interactionGraph.vertex_label(v) for v in clique.vertices()]
    for m in movers:
      states.append(m.CoarsePositions())
      numStates.append(len(states[-1].positions))
    for curStateValues in _generateAllStates(numStates):

      # Set all movers to match the state list.
      # @todo Optimize this so that it only changes states that differed from last time.
      for i in range(len(movers)):
        self._setMoverState(states[i], curStateValues[i])
        self._coarseLocations[movers[i]] = curStateValues[i]

      # Compute the score over all atoms in all Movers and see if it is the best.  If so,
      # update the best.
      score = 0
      for i in range(len(movers)):
        for a in states[i].atoms:
          score += self._scoreAtom(a)
      self._infoString += _VerboseCheck(5,"Score is {:.2f} at {}\n".format(score, curStateValues))
      if score > bestScore or bestState is None:
        self._infoString += _VerboseCheck(4,"New best score is {:.2f} at {}\n".format(score, curStateValues))
        bestScore = score
        bestState = curStateValues[:]

    # Put each Mover into its best state and compute its high-score value.
    # Compute the best individual scores for these Movers for use in later fine-motion
    # processing.  Return the total score.
    ret = 0.0
    for i,m in enumerate(movers):
      self._setMoverState(states[i], bestState[i])
      self._coarseLocations[movers[i]] = bestState[i]
      self._highScores[m] = 0
      for a in states[i].atoms:
        self._highScores[m] += self._scoreAtom(a)
      self._infoString += _VerboseCheck(3,"Setting Mover in clique to coarse orientation {}".format(bestState[i])+
        ", max score = {:.2f}\n".format(self._highScores[m]))
      ret += self._highScores[m]
    return ret

class _CliqueOptimizer(_BruteForceOptimizer):
  def __init__(self, addFlipMovers, model, modelIndex, altID,
                bondedNeighborDepth,
                probeRadius,
                useNeutronDistances,
                probeDensity,
                minOccupancy,
                preferenceMagnitude
              ):
    """Constructor for _CliqueOptimizer.  This uses a recursive algorithm to break down the total
    clique into sets of smaller cliques.  It looks for a vertex cut in the Clique it is called with
    that will separate the remaining vertices into two more more connected subcomponents.  It then tests
    each combined state of the set of Movers in the vertex cut to find the one with the best overall
    maximum score.  For each state, it first recursively optimizes all of the connected subcomponents
    and then (with each of the subcomponents in its optimal state) computes the score for the Movers in
    the vertex cut.  Recursion terminates when there are two or fewer Movers in the Clique or when no
    vertex cut can be found; the parent-class Clique solver is used in these cases.
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
    hierarchy.
    :param altID: The conformer alternate location specifier to use.  The value "" will
    cause it to run on the first conformer found in each model.  If this is set to None,
    optimization will be run sequentially for every conformer in the model, starting with
    the last and ending with the first.  This will leave the initial conformer's values as the
    final location for atoms that are not inside a conformer or are in the first conformer.
    :param bondedNeighborDepth: How many hops to ignore bonding when doing Probe calculations.
    A depth of 3 will ignore my bonded neighbors and their bonded
    neighbors and their bonded neighbors.
    :param probeRadius: Radius of the probe to be used in Probe calculations (Angstroms).
    :param useNeutronDistances: False will use X-ray/electron cloud distances.  If set to
    True, it will use neutron (nuclear) distances.  This must be set consistently with the
    values used to generate the hydrogens and to run PDB interpretation.
    :param probeDensity: How many dots per sq Angstroms in VDW calculations.
    :param minOccupancy: Minimum occupancy for an atom to be considered in the Probe score.
    :param preferenceMagnitude: Multiplier for the preference energies expressed
    by some Movers for particular orientations.
    """
    super(_CliqueOptimizer, self).__init__(addFlipMovers, model, modelIndex = modelIndex, altID = altID,
                bondedNeighborDepth = bondedNeighborDepth,
                probeRadius = probeRadius, useNeutronDistances = useNeutronDistances, probeDensity = probeDensity,
                minOccupancy = minOccupancy, preferenceMagnitude = preferenceMagnitude)

  def _optimizeCliqueCoarse(self, clique):
    """
    Looks for a vertex cut in the Clique that will separate the remaining vertices into two more
    more connected subcomponents.  Test each combined state of the set of Movers in the vertex cut
    to find the one with the best overall maximum score.
      For each state, recursively optimize all of the connected subcomponents and then (with each
    of the subcomponents in its optimal state) compute the score for the Movers in the vertex cut.
      Recursion terminates when there are two or fewer Movers in the Clique or when no
    vertex cut can be found; the parent-class Clique solver is used in these cases.
    :param clique: Boost graph whose vertices contain all of the Movers to be optimized
    and whose edges contain a description of all pairwise interections between them.
    :return: the score for the Movers in their optimal state.
    :side_effect: self._setMoverState() is called to put the Movers into the best combined state.
    :side_effect: self._coarseLocations is set to the Mover's best state.
    :side_effect: self._highScores is set to the individual score for each of the Movers.
    """
    self._infoString += _VerboseCheck(3,"Optimizing clique of size {} using recursion\n".format(len(list(clique.vertices()))))

    # If we've gotten down to a clique of size 2, we terminate recursion and call our parent's method
    # because we can never split this into two connected components.
    if len(list(clique.vertices())) <= 2:
      self._infoString += _VerboseCheck(3,"Recursion terminated at clique of size {}\n".format(len(list(clique.vertices()))))
      ret = super(_CliqueOptimizer, self)._optimizeCliqueCoarse(clique)
      return ret

    # Find all of the Movers in the clique so that we know which ones to capture the state for.
    movers = [self._interactionGraph.vertex_label(v) for v in clique.vertices()]
    states = {}         # Coarse position state return for each Mover in the entire Clique
    for m in movers:
      states[m]= m.CoarsePositions()

    # Find a vertex cut for the graph we were given.  Run through all states of the vertex cut
    # and for each recursively find the best score for all of the connected components, followed by the score
    # for the vertex cut.  Keep track of the best state and score across all of them and set back
    # to that at the end.  If we have no Movers in the vertex cut, none was found so we don't recur.
    cutMovers, cutGraph = _vertexCut(clique)
    if len(cutMovers) > 0:
      self._infoString += _VerboseCheck(3,"Found vertex cut of size {}\n".format(len(cutMovers)))

      score = 0.0
      bestState = None
      bestScore = -1e100
      numStates = []      # Number of positions for each Mover
      for m in cutMovers:
        numStates.append(len(states[m].positions))
      for curStateValues in _generateAllStates(numStates):
        # Set all cutMovers to match the state list.
        # @todo Optimize this so that it only changes states that differed from last time.
        for i,m in enumerate(cutMovers):
          self._setMoverState(states[m], curStateValues[i])
          self._coarseLocations[m] = curStateValues[i]

        # Recursively compute the best score across all connected components in the cutGraph.
        # This will leave each subgraph in its best state for this set of cutMovers states.
        score = 0
        components = cca.connected_components( graph = cutGraph )
        for g in components:
          subMovers = [cutGraph.vertex_label(i) for i in g]
          subset = _subsetGraph(cutGraph, subMovers)
          score += self._optimizeCliqueCoarse(subset)

        # Add the score over all atoms in the vertex-cut Movers and see if it is the best.  If so,
        # update the best.
        for m in cutMovers:
          for a in states[m].atoms:
            score += self._scoreAtom(a)
        self._infoString += _VerboseCheck(5,"Score is {:.2f} at {}\n".format(score, curStateValues))
        if score > bestScore or bestState is None:
          self._infoString += _VerboseCheck(4,"New best score is {:.2f} at {}\n".format(score, curStateValues))
          bestScore = score
          # Get the current state for all Movers in the Clique, not just the vertex-cut Movers
          bestState = [self._coarseLocations[m] for m in movers]

      # Put each Mover in the entire Clique into its best state and compute its high-score value.
      # Compute the best individual scores for these Movers for use in later fine-motion
      # processing.  Return the total score.
      ret = 0.0
      for i,m in enumerate(movers):
        self._setMoverState(states[m], bestState[i])
        self._coarseLocations[m] = bestState[i]
        self._highScores[m] = 0
        for a in states[m].atoms:
          self._highScores[m] += self._scoreAtom(a)
        self._infoString += _VerboseCheck(3,"Setting Mover in clique to coarse orientation {}".format(bestState[i])+
          ", max score = {:.2f}\n".format(self._highScores[m]))
        ret += self._highScores[m]
      return ret

    # Give up and use our parent's method.
    self._infoString += _VerboseCheck(3,"No vertex cut for clique of size {}, calling parent\n".format(len(list(clique.vertices()))))
    ret = super(_CliqueOptimizer, self)._optimizeCliqueCoarse(clique)
    return ret

class FastOptimizer(_CliqueOptimizer):
  def __init__(self, addFlipMovers, model, modelIndex = 0, altID = None,
                bondedNeighborDepth = 3,
                probeRadius = 0.25,
                useNeutronDistances = False,
                probeDensity = 16.0,
                minOccupancy = 0.02,
                preferenceMagnitude = 1.0
              ):
    """Constructor for FastOptimizer.  This uses the same algorithm as the
    parent-class but first constructs a cache for every atom in every Mover
    of all the Movers whose positions can affect its answer.  The _scoreAtom() method
    is overridden to use this cached value in clique optimization (but not in singleton
    or fine optimization) when it has already been computed for a given configuration
    of Movers.
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
    hierarchy.
    :param altID: The conformer alternate location specifier to use.  The value "" will
    cause it to run on the first conformer found in each model.  If this is set to None
    (the default), optimization will be run sequentially for every conformer in the model, starting with
    the last and ending with the first.  This will leave the initial conformer's values as the
    final location for atoms that are not inside a conformer or are in the first conformer.
    :param bondedNeighborDepth: How many hops to ignore bonding when doing Probe calculations.
    The default is to ignore interactions to a depth of 3 (my bonded neighbors and their bonded
    neighbors and their bonded neighbors).
    :param probeRadius: Radius of the probe to be used in Probe calculations (Angstroms).
    :param useNeutronDistances: Defaults to using X-ray/electron cloud distances.  If set to
    True, it will use neutron (nuclear) distances.  This must be set consistently with the
    values used to generate the hydrogens and to run PDB interpretation.
    :param probeDensity: How many dots per sq Angstroms in VDW calculations.
    :param minOccupancy: Minimum occupancy for an atom to be considered in the Probe score.
    :param preferenceMagnitude: Multiplier for the preference energies expressed
    by some Movers for particular orientations.
    """
    # Set a flag that will be used by the overridden _scoreAtom() method to determine whether it
    # should be doing caching or not.  Initially, it should not be -- it should only be doing this
    # when we're inside a Clique optimization in the overridden _optimizeCliqueCoarse() method.
    # We do this before constructing the parent class because it is needed by the functions that
    # it calls.
    self._doScoreCaching = False

    super(FastOptimizer, self).__init__(addFlipMovers, model, modelIndex = modelIndex, altID = altID,
                bondedNeighborDepth = bondedNeighborDepth,
                probeRadius = probeRadius, useNeutronDistances = useNeutronDistances, probeDensity = probeDensity,
                minOccupancy = minOccupancy, preferenceMagnitude = preferenceMagnitude)

  def _initializeAlternate(self):
    # Ensure that we have a per-atom _scoreCache dictionary that will store already-computed
    # results for a given atom based on the configurations of the Movers that can affect its
    # results.  The entries will be empty to start with and will be filled in as they are computed.
    # We build entries for all atoms, even those not in the Movers to avoid having to traverse
    # the Movers.
    # This structure is a dictionary (looked up by atom) of dictionaries (looked up by tuple)
    # of values (scores).
    self._scoreCache = {}
    for a in self._atoms:
      self._scoreCache[a] = {}
    return

  def _scoreAtom(self, atom):

    if self._doScoreCaching:
      # Construct a tuple that holds the entries for the coarse position of all of the Movers
      # that this atom depends on, using the _atomMoverSets to determine which ones to look up.
      # See if this result is already in the dictionary for that atom.  If so, use it.  If not,
      # compute and store it and then return that value.
      state = tuple([self._coarseLocations[m] for m in self._atomMoverSets[atom]])
      try:
        return self._scoreCache[atom][state]
      except Exception:
        self._scoreCache[atom][state] = super(FastOptimizer, self)._scoreAtom(atom)
        return self._scoreCache[atom][state]
    else:
      return super(FastOptimizer, self)._scoreAtom(atom)

  def _optimizeCliqueCoarse(self, clique):
    """
    The FastOptimizer class generates a per-atom score cache object and uses it along
    with an overridden _scoreAtom() method to avoid recomputing scores for atoms where
    there has been no change in any of the Movers they depend on.
    It wraps the parent-class method after setting things up to use the cache, and
    then turns off the cache before returning.
    :param clique: Boost graph whose vertices contain all of the Movers to be optimized
    and whose edges contain a description of all pairwise interections between them.
    :return: the score for the Movers in their optimal state.
    :side_effect: self._setMoverState() is called to put the Movers into the best combined state.
    :side_effect: self._coarseLocations is set to the Mover's best state.
    :side_effect: self._highScores is set to the individual score for each of the Movers.
    """
    self._infoString += _VerboseCheck(3,"Optimizing clique of size {} using atom-score cache\n".format(len(list(clique.vertices()))))

    # Call the parent-class optimizer, turning on and off the cache behavior before
    # and after.
    self._doScoreCaching = True
    ret = super(FastOptimizer, self)._optimizeCliqueCoarse(clique)
    self._doScoreCaching = False
    return ret

##################################################################################
# Placement

class _PlaceMoversReturn(object):
  # Return type from PlaceMovers() call.  List of movers and then an information string
  # that may contain information the user would like to know (where Movers were placed,
  # failed Mover placements due to missing Hydrogens, etc.).  Also returns a list of
  # atoms that should be deleted as a result of situations determined during the
  # placement.  Also returns a dictionary mapping from Movers to strings describing each Mover.
  # @todo Make this a method of the base Optimizer so we're not passing things in and
  # out and immediately reassigning them to member variables.
  def __init__(self, moverList, infoString, deleteAtoms, moverInfo):
    self.moverList = moverList
    self.infoString = infoString
    self.deleteAtoms = deleteAtoms
    self.moverInfo = moverInfo

def _PlaceMovers(atoms, rotatableHydrogenIDs, bondedNeighborLists, hParameters, spatialQuery,
                  extraAtomInfo, maxVDWRadius, addFlipMovers):
  """Produce a list of Movers for atoms in a pdb.hierarchy.conformer that has added Hydrogens.
  :param atoms: flex array of atoms to search.  This must have all Hydrogens needed by the
  Movers present in the structure already.
  :param rotateableHydrogenIDs: List of sequence IDs for single hydrogens that are rotatable.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  :param hParameters: List indexed by sequence ID that stores the riding
  coefficients for hydrogens that have associated dihedral angles.  This can be
  obtained by calling model.setup_riding_h_manager() and then model.get_riding_h_manager().
  :param spatialQuery: Probe.SpatialQuery structure to rapidly determine which atoms
  are within a specified distance of a location.
  :param extraAtomInfo: Probe.ExtraAtomInfo mapper that provides radius and other
  information about atoms beyond what is in the pdb.hierarchy.  Used here to determine
  which atoms may be acceptors.
  :param maxVDWRadius: Maximum VdW radius of the atoms.
  :param addFlipMovers: Do we add flip Movers along with other types?
  :return: _PlaceMoversReturn giving the list of Movers found in the conformation and
  an error string that is empty if no errors are found during the process and which
  has a printable message in case one or more errors are found.
  """

  # List of Movers to return
  movers = []

  # List of atoms to delete
  deleteAtoms = []

  # Information string to return, filled in when we cannot place a Mover.
  infoString = ""

  # Dictionary mapping from Mover to information about the Mover, so that we can keep track of
  # everything from its placement to its final state in the same string.
  moverInfo = {}

  # The radius of a Polar Hydrogen (one that forms Hydrogen bonds)
  # @todo Can we look this up from inside CCTBX?
  polarHydrogenRadius = 1.05

  # For each atom, check to see if it is the main atom from a Mover.  If so, attempt to
  # construct the appropriate Mover.  The construction may fail because not all required
  # atoms are present (especially the Hydrogens).  If the construction succeeds, add the
  # Mover to the list of those to return.  When they fail, add the output to the information
  # string.

  for a in atoms:
    # Find the stripped upper-case atom and residue names and the residue ID,
    # and an identifier string
    aName = a.name.strip().upper()
    resName = a.parent().resname.strip().upper()
    resNameAndID = _ResNameAndID(a)

    # See if we should construct a MoverSingleHydrogenRotator here.
    # @todo This is placing on atoms C and N in CYS 352; and on HE2 and CA in SER 500 of 4z4d
    if a.i_seq in rotatableHydrogenIDs:
      try:
        # Skip Hydrogens that are not bound to an atom that only has a single other
        # atom bound to it -- we will handle those in other cases.
        neighbor = bondedNeighborLists[a][0]
        if len(bondedNeighborLists[neighbor]) == 2:

          # Construct a list of nearby atoms that are potential acceptors.
          potentialAcceptors = []

          # Get the list of nearby atoms.  The center of the search is the atom that
          # the Hydrogen is bound to and its radius is 4 (these values are pulled from
          # the Reduce C++ code).
          maxDist = 4.0
          nearby = spatialQuery.neighbors(neighbor.xyz, extraAtomInfo.getMappingFor(neighbor).vdwRadius, maxDist)

          # Check each nearby atom to see if it distance from the neighbor is within
          # the sum of the hydrogen-bond length of the neighbor atom, the radius of
          # a polar Hydrogen, and the radius of the nearby atom, indicating potential
          # overlap.
          # O-H & N-H bond len == 1.0, S-H bond len == 1.3
          XHbondlen = 1.0
          if neighbor.element == "S":
            XHbondlen = 1.3
          candidates = []
          for n in nearby:
            d = (Helpers.rvec3(neighbor.xyz) - Helpers.rvec3(n.xyz)).length()
            if d <= XHbondlen + extraAtomInfo.getMappingFor(n).vdwRadius + polarHydrogenRadius:
              candidates.append(n)

          # See if each nearby atom is a potential acceptor or a flip partner from
          # Histidine or NH2 flips (which may be moved into that atom's position during a flip).
          # We check the partner (N's for GLN And ASN, C's for HIS) because if it is the O or
          # N atom, we will already be checking it as an acceptor now.
          for c in candidates:
            cName = c.name.strip().upper()
            resName = c.parent().resname.strip().upper()
            flipPartner = (
              (cName == 'ND2' and resName == 'ASN') or
              (cName == 'NE2' and resName == 'GLN') or
              (cName == 'CE1' and resName == 'HIS') or (cName == 'CD2' and resName == 'HIS') )
            acceptor = extraAtomInfo.getMappingFor(c).isAcceptor
            if acceptor or flipPartner:
              potentialAcceptors.append(c)

          movers.append(Movers.MoverSingleHydrogenRotator(a, bondedNeighborLists, hParameters, potentialAcceptors))
          infoString += _VerboseCheck(1,"Added MoverSingleHydrogenRotator "+str(len(movers))+" to "+resNameAndID+" "+aName+
            " with "+str(len(potentialAcceptors))+" potential nearby acceptors\n")
          moverInfo[movers[-1]] = "SingleHydrogenRotator at "+resNameAndID+" "+aName;
      except Exception as e:
        infoString += _VerboseCheck(0,"Could not add MoverSingleHydrogenRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")

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
          movers.append(Movers.MoverNH3Rotator(a, bondedNeighborLists, hParameters))
          infoString += _VerboseCheck(1,"Added MoverNH3Rotator "+str(len(movers))+" to "+resNameAndID+"\n")
          moverInfo[movers[-1]] = "NH3Rotator at "+resNameAndID+" "+aName;
        except Exception as e:
          infoString += _VerboseCheck(0,"Could not add MoverNH3Rotator to "+resNameAndID+": "+str(e)+"\n")

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
        # so that they Hydrogens will be staggered but do not add it to those to be optimized.
        if len(bondedNeighborLists[neighbor]) == 3:
          try:
            movers.append(Movers.MoverAromaticMethylRotator(a, bondedNeighborLists, hParameters))
            infoString += _VerboseCheck(1,"Added MoverAromaticMethylRotator "+str(len(movers))+" to "+resNameAndID+" "+aName+"\n")
            moverInfo[movers[-1]] = "AromaticMethylRotator at "+resNameAndID+" "+aName;
          except Exception as e:
            infoString += _VerboseCheck(0,"Could not add MoverAromaticMethylRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")
        else:
          try:
            ignored = Movers.MoverTetrahedralMethylRotator(a, bondedNeighborLists, hParameters)
            infoString += _VerboseCheck(1,"Used MoverTetrahedralMethylRotator to stagger "+resNameAndID+" "+aName+"\n")
          except Exception as e:
            infoString += _VerboseCheck(0,"Could not add MoverTetrahedralMethylRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n")

    # See if we should insert a MoverAmideFlip here.
    # @todo Is there a more general way than looking for specific names?
    if addFlipMovers and ((aName == 'XD2' and resName == 'ASX') or (aName == 'XE2' and resName == 'GLX')):
      infoString += _VerboseCheck(1,"Not attempting to adjust "+resNameAndID+" "+aName+"\n")
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
        myRad = extraAtomInfo.getMappingFor(a).vdwRadius
        minDist = myRad
        maxDist = 0.25 + myRad + maxVDWRadius
        neighbors = spatialQuery.neighbors(oxygen.xyz, minDist, maxDist)
        for n in neighbors:
          if n.element_is_positive_ion():
            dist = (Helpers.rvec3(oxygen.xyz) - Helpers.rvec3(n.xyz)).length()
            expected = myRad + extraAtomInfo.getMappingFor(n).vdwRadius
            infoString += _VerboseCheck(5,'Checking AmideFlip '+str(i)+' against '+n.name.strip()+' at '+str(n.xyz)+' from '+str(pos)+
              ' dist = '+str(dist)+', expected = '+str(expected)+'; N rad = '+str(myRad)+
              ', '+n.name.strip()+' rad = '+str(extraAtomInfo.getMappingFor(n).vdwRadius)+'\n')
            # @todo Why are we using -0.65 here and -0.55 for Histidine?
            if dist >= (expected - 0.65) and dist <= (expected + 0.25):
              foundIon = True

      if not foundIon:
        try:
          movers.append(Movers.MoverAmideFlip(a, "CA", bondedNeighborLists))
          infoString += _VerboseCheck(1,"Added MoverAmideFlip "+str(len(movers))+" to "+resNameAndID+"\n")
          moverInfo[movers[-1]] = "AmideFlip at "+resNameAndID+" "+aName;
        except Exception as e:
          infoString += _VerboseCheck(0,"Could not add MoverAmideFlip to "+resNameAndID+": "+str(e)+"\n")

    # See if we should insert a MoverHisFlip here.
    # @todo Is there a more general way than looking for specific names?
    if aName == 'NE2' and resName == 'HIS':
      try:
        # Get a potential Mover and test both of its Nitrogens in the original and flipped
        # locations.  If one or both of them are near enough to be ionically bound to an
        # ion, then we remove the Hydrogen(s) and lock the Histidine at that orientation
        # rather than inserting the Mover into the list of those to be optimized.
        # @todo Consider checking both configurations to see if either one has two bonds.
        hist = Movers.MoverHisFlip(a, bondedNeighborLists, extraAtomInfo)

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
        myRad = extraAtomInfo.getMappingFor(a).vdwRadius
        minDist = myRad
        maxDist = 0.25 + myRad + maxVDWRadius
        bondedConfig = None
        for i,pos in enumerate([ne2Orig, nd1Orig, ne2Flip, nd1Flip]):
          neighbors = spatialQuery.neighbors(pos, minDist, maxDist)
          for n in neighbors:
            if n.element_is_positive_ion():
              dist = (Helpers.rvec3(pos) - Helpers.rvec3(n.xyz)).length()
              expected = myRad + extraAtomInfo.getMappingFor(n).vdwRadius
              infoString += _VerboseCheck(5,'Checking Histidine '+str(i)+' against '+n.name.strip()+' at '+str(n.xyz)+' from '+str(pos)+
                ' dist = '+str(dist)+', expected = '+str(expected)+'; N rad = '+str(myRad)+
                ', '+n.name.strip()+' rad = '+str(extraAtomInfo.getMappingFor(n).vdwRadius)+'\n')
              if dist >= (expected - 0.55) and dist <= (expected + 0.25):
                # The first two elements in the enumeration come from Histidine configuration 0
                # and the second two from configuration 4, so we figure out which one we should
                # be in.
                bondedConfig = (i // 2) * 4
                infoString += _VerboseCheck(5,'Found ionic bond in coarse configuration '+str(bondedConfig)+'\n')
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
              extraAtomInfo.setMappingFor(a, fixUp.extraInfos[i])

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
            myRad = extraAtomInfo.getMappingFor(nitro).vdwRadius
            minDist = myRad
            maxDist = 0.25 + myRad + maxVDWRadius
            neighbors = spatialQuery.neighbors(coarseNitroPos, minDist, maxDist)
            for n in neighbors:
              if n.element_is_positive_ion():
                dist = (Helpers.rvec3(coarseNitroPos) - Helpers.rvec3(n.xyz)).length()
                expected = myRad + extraAtomInfo.getMappingFor(n).vdwRadius
                if dist >= (expected - 0.55) and dist <= (expected + 0.25):
                  result += _VerboseCheck(1,'Not adding Hydrogen to '+resNameAndID+nitro.name+' and marking as an acceptor '+
                    '(ionic bond to '+n.name.strip()+')\n')
                  extra = extraAtomInfo.getMappingFor(nitro)
                  extra.isAcceptor = True
                  extraAtomInfo.setMappingFor(nitro, extra)
                  deleteAtoms.append(hydro)
                  break
            return result

          infoString += _modifyIfNeeded(fixUp.atoms[0], coarsePositions[0], fixUp.atoms[1])
          infoString += _modifyIfNeeded(fixUp.atoms[4], coarsePositions[4], fixUp.atoms[5])

          infoString += _VerboseCheck(1,"Set MoverHisFlip on "+resNameAndID+" to state "+str(bondedConfig)+"\n")
        elif addFlipMovers: # Add a Histidine flip Mover if we're adding flip Movers
          movers.append(hist)
          infoString += _VerboseCheck(1,"Added MoverHisFlip "+str(len(movers))+" to "+resNameAndID+"\n")
          moverInfo[movers[-1]] = "HisFlip at "+resNameAndID+" "+aName;
      except Exception as e:
        infoString += _VerboseCheck(0,"Could not add MoverHisFlip to "+resNameAndID+": "+str(e)+"\n")

  return _PlaceMoversReturn(movers, infoString, deleteAtoms, moverInfo)

##################################################################################
# Private helper functions.

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

def _vertexCut(clique):
  """
  Return a "vertex cut" of clique such that removal of a subset of the vertices
  causes the remainder of the vertices to become disconnected.
  :param clique: Boost graph of Movers that describes a Clique (or a subset of a Clique) and
  includes only vertices and edges within the subset for which a vertex cut is to be found.
  This must not be self._interactionGraph() or any other graph that already has multiple
  connected components.
  :return: A tuple: (1) List of Movers that were removed.  If this list is empty, no vertex
  cut was found. (2) New graph with the vertices and edges associated with the removed Movers
  removed from the graph.  This graph will have at least two connected components.  If no
  vertex cut is found, this will be a copy of the original graph.
  """

  # Check all vertex cut sizes from 1 to 2 less than the number of vertices (we must)
  # have at least 2 vertices left to have a disconnected graph).
  for n in range(1, len(list(clique.vertices()))+1-2):
    # Iterate over all sets of vertices of size n that might be removed
    for removed in itertools.combinations(clique.vertices(),n):
      movers = [clique.vertex_label(v) for v in removed]
      # Find the list of remaining vertices, which is all but the ones to be removed.
      remain = set(clique.vertex_label(v) for v in clique.vertices())
      for m in movers:
        remain.remove(m)
      # Make a copy of the clique with those vertices removed.
      newGraph = _subsetGraph(clique, remain)

      # If the graph has 2 or more connected components, we've found our answer.
      components = cca.connected_components( graph = newGraph )
      if len(components) > 1:
        return movers, newGraph

  # We didn't find an answer.  Return a complete copy of the graph and an empty set of movers.
  movers = []
  newGraph = _subsetGraph(clique, [clique.vertex_label(v) for v in clique.vertices()])

  return movers, newGraph

def _generateAllStates(numStates):
  """ This is a generator function that will yield all combinations of states within the
  counts specified in the state vector.  For example, if numStates is [2, 2] then this will
  produce the following outputs, one at a time: [0, 0], [1, 0], [0, 1], [1, 1].
  :param numStates: List of integers specifying the number of states for each element.
  :return: Yields all combinations of all elements in the state vector.
  """

  # Cycle the curStateValues[] list through all possible states for each element,
  # incremementing each until it rolls over and then jumping up to the next.
  # This is similar to doing +1 arithmetic with carry on a multi-digit number.
  curStateValues = [0] * len(numStates) # Cycles through available counts, one per element
  curState = 0

  # Increment the state.  We do this by increasing the current element until it reaches its
  # number of values, then we bump it and all of its neighbors to the right back to 0 and, if
  # we're not at the left end, bump the next one up.  If we are at the left end, we're done.
  # When done, leave all states at 0.
  while True:
    # Yield the current value
    yield curStateValues

    # Go to the next value, if there is one.
    curStateValues[curState] += 1
    rippled = False
    while curStateValues[curState] == numStates[curState]:  # Ripple to the left if we overflow more than one
      # Clear all of the values to the right, and ours, because we're rolling over
      for i in range(curState+1):
        curStateValues[i] = 0
      # If we're the left-most state, we're done
      if curState+1 >= len(numStates):
        return
      else:
        curState += 1
        curStateValues[curState] += 1
        rippled = True
    # If we rippled, bump back to the right-most column and start counting there again in
    # the next iteration.
    if rippled:
      curState = 0

##################################################################################
# Test function and associated data and helpers to verify that all functions behave properly.

def Test(inFileName = None, dumpAtoms = False):
  """Test function for all functions provided above.
  :param inFileName: Name of a PDB or CIF file to load (default makes a small molecule)
  :return: Empty string on success, string describing the problem on failure.
  """

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
  class philLike:
    def __init__(self, useImplicitHydrogenDistances = False):
      self.implicit_hydrogens = useImplicitHydrogenDistances
      self.set_polar_hydrogen_radius = True
  global probePhil
  probePhil = philLike(False)
  ret = Helpers.getExtraAtomInfo(model = model, bondedNeighborLists = bondedNeighborLists,
      useNeutronDistances=False,probePhil=probePhil)
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

  # Get the spatial-query information needed to quickly determine which atoms are nearby
  sq = probeExt.SpatialQuery(atoms)

  # Place the movers, which should be none because the Histidine flip
  # will be constrained by the ionic bonds.
  ret = _PlaceMovers(atoms, model.rotatable_hd_selection(iselection=True),
                     bondedNeighborLists, None, sq, extra, maxVDWRad, True)
  movers = ret.moverList
  if len(movers) != 0:
    return "Optimizers.Test(): Incorrect number of Movers for 1xso Histidine test: " + str(len(movers))

  # Make sure that the two ring Nitrogens have been marked as acceptors.
  # Make sure that the two hydrogens have been marked for deletion.
  for a in model.get_hierarchy().models()[0].atoms():
    name = a.name.strip()
    if name in ["ND1", "NE2"]:
      if not extra.getMappingFor(a).isAcceptor:
        return 'Optimizers.Test(): '+name+' in 1xso Histidine test was not an acceptor'
    if name in ["HD1", "HE2"]:
      if not a in ret.deleteAtoms:
        return 'Optimizers.Test(): '+name+' in 1xso Histidine test was not set for deletion'

  #========================================================================
  # Check a case where a HisFlip would be flipped and locked down and have its Hydrogen removed.
  # @todo

  #========================================================================
  # Check a case where an AmideFlip would be locked down and have its Hydrogen removed.
  # @todo

  #========================================================================
  # Check that the occupancy and B-factor cut-offs for water Oxygens are causing them
  # to be ignored in the calculations by putting an atom in the way and making sure it
  # is ignored.
  # @todo

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

  # Run a fast optimizer on the model.
  print('Testing FastOptimizer')
  opt = FastOptimizer(True, model,probeRadius=0.25)

  # Write debugging output if we've been asked to
  if dumpAtoms:
    f = open("deleteme.pdb","w")
    f.write(model.model_as_pdb())
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
