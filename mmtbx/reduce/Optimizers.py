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

import argparse

import Movers
import InteractionGraph

from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
import mmtbx
from scitbx.array_family import flex

from mmtbx.probe import AtomTypes
from mmtbx.probe import Helpers
import mmtbx_probe_ext as probeExt

# To enable addition of Hydrogens
# @todo See if we can remove the shift and box once reduce_hydrogen is complete
from cctbx.maptbx.box import shift_and_box_model
from mmtbx.hydrogens import reduce_hydrogen

##################################################################################
# This is a set of functions that implement placement and optimization of
# Reduce's "Movers".

class PlaceMoverstReturn:
  # Return type from PlaceMovers() call.  List of movers and then an information string
  # that may contain information the user would like to know (where Movers were placed,
  # failed Mover placements due to missing Hydrogens, etc.).
  def __init__(self, moverList, infoString):
    self.moverList = moverList
    self.infoString = infoString

def PlaceMovers(atoms, rotatableHydrogenIDs, bondedNeighborLists, spatialQuery, extraAtomInfo):
  """Produce a list of Movers for atoms in a pdb.hierarchy.conformer that has added Hydrogens.
  :param atoms: flex array of atoms to search.  This must have all Hydrogens needed by the
  Movers present in the structure already.
  :param rotateableHydrogenIDs: List of sequence IDs for single hydrogens that are rotatable.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  :param spatialQuery: Probe.SpatialQuery structure to rapidly determine which atoms
  are within a specified distance of a location.
  :param extraAtomInfo: Probe.ExtraAtomInfo structure that provides radius and other
  information about atoms beyond what is in the pdb.hierarchy.  Used here to determine
  which atoms may be acceptors.
  :returns PlaceMoverstReturn giving the list of Movers found in the conformation and
  an error string that is empty if no errors are found during the process and which
  has a printable message in case one or more errors are found.
  """

  # List of Movers to return
  movers = []

  # Information string to return, filled in when we cannot place a Mover.
  infoString = ""

  # The radius of a Polar Hydrogen (one that forms Hydrogen bonds)
  # @todo Can we look this up from inside CCTBX?
  polarHydrogenRadius = 1.05

  # For each atom, check to see if it is the main atom from a Mover.  If so, attempt to
  # construct the appropriate Mover.  The construction may fail because not all required
  # atoms are present (especially the Hydrogens).  If the construction succeeds, add the
  # Mover to the list of those to return.  When they fail, add the output to the information
  # string.
  print('XXX Number of atoms =',len(atoms))

  for a in atoms:
    # Find the stripped upper-case atom and residue names and the residue ID,
    # and an identifier string
    aName = a.name.strip().upper()
    resName = a.parent().resname.strip().upper()
    resID = str(a.parent().parent().resseq_as_int())
    resNameAndID = resName+" "+resID

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
          nearby = spatialQuery.neighbors(neighbor.xyz, extraAtomInfo[neighbor.i_seq].vdwRadius, maxDist)

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
            d = (Movers._rvec3(neighbor.xyz) - Movers._rvec3(n.xyz)).length()
            if d <= XHbondlen + extraAtomInfo[n.i_seq].vdwRadius + polarHydrogenRadius:
              candidates.append(n)

          # See if each nearby atom is a potential acceptor or a flip partner from
          # Histidine or NH2 flips (which may be moved into that atom's position during a flip).
          # We check the partner (N's for GLN And ASN, C's for HIS) because if it is the O or
          # N atom, we will already be checking it as an acceptor now.
          # @todo Ensure that Hydrogenate does not change the acceptor state of the HIS Nitrogens when it adds
          # Hydrogen to them so that they will show up here.
          for c in candidates:
            aName = c.name.strip().upper()
            resName = c.parent().resname.strip().upper()
            flipPartner = ( (aName == 'ND2' and resName == 'ASN') or (aName == 'NE2' and resName == 'GLN') or
              (aName == 'CE1' and resName == 'HIS') or (aName == 'CD2' and resName == 'HIS') )
            acceptor = extraAtomInfo[c.i_seq].isDonor
            if acceptor or flipPartner:
              potentialAcceptors.append(c)

          movers.append(Movers.MoverSingleHydrogenRotator(a, bondedNeighborLists, potentialAcceptors))
          infoString += ("Added MoverSingleHydrogenRotator to "+resNameAndID+" "+aName+
            " with "+str(len(potentialAcceptors))+" potential nearby acceptors\n")
      except Exception as e:
        infoString += "Could not add MoverSingleHydrogenRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n"

    # See if we should construct a MoverNH3Rotator here.
    # Find any Nitrogen that has four total bonded neighbors, three of which are Hydrogens.
    # @todo Is this a valid way to search for them?
    if a.element == 'N' and len(bondedNeighborLists[a]) == 4:
      numH = 0
      for n in bondedNeighborLists[a]:
        if n.element == "H":
          numH += 1
      if numH == 3:
        try:
          movers.append(Movers.MoverNH3Rotator(a, bondedNeighborLists))
          infoString += "Added MoverNH3Rotator to "+resNameAndID+"\n"
        except Exception as e:
          infoString += "Could not add MoverNH3Rotator to "+resNameAndID+": "+str(e)+"\n"

    # See if we should construct a MoverAromaticMethylRotator or MoverTetrahedralMethylRotator here.
    # Find any Carbon that has four total bonded neighbors, three of which are Hydrogens.
    # @todo Is this a valid way to search for them?
    if a.element == 'C' and len(bondedNeighborLists[a]) == 4:
      numH = 0
      neighbor = None
      for n in bondedNeighborLists[a]:
        if n.element == "H":
          numH += 1
        else:
          neighbor = n
      if numH == 3:
        # See if the Carbon's other neighbor is attached to two other atoms (3 total).  If so,
        # then insert a MoverAromaticMethylRotator and if not, generate a MoverTetrahedralMethylRotator
        # so that they Hydrogens will be staggered but do not add it to those to be optimized.
        if len(bondedNeighborLists[neighbor]) == 3:
          try:
            movers.append(Movers.MoverAromaticMethylRotator(a, bondedNeighborLists))
            infoString += "Added MoverAromaticMethylRotator to "+resNameAndID+" "+aName+"\n"
          except Exception as e:
            infoString += "Could not add MoverAromaticMethylRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n"
        else:
          try:
            ignored = Movers.MoverTetrahedralMethylRotator(a, bondedNeighborLists)
            infoString += "Used MoverTetrahedralMethylRotator to stagger "+resNameAndID+" "+aName+"\n"
          except Exception as e:
            infoString += "Could not add MoverTetrahedralMethylRotator to "+resNameAndID+" "+aName+": "+str(e)+"\n"

    # See if we should insert a MoverNH2Flip here.
    # @todo Is there a more general way than looking for specific names?
    if (aName == 'ND2' and resName == 'ASN') or (aName == 'NE2' and resName == 'GLN'):
      try:
        movers.append(Movers.MoverNH2Flip(a, "CA", bondedNeighborLists))
        infoString += "Added MoverNH2Flip to "+resNameAndID+"\n"
      except Exception as e:
        infoString += "Could not add MoverNH2Flip to "+resNameAndID+": "+str(e)+"\n"

    # See if we should insert a MoverHistidineFlip here.
    if aName == 'NE2' and resName == 'HIS':
      try:
        movers.append(Movers.MoverHistidineFlip(a, bondedNeighborLists))
        infoString += "Added MoverHistidineFlip to "+resNameAndID+"\n"
      except Exception as e:
        infoString += "Could not add MoverHistidineFlip to "+resNameAndID+": "+str(e)+"\n"

  return PlaceMoverstReturn(movers, infoString)

##################################################################################
# Helper functions

def AlternatesInModel(model):
  """Returns a list of altloc names of all conformers in all chains.
  :returns Set of strings.  The set is will include only the empty string if no
  chains in the model have alternates.  It has an additional entry for every altloc
  found in every chain.
  """
  ret = set()
  for c in model.chains():
    for alt in c.conformers():
      ret.add(alt.altloc)
  return ret

def GetAtomsForConformer(model, conf):
  """Returns a list of atoms in the named conformer.  It also includes atoms from
  the first conformation in chains that do not have this conformation, to provide a
  complete model.  For a model whose chains only have one conformer, it will return
  all of the atoms.
  :param model: pdb.hierarchy.model to search for atoms.
  :param conf: String name of the conformation to find.  Can be "".
  :returns List of atoms consistent with the specified conformer.
  """
  ret = []
  for ch in model.chains():
    confs = ch.conformers()
    which = 0
    for i in range(1,len(confs)):
      if confs[i].altloc == conf:
        which = i
        break
    ret.extend(ch.atoms())
  return ret

##################################################################################
# Test function to verify that all functions behave properly.

def Test(inFileName = None):
  """Test function for all functions provided above.
  :param inFileName: Name of a PDB or CIF file to load (default makes a small molecule)
  :returns Empty string on success, string describing the problem on failure.
  """

  #========================================================================
  # Test the AlternatesInModel() function.
  # @todo

  #========================================================================
  # Test the GetAtomsForConformer() function.
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

  # Make sure we have a valid unit cell.  Do this before we interpret the model because
  # it clobbers information that we need from interpretation.
  model = shift_and_box_model(model = model)

  print('Interpreting model')
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  model.set_pdb_interpretation_params(params = p)
  model.process_input_model(make_restraints=True) # make restraints

  # Add Hydrogens to the model
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()
  model.process_input_model(make_restraints=True) # make restraints

  # Get the first model in the hierarchy.
  firstModel = model.get_hierarchy().models()[0]

  # Get the list of alternate conformation names present in all chains for this model.
  alts = AlternatesInModel(firstModel)

  # Get the atoms from the first conformer in the first model (the empty string is the name
  # of the first conformation in the model; if there is no empty conformation, then it will
  # pick the first available conformation for each atom group.
  atoms = GetAtomsForConformer(firstModel, "")

  ################################################################################
  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  ################################################################################
  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = Helpers.getBondedNeighborLists(atoms, bondProxies)

  ################################################################################
  # Get the spatial-query information needed to quickly determine which atoms are nearby
  sq = probeExt.SpatialQuery(atoms)

  ################################################################################
  # Get the probeExt.ExtraAtomInfo needed to determine which atoms are potential acceptors.
  ret = Helpers.getExtraAtomInfo(model)
  extra = ret.extraAtomInfo
  infoString += ret.warnings

  ################################################################################
  # Get the list of Movers
  ret = PlaceMovers(atoms, model.rotatable_hd_selection(iselection=True), bondedNeighborLists, sq, extra)
  infoString += ret.infoString
  movers = ret.moverList
  print('XXX info:\n'+infoString)
  print('XXX Found',len(movers),'Movers')

  # @todo
  return ""

##################################################################################
# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  parser = argparse.ArgumentParser(description='Test mmtbx.reduce.Optimizers.')
  parser.add_argument('inputFile', nargs='?', default="")
  args = parser.parse_args()

  ret = Test(args.inputFile)
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
