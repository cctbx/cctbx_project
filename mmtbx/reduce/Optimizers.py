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

# To enable addition of Hydrogens
# @todo 
from cctbx.maptbx.box import shift_and_box_model
from mmtbx.hydrogens import reduce_hydrogen

##################################################################################
# This is a set of functions that implement placement and optimization of
# Reduce's "Movers".

class PlacementReturn:
  # Return type from PlaceMovers() call.  List of movers and then an information string
  # that may contain information the user would like to know (where Movers were placed,
  # failed Mover placements due to missing Hydrogens, etc.).
  def __init__(self, moverList, infoString):
    self.moverList = moverList
    self.infoString = infoString

def PlaceMovers(atoms, bondedNeighborLists):
  """Produce a list of Movers for atoms in a pdb.hierarchy.conformer that has added Hydrogens.
  :param atoms: flex array of atoms to search.  This must have all Hydrogens needed by the
  Movers present in the structure already.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling getBondedNeighborLists().
  :returns PlacementReturn giving the list of Movers found in the conformation and
  an error string that is empty if no errors are found during the process and which
  has a printable message in case one or more errors are found.
  """

  # List of Movers to return
  movers = []

  # Information string to return, filled in when we cannot place a Mover.
  infoString = ""

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
    # @todo

    # See if we should construct a MoverNH3Rotator here.
    # @todo

    # See if we should construct a MoverAromaticMethylRotator here.
    # @todo

    # See if we should construct a MoverTetrahedralMethylRotator here so that we
    # ensure that the Hydrogens are staggered but don't add it to those that are
    # optimized.
    # @todo

    # See if we should insert a MoverNH2Flip here.
    # @todo Is there a more general way than looking for specific names?
    if (aName == 'ND2' and resName == 'ASN') or (aName == 'NE2' and resName == 'GLN'):
      # Find the alpha carbon associated with this residue
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

  return PlacementReturn(movers, infoString)

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

  print('Interpreting model')
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  model.set_pdb_interpretation_params(params = p)
  model.process_input_model(make_restraints=True) # make restraints

  # Add Hydrogens to the model
  model = shift_and_box_model(model = model)
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Get the first model in the hierarchy.
  firstModel = model.get_hierarchy().models()[0]

  # Get the list of alternate conformation names present in all chains for this model.
  alts = AlternatesInModel(firstModel)

  # Get the atoms from the first conformer in the first model.
  atoms = GetAtomsForConformer(firstModel, "")

  # Get the Cartesian positions of all of the atoms we're considering for this alternate
  # conformation.
  carts = flex.vec3_double()
  for a in atoms:
    carts.append(a.xyz)

  # Get the bond proxies for the atoms in the model and conformation we're using and
  # use them to determine the bonded neighbor lists.
  bondProxies = model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
  bondedNeighborLists = Movers.getBondedNeighborLists(atoms, bondProxies)

  # Get the list of Movers
  ret = PlaceMovers(atoms, bondedNeighborLists)
  movers = ret.moverList
  print('XXX info:\n'+ret.infoString)
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
