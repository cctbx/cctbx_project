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

import sys
import iotbx.map_model_manager
import iotbx.data_manager
import cctbx.maptbx.box
import mmtbx
import mmtbx_probe_ext as probe
import mmtbx.probe.AtomTypes

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
    the atoms that are bonded to it.
  """
  atomDict = {}
  for a in atoms:
    atomDict[a.i_seq] = a
  bondedNeighbors = {}
  for a in atoms:
    bondedNeighbors[a] = []
  for bp in bondProxies:
    bondedNeighbors[atomDict[bp.i_seqs[0]]].append(atomDict[bp.i_seqs[1]])
    bondedNeighbors[atomDict[bp.i_seqs[1]]].append(atomDict[bp.i_seqs[0]])
  return bondedNeighbors

class getExtraAtomInfoReturn:
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

def getExtraAtomInfo(model, useNeutronDistances = False):
  """
    Helper function to provide a mapper for ExtraAtomInfo needed by Probe when scoring
    models.  It first tries to find the information in CCTBX.  If it cannot, it looks
    the information up using the original C-code Probe tables and algorithms.
    :param model: Map Model Manager's Model containing all of the atoms to be described.
    PDB interpretation must have been done on the model, perhaps by calling
    model.process_input_model(make_restraints=True), with useNeutronDistances matching
    the parameter to this function.
    :param useNeutronDistances: Default is to use x-ray distances, but setting this to
    True uses neutron distances instead.
    Can be obtained by calling iotbx.map_model_manager.map_model_manager().model().
    :returns a ExtraAtomInfoMap with an entry for every atom in the model suitable for
    passing to the scoring functions.
  """

  warnings = ""

  # Construct the AtomTypes object we're going to use, telling it whether to use neutron distances.
  at = mmtbx.probe.AtomTypes.AtomTypes(useNeutronDistances)

  # Traverse the hierarchy and look up the extra data to be filled in.
  extras = probe.ExtraAtomInfoMap([],[])
  mon_lib_srv = model.get_mon_lib_srv()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  ph = model.get_hierarchy()
  for m in ph.models():
    for chain in m.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
                residue_name=ag.resname, atom_names=ag.atoms().extract_name())
          atom_dict = md.atom_dict()

          for a in ag.atoms():
            extra = probe.ExtraAtomInfo()
            # See if we can find out about its Hydrogen-bonding status from the
            # model.  If so, we fill it and the vdwRadius information from
            # CCTBX.
            try:
              hb_type = model.get_specific_h_bond_type(a)
              if isinstance(hb_type, str):
                if hb_type == "A" or hb_type == "B":
                  extra.isAcceptor = True
                if hb_type == "D" or hb_type == "B":
                  extra.isDonor = True
                # @todo How to tell this code whether or not to use neutron distances?
                extra.vdwRadius = model.get_specific_vdw_radii(a)
                continue

              # Did not find the information from CCTBX, so look it up using
              # the original Probe approach by dropping through to below
              else:
                warnings += "Could not find "+a.name+" in CCTBX, using Probe tables\n"
            except:
              # Warn and drop through to below.
              warnings += ("Could not look up "+a.name+" in CCTBX "+
                "(perhaps interpretation was not run on the model?), using Probe tables\n")

            # Did not find what we were looking for in CCTBX, so drop through to Probe
            extra, warn = at.FindProbeExtraAtomInfo(a)
            if len(warn) > 0:
              warnings += "  Probe says: "+warn+"\n"
            extras.setMappingFor(a, extra)

  return getExtraAtomInfoReturn(extras, warnings)

def Test(inFileName = None):

  #========================================================================
  # Generate an example data model with a small molecule in it or else read
  # from the specified file.
  if inFileName is not None and len(inFileName) > 0:
    # Read a model from a file using the DataManager
    print('Reading model from',inFileName)
    dm = iotbx.data_manager.DataManager()
    dm.process_model_file(inFileName)
    model = dm.get_model(inFileName)
  else:
    # Generate a small-molecule model using the map model manager
    print('Generating model')
    # get an initialized instance of the map_model_manager
    mmm=iotbx.map_model_manager.map_model_manager()
    mmm.generate_map()     #   get a model from a generated small library model and calculate a map for it
    model = mmm.model()    #   get the model

  # Fix up bogus unit cell when it occurs by checking crystal symmetry.
  cs =model.crystal_symmetry()
  if (cs is None) or (cs.unit_cell() is None):
    model = cctbx.maptbx.box.shift_and_box_model(model = model)

  # Run PDB interpretation on the model to fill in the required CCTBX information.
  print('Interpreting model')
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = False
  model.set_pdb_interpretation_params(params = p)
  model.process_input_model(make_restraints=True) # make restraints

  print('Getting extraAtomInfo')
  ret = getExtraAtomInfo(model)
  if len(ret.warnings) > 0:
    print('Warnings returned by getExtraAtomInfo():\n'+ret.warnings)

  # Run unit test on getBondedNeighborLists().
  # @todo

  # Run unit test on getExtraAtomInfo().
  # @todo

  return ""

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB/CIF file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    fileName = sys.argv[i]

  ret = Test(fileName)
  if len(ret) == 0:
    print('Success!')

  assert (len(ret) == 0)
