##################################################################################
# This is a test program to validate that the Python wrapping of Probe worked.
#

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
import argparse

import boost_adaptbx.boost.python as bp
bp.import_ext("mmtbx_probe_ext")
import mmtbx_probe_ext as probeext
import mmtbx_probe_ext as probe

from mmtbx.probe import AtomTypes

from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
from cctbx.maptbx.box import shift_and_box_model
import mmtbx
import mmtbx.probe
import mmtbx.probe.Helpers

# To enable addition of Hydrogens
# @todo See if we can remove the shift and box once reduce_hydrogen is complete
from cctbx.maptbx.box import shift_and_box_model
from mmtbx.hydrogens import reduce_hydrogen

def RunProbeVsCCTBXTest(inFileName, useNeutronDistances = False):

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

  # Fix up bogus unit cell when it occurs by checking crystal symmetry.
  cs = model.crystal_symmetry()
  if (cs is None) or (cs.unit_cell() is None):
    model = shift_and_box_model(model = model)

  # Add Hydrogens to the model
  print('Adding Hydrogens')
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model)
  reduce_add_h_obj.run()
  model = reduce_add_h_obj.get_model()

  # Run PDB interpretation on the model
  # @todo Not sure this does anything... setting Neutron distances does not change radii?
  print('Interpreting model')
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = useNeutronDistances
  model.process(make_restraints=True, pdb_interpretation_params=p) # make restraints

  # Construct the AtomTypes object we're going to use, telling it whether to use neutron distances
  # or not.
  at = AtomTypes.AtomTypes(useNeutronDistances)

  # Traverse the hierarchy and look up the extra data from Probe and CCTBX.
  # print information on differences, keeping track of how many atoms differed.
  mon_lib_srv = model.get_mon_lib_srv()
  ener_lib = mmtbx.monomer_library.server.ener_lib(use_neutron_distances = useNeutronDistances)
  ph = model.get_hierarchy()
  diffs = 0
  count = 0
  for m in ph.models():
    for chain in m.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          for a in ag.atoms():

            count += 1

            # Look up in CCTBX
            ccei = probe.ExtraAtomInfo()
            if mmtbx.probe.Helpers.isMetallic(a):
              ccei.vdwRadius = model.get_specific_ion_radius(a.i_seq)
            else:
              ccei.vdwRadius = model.get_specific_vdw_radius(a.i_seq)
            hb_type = model.get_specific_h_bond_type(a.i_seq)
            if hb_type == "A" or hb_type == "B":  # B is for Both
              ccei.isAcceptor = True
            if hb_type == "D" or hb_type == "B":  # B is for Both
              ccei.isDonor = True

            # Look up in Probe
            ei, warn = at.FindProbeExtraAtomInfo(a)
            if len(warn) != 0:
              print(warn)

            # Compare the results.
            ''' @todo Donor and Acceptor info relies on information other than the tables.'''
            if ei.isDonor != ccei.isDonor:
              print(ag.resname,a.name,'Mismatched Donor status: Probe thinks',ei.isDonor)
              diffs += 1
            if ei.isAcceptor != ccei.isAcceptor:
              print(ag.resname,a.name,'Mismatched Acceptor status: Probe thinks',ei.isAcceptor)
              diffs += 1
            if abs(ei.vdwRadius - ccei.vdwRadius) > 1.0E-4:
              print(ag.resname,a.name,'Mismatched radii: Probe thinks',ei.vdwRadius,'CCTBX thinks',ccei.vdwRadius)
              diffs += 1

  if diffs == 0:
    return ""
  else:
    print("Found "+str(diffs)+" differences between Probe and CCTBX in "+str(count)+" atoms")
    return "Found "+str(diffs)+" differences between Probe and CCTBX "+str(count)+" atoms"

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  parser = argparse.ArgumentParser(description='Compare Probe vs. CCTBX atom information')
  parser.add_argument('-n','--neutronDistances', default=False, action='store_true',
                      help='Use neutron distances (default X-ray electron cloud distances)')
  parser.add_argument('inputFile', nargs='?', default="")
  args = parser.parse_args()

  ret = RunProbeVsCCTBXTest(args.inputFile, args.neutronDistances)
  if len(ret) == 0:
    print('Success!')

  assert (len(ret) == 0)
