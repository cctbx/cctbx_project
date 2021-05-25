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

import sys
import argparse

import mmtbx_probe_ext as probe
import AtomTypes

from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
from cctbx.maptbx.box import shift_and_box_model
import mmtbx

def BestMatch(name, d):
  '''Find the best match for the name in the dictionary keys.  It must match the
  first character of the name and then pick the one with the maximum number of
  matched characters.
  :param name: Name of the atom to match.
  :param d: atom_dict dictionary mapping atom names to type and energy.
  :returns Best match on success, raises KeyError on no match.
  '''

  # Get a list of names from the dictionary that match the first character.
  keys = d.keys()
  matches = []
  for k in keys:
    if k[0] == name[0]:
      matches.append(k)
  if len(matches) == 0:
    raise KeyError('Could not find atom name with matching first character from '+name)

  # Keep track of the match that has a maximum number of matching characters.
  maxCount = 0
  bestMatch = ""
  for m in matches:
    count = 0
    for c in name:
      if c in m:
        count += 1
    if count > maxCount:
      maxCount = count
      bestMatch = m

  return bestMatch

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
  cs =model.crystal_symmetry()
  if (cs is None) or (cs.unit_cell() is None):
    model = shift_and_box_model(model = model)

  # Run PDB interpretation on the model
  # @todo Not sure this does anything... setting Neutron distances does not change radii?
  print('Interpreting model')
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = useNeutronDistances
  model.set_pdb_interpretation_params(params = p)
  model.process_input_model(make_restraints=True) # make restraints

  # Construct the AtomTypes object we're going to use, telling it whether to use neutron distances
  # or not.
  class FakePhil:
    pass
  fakePhil = FakePhil()
  fakePhil.useNeutronDistances = useNeutronDistances
  at = AtomTypes.AtomTypes(fakePhil)

  # Traverse the hierarchy and look up the extra data from Probe and CCTBX.
  # print information on differences, keeping track of how many atoms differed.
  mon_lib_srv = model.get_mon_lib_srv()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  ph = model.get_hierarchy()
  diffs = 0
  count = 0
  for m in ph.models():
    for chain in m.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
                residue_name=ag.resname, atom_names=ag.atoms().extract_name())
          atom_dict = md.atom_dict()

          for a in ag.atoms():

            count += 1

            # Look up in CCTBX
            name = a.name.strip()
            try:
              te = atom_dict[name].type_energy
            except KeyError:
              match = BestMatch(name, atom_dict)
              print('Warning: Could not find entry in CCTBX for',name,'replacing with',match, flush=True)
              te = atom_dict[match].type_energy
            ccei = probe.ExtraAtomInfo()
            ccei.vdwRadius = ener_lib.lib_atom[te].vdw_radius
            hb_type = ener_lib.lib_atom[te].hb_type
            if hb_type == "A":
              ccei.isAcceptor = True
            if hb_type == "D":
              ccei.isDonor = True

            # Look up in Probe
            ei, warn = at.FindProbeExtraAtomInfo(a)
            if len(warn) != 0:
              print(warn)

            # Compare the results.
            ''' @todo Donor and Acceptor info relies on information other than the tables.
            if ei.isDonor != ccei.isDonor:
              print(ag.resname,a.name,'Mismatched Donor status: Probe thinks',ei.isDonor)
              diffs += 1
            if ei.isAcceptor != ccei.isAcceptor:
              print(ag.resname,a.name,'Mismatched Acceptor status: Probe thinks',ei.isAcceptor)
              diffs += 1
            '''
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
