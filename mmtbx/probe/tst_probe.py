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
import mmtbx_probe_ext as probe
from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
from cctbx.maptbx.box import shift_and_box_model
import mmtbx

import Helpers

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

def RunProbeTests(inFileName):

  #========================================================================
  # Make sure we can get at the DotSphere objects and their methods
  print('Making DotSphere from cache and getting its dots')
  cache = probe.DotSphereCache(10)
  sphere1 = cache.get_sphere(1)
  dots = sphere1.dots()
  print(' Found',len(dots),'dots')
  print('First dot is at',dots[0][0],dots[0][1],dots[0][2])

  #========================================================================
  # Make sure we can fill in an ExtraAtomInfoList and pass it to scoring
  # Generate an example data model with a small molecule in it
  print('Generating model')
  if inFileName is not None and len(inFileName) > 0:
    # Read a model from a file using the DataManager
    dm = DataManager()
    dm.process_model_file(inFileName)
    model = dm.get_model(inFileName)
  else:
    # Generate a small-molecule model using the map model manager
    mmm=map_model_manager()         #   get an initialized instance of the map_model_manager
    mmm.generate_map()              #   get a model from a generated small library model and calculate a map for it
    model = mmm.model()             #   get the model

  # Fix up bogus unit cell when it occurs by checking crystal symmetry.
  cs =model.crystal_symmetry()
  if (cs is None) or (cs.unit_cell() is None):
    model = shift_and_box_model(model = model)

  # Get the list of all atoms in the model
  atoms = model.get_atoms()

  # Get the bonding information we'll need to exclude our bonded neighbors.
  try:
    p = mmtbx.model.manager.get_default_pdb_interpretation_params()
    model.set_pdb_interpretation_params(params = p)
    model.process_input_model(make_restraints=True) # make restraints
    geometry = model.get_restraints_manager().geometry
    sites_cart = model.get_sites_cart() # cartesian coordinates
    bond_proxies_simple, asu = \
        geometry.get_all_bond_proxies(sites_cart = sites_cart)
  except Exception as e:
    return "Could not get bonding information for input file: " + str(e)
  bondedNeighbors = Helpers.getBondedNeighborLists(atoms, bond_proxies_simple)

  # Traverse the hierarchy and look up the extra data to be filled in.
  print('Filling in extra atom information needed for probe score')
  ret = Helpers.getExtraAtomInfo(model)
  extra = ret.extraAtomInfo
  if len(ret.warnings) > 0:
    print('Warnings returned by getExtraAtomInfo():\n'+ret.warnings)

  # Construct a SpatialQuery and fill in the atoms.  Ensure that we can make a
  # query within 1000 Angstroms of the origin.
  sq = probe.SpatialQuery(atoms)
  nb = sq.neighbors((0,0,0), 0, 1000)
  print('Found this many atoms within 1000A of the origin:', len(nb))

  # Construct a DotScorer object.
  # Find the radius of each atom in the structure and construct dot spheres for
  # them. Find the atoms that are bonded to them and add them to an excluded list.
  # Then compute the score for each of them and report the summed score over the
  # whole molecule the way that Reduce will.
  ds = probe.DotScorer(extra)
  total = 0
  badBumpTotal = 0
  for a in atoms:
    rad = extra[a.i_seq].vdwRadius
    sphere = cache.get_sphere(rad)

    # Excluded atoms that are bonded to me or to one of my neightbors.
    # It has the side effect of excluding myself if I have any neighbors.
    # Construct as a set to avoid duplicates.
    exclude = set()
    for n in bondedNeighbors[a]:
      exclude.add(n)
      for n2 in bondedNeighbors[n]:
        exclude.add(n2)
    exclude = list(exclude)

    dots = sphere.dots()
    res = ds.score_dots(a, 1.0, sq, rad*3, 0.25, exclude, sphere.dots(), sphere.density(), False)
    total += res.totalScore()
    if res.hasBadBump:
      badBumpTotal += 1
  print('Summed probe score for molecule =',total,'with',badBumpTotal,'bad bumps')

  # Test calling the single-dot checking code as will be used by Probe to make sure
  # all of the Python linkage is working
  dotOffset = [1, 0, 0]
  check = ds.check_dot(atoms[0], dotOffset, 1, atoms, [atoms[0]])
  print ('Check dot overlap type:',check.overlapType)

  #========================================================================
  # Call the test functions for all of our files.

  print('Testing DotSphere objects')
  ret = probe.DotSpheres_test()
  if len(ret) > 0:
    return "DotSpheres_test() failed: " + ret

  print('Testing SpatialQuery objects')
  ret = probe.SpatialQuery_test()
  if len(ret) > 0:
    return "SpatialQuery_test() failed: " + ret

  print('Testing Scoring objects')
  ret = probe.Scoring_test()
  if len(ret) > 0:
    return "Scoring_test() failed: " + ret

  return ret

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    fileName = sys.argv[i]

  # Do not print anything on success because it may be counted as a warning.
  ret = RunProbeTests(fileName)
  if len(ret) == 0:
    print('OK')

  assert (len(ret) == 0), "Failure: " + ret
