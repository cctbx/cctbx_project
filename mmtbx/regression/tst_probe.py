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

from __future__ import absolute_import, division, print_function
from libtbx.utils import format_cpu_times
import sys, os, subprocess, tempfile, platform
from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
from cctbx.maptbx.box import shift_and_box_model
import mmtbx
import libtbx.load_env

import boost_adaptbx.boost.python as bp
bp.import_ext("mmtbx_probe_ext")
import mmtbx_probe_ext as probeext

from mmtbx.probe import Helpers, AtomTypes

def RunProbeTests(inFileName):

  #========================================================================
  # Call the test functions for the libraries we test.

  ret = probeext.DotSpheres_test()
  assert len(ret) == 0, "DotSpheres_test() failed: " + ret

  ret = probeext.SpatialQuery_test()
  assert len(ret) == 0, "SpatialQuery_test() failed: " + ret

  ret = probeext.Scoring_test()
  assert len(ret) == 0, "Scoring_test() failed: " + ret

  AtomTypes.Test()
  Helpers.Test()

  #========================================================================
  # Now ensure that we can use the C++-wrapped classes as intended to make sure
  # that the wrapping code or parameters have not changed.

  #========================================================================
  # Make sure we can get at the DotSphere objects and their methods
  cache = probeext.DotSphereCache(10)
  sphere1 = cache.get_sphere(1)
  dots = sphere1.dots()

  #========================================================================
  # Make sure we can fill in an ExtraAtomInfoList and pass it to scoring
  # Generate an example data model with a small molecule in it unless we've
  # been given a file name to open.
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
  cs = model.crystal_symmetry()
  if (cs is None) or (cs.unit_cell() is None):
    model = shift_and_box_model(model = model)

  # Get the list of all atoms in the model
  atoms = model.get_atoms()

  # Get the bonding information we'll need to exclude our bonded neighbors.
  try:
    p = mmtbx.model.manager.get_default_pdb_interpretation_params()
    model.process(make_restraints=True, pdb_interpretation_params=p) # make restraints
    geometry = model.get_restraints_manager().geometry
    sites_cart = model.get_sites_cart() # cartesian coordinates
    bond_proxies_simple, asu = \
        geometry.get_all_bond_proxies(sites_cart = sites_cart)
  except Exception as e:
    raise Exception("Could not get bonding information for input file: " + str(e))
  bondedNeighbors = Helpers.getBondedNeighborLists(atoms, bond_proxies_simple)

  # Traverse the hierarchy and look up the extra data to be filled in.
  class philLike:
    def __init__(self, useImplicitHydrogenDistances = False):
      self.implicit_hydrogens = useImplicitHydrogenDistances
      self.set_polar_hydrogen_radius = True
  ret = Helpers.getExtraAtomInfo(model,bondedNeighbors,
          useNeutronDistances=False,probePhil=philLike(False))
  extra = ret.extraAtomInfo

  # Construct a SpatialQuery and fill in the atoms.  Ensure that we can make a
  # query within 1000 Angstroms of the origin.
  sq = probeext.SpatialQuery(atoms)
  nb = sq.neighbors((0,0,0), 0, 1000)

  # Construct a DotScorer object.
  # Find the radius of each atom in the structure and construct dot spheres for
  # them. Find the atoms that are bonded to them and add them to an excluded list.
  # Then compute the score for each of them and report the summed score over the
  # whole molecule the way that Reduce will.
  ds = probeext.DotScorer(extra)
  total = 0
  badBumpTotal = 0
  for a in atoms:
    rad = extra.getMappingFor(a).vdwRadius
    assert rad > 0, "Invalid radius for atom look-up: "+a.name+" rad = "+str(rad)
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
    res = ds.score_dots(a, 1.0, sq, rad*3, 0.25, exclude, sphere.dots(), sphere.density(), False, False)
    total += res.totalScore()
    if res.hasBadBump:
      badBumpTotal += 1

  # Test calling the single-dot checking code as will be used by Probe to make sure
  # all of the Python linkage is working
  dotOffset = [1, 0, 0]
  check = ds.check_dot(atoms[0], dotOffset, 1, atoms, [atoms[0]])
  overlapType = check.overlapType

  # Test calling the interaction_type method to be sure Python linkage is working
  interactionType = ds.interaction_type(check.overlapType, check.gap)

  #========================================================================
  # Regression test a Probe2 run against a snippet of a file, comparing the output
  # to the output generated by a previous version of the program.  If there are
  # differences, report that this is the case and recommend verifying that the
  # differences are intentional and replacing the stored output.
  data_dir = libtbx.env.under_dist(
    module_name = "mmtbx",
    path = os.path.join("regression","pdbs"),
    test = os.path.isdir)
  model_file = os.path.join(data_dir,'Fe_1brf_snip_reduced.pdb')
  kin_dir = libtbx.env.under_dist(
    module_name = "mmtbx",
    path = os.path.join("regression","kins"),
    test = os.path.isdir)
  kin_file = os.path.join(kin_dir,'Fe_1brf_snip_reduced.kin')
  temp_file = os.path.join(tempfile._get_default_tempdir(),
    next(tempfile._get_candidate_names())+".kin" )
  try:
    my_env = os.environ
    exe_name = 'mmtbx.probe2'
    if platform.system() == 'Windows':
      exe_name += '.bat'
    if subprocess.check_call([exe_name
                              ,'source_selection="all"'
                              ,'approach=self'
                              ,'output.separate_worse_clashes=True'
                              ,'output.file_name='+temp_file
                              ,'include_mainchain_mainchain=True'
                              ,'output.add_kinemage_keyword=True'
                              ,model_file
                            ], env = my_env
                             , stdout = subprocess.DEVNULL
                             , stderr = subprocess.DEVNULL):
      raise Exception("Call to subprocess to regression test had nonzero return")
  except Exception as e:
    raise Exception("Could not call subprocess to do regression test: "+str(e))
  with open(temp_file) as ft:
    ft_text = ft.readlines()
  with open(kin_file) as fk:
    fk_text = fk.readlines()
  instructions = ("Use KiNG or another program to see what changed and then determine if the "+
      "differences are expected.  If so, replace "+kin_file+" with the new file.")
  if len(ft_text) != len(fk_text):
    raise Exception("Different number of lines in "+temp_file+" and "+kin_file+instructions)
  for i in range(3,len(ft_text)):
    if ft_text[i] != fk_text[i]:
      print('Line',i,'from each file:')
      print(ft_text[i])
      print(fk_text[i])
      raise Exception("Line "+str(i)+" in "+temp_file+" and "+kin_file+" differ.  "+instructions)

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    fileName = sys.argv[i]

  RunProbeTests(fileName)
  print(format_cpu_times())
  print('OK')
