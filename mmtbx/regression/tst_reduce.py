##################################################################################
# This is a test program to validate that the Python wrapping of Reduce worked.
#

#                Copyright 2021-2023  Richardson Lab at Duke University
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
from mmtbx.reduce import Movers
from mmtbx.reduce import InteractionGraph
from mmtbx.reduce import Optimizers
from mmtbx.programs import reduce2
from iotbx.cli_parser import run_program
from six.moves import cStringIO as StringIO
from libtbx.utils import Sorry, Usage, null_out
import libtbx.load_env
import os.path as op
import os
import math

def RunReduceTests():

  #========================================================================
  # Call the test functions for all of our files.

  print('Testing Movers objects')
  ret = Movers.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0)

  print('Testing InteractionGraph objects')
  ret = InteractionGraph.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0)

  print('Testing Optimizers')
  ret = Optimizers.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0)

  return ret

# Each test case has a name, a raw PDB file, a chain ID, a residue ID, an atom name,
# a list of possible locations for the atom (for sets of 3 hydrogens, any 60-degree
# rotation is equivalent in score), and a maximum distance threshold.
testCases = [
  ["7c31_single_hydrogen_rotator",
   """\
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
""",
   "A",
   5,
   "HG",
   [ (-33.788, 46.891, 0.809)
   ],
   0.1,
  ]
]

def RunRegressionTests():
  for tc in testCases:
    name = tc[0]
    pdb_raw = tc[1]
    chain = tc[2]
    resID = tc[3]
    atomName = tc[4]
    positions = tc[5]
    maxDist = tc[6]

    print('Testing regression on', name)

    pdb_file = "./deleteme.pdb"
    with open(pdb_file, "w") as f:
      f.write(pdb_raw)
    if (op.exists("./deletemeFH.pdb")):
      os.remove("./deletemeFH.pdb")
    if (op.exists("./deleteme_description.txt")):
      os.remove("./deleteme_description.txt")
    out = StringIO()
    try:
      # Run the program
      results = run_program(program_class=reduce2.Program,
        logger=out,
        args=[pdb_file, "add_flip_movers=True", "output.description_file=./deleteme_description.txt"])
      # Check the position of the atom to see if it is near enough to one of the expected locations.
      found = False
      for c in results.model.chains():
        if c.id == chain:
          for rg in c.residue_groups():
            if rg.resseq_as_int() == resID:
              for atom in c.atoms():
                if atom.name.strip().upper() == atomName:
                  found = True
                  loc = atom.xyz
                  closeEnough = False
                  for pos in positions:
                    dist = math.sqrt( (loc[0]-pos[0])*(loc[0]-pos[0]) +
                      (loc[1]-pos[1])*(loc[1]-pos[1]) +
                      (loc[2]-pos[2])*(loc[2]-pos[2])
                    )
                    if dist <= maxDist:
                      closeEnough = True
                  if not closeEnough:
                    return "Atom "+chain+" "+str(resID)+" "+atomName+" in "+name+" too far from expected locations"
      if not found:
        return "Did not find atom "+atomName+" in chain "+chain+" residue "+str(resID)+" of "+name
    except Exception as e:
      return "Exception when running reduce2: "+str(e)
    if (op.exists(pdb_file)):
      os.remove(pdb_file)
    if (op.exists("./deletemeFH.pdb")):
      os.remove("./deletemeFH.pdb")
    if (op.exists("./deleteme_description.txt")):
      os.remove("./deleteme_description.txt")

  #========================================================================
  # Regression test a Reduce2 run against 1xso, comparing flips
  # to the output generated by a previous version of the program.  If there are
  # differences, report that this is the case and recommend verifying that the
  # differences are intentional and replacing the stored output.
  data_dir = libtbx.env.under_dist(
    module_name = "mmtbx",
    path = os.path.join("regression","pdbs"),
    test = os.path.isdir)
  pdb_file = os.path.join(data_dir,'1xso.pdb')
  if (op.exists("./deletemeFH.pdb")):
    os.remove("./deletemeFH.pdb")
  if (op.exists("./deleteme_description.txt")):
    os.remove("./deleteme_description.txt")
  out = StringIO()
  try:
    # Run the program
    results = run_program(program_class=reduce2.Program,
      logger=out,
      args=[pdb_file, "add_flip_movers=True",
            "output.file_name=./deletemeFH.pdb",
            "output.description_file=./deleteme_description.txt",
            "alt_id=a", "output.overwrite=True"])
    # Parse the description file to see whether the flips are as expected.
    expected = ['Unflipped', 'Unflipped', 'Unflipped',
                'Unflipped',
                'Unflipped', 'Unflipped', 'Flipped', 'Unflipped',
                'Unflipped', 'Flipped', 'Unflipped', 'Flipped',
                'Unflipped',
                'Unflipped', 'Unflipped', 'Unflipped', 'Flipped', 'Unflipped', 'Unflipped', 'Unflipped',
                'Flipped',
                'Unflipped', 'Unflipped', 'Unflipped',
                'Flipped', 'Unflipped', 'Unflipped', 'Unflipped']
    with open('./deleteme_description.txt', 'r') as f:
      lines = f.readlines()
      which = 0
      for line in lines:
        if 'lipped' in line:
          if not expected[which] in line:
           return "Flips in 1xso don't match expected flips for flip {}".format(which)
          which += 1

  except Exception as e:
    return "Exception when running reduce2: "+str(e)
  if (op.exists("./deletemeFH.pdb")):
    os.remove("./deletemeFH.pdb")
  if (op.exists("./deleteme_description.txt")):
    os.remove("./deleteme_description.txt")

  return ""

if __name__ == '__main__':

  ret = RunReduceTests()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)

  ret = RunRegressionTests()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
