import sys, os, inspect

#        Copyright 2021  Richardson Lab at Duke University
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

from suitenamedefs import Residue, Suite, globals
from iotbx.data_manager import DataManager    #   Load in the DataManager
from mmtbx.validation import utils
from cctbx import geometry_restraints
from collections import defaultdict
from myangle import getResidueDihedrals, residueString


def main(options, outFile=None):
  from mmtbx.suitename.suitename import compute, write, finalStats, clearStats
  setOptions(options)
  import suiteninput  # AFTER setOptions

  if not outFile:
    outFile = sys.stdout
  inFile = options.infile
  model = loadModel(inFile)


  residues = getResidueDihedrals(model, options.altid, 
                                 name=os.path.splitext(inFile)[0])
  ### to print mp_geo-like output:
  # for r in residues:
  #   print(residueString(r))
                                   
  if len(residues) == 0:
      sys.stderr.write("read no residues: perhaps wrong alternate code\n")
      sys.exit(1)
  suiteList = suiteninput.buildSuites(residues)
  suiteList = suiteList[:-1]
  
  suiteList = compute(suiteList)
  finalStats()
  write(outFile, suiteList)
  clearStats()


def setOptions(optionsIn):
  from mmtbx.suitename.suitename import loadOptions
  global options
  options = optionsIn
  globals.options = options
  loadOptions(optionsIn)


def loadModel(filename):
  dm = DataManager()             #   Initialize the DataManager and call it dm
  dm.set_overwrite(True)         #   tell the DataManager to overwrite files with the same name
  #print("Reading file")
  model = dm.get_model(filename)
  #print("Processing model")
  #model.process_input_model(make_restraints=True)
  # removed because Restraints Manager will not operate 
  # on unfamiliar residues  KPB 6/10/2021
  return model


def testResidues(model):
  #print("computing dihedrals")
  residues = getResidueDihedrals(model)
  for r in residues:
    print(r.pointIDs, " : ", r.angle)

