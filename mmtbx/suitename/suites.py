from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import  with_statement, print_function
import sys, os

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

from suitenamedefs import globals
from iotbx.data_manager import DataManager    #   Load in the DataManager
from libtbx import phil
from libtbx.utils import Sorry
# from mmtbx.validation import utils
# from cctbx import geometry_restraints
# from collections import defaultdict
from diangle import getResidueDihedrals


# IMPORT TO EXPORT:
from mmtbx.suitename.suitename import compute, write, \
    finalStats, clearStats


# The following are the options available, in Phil format,
# for human and computer comprehension.
philOptions = """
  suitename {
    # input
      infile=""
        .type=str
        .help="the file to process"
      anglefields = 9
        .type=int
        .help="number of angle fields provided, for textual input only"
      pointidfields = 7
        .type=int
        .help="number of point id fields before the angle fields"
      ptid=0
        .type=int
        .help="number of point id fields before the angle fields"
      residuein=false
        .type=bool
        .help="expect dangle format giving residues"
      suitein=false
        .type=bool
        .help="expect kinemage format giving suites directly"
    # output
      string=False
        .type=bool
        .help="output in string format, 3 characters per suite"
      kinemage=False
        .type=bool
        .help="output in kinemage format, useful for visualization"
      report=true
        .type=bool
        .help="output as a report, giving statistical details"
      chart=False
        .type=bool
        .help="modifier to standard report, output without statistical summary"
      nosequence = False
        .type=bool
        .help="modifier to string format, do not include base letters"
      causes=False
        .type=bool
        .help="output extra details concerning the causes of each assignment made"
      test=False
        .type=bool
        .help="display a lot of additional information about program internals"
    # compute
      satellites=False
        .type=bool
        .help="use the special satelliteWidths values for satellites"
      nowannabe=False
        .type=bool
        .help="do not consider 'wannabe' clusters"
      noinc=False
        .type=bool
        .help="do not display incomplete suites"
      etatheta=False
        .type=bool
      altid="A"
        .type=str
        .help="which alternate conformer to use (A, B, etc)"
      altidfield = 6
        .type=int
        .help="which field (1-based) gives the alternate conformer code"
      version=false
        .type=bool
        .help="give the version number of suite name"
    # deprecated and automatically true:
      oneline=false
        .type=bool
    }
"""

def main(options, outFile=None, errorFile=None):
  """The main track for handling PDB and CIF input formats, which will involve
     parsing the model hierarchy to get the dihedral angles for ourselves"""
  setOptions(options)
  import suiteninput  # must be AFTER setOptions

  if not outFile: outFile = sys.stdout
  if not errorFile: errorFile = sys.stderr
  inFile = options.infile
  model = loadModel(inFile)

  residues = getResidueDihedrals(model, options.altid,
                                 name=os.path.splitext(inFile)[0],
                                 errorFile=errorFile)
  ### to print mp_geo-like output:
  # for r in residues:
  #   print(residueString(r))
  # useful for seeing what suites were generated

  if residues is not None and len(residues) > 0:
    suiteList = suiteninput.buildSuites(residues)
    suiteList = suiteList[:-1]
    suiteList = compute(suiteList)
    finalStats()
    write(outFile, suiteList)
    clearStats()


def parseOptions(optionString, errorFile=None):
  """ Use optionString to modify the defaults given in philOptions above.
  Returns a Python object that has an attribute for every option listed
  in philOptions.  Example: "chart=true noinc=true causes=true"
  The values in optionString are case insensitive.
  """
  opt2 = """  # use this for more complex option types e.g. multiples
  suitename {
    report=true
    chart=true
}  """
  # user_phil = phil.parse(opt2)

  master_phil = phil.parse(philOptions)
  interp = master_phil.command_line_argument_interpreter()
  optionList = optionString.split()
  try:
    user_phil = interp.process(args=optionList)
  except Sorry as e:
    if errorFile is None:  errorFile = sys.stderr
    print(e, file=errorFile)

  working_phil = master_phil.fetch(sources=user_phil)
  full_options = working_phil.extract()
  return full_options.suitename


def setOptions(optionsIn):
  """optionsIn may be the result of parseOptions above
     or the result of an argparse parse_args operation"""
  from mmtbx.suitename.suitename import loadOptions
  globals.options = optionsIn
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

