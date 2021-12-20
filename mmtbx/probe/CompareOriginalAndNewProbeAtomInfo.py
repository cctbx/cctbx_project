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

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import

import subprocess
import os,sys
import argparse
from mmtbx.probe import Helpers

##################################################################################
# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  parser = argparse.ArgumentParser(description='Test mmtbx.reduce.Optimizers.')
  parser.add_argument("--distanceThreshold", type=float, help="Same atom must be >= this distance between files to report", default=0.02)
  parser.add_argument('inputFile', type=str, help="PDB formatted model file to compare")
  args = parser.parse_args()

  #==============================================================
  # Check out and build the original probe program in a subdirectory.
  print('Cloning original Probe repository')
  try:
    process = subprocess.Popen(["git","clone","https://github.com/rlabduke/probe"],
     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
  except Exception as e:
    print('Error during clone; git must be installed:',str(e))
    sys.exit(1)

  print('Building original Probe')
  cwd = os.getcwd();
  os.chdir("./probe")
  try:
    process = subprocess.Popen(["make"],
      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
  except Exception as e:
    print("Exception running make (note; this will not work on Windows):",str(e))
    os.chdir(cwd)
    sys.exit(2)
  os.chdir(cwd)
  if len(stderr) != 0:
    print('Error running make (note; this will not work on Windows):',stderr)
    sys.exit(3)

  #==============================================================
  # Run the original Probe on the input file.
  print('Running old and new Probe on',args.inputFile)
  with open(os.path.basename(args.inputFile)+".orig.out","w") as fo:
    with open(os.path.basename(args.inputFile)+".orig.err","w") as fe:
      process = subprocess.Popen(["./probe/probe","-quiet","-kin","-mc","-self","all",
        "-count","-sepworse","-dumpatoms",os.path.basename(args.inputFile)+".orig.dump",args.inputFile],
        stdout=fo, stderr=fe)
      process.wait()

  with open(os.path.basename(args.inputFile)+".new.out","w") as fo:
    with open(os.path.basename(args.inputFile)+".new.err","w") as fe:
      process = subprocess.Popen(["mmtbx.probe2",'source_selection="all"','record_added_hydrogens=False',
        'approach=self','count_dots=True','output.separate_worse_clashes=True',
        'output.file_name='+os.path.basename(args.inputFile)+".new.out",
        'output.dump_file_name='+os.path.basename(args.inputFile)+".new.dump",
        args.inputFile],
        stdout=fo, stderr=fe)
      process.wait()

  #==============================================================
  # Compare the two dump files.  Store the results in a comparison file and say
  # to look there.
  compare = Helpers.compareAtomInfoFiles(os.path.basename(args.inputFile)+".orig.dump",
    os.path.basename(args.inputFile)+".new.dump",
    distanceThreshold=args.distanceThreshold)
  compareFileName = "./"+os.path.basename(args.inputFile)+".compare"
  with open(compareFileName, "w") as f:
    f.write(compare)

  print('Look in',compareFileName,'for differences')
