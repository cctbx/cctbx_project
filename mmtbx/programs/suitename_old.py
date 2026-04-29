"""Identify suitenames for nucleic acid residues (old version)"""
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

from __future__ import division
from __future__ import nested_scopes, generators, absolute_import
from __future__ import with_statement, print_function
import sys, os

from mmtbx.suitename import dualparse, suites
from mmtbx.suitename.suitename import main, version

from iotbx.cli_parser import CCTBXParser
from libtbx.program_template import ProgramTemplate
from libtbx.utils import multi_out, show_total_time

#import libtbx.load_env
#from libtbx.utils import Usage

import os, sys

def run(args):
  "The main program, if run from CCTBX / PHENIX."
  logger = multi_out()
  logger.register('stderr', sys.stderr)
  logger2 = multi_out()
  logger2.register('stdout', sys.stdout)
  logger3 = multi_out()
  # TO DIAGNOSE OPTIONS TROUBLES:
  # logger3.register('verbiage', open("suitename.stderr.log", "w"))

  parser = dualparse.parseArgs(Program, logger3)
  working_phil = parser.working_phil
  options = working_phil.extract().suitename

  # now we call into the core of suitename itself
  if options.version:
      print(version, file=logger2)
      return
  if options.infile == "" or options.infile =="-" or options.residuein or options.suitein:
      # let the core figure out the input
      main(optionsIn=options, outFile=logger2, errorFile=logger)
  else:
      type, ext = analyzeFileType(options.infile)
      if type=="":
        logger.write("File extension "+str(ext)+" not recognized\n")
        return
      if type == "pdb":
          suites.main(options=options, outFile=logger2, errorFile=logger)
      else:
        # help the core figure out the input file type
        if type == "kinemage":
          options.suitein = True
        elif type == "dangle":
          options.residuein = True
        main(optionsIn=options, outFile=logger2, errorFile=logger)


extensionList={
    "pdb": "pdb",
    "cif": "pdb",
    "kin": "kinemage",
    "dangle": "dangle",
    "suitegeom": "dangle",
}

def analyzeFileType(filename):
    base, extension = os.path.splitext(filename)
    extension = extension[1:]
    type = extensionList.get(extension, "")
    return type, extension


class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  program_name = "suitename"
  description="""
   < insert help text here>
""" % locals()

  # The substructure below has been commented out because it would make it more
  # difficult to integrate with suitename's ability to run independently
  # of CCTBX
  master_phil_str = """
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
          .help="display a lat of additional information about program internals"
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
        version=False
          .type=bool
          .help="give the version number of suite name"
      # deprecated
        oneline=False
          .type=bool
      }
"""

# might add:
        # altidval="A"
        #   .type=str
        #   .help="which alternate conformer to use (A, B, etc)"
#
  datatypes = ['model', 'phil']  # also
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    pass

  def run(self, args):
    pass


  # end of class Program


# def rest_of_old_run(args):  -- an intermediate state, outdated technique
#   namespace, others = parser.parse_args(sys.argv[1:])
#   # whatever the old fashioned parser won't use becomes part of <others>
#   # and will be given to the phil parser
#
#   parser = CCTBXParser(
#     program_class=Program,
#     logger=logger)
#   namespace = parser.parse_args(others)
#
#   # start program
#   print('Starting job', file=logger)
#   print('='*79, file=logger)
#   phil1 = parser.working_phil.extract()
#   args2 = parseCommandLine()
#   task = Program(
#     parser.data_manager, phil1, logger=logger2)
#   main()

#=============================================================================
# def run(args):
#
#   # create parser
#   logger = multi_out()
#   logger.register('stderr', sys.stderr)
#   logger2 = multi_out()
#   logger2.register('stdout', sys.stdout)
#
#   parser = CCTBXParser(
#     program_class=cablam.Program,
#     logger=logger)
#   namespace = parser.parse_args(sys.argv[1:])
#
#   # start program
#   print('Starting job', file=logger)
#   print('='*79, file=logger)
#   task = cablam.Program(
#     parser.data_manager, parser.working_phil.extract(), logger=logger2)
#
#   # validate inputs
#   task.validate()
#
#   # run program
#   task.run()
#
#   # stop timer
#   print('', file=logger)
#   print('='*79, file=logger)
#   print('Job complete', file=logger)
#   show_total_time(out=logger)

# =============================================================================
if __name__ == '__main__':
  run(sys.argv[1:])
