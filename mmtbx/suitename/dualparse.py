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

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import  with_statement, print_function
import argparse
from suiteninit import buildParser

from iotbx.cli_parser import CCTBXParser
import sys


# suitename first
def parseArgs(programClass, logger):
  for i, arg in enumerate(sys.argv):  # make case insensitive
    if arg[0]=="-":
      sys.argv[i] = arg.lower()
  parser = argparse.ArgumentParser()
  buildParser(parser)  # set up the arguments for suite name
  args, others = parser.parse_known_args()
  args2 = argsToPhilForm(args, others)
  args2.extend(others)

  # create the cctbx aspect of the parser
  parser2 = CCTBXParser(
    program_class = programClass,
    # custom_process_arguments = custom_args_proc,
    logger=logger
    )
  # do it
  namespace = parser2.parse_args(args2)
  extracted = parser2.working_phil.extract()
  # for name, value in extracted.suitename.__dict__.items():
    # print(name, ":", value)
  return parser2


class TableItem(object):
    main_name = ""
    alias_name = ""
    main_default = "false"
    alias_default = "false"

    def __init__(self, name1, name2, default1="", default2=""):
        self.main_name = name1
        self.alias_name = name2
        self.main_default = default1
        self.alias_default = default2


aliasList=(
('residuein', 'residuesin', False, False) ,
('suitein', 'suitesin', False, False) ,
('altid', 'altidval', 'A', 'A') ,
('pointidfields', 'ptid', 7, 0) ,
('etatheta', 'thetaeta', False, False),
(None, 'angles', None, 9),  # deprecated items
(None, 'resAngles', None, 6)
)

def argsToPhilForm(namespace, others):
  """
  Converting the parseArgs data to phil input format:
  The raw namespace as an entry for every argument that could possibly appear,
  including some that are aliases for others, some that are deprecated that
  we will not support in this context. We determine that an alias has been used
  if its value is different from its default. If so, the alias overrides the
  main item. In any case, only one of a main/alias pair will be passed on.

  The output of this function is a list of command line arguments suitable
  for phil.
  """
  aliasTable = [TableItem(*entry) for entry in aliasList]
  aliases = [item.alias_name for item in aliasTable]
  aliasDict = {aliases[i]: aliasTable[i] for i in range(len(aliasTable))}
  ###### res = {test_keys[i]: test_values[i] for i in range(len(test_keys))}
  dictionary = vars(namespace)
  for  key,value in dictionary.items(): # vars turns namespace to dict
    if key in aliasDict:  # if key is an alias or deprecated
      item = aliasDict[key]
      if item.main_name is None:
        pass
      if value != item.alias_default:
        # If the caller has actually specified a value for the alias,
        # lift the alias's value to the main entry
        dictionary[item.main_name] = value
      dictionary[key] = None
      # -- all aliases and deprecated items are suppressed from the dictionary

  argsOut = []
  for  key,value in dictionary.items(): # vars turns namespace to dict
    if value is not None and value != "":
      argsOut.append(str(key)+"="+str(value))
  return argsOut

# CCTBX first
def parseArgs1(Program, logger):
  # create the cc tbx aspect of the parser
  parser = CCTBXParser(
    program_class = Program,
    # custom_process_arguments = custom_args_proc,
    logger=logger)

  # add the old suitename aspect of the parser
  buildParser(parser)
  # do it
  namespace = parser.parse_args(sys.argv[1:])
  return parser.working_phil

