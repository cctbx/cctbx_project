from __future__ import division
from __future__ import print_function

# XXX most of this code is unused when run from the command line, but the
# PHENIX GUI includes a simple frontend that uses the phil interface.

import libtbx.phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
#from libtbx import easy_run
import sys
import os

from mmtbx import conformation_dependent_library
from mmtbx.conformation_dependent_library import cdl_utils
from mmtbx.conformation_dependent_library.cdl_database import cdl_database

master_phil = libtbx.phil.parse("""
cdl_lookup
  .caption = CDL
{
  residue_names = None
    .type = str
  residue_group_class = None
    .type = choice
  phi_psi_angles = None
    .type = str
}""")

def run2(args=(), out=sys.stdout):
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="cdl_lookup")
  phils = []
  phil_args = []
  pdbs = []
  for arg in args:
    if os.path.isfile(arg):
      if iotbx.pdb.is_pdb_file(arg):
        pdbs.append(arg)
        continue
      try :
        file_phil = phil.parse(file_name=arg)
      except RuntimeError :
        pass
      else :
        phils.append(file_phil)
    else :
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  working_phil.show()
  working_params = working_phil.extract()

  if working_params.cdl_lookup.residue_group_class is None:
    working_params.cdl_lookup.residue_group_class = cdl_utils.get_res_type_group(
      *tuple(working_params.cdl_lookup.residue_names.split(",")[1:])
      )
    print("\nPeptide triplet class : %s" % working_params.cdl_lookup.residue_group_class)
  key = working_params.cdl_lookup.phi_psi_angles.split(",")
  key[0] = int(key[0])
  key[1] = int(key[1])
  key = tuple(key)
  restraints_values = cdl_database[working_params.cdl_lookup.residue_group_class][key]
  outl = conformation_dependent_library.restraints_show(restraints_values)
  print("\nCDL values\n%s" % outl)
  return restraints_values

def validate_params(params):
  if (params.fetch_pdb.pdb_ids is None) or (len(params.fetch_pdb.pdb_ids)==0):
    raise Sorry("No PDB IDs specified!")
  return True

if __name__ == "__main__" :
  run2(sys.argv[1:])
