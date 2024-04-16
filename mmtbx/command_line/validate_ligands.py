from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.validate_ligands
from iotbx.cli_parser import run_program
from mmtbx.programs import validate_ligands

if __name__ == '__main__':
  run_program(program_class=validate_ligands.Program)

#old stuff

#from __future__ import absolute_import, division, print_function
#from libtbx.str_utils import make_sub_header
#from libtbx.utils import Sorry
#import os
#import sys
#
#def master_phil():
#  from mmtbx.command_line import generate_master_phil_with_inputs
#  return generate_master_phil_with_inputs(
#    enable_automatic_twin_detection=True,
#    phil_string="""
#ligand_code = None
#  .type = str
#  .multiple = True
#reference_structure = None
#  .type = path
#only_segid = None
#  .type = str
#verbose = False
#  .type = bool
#""")
#
#def run(args, out=sys.stdout):
#  usage_string = """\
#mmtbx.validate_ligands model.pdb data.mtz LIGAND_CODE [...]
#
#Print out basic statistics for residue(s) with the given code(s), including
#electron density values/CC.
#"""
#  import mmtbx.validation.ligands
#  import mmtbx.command_line
#  args_ = []
#  for arg in args :
#    if (len(arg) == 3) and arg.isalnum() and (not os.path.exists(arg)):
#      args_.append("ligand_code=%s" % arg)
#    else :
#      args_.append(arg)
#  cmdline = mmtbx.command_line.load_model_and_data(
#    args=args_,
#    master_phil=master_phil(),
#    process_pdb_file=False,
#    usage_string=usage_string)
#  params = cmdline.params
#  if (params.ligand_code is None) or (len(params.ligand_code) == 0):
#    raise Sorry("Ligand code required!")
#  make_sub_header("Validating ligands", out=out)
#  for ligand_code in params.ligand_code :
#    validations = mmtbx.validation.ligands.validate_ligands(
#      pdb_hierarchy=cmdline.pdb_hierarchy,
#      fmodel=cmdline.fmodel,
#      ligand_code=ligand_code,
#      reference_structure=params.reference_structure,
#      only_segid=params.only_segid)
#    if (validations is None):
#      raise Sorry("No ligands named '%s' found." % ligand_code)
#    mmtbx.validation.ligands.show_validation_results(validations=validations,
#      out=out,
#      verbose=params.verbose)
#
#if (__name__ == "__main__"):
#  run(sys.argv[1:])
