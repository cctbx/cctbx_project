
from __future__ import division
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry, Usage
import os
import sys

master_phil_str = """
include scope mmtbx.utils.cmdline_input_phil_str
ligand_code = None
  .type = str
  .multiple = True
reference_structure = None
  .type = path
only_segid = None
  .type = str
verbose = False
  .type = bool
"""

def run (args, out=sys.stdout) :
  if (len(args) == 0) or ("--help" in args) :
    raise Usage("""\
mmtbx.validate_ligands model.pdb data.mtz LIGAND_CODE [...]

Print out basic statistics for residue(s) with the given code(s), including
electron density values/CC.
""")
  import mmtbx.validation.ligands
  import mmtbx.utils
  import iotbx.phil
  master_phil = iotbx.phil.parse(master_phil_str, process_includes=True)
  args_ = []
  for arg in args :
    if (len(arg) == 3) and arg.isalnum() and (not os.path.exists(arg)) :
      args_.append("ligand_code=%s" % arg)
    else :
      args_.append(arg)
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args_,
    master_phil=master_phil,
    process_pdb_file=False)
  params = cmdline.params
  if (params.ligand_code is None) or (len(params.ligand_code) == 0) :
    raise Sorry("Ligand code required!")
  make_sub_header("Validating ligands", out=out)
  for ligand_code in params.ligand_code :
    validations = mmtbx.validation.ligands.validate_ligands(
      pdb_hierarchy=cmdline.pdb_hierarchy,
      fmodel=cmdline.fmodel,
      ligand_code=ligand_code,
      reference_structure=params.reference_structure,
      only_segid=params.only_segid)
    if (validations is None) :
      raise Sorry("No ligands named '%s' found." % ligand_code)
    mmtbx.validation.ligands.show_validation_results(validations=validations,
      out=out,
      verbose=params.verbose)

if (__name__ == "__main__") :
  run(sys.argv[1:])
