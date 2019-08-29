
from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_sub_header
import sys

def run(args, out=sys.stdout):
  from mmtbx.validation import waters
  import mmtbx.command_line
  master_phil = mmtbx.command_line.generate_master_phil_with_inputs("")
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil,
    process_pdb_file=False,
    out=out)
  result = waters.waters(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    xray_structure=cmdline.xray_structure,
    fmodel=cmdline.fmodel,
    collect_all=True)
  make_sub_header("Solvent analysis", out=out)
  result.show(out=out)
  return result

if (__name__ == "__main__"):
  run(sys.argv[1:])
