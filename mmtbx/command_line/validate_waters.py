
from __future__ import division
from libtbx.str_utils import make_sub_header
import libtbx.phil
import sys

master_phil = libtbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
""", process_includes=True)

def run (args, out=sys.stdout) :
  from mmtbx.validation import waters
  import mmtbx.utils
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
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

if (__name__ == "__main__") :
  run(sys.argv[1:])
