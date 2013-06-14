
from __future__ import division
import libtbx.phil
import sys

master_phil = libtbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
""", process_includes=True)

def run (args, out=sys.stdout) :
  import mmtbx.utils
  from mmtbx.validation.waters import analyze_waters
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    process_pdb_file=False,
    out=out)
  waters = analyze_waters(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    xray_structure=cmdline.xray_structure,
    fmodel=cmdline.fmodel)
  return waters

if (__name__ == "__main__") :
  run(sys.argv[1:])
