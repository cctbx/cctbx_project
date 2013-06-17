
from __future__ import division
import iotbx.phil
from libtbx import easy_pickle
import sys

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
output_file = fmodel.pkl
  .type = path
""", process_includes=True)

def run (args, out=sys.stdout) :
  import mmtbx.utils
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    update_f_part1_for="refinement",
    args=args,
    master_phil=master_phil,
    process_pdb_file=False,
    out=out,
    usage_string="mmtbx.fmodel_simple model.pdb data.mtz [options]")
  fmodel = cmdline.fmodel
  fmodel_info = fmodel.info()
  fmodel_info.show_rfactors_targets_scales_overall(out=out)
  easy_pickle.dump(cmdline.params.output_file, fmodel)
  print >> out, "Wrote fmodel to %s" % cmdline.params.output_file
  return fmodel

if (__name__ == "__main__") :
  run(sys.argv[1:])
