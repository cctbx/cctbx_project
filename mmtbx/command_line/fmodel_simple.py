
from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
import sys

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    phil_string="""\
output_file = fmodel.pkl
  .type = path
""")

def run(args, out=sys.stdout):
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    process_pdb_file=False,
    out=out,
    usage_string="mmtbx.fmodel_simple model.pdb data.mtz [options]")
  fmodel = cmdline.fmodel
  fmodel_info = fmodel.info()
  fmodel_info.show_rfactors_targets_scales_overall(out=out)
  easy_pickle.dump(cmdline.params.output_file, fmodel)
  print("Wrote fmodel to %s" % cmdline.params.output_file, file=out)
  return fmodel

if (__name__ == "__main__"):
  run(sys.argv[1:])
