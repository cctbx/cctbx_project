
from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import os.path
import sys

master_phil = """
  model = None
    .type = path
  include scope mmtbx.refinement.select_best_starting_model.strip_model_params
  preserve_remarks = False
    .type = bool
  preserve_symmetry = True
    .type = bool
  output_file = None
    .type = bool
"""

def run(args, out=sys.stdout):
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil,
    pdb_file_def="model",
    usage_string="""\
mmtbx.strip_model model.pdb [options]

Prepare a model for use in MR or isomorphous substitution.  See also
phenix.pdbtools, which overlaps considerably in functionality.
""")
  params = cmdline.work.extract()
  if (params.model is None):
    raise Sorry("No model specified.")
  from mmtbx.refinement import select_best_starting_model
  output_file = params.output_file
  if (output_file is None):
    output_file = os.path.basename(os.path.splitext(params.model)[0]) + \
      "_modified.pdb"
  select_best_starting_model.strip_model(
    file_name=params.model,
    params=params,
    preserve_remarks=params.preserve_remarks,
    preserve_symmetry=params.preserve_symmetry,
    output_file=output_file,
    log=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
