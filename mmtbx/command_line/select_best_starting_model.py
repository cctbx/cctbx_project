
from __future__ import absolute_import, division, print_function
import libtbx.phil
import sys
from iotbx import extract_xtal_data

master_phil_str = """
input {
  model = None
    .type = path
    .multiple = True
  include scope iotbx.extract_xtal_data.xray_data_str
  skip_twin_detection = False
    .type = bool
}
include scope mmtbx.refinement.select_best_starting_model.master_phil
nproc = Auto
  .type = int
output {
  write_files = False
    .type = bool
  model_file_name = best_model.pdb
    .type = path
  data_file_name = best_model_data.mtz
    .type = path
}
"""

def master_params():
  return libtbx.phil.parse(master_phil_str, process_includes=True)

def run(args, external_params=None, out=sys.stdout):
  from mmtbx.refinement import select_best_starting_model
  from iotbx.file_reader import any_file
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_params(),
    pdb_file_def="input.model",
    reflection_file_def="input.xray_data.file_name",
    usage_string="""\
mmtbx.select_best_starting_model data.mtz model1.pdb model2.pdb [...]

Given experimental data and a set of models, determines whether any of the
models can be refined directly (i.e. if isomorphous), and picks the one
with the best R-free below a specified cutoff.  Will optionally perform
rigid-body refinement on suitable models if requested.""")
  params = cmdline.work.extract()
  validate_params(params)
  hkl_in = any_file(params.input.xray_data.file_name)
  data_and_flags = extract_xtal_data.run(
    reflection_file_server=hkl_in.file_server,
    parameters=params.input.xray_data)
  model_data = []
  for file_name in params.input.model :
    model_in = cmdline.get_file(
      file_name=file_name,
      force_type="pdb").file_object
    pdb_hierarchy = model_in.hierarchy
    xray_structure = model_in.xray_structure_simple()
    model_data.append((pdb_hierarchy, xray_structure))
  # we can optionally pass a parameter block from elsewhere (e.g. phenix ligand
  # pipeline) and just use this run() method to load files
  params_ = external_params
  if (params_ is None):
    params_ = params
  result = select_best_starting_model.select_model(
    model_names=params.input.model,
    model_data=model_data,
    f_obs=data_and_flags.f_obs,
    r_free_flags=data_and_flags.r_free_flags,
    params=params_,
    nproc=params.nproc,
    skip_twin_detection=params.input.skip_twin_detection,
    log=out)
  if result.success() and params.output.write_files :
    result.save_best_model(file_name=params.output.model_file_name)
    print("", file=out)
    print("Wrote best model to %s" % params.output.model_file_name, file=out)
    result.save_updated_data(file_name=params.output.data_file_name)
    print("Wrote updated data to %s" % params.output.data_file_name, file=out)
  return result

def validate_params(params):
  if (len(params.input.model) == 0):
    raise Sorry("No models specified.")
  elif (params.input.xray_data.file_name is None):
    raise Sorry("No data file specified.")
  return True

if (__name__ == "__main__"):
  run(sys.argv[1:])
