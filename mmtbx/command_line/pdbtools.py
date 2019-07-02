from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdbtools

import mmtbx.pdbtools
import sys, os
from libtbx import runtime_utils
import mmtbx.utils
import mmtbx.model
import iotbx.pdb
from iotbx.pdb import combine_unique_pdb_files
from iotbx.phil import process_command_line_with_files
from scitbx.array_family import flex
from libtbx import group_args
from libtbx.utils import check_if_output_directory_exists

master_params_str = """\
model_file_name = None
  .type = path
  .multiple = True
  .help = Model file name
  .short_caption = Input model
  .style = bold input_file file_type:pdb file_type_default
output
  .help = Write out file with modified model (file name is defined in \
          write_modified)
  .style = noauto
{
  file_name=None
    .type=path
    .input_size=400
    .short_caption = Output model file
    .help = Default is the original file name with the file extension \
            replaced by "_modified.pdb".
    .style = bold new_file file_type:pdb
  format = *pdb *mmcif
    .type = choice(multi=True)
    .help = Choose the output format of coordinate file (PDB or mmCIF)
}
include scope mmtbx.pdbtools.master_params
include scope libtbx.phil.interface.tracking_params
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def broadcast(m, log):
  print("-"*79, file=log)
  print(m, file=log)
  print("*"*len(m), file=log)

def get_inputs(args, log, master_params):
  """
  Eventually, this will be centralized.
  """
  cmdline = process_command_line_with_files(
    args         = args,
    master_phil  = master_params,
    pdb_file_def = 'model_file_name'
  )
  params = cmdline.work.extract()
  # Model
  file_names = params.model_file_name
  pdb_combined = combine_unique_pdb_files(file_names = file_names)
  pdb_inp = iotbx.pdb.input(
    source_info = None,
    lines       = flex.std_string(pdb_combined.raw_records),
    raise_sorry_if_format_error = True)
  # Crystal symmetry
  fake_crystal_symmetry = False
  crystal_symmetry = pdb_inp.crystal_symmetry()
  if(crystal_symmetry is None or
     crystal_symmetry.is_empty() or
     crystal_symmetry.is_nonsence()):
    fake_crystal_symmetry = True
    from cctbx import uctbx
    crystal_symmetry = \
      uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
        sites_cart=pdb_inp.atoms().extract_xyz(),
        buffer_layer=5).crystal_symmetry()
  #
  model = mmtbx.model.manager(
      model_input = pdb_inp,
      crystal_symmetry = crystal_symmetry)
  #
  return group_args(
    params                = params,
    pdb_file_names        = file_names,
    model                 = model,
    fake_crystal_symmetry = fake_crystal_symmetry)

def run(args, out=sys.stdout, replace_stderr=True):
  log = mmtbx.utils.set_log(args, out=out, replace_stderr=replace_stderr)
  broadcast(m="phenix.pdbtools tools for PDB model manipulations.", log=log)
  if(len(args)==0):
    master_params().show(out = log)
    return None
  inputs = get_inputs(args=args, log=log, master_params=master_params())
  ### get i/o file names
  ofn = inputs.params.output.file_name
  ifn = inputs.pdb_file_names
  if(ofn is None):
    if(len(ifn)==1): ofn = os.path.basename(ifn[0]) + "_modified."
    elif(len(ifn)>1): ofn = os.path.basename(ifn[0]) + "_et_al_modified."
    else:
      pdbout = os.path.basename(inputs.pdb_file_names[0])
      ofn = pdbout+"_modified."
  else:
    if ofn[-4:] == ".pdb" or ofn[-4:] == ".cif":
      ofn = ofn[:-3]
  output_pdb_name = ofn + 'pdb'
  output_cif_name = ofn + 'cif'
  # Show parameters
  broadcast(m="Complete set of parameters", log=log)
  master_params().format(inputs.params).show(out = log)
  print(file=log)
  # Run calcs
  broadcast(m="Performing manipulations", log=log)
  task_obj = mmtbx.pdbtools.modify(
    model          = inputs.model,
    params         = inputs.params.modify,
    log            = log)
  results = task_obj.get_results()
  #
  broadcast(m="Writing output model", log=log)
  if "pdb" in inputs.params.output.format:
    print("Output model file name: ", output_pdb_name, file=log)
    with open(output_pdb_name, 'w') as f:
      f.write(inputs.model.model_as_pdb(output_cs= not inputs.fake_crystal_symmetry))
  if "mmcif" in inputs.params.output.format:
    print("Output model file name: ", output_cif_name, file=log)
    with open(output_cif_name, 'w') as f:
      f.write(inputs.model.model_as_mmcif(output_cs= not inputs.fake_crystal_symmetry))
  broadcast(m="All done.", log=log)
  return None

def validate_params(params, callback=None):
  if (params.model_file_name is None):
    raise Sorry("No PDB file(s) specified.")
  elif (params.output.file_name is not None):
    if (os.path.isdir(params.output.file_name)):
      raise Sorry("The specified output file is a currently existing "+
                  "directory.")
    check_if_output_directory_exists(params.output.file_name)
  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    return run(args=list(self.args), out=sys.stdout)

if (__name__ == "__main__"):
  assert run(args=sys.argv[1:]) is None # assert here is intentional
