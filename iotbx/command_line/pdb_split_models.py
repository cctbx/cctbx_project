"""separate a multi-model PDB file (such as an \
    NMR ensemble) into individual files for each model.  The output files \
    will be named similarly to the input file but ending in _1.pdb, _2.pdb, \
    etc.
"""

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb.split_models

from libtbx.utils import Sorry, Usage, null_out
import os
import sys

master_phil = """
split_models
  .short_caption = Split multi-model PDB file
  .caption = This utility will separate a multi-model PDB file (such as an \
    NMR ensemble) into individual files for each model.  The output files \
    will be named similarly to the input file but ending in _1.pdb, _2.pdb, \
    etc.
  .style = auto_align box caption_img:icons/custom/iotbx.pdb.join_fragment_files.png
{
  file_name = None
    .type = path
    .short_caption = PDB file
    .style = file_type:pdb input_file bold
  output_dir = None
    .type = path
    .short_caption = Output directory
    .style = directory default_cwd bold
}
"""

def run(args=(), params=None, out=None):
  if (out is None):
    out = sys.stdout
  if (params is None):
    if (len(args) == 0):
      raise Usage("""
iotbx.pdb.split_models ensemble.pdb [output_dir=/path/...]

Splits a multi-model PDB file into separate files for each model.
""")
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil_string=master_phil,
      pdb_file_def="split_models.file_name",
      directory_def="split_models.output_dir")
    params = cmdline.work.extract()
    validate_params(params)
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(params.split_models.file_name)
  hierarchy = pdb_in.construct_hierarchy()
  if (len(hierarchy.models()) <= 1):
    raise Sorry("The PDB file %s already has a single model." %
      params.split_models.file_name)
  pdb_rel_path = os.path.basename(params.split_models.file_name)
  if (pdb_rel_path.endswith(".gz")):
    pdb_rel_path = pdb_rel_path[:-3]
  elif (pdb_rel_path.endswith(".Z")):
    pdb_rel_path = pdb_rel_path[:-2]
  base_name = os.path.splitext(pdb_rel_path)[0]
  if (params.split_models.output_dir is None):
    params.split_models.output_dir = os.getcwd()
  output_base = os.path.join(params.split_models.output_dir, base_name)
  return split_models(
    hierarchy=hierarchy,
    crystal_symmetry=pdb_in.crystal_symmetry(),
    output_base=output_base,
    original_file=params.split_models.file_name,
    log=out)

def split_models(hierarchy,
                  crystal_symmetry,
                  output_base,
                  original_file=None,
                  log=None):
  if (log is None) : log = null_out()
  import iotbx.pdb.hierarchy
  n_models = len(hierarchy.models())
  file_names = []
  for k, model in enumerate(hierarchy.models()):
    k += 1
    new_hierarchy = iotbx.pdb.hierarchy.root()
    new_hierarchy.append_model(model.detached_copy())
    if (model.id == ""):
      model_id = str(k)
    else :
      model_id = model.id.strip()
    output_file = "%s_%s.pdb" % (output_base, model_id)
    final_fname = new_hierarchy.write_pdb_or_mmcif_file(
        target_filename=output_file,
        crystal_symmetry=crystal_symmetry)
    file_names.append(final_fname)
    print("Wrote %s" % final_fname, file=log)
  return file_names

def validate_params(params):
  if (params.split_models.file_name is None):
    raise Sorry("Please specify a PDB file to split!")
  elif (not os.path.isfile(params.split_models.file_name)):
    raise Sorry("The PDB file '%s' does not exist or is not a file." %
      params.split_models.file_name)
  if (params.split_models.output_dir is not None):
    if (not os.path.isdir(params.split_models.output_dir)):
      raise Sorry(("The specified output directory '%s' does not exist or is "+
        "not a directory.") % params.split_models.output_dir)
  return True

if (__name__ == "__main__"):
  run(sys.argv[1:])
