"""Example of how to split a model by chains into new models"""
from __future__ import absolute_import, division, print_function

import libtbx.phil
from libtbx.utils import Sorry, Usage
import iotbx.pdb
from iotbx.pdb.hierarchy import new_hierarchy_from_chain

import os
import sys

master_phil = libtbx.phil.parse("""
split_chains
{
  pdb_file = None
    .type = path
    .short_caption = PDB file
    .style = file_type:pdb input_file
  output_dir = None
    .type = path
    .short_caption = Output directory
    .style = directory
  output_base = None
    .type = str
    .short_caption = Base file name
  exclude_heteroatoms = True
    .type = bool
  preserve_symmetry = True
    .type = bool
}
""")

def run(args=(), params=None, out=None):
  if (out is None) : out = sys.stdout
  if (params is None):
    if (len(args) == 0):
      raise Usage("pdb_split_chains.py model.pdb")
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      pdb_file_def="split_chains.pdb_file")
    params = cmdline.work.extract()
  validate_params(params)
  params = params.split_chains
  if (params.output_dir is None):
    params.output_dir = os.getcwd()
  if (params.output_base is None):
    params.output_base = os.path.basename(os.path.splitext(params.pdb_file)[0])
  # from iotbx import file_reader
  pdb_in = iotbx.pdb.input(params.pdb_file)
  cs = None
  if params.preserve_symmetry:
    cs = pdb_in.crystal_symmetry()
  hierarchy = pdb_in.construct_hierarchy()
  if (len(hierarchy.models()) > 1):
    raise Sorry("Multi-model PDB files are not supported.  You can use "+
      "iotbx.pdb.split_models to break the structure into individual model "+
      "files.")
  id_counts = {}
  outputs = []
  for chain in hierarchy.models()[0].chains():
    if (params.exclude_heteroatoms):
      if not chain.is_protein() and not chain.is_na():
        continue
      elif (len(chain.residue_groups()) == 1):
        continue
    id = chain.id
    if (id == " "):
      id = "_"
    if (id_counts.get(id, 0) > 0):
      suffix = "%s-%d" % (id, id_counts[id] + 1)
    else :
      suffix = id
    if (not id in id_counts):
      id_counts[id] = 0
    id_counts[id] += 1
    output_file = os.path.join(params.output_dir, "%s_%s.pdb" %
      (params.output_base, suffix))
    new_hierarchy = new_hierarchy_from_chain(chain)
    ofname = new_hierarchy.write_pdb_or_mmcif_file(target_filename=output_file, crystal_symmetry=cs)
    outputs.append(ofname)
    print("Wrote chain '%s' to %s" % (chain.id, ofname), file=out)
  return outputs

def validate_params(params):
  if (params.split_chains.pdb_file is None):
    raise Sorry("PDB file not defined!")

if (__name__ == "__main__"):
  run(sys.argv[1:])
