
import libtbx.phil
from libtbx.utils import Sorry, Usage
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

def run (args=(), params=None, out=None) :
  if (out is None) : out = sys.stdout
  if (params is None) :
    if (len(args) == 0) :
      raise Usage("pdb_split_chains.py model.pdb")
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      pdb_file_def="split_chains.pdb_file")
    params = cmdline.work.extract()
  validate_params(params)
  params = params.split_chains
  if (params.output_dir is None) :
    params.output_dir = os.getcwd()
  if (params.output_base is None) :
    params.output_base = os.path.basename(os.path.splitext(params.pdb_file)[0])
  from iotbx import file_reader
  from iotbx.pdb.hierarchy import new_hierarchy_from_chain
  pdb_in = file_reader.any_file(params.pdb_file, force_type="pdb")
  pdb_in.check_file_type("pdb")
  symm = pdb_in.file_object.crystallographic_section()
  hierarchy = pdb_in.file_object.construct_hierarchy()
  if (len(hierarchy.models()) > 1) :
    raise Sorry("Multi-model PDB files are not supported.  You can use "+
      "iotbx.pdb.split_models to break the structure into individual model "+
      "files.")
  id_counts = {}
  outputs = []
  for chain in hierarchy.models()[0].chains() :
    main_conf = chain.conformers()[0]
    if (params.exclude_heteroatoms) :
      if ((not main_conf.is_protein()) and (not main_conf.is_na())) :
        continue
      elif (len(chain.residue_groups()) == 1) :
        continue
    id = chain.id
    if (id == " ") :
      id = "_"
    if (id_counts.get(id, 0) > 0) :
      suffix = "%s-%d" % (id, id_counts[id] + 1)
    else :
      suffix = id
    if (not id in id_counts) :
      id_counts[id] = 0
    id_counts[id] += 1
    output_file = os.path.join(params.output_dir, "%s_%s.pdb" %
      (params.output_base, suffix))
    f = open(output_file, "w")
    if (params.preserve_symmetry) :
      f.write("\n".join(symm))
      f.write("\n")
    new_hierarchy = new_hierarchy_from_chain(chain)
    f.write(new_hierarchy.as_pdb_string())
    f.close()
    outputs.append(output_file)
    print >> out, "Wrote chain '%s' to %s" % (chain.id, output_file)
  return outputs

def validate_params (params) :
  if (params.split_chains.pdb_file is None) :
    raise Sorry("PDB file not defined!")

if (__name__ == "__main__") :
  run(sys.argv[1:])
