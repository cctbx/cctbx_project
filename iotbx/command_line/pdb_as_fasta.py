"""Create fasta sequence file from a model file (PDB or mmCIF format)"""

from __future__ import absolute_import, division, print_function
import libtbx.phil
from libtbx.utils import Sorry
import sys, os, re

def run(args=(), params=None, out=sys.stdout):
  import iotbx.pdb
  import iotbx.phil
  if (params is None):
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil_string=master_phil_str_template % "None",
      pdb_file_def="pdb_as_fasta.file_name")
    params = cmdline.work.extract()
    validate_params(params)
  params = params.pdb_as_fasta
  fasta_seqs = []
  for pdb_file in params.file_name :
    name_base = ""
    if (len(params.file_name) > 1):
      name_base = "%s " % os.path.splitext(os.path.basename(pdb_file))[0]
    assert os.path.isfile(pdb_file)
    pdb_obj = iotbx.pdb.input(file_name=pdb_file)
    hierarchy = pdb_obj.construct_hierarchy()
    for model in hierarchy.models():
      for chain in model.chains():
        if (chain.is_protein() or chain.is_na(min_content=0.5)):
          if (params.pad_missing_residues):
            seq = chain.as_padded_sequence(
              skip_insertions=(not params.include_insertion_residues))
          elif (not params.include_insertion_residues):
            seq = chain.as_padded_sequence(skip_insertions=True,
              pad=False)
          else :
            seq = "".join(chain.as_sequence())
          if (params.ignore_missing_residues_at_start):
            seq = re.sub("^X*", "", seq)
          seq_lines = []
          k = 0
          while (k < len(seq)):
            if ((k % 70) == 0):
              seq_lines.append("")
            seq_lines[-1] += seq[k]
            k += 1
          seq = "\n".join(seq_lines)
          fasta_seqs.append(">%schain '%2s'\n%s" % (name_base,chain.id, seq))
  if (params.output_file is not None):
    f = open(params.output_file, "w")
    f.write("\n".join(fasta_seqs))
    f.close()
  else :
    print("\n".join(fasta_seqs))
  return params.output_file

def validate_params(params):
  if (len(params.pdb_as_fasta.file_name) == 0):
    raise Sorry("No PDB files defined!")

master_phil_str_template = """
pdb_as_fasta
  .short_caption = Extract sequence from PDB file(s)
  .caption = Only protein and nucleic acid chains are supported by this \
    program.  Except for selenomethionine (MSE), non-standard residues and \
    bases will be replaced by 'X' and 'N', respectively.
  .style = caption_img:icons/custom/phenix.pdbtools.png
{
  file_name = None
    .type = path
    .multiple = True
    .short_caption = PDB file
    .style = bold noauto file_type:pdb
  output_file = %s
    .type = path
    .style = bold noauto file_type:seq new_file
  pad_missing_residues = True
    .type = bool
    .short_caption = \"\"\"Use 'X' in place of missing residues\"\"\"
    .style = bold noauto
  ignore_missing_residues_at_start = False
    .type = bool
    .style = bold noauto
  include_insertion_residues = True
    .type = bool
    .short_caption = Include insertion residues
    .style = bold noauto
}
"""

# XXX for phenix gui
master_phil = libtbx.phil.parse(master_phil_str_template % "pdb_sequences.fa")

if (__name__ == "__main__"):
  run(sys.argv[1:])
