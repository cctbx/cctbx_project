import libtbx.phil
from libtbx.utils import Sorry
import sys, os, re

def run(args=(), params=None, out=sys.stdout):
  import iotbx.pdb
  if (params is None) :
    for arg in args:
      assert os.path.isfile(arg)
      pdb_obj = iotbx.pdb.input(file_name=arg)
      hierarchy = pdb_obj.construct_hierarchy()
      for model in hierarchy.models():
        for chain in model.chains():
          for conformer in chain.conformers():
            f = conformer.format_fasta()
            if (f is not None):
              print >> out, "\n".join(f)
  else :
    params = params.pdb_as_fasta
    fasta_seqs = []
    for pdb_file in params.file_name :
      name_base = ""
      if (len(params.file_name) > 1) :
        name_base = "%s " % os.path.splitext(os.path.basename(pdb_file))[0]
      assert os.path.isfile(pdb_file)
      pdb_obj = iotbx.pdb.input(file_name=pdb_file)
      hierarchy = pdb_obj.construct_hierarchy()
      for model in hierarchy.models() :
        for chain in model.chains() :
          conformer = chain.conformers()[0]
          if (conformer.is_protein() or conformer.is_na(min_content=0.5)) :
            if (params.pad_missing_residues) :
              seq = conformer.as_padded_sequence(
                skip_insertions=(not params.include_insertion_residues))
            elif (not params.include_insertion_residues) :
              seq = conformer.as_padded_sequence(skip_insertions=True,
                pad=False)
            else :
              seq = "".join(conformer.as_sequence())
            if (params.ignore_missing_residues_at_start) :
              seq = re.sub("^X*", "", seq)
            seq_lines = []
            k = 0
            while (k < len(seq)) :
              if ((k % 70) == 0) :
                seq_lines.append("")
              seq_lines[-1] += seq[k]
              k += 1
            seq = "\n".join(seq_lines)
            fasta_seqs.append(">%schain '%2s'\n%s" % (name_base,chain.id, seq))
    f = open(params.output_file, "w")
    f.write("\n".join(fasta_seqs))
    f.close()
    return params.output_file

def validate_params (params) :
  if (len(params.pdb_as_fasta.file_name) == 0) :
    raise Sorry("No PDB files defined!")

# XXX for phenix gui
master_phil = libtbx.phil.parse("""
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
  output_file = pdb_sequences.fa
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
""")

if (__name__ == "__main__"):
  run(sys.argv[1:])
