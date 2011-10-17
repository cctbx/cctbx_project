
import libtbx.phil
from libtbx.utils import Sorry
import os
import sys

master_phil = libtbx.phil.parse("""
blast_pdb
  .caption = This program will run a BLAST search on the NCBI's web servers. \
    You may use any format sequence file, but only a single sequence may be \
    searched at a time.
  .short_caption = NCBI BLAST search of PDB
  .style = box auto_align
{
  file_name = None
    .type = path
    .style = bold input_file file_type:seq
  output_file = None
    .type = path
    .style = bold new_file
  blast_type = *blastp blastn
    .type = choice
    .caption = Protein_(blastp) Nucleotide_(blastn)
    .short_caption = Search type
}""")

def run (args=(), params=None, out=None) :
  if (out is None) :
    out = sys.stdout
  if (params is None) :
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      seq_file_def="blast_pdb.file_name")
    params = cmdline.work.extract()
  validate_params(params)
  params = params.blast_pdb
  from iotbx.pdb.fetch import get_ncbi_pdb_blast
  from iotbx.file_reader import any_file
  seq_file = any_file(params.file_name, force_type="seq")
  seq_file.check_file_type("seq")
  seq_objects = seq_file.file_object
  if (len(seq_objects) == 0) :
    raise Sorry("Empty sequence file!")
  elif (len(seq_objects) > 1) :
    print >> out, "WARNING: multiple sequences provided; searching only the 1st"
  sequence = seq_objects[0].sequence
  if (params.output_file is None) :
    params.output_file = "blast.xml"
  blast_out = get_ncbi_pdb_blast(sequence,
    file_name=params.output_file,
    blast_type=params.blast_type)
  print >> out, "Wrote results to %s" % params.output_file
  return True

def validate_params (params) :
  if (params.blast_pdb.file_name is None) :
    raise Sorry("A sequence file is required as input.")
  elif (not os.path.isfile(params.blast_pdb.file_name)) :
    raise Sorry("%s is not a file." % params.blast_pdb.file_name)
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
