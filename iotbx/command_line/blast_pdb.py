
import libtbx.phil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import cStringIO
import os
import sys

master_phil = libtbx.phil.parse("""
blast_pdb
  .caption = This program will run a BLAST search on the NCBI's web servers. \
    You may use any format sequence file, but only a single sequence may be \
    searched at a time.  (Please limit your use of this service, as it is a \
    shared public resource!)
  .short_caption = NCBI BLAST search of PDB
  .style = box auto_align caption_img:icons/custom/pdb_import.png
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
  expect = 0.01
    .type = float
    .short_caption = E-value cutoff
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
    blast_type=params.blast_type,
    expect=params.expect)
  print >> out, "Wrote results to %s" % params.output_file
  results = summarize_blast_output(blast_out)
  if (len(args) != 0) : # command-line mode
    print >> out, ""
    print >> out, "%d matching structures" % len(results)
    print >> out, ""
    print >> out, "ID    Chain     evalue  length  %ident    %pos"
    print >> out, "-" * 46
    for result in results :
      result.show(out)
  if (len(results) > 0) :
    return sequence, os.path.abspath(params.output_file)
  else :
    return sequence, None

class blast_hit (object) :
  def __init__ (self, hit_num, pdb_id, chain_id, evalue, length, identity,
      positives) :
    adopt_init_args(self, locals())

  def show (self, out=None) :
    if (out is None) : out = sys.stdout
    print >> out, "%3s  %1s   %12g  %6d  %6.2f  %6.2f" % (self.pdb_id,
      self.chain_id, self.evalue, self.length, self.identity, self.positives)

# XXX this should probably use Bio.Blast.NCBIXML
def summarize_blast_output (blast_out=None, blast_file=None) :
  """
  Parse NCBI BLAST XML output and convert to a list of simple summary
  objects.  Note that this is very specific to searching the PDB, and returns
  incomplete information (suitable for summarizing in a flat table).
  """
  assert ([blast_out, blast_file].count(None) == 1)
  from Bio.Blast import NCBIXML
  if (blast_out is not None) :
    blast_in = cStringIO.StringIO(blast_out)
  else :
    assert os.path.isfile(blast_file)
    blast_in = open(blast_file)
  parsed = NCBIXML.parse(blast_in)
  blast = parsed.next()
  if (len(blast.alignments) == 0) :
    raise Sorry("No matching sequences!")
  results = []
  for i, hit in enumerate(blast.alignments) :
    pdb_chain_id = str(hit.accession)
    pdb_id, chain_id = pdb_chain_id.split("_")
    hsp = hit.hsps[0]
    assert (hsp.align_length > 0)
    summary = blast_hit(
      hit_num=i+1,
      pdb_id=pdb_id,
      chain_id=chain_id,
      evalue=hsp.expect,
      length=hsp.align_length,
      identity=100*hsp.identities/hsp.align_length,
      positives=100*hsp.positives/hsp.align_length)
    results.append(summary)
  return results

def validate_params (params) :
  if (params.blast_pdb.file_name is None) :
    raise Sorry("A sequence file is required as input.")
  elif (not os.path.isfile(params.blast_pdb.file_name)) :
    raise Sorry("%s is not a file." % params.blast_pdb.file_name)
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
