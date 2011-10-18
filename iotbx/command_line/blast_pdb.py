
import libtbx.phil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
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
  if (len(args) != 0) : # command-line mode
    results = summarize_blast_output(blast_out)
    print >> out, ""
    print >> out, "%d matching structures" % len(results)
    print >> out, ""
    print >> out, "ID    Chain     evalue  length  %ident    %pos"
    print >> out, "-" * 46
    for result in results :
      result.show(out)
  return os.path.abspath(params.output_file)

class blast_hit (object) :
  def __init__ (self, hit_num, pdb_id, chain_id, evalue, length, identity,
      positives) :
    adopt_init_args(self, locals())

  def show (self, out=None) :
    if (out is None) : out = sys.stdout
    print >> out, "%3s  %1s   %12g  %6d  %6.2f  %6.2f" % (self.pdb_id,
      self.chain_id, self.evalue, self.length, self.identity, self.positives)

def summarize_blast_output (blast_out=None, blast_file=None) :
  """
  Parse NCBI BLAST XML output and convert to a list of simple summary
  objects.  Note that this is very specific to searching the PDB, and returns
  incomplete information (suitable for summarizing in a flat table).
  """
  assert ([blast_out, blast_file].count(None) == 1)
  from xml.dom import minidom
  if (blast_out is None) :
    assert os.path.isfile(blast_file)
    blast_out = open(blast_file).read()
  def getText (node) :
    return str(node.childNodes[0].data)
  b = minidom.parseString(blast_out)
  hits = b.getElementsByTagName("Hit")
  results = []
  for hit in hits :
    hit_num = int(getText(hit.getElementsByTagName("Hit_num")[0]))
    pdb_chain_id = getText(hit.getElementsByTagName("Hit_accession")[0])
    pdb_id, chain_id = pdb_chain_id.split("_")
    hsp = hit.getElementsByTagName("Hsp")[0]
    evalue = float(getText(hsp.getElementsByTagName("Hsp_evalue")[0]))
    assert (evalue >= 0)
    aln_len = int(getText(hsp.getElementsByTagName("Hsp_align-len")[0]))
    assert (aln_len > 0)
    ident_ = int(getText(hsp.getElementsByTagName("Hsp_identity")[0]))
    pos_ = int(getText(hsp.getElementsByTagName("Hsp_positive")[0]))
    summary = blast_hit(
      hit_num=hit_num,
      pdb_id=pdb_id,
      chain_id=chain_id,
      evalue=evalue,
      length=aln_len,
      identity=100*ident_/aln_len,
      positives=100*pos_/aln_len)
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
