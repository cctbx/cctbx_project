"""Tools to run BLAST searches"""

from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
from libtbx import slots_getstate_setstate_default_initializer
import libtbx.utils
from six.moves import urllib
from six.moves import cStringIO as StringIO
import time, os

def get_ncbi_pdb_blast(sequence, file_name=None, blast_type="blastp",
    expect=0.01):
  """
  Run BLAST against the PDB sequence database on the NCBI's web server.
  Basically just a frontend to a BioPython module.

  :param sequence: plaintext amino-acid or nucleotide sequence
  :param file_name: optional output file name
  :param blast_type: program to run ("blastp" or "blastn")
  :param expect: BLAST e-value cutoff

  :returns: XML string
  """
  assert (blast_type in ["blastp", "blastn"])
  if (sequence[-1] == '*'):
    sequence = sequence[:-1]
  if (not sequence.isalpha()):
    raise Sorry("The sequence contains non-alphabetical characters; in "+
      "addition to A-Z, only an asterisk denoting a stop codon is permitted.")
  assert (expect >= 0)
  try :
    from Bio.Blast import NCBIWWW
  except ImportError :
    raise Sorry("You need to have BioPython installed to use this function.")
  # FIXME will this use the HTTP proxy if defined?
  blast = NCBIWWW.qblast(blast_type, "pdb", sequence, expect=expect)
  blast_out = blast.read()
  if (file_name is not None):
    f = open(file_name, "w")
    f.write(blast_out)
    f.close()
  return blast_out

class blast_hit(slots_getstate_setstate_default_initializer):
  __slots__ = [ "hit_num", "pdb_id", "chain_id", "evalue", "length",
      "identity", "positives", "all_ids", ]

  def show(self, out=None):
    if (out is None) : out = sys.stdout
    print("%3s  %1s   %12g  %6d  %6.2f  %6.2f  %4d" % (self.pdb_id,
      self.chain_id, self.evalue, self.length, self.identity, self.positives,
      len(self.all_ids)), file=out)

def summarize_blast_output(blast_out=None, blast_file=None,
    min_identity=None, expect=None, stop_if_no_alignment=True):
  """
  Parse NCBI BLAST XML output and convert to a list of simple summary
  objects.  Note that this is very specific to searching the PDB, and returns
  incomplete information (suitable for summarizing in a flat table).
  """
  assert ([blast_out, blast_file].count(None) == 1)
  from Bio.Blast import NCBIXML
  import iotbx.pdb.fetch
  if (blast_out is not None):
    blast_in = StringIO(blast_out)
  else :
    assert os.path.isfile(blast_file)
    blast_in = open(blast_file)
  parsed = NCBIXML.parse(blast_in)
  blast = next(parsed)
  if (len(blast.alignments) == 0):
    if stop_if_no_alignment:
      raise Sorry("No matching sequences!")
    else: return list()
  results = []
  for i_hit, hit in enumerate(blast.alignments):
    pdb_chain_id = str(hit.accession)
    #hit.accession may only have pdb_id, e.g. 1EMB
    if len(pdb_chain_id.split("_")) > 1:
      pdb_id, chain_id = pdb_chain_id.split("_")
    else:
      pdb_id = pdb_chain_id
      chain_id = None
    #
    hsp = hit.hsps[0]
    assert (hsp.align_length > 0)
    identity = 100 * hsp.identities / hsp.align_length
    if (min_identity is not None) and (identity < min_identity):
      continue
    # XXX this is really appalling, but the NCBI groups together identical
    # sequences in its BLAST output, so I need to parse the accession code
    # strings to extract the individual PDB IDs
    hit_def_fields = hit.hit_def.split("|")
    all_ids = []
    all_ids.append([pdb_id,chain_id])
    for i_field, field in enumerate(hit_def_fields):
      if (field == "pdb") and (i_field < len(hit_def_fields) -1):
        next_pdb_id = hit_def_fields[i_field + 1]
        if "Chain" in hit_def_fields[i_field + 2]:
          next_chain_id = hit_def_fields[i_field + 2].split()[0]
        else:
          next_chain_id = None
        if (iotbx.pdb.fetch.valid_pdb_id(next_pdb_id)):
          all_ids.append([next_pdb_id,next_chain_id])
    summary = blast_hit(
      hit_num=i_hit+1,
      pdb_id=pdb_id,
      chain_id=chain_id,
      evalue=hsp.expect,
      length=hsp.align_length,
      identity=identity,
      positives=100*hsp.positives/hsp.align_length,
      hsp = hsp,
      all_ids=all_ids)
    results.append(summary)
  return results

def get_all_matching_pdb_ids(sequence, min_identity):
  blast_out = get_ncbi_pdb_blast(sequence=sequence)
  results = summarize_blast_output(blast_out=blast_out,
    min_identity=min_identity)
  ids = []
  for result in results :
    ids.extend([ pdb_chain_id[0].lower() for pdb_chain_id in result.all_ids ])
  return ids

def get_ebi_pdb_wublast(sequence, email, file_name=None, blast_type="blastp",
    sequence_type="protein", exp="1e-3"):
  """
  Run WU-BLAST against the PDB sequence database on the EBI's web server.
  Somewhat more complicated than the NCBI BLAST service, because of two-step
  process (submission and retrieval).  An email address is required to submit
  a job.

  :param sequence: plaintext amino-acid or nucleotide sequence
  :param email: user email address
  :param file_name: optional output file name
  :param blast_type: program to run ("blastp" or "blastn")
  :param sequence_type: currently ignored
  :param exp: BLAST e-value cutoff

  :returns: XML string
  """
  assert (email is not None)
  url = "http://www.ebi.ac.uk/Tools/services/rest/wublast/run/"
  params = urllib.parse.urlencode({
    'sequence': sequence,
    'program' : program,
    'email'   : email,
    'exp'     : exp,
    'database': 'pdb',
    'stype'   : 'protein',
  })
  job_id = libtbx.utils.urlopen(url, params).read()
  while (True):
    time.sleep(1)
    url = "http://www.ebi.ac.uk/Tools/services/rest/wublast/status/%s" % job_id
    status = libtbx.utils.urlopen(url).read()
    if (status == "RUNNING"):
      continue
    elif (status == "FINISHED"):
      url = "http://www.ebi.ac.uk/Tools/services/rest/wublast/result/%s/xml" %\
        job_id
      result = libtbx.utils.urlopen(url).read()
      return result
    elif (status == "ERROR"):
      raise RuntimeError("The EBI server reported an error.")
    elif (status == "FAILURE"):
      raise Sorry("Search failed!")
    elif (status == "NOT_FOUND"):
      raise RuntimeError("The EBI server can't find the job!")
    else :
      raise RuntimeError("Unknown status %s" % status)
