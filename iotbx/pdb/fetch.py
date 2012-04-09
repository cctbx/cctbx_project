# TODO other PDB sites?
#
# See RCSB documentation at:
# http://www.rcsb.org/pdb/static.do?p=download/http/index.html
#
# File format  Compression   Example URL
# PDB  uncompressed http://www.rcsb.org/pdb/files/2vz8.pdb
# CIF  uncompressed http://www.rcsb.org/pdb/files/2vz8.cif
# XML  uncompressed http://www.rcsb.org/pdb/files/2vz8.xml
# Data uncompressed http://www.rcsb.org/pdb/files/2vz8-sf.cif
# CIF  uncompressed http://www.rcsb.org/pdb/files/ligand/ATP.cif
#
# For the PDBe:
# http://www.ebi.ac.uk/pdbe-srv/view/files/2vz8.ent
# http://www.ebi.ac.uk/pdbe-srv/view/files/r2vz8sf.ent

from libtbx.utils import Sorry
import libtbx.load_env
import urllib2
import urllib
import time
import re
import os

def validate_pdb_id (id) :
  if (len(id) != 4) or (not re.match("[1-9]{1}[a-zA-Z0-9]{3}", id)) :
    raise RuntimeError(("Invalid PDB ID '%s'.  IDs must be exactly four "+
      "alphanumeric characters, starting with a number from 1-9.") % id)

def validate_pdb_ids (id_list) :
  for id in id_list :
    try :
      validate_pdb_id(id)
    except RuntimeError, e :
      raise Sorry(str(e))

def fetch (id, data_type="pdb", format="pdb", mirror="rcsb") :
  assert data_type in ["pdb", "xray", "fasta"]
  assert format in ["cif", "pdb", "xml"]
  assert mirror in ["rcsb", "pdbe"]
  validate_pdb_id(id)

  id = id.lower()
  if (mirror == "rcsb") :
    url_base = "http://www.rcsb.org/pdb/files/"
    pdb_ext = ".pdb"
    sf_prefix = ""
    sf_ext = "-sf.cif"
  elif (mirror == "pdbe") :
    url_base = "http://www.ebi.ac.uk/pdbe-srv/view/files/"
    pdb_ext = ".ent"
    sf_prefix = "r"
    sf_ext = "sf.ent"
  if data_type == "fasta" :
    # XXX the RCSB doesn't appear to have a simple URL for FASTA files
    url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=FASTA&compression=NO&structureId=%s" % id
    try :
      data = urllib2.urlopen(url)
    except urllib2.HTTPError, e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download sequence for %s." % id)
      else :
        raise
  elif data_type == "xray" :
    url = url_base + sf_prefix + id + sf_ext
    try :
      data = urllib2.urlopen(url)
    except urllib2.HTTPError, e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download structure factors for %s." % id)
      else :
        raise
  else :
    if format == "pdb" :
      url = url_base + id + pdb_ext
    else :
      url = url_base + id + "." + format
    try :
      data = urllib2.urlopen(url)
    except urllib2.HTTPError, e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download model for %s." % id)
      else :
        raise
  return data

def get_pdb (id, data_type, mirror, log, quiet=False) :
  try :
    data = fetch(id, data_type, mirror=mirror)
  except RuntimeError, e :
    raise Sorry(str(e))
  file_name = None
  if data_type == "xray" :
    file_name = os.path.join(os.getcwd(), "%s-sf.cif" % id)
    open(file_name, "w").write(data.read())
    if not quiet :
      print >> log, "Structure factors saved to %s" % file_name
  elif data_type == "fasta" :
    file_name = os.path.join(os.getcwd(), "%s.fa" % id)
    open(file_name, "w").write(data.read())
    if not quiet :
      print >> log, "Sequence saved to %s" % file_name
  else :
    file_name = os.path.join(os.getcwd(), "%s.pdb" % id)
    open(file_name, "w").write(data.read())
    if not quiet :
      print >> log, "Model saved to %s" % file_name
  return file_name

def get_ncbi_pdb_blast (sequence, file_name=None, blast_type="blastp",
    expect=0.01) :
  assert (blast_type in ["blastp", "blastn"])
  if (sequence[-1] == '*') :
    sequence = sequence[:-1]
  if (not sequence.isalpha()) :
    raise Sorry("The sequence contains non-alphabetical characters; in "+
      "addition to A-Z, only an asterisk denoting a stop codon is permitted.")
  assert (expect >= 0)
  try :
    from Bio.Blast import NCBIWWW
  except ImportError :
    raise Sorry("You need to have BioPython installed to use this function.")
  blast = NCBIWWW.qblast(blast_type, "pdb", sequence)
  blast_out = blast.read()
  if (file_name is not None) :
    f = open(file_name, "w")
    f.write(blast_out)
    f.close()
  return blast_out

def get_ebi_pdb_wublast (sequence, email, file_name=None, blast_type="blastp",
    sequence_type="protein", exp="1e-3") :
  assert (email is not None)
  url = "http://www.ebi.ac.uk/Tools/services/rest/wublast/run/"
  params = urllib.urlencode({
    'sequence': sequence,
    'program' : program,
    'email'   : email,
    'exp'     : exp,
    'database': 'pdb',
    'stype'   : 'protein',
  })
  job_id = urllib.urlopen(url, params).read()
  while (True) :
    time.sleep(1)
    url = "http://www.ebi.ac.uk/Tools/services/rest/wublast/status/%s" % job_id
    status = urllib.urlopen(url).read()
    if (status == "RUNNING") :
      continue
    elif (status == "FINISHED") :
      url = "http://www.ebi.ac.uk/Tools/services/rest/wublast/result/%s/xml" %\
        job_id
      result = urllib.urlopen(url).read()
      return result
    elif (status == "ERROR") :
      raise RuntimeError("The EBI server reported an error.")
    elif (status == "FAILURE") :
      raise Sorry("Search failed!")
    elif (status == "NOT_FOUND") :
      raise RuntimeError("The EBI server can't find the job!")
    else :
      raise RuntimeError("Unknown status %s" % status)

def get_chemical_components_cif (code, return_none_if_already_present=False) :
  assert (code is not None)
  if (len(code) == 0) or (len(code) > 3) :
    raise Sorry("PDB residue codes must be at least 1 but no more than 3 "+
      "characters.")
  first_char = code[0].lower()
  code = code.upper()
  chem_comp_cif = libtbx.env.find_in_repositories(
    relative_path="chem_data/chemical_components/%s/data_%s.cif" % (first_char,
      code),
    test=os.path.isfile)
  if (chem_comp_cif is None) :
    url = "http://www.rcsb.org/pdb/files/ligand/%s.cif" % code
    try :
      data = urllib2.urlopen(url)
    except urllib2.HTTPError, e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download sequence for %s." % id)
      else :
        raise
    else :
      file_name = "%s.cif" % code
      open(file_name, "w").write(data.read())
      return file_name
  elif (not return_none_if_already_present) :
    return chem_comp_cif
  return None
