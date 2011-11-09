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
#
# For the PDBe:
# http://www.ebi.ac.uk/pdbe-srv/view/files/2vz8.ent
# http://www.ebi.ac.uk/pdbe-srv/view/files/r2vz8sf.ent

from libtbx.utils import Sorry, Usage
from libtbx import easy_run
import urllib2
import urllib
import time
import re
import os
import sys

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

def run (args, log=sys.stdout) :
  if len(args) < 1 :
    raise Usage("phenix.fetch_pdb [-x|-f|--all] [--mtz] [-q] ID1 [ID2, ...]")
  quiet = False
  convert_to_mtz = False
  data_type = "pdb"
  ids = []
  for arg in args :
    if (arg == "--all") :
      data_type = "all"
    elif (arg == "-x") :
      data_type = "xray"
    elif (arg == "-f") :
      data_type = "fasta"
    elif (arg == "-q") :
      quiet = True
    elif (arg == "--mtz") :
      convert_to_mtz = True
    elif (arg.startswith("--mirror=")) :
      mirror = arg.split("=")[1]
      if (not mirror in ["rcsb", "pdbe"]) :
        raise Sorry("Unrecognized mirror site '%s' (choices: rcsb, pdbe)" %
          mirror)
    else :
      ids.append(arg)
  if (len(ids) == 0) :
    raise Sorry("No PDB IDs specified.")
  mirror = "rcsb"
  if (data_type != "all") :
    files = []
    for id in ids :
      files.append(get_pdb(id, data_type, mirror, log))
    if (len(files) == 1) :
      return files[0]
    return files
  else :
    files = []
    for id in ids :
      for data_type_ in ["pdb", "fasta", "xray"] :
        files.append(get_pdb(id, data_type_, mirror, log))
      if (convert_to_mtz) :
        easy_run.call("phenix.cif_as_mtz %s-sf.cif --symmetry=%s.pdb --merge" %
          (id,id))
        if os.path.isfile("%s-sf.mtz" % id) :
          os.rename("%s-sf.mtz" % id, "%s.mtz" % id)
          os.remove("%s-sf.cif" % id)
        files[-1] = os.path.abspath("%s.mtz" % id)
    return files

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
