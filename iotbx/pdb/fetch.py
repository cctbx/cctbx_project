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

import sys, os, re
import urllib2
from libtbx.utils import Sorry, Usage

def fetch (id, data_type="pdb", format="pdb") :
  assert data_type in ["pdb", "xray", "fasta"]
  assert format in ["cif", "pdb", "xml"]
  if (len(id) != 4) or (not re.match("[1-9]{1}[a-zA-Z0-9]{3}", id)) :
    raise RuntimeError(("Invalid PDB ID '%s'.  IDs must be exactly four "+
      "alphanumeric characters, starting with a number from 1-9.") % id)

  id = id.lower()
  url_base = "http://www.rcsb.org/pdb/files/"
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
    url = url_base + id + "-sf.cif"
    try :
      data = urllib2.urlopen(url)
    except urllib2.HTTPError, e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download structure factors for %s." % id)
      else :
        raise
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

def run (args, log=sys.stdout) :
  if len(args) < 1 :
    raise Usage("phenix.fetch_pdb [-x|-f] [-q] ID")
  quiet = False
  data_type = "pdb"
  if "-x" in args :
    data_type = "xray"
  elif "-f" in args :
    data_type = "fasta"
  if "-q" in args :
    quiet = True
  id = args[-1]
  try :
    data = fetch(id, data_type)
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
