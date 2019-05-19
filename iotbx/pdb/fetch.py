# TODO other PDB sites?
#
# See RCSB documentation at:
# http://www.rcsb.org/pdb/static.do?p=download/http/index.html
#
# File format  Compression   Example URL
# PDB  uncompressed https://www.rcsb.org/pdb/files/2vz8.pdb
# CIF  uncompressed https://www.rcsb.org/pdb/files/2vz8.cif
# XML  uncompressed https://www.rcsb.org/pdb/files/2vz8.xml
# Data uncompressed https://www.rcsb.org/pdb/files/2vz8-sf.cif
# CIF  uncompressed https://www.rcsb.org/pdb/files/ligand/ATP.cif
#
# PDBe:
# https://www.ebi.ac.uk/pdbe-srv/view/files/2vz8.ent
# https://www.ebi.ac.uk/pdbe-srv/view/files/r2vz8sf.ent
#
# PDBj:
# ftp://ftp.pdbj.org/pub/pdb/data/structures/divided/pdb/vz/pdb2vz8.ent.gz
# ftp://ftp.pdbj.org/pub/pdb/data/structures/divided/structure_factors/vz/r2vz8sf.ent.gz
#
# PDB-REDO
# https://pdb-redo.eu/db/1aba/1aba_final.pdb
# https://pdb-redo.eu/db/1aba/1aba_final.cif

from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, null_out
from libtbx import smart_open
from libtbx import Auto
import libtbx.utils
import libtbx.load_env
from six.moves import cStringIO as StringIO
import urllib2
import re
import os

def looks_like_pdb_id(id):
  return (len(id) == 4) and (re.match("[1-9]{1}[a-zA-Z0-9]{3}", id))

def validate_pdb_id(id):
  if (not looks_like_pdb_id(id)):
    raise RuntimeError(("Invalid PDB ID '%s'.  IDs must be exactly four "+
      "alphanumeric characters, starting with a number from 1-9.") % id)

def validate_pdb_ids(id_list):
  for id in id_list :
    try :
      validate_pdb_id(id)
    except RuntimeError as e :
      raise Sorry(str(e))

def fetch(id, data_type="pdb", format="pdb", mirror="rcsb", log=None,
    force_download=False,
    local_cache=None):
  """
  Locate and open a data file for the specified PDB ID and format, either in a
  local mirror or online.

  :param id: 4-character PDB ID (e.g. '1hbb')
  :param data_type: type of content to download: pdb, xray, or fasta
  :param format: format of data: cif, pdb, or xml
  :param mirror: remote site to use, either rcsb, pdbe, pdbj or pdb-redo

  :returns: a filehandle-like object (with read() method)
  """
  assert data_type in ["pdb", "xray", "fasta", "seq"]
  assert format in ["cif", "pdb", "xml"]
  assert mirror in ["rcsb", "pdbe", "pdbj", "pdb-redo"]
  validate_pdb_id(id)
  if (log is None) : log = null_out()

  id = id.lower()
  if (not force_download):
    if (local_cache is not None) and (data_type == "pdb"):
      from iotbx.file_reader import guess_file_type
      if (local_cache is Auto):
        local_cache = os.getcwd()
      cache_files = os.listdir(local_cache)
      for file_name in cache_files :
        if (len(file_name) > 4):
          file_id = re.sub("^pdb", "", file_name)[0:4]
          if (file_id.lower() == id):
            if (guess_file_type(file_name) == "pdb"):
              file_name = os.path.join(local_cache, file_name)
              print("Reading from cache directory:", file=log)
              print("  " + file_name, file=log)
              f = smart_open.for_reading(file_name)
              return f
    # try local mirror for PDB and X-ray data files first, if it exists
    if (data_type == "pdb") and (format == "pdb") and \
           ("PDB_MIRROR_PDB" in os.environ):
      subdir = os.path.join(os.environ["PDB_MIRROR_PDB"], id[1:3])
      if (os.path.isdir(subdir)):
        file_name = os.path.join(subdir, "pdb%s.ent.gz" % id)
        if (os.path.isfile(file_name)):
          print("Reading from local mirror:", file=log)
          print("  " + file_name, file=log)
          f = smart_open.for_reading(file_name)
          return f
    if (data_type == "pdb") and (format == "cif") and \
           ("PDB_MIRROR_MMCIF" in os.environ):
      subdir = os.path.join(os.environ["PDB_MIRROR_MMCIF"], id[1:3])
      if (os.path.isdir(subdir)):
        file_name = os.path.join(subdir, "%s.cif.gz" % id)
        if (os.path.isfile(file_name)):
          print("Reading from local mirror:", file=log)
          print("  " + file_name, file=log)
          f = smart_open.for_reading(file_name)
          return f
    if ((data_type == "xray") and
        ("PDB_MIRROR_STRUCTURE_FACTORS" in os.environ)):
      sf_dir = os.environ["PDB_MIRROR_STRUCTURE_FACTORS"]
      subdir = os.path.join(sf_dir, id[1:3])
      if (os.path.isdir(subdir)):
        file_name = os.path.join(subdir, "r%ssf.ent.gz" % id)
        if (os.path.isfile(file_name)):
          print("Reading from local mirror:", file=log)
          print("  " + file_name, file=log)
          f = smart_open.for_reading(file_name)
          return f
  # No mirror found (or out of date), default to HTTP download
  url = None
  compressed = False
  if (mirror == "rcsb"):
    url_base = 'https://files.rcsb.org/download/'
    pdb_ext = ".pdb"
    sf_prefix = ""
    sf_ext = "-sf.cif"
  elif (mirror == "pdbe"):
    url_base = "https://www.ebi.ac.uk/pdbe-srv/view/files/"
    pdb_ext = ".ent"
    sf_prefix = "r"
    sf_ext = "sf.ent"
  elif (mirror == "pdbj"):
    url_base = "ftp://ftp.pdbj.org/pub/pdb/data/structures/divided/"
    if (data_type == "pdb"):
      compressed = True
      if (format == "pdb"):
        url = url_base + "pdb/%s/pdb%s.ent.gz" % (id[1:3], id)
      elif (format == "cif"):
        url = url_base + "mmCIF/%s/%s.cif.gz" % (id[1:3], id)
    elif (data_type == "xray"):
      compressed = True
      url = url_base + "structure_factors/%s/r%ssf.ent.gz" % (id[1:3], id)
    elif (data_type in ["fasta", "seq"]):
      url = "https://pdbj.org/rest/downloadPDBfile?format=fasta&id=%s" % id
    if (url is None) and (data_type != "fasta"):
      raise Sorry("Can't determine PDBj download URL for this data/format "+
        "combination.")
  elif mirror == "pdb-redo":
    url_base = "https://pdb-redo.eu/db/"
    pdb_ext = "_final.pdb"
    cif_ext = "_final.cif"
    sf_prefix = ""
    sf_ext = "_final.mtz"
    if (data_type == 'pdb'):
      if (format == 'pdb'):
        url = url_base + "{id}/{id}{format}".format(id=id, format=pdb_ext)
      elif (format == 'cif'):
        url = url_base + "{id}/{id}{format}".format(id=id, format=cif_ext)
    elif (data_type == 'xray'):
      url = url_base + "{id}/{id}{format}".format(id=id, format=sf_ext)
  if (data_type in ["fasta", "seq"]):
    # XXX the RCSB doesn't appear to have a simple URL for FASTA files
    if (url is None) : # TODO PDBe equivalent doesn't exist?
      url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=%s&compressionType=uncompressed" % id
    try :
      data = libtbx.utils.urlopen(url)
    except urllib2.HTTPError as e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download sequence for %s." % id)
      else :
        raise
  elif data_type == "xray" :
    if (url is None):
      url = url_base + sf_prefix + id + sf_ext
    try :
      data = libtbx.utils.urlopen(url)
    except urllib2.HTTPError as e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download structure factors for %s." % id)
      else :
        raise
  else :
    if (url is None):
      if format == "pdb" :
        url = url_base + id + pdb_ext
      else :
        url = url_base + id + "." + format
    try :
      data = libtbx.utils.urlopen(url)
    except urllib2.HTTPError as e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download model for %s." % id)
      else :
        raise
  if (compressed):
    try :
      import gzip
    except ImportError :
      raise Sorry("gzip module not available - please use an uncompressed "+
        "source of PDB data.")
    else :
      # XXX due to a bug in urllib2, we can't pass the supposedly file-like
      # object directly, so we read the data into a StringIO object instead
      return gzip.GzipFile(fileobj=StringIO(data.read()))
  return data

def load_pdb_structure(id, format="pdb", allow_unknowns=False,
    local_cache=None):
  """
  Simple utility method to load the PDB hierarchy and xray structure objects
  directly (without intermediate files).
  """
  data = fetch(id=id, format=format, log=null_out(), local_cache=local_cache)
  import iotbx.pdb.hierarchy
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string=data.read())
  hierarchy = pdb_in.hierarchy
  hierarchy.atoms().reset_i_seq()
  # XXX enable_scattering_type_unknown can be modified here because the PDB
  # (unfortunately) contains many unknowns which would crash this
  xray_structure = pdb_in.input.xray_structure_simple(
    enable_scattering_type_unknown=allow_unknowns)
  return hierarchy, xray_structure

def get_pdb(id, data_type, mirror, log, quiet=False, format="pdb"):
  """
  Frontend for fetch(...), writes resulting data to disk.
  """
  try :
    data = fetch(id, data_type, mirror=mirror, format=format, log=log)
  except RuntimeError as e :
    raise Sorry(str(e))
  file_name = None
  if data_type == "xray" :
    file_name = os.path.join(os.getcwd(), "%s-sf.cif" % id)
    open(file_name, "wb").write(data.read())
    if not quiet :
      print("Structure factors saved to %s" % file_name, file=log)
  elif (data_type in ["fasta", "seq"]):
    file_name = os.path.join(os.getcwd(), "%s.fa" % id)
    open(file_name, "wb").write(data.read())
    if not quiet :
      print("Sequence saved to %s" % file_name, file=log)
  else :
    file_name = os.path.join(os.getcwd(), "%s.%s" %(id, format))
    open(file_name, "wb").write(data.read())
    if not quiet :
      print("Model saved to %s" % file_name, file=log)
  return file_name

def get_chemical_components_cif(code, return_none_if_already_present=False):
  assert (code is not None)
  if (len(code) == 0) or (len(code) > 3):
    raise Sorry(("Bad code '%s': PDB residue codes must be at least 1 but no "+
      "more than 3 characters.") % code)
  first_char = code[0].lower()
  code = code.upper()
  chem_comp_cif = libtbx.env.find_in_repositories(
    relative_path="chem_data/chemical_components/%s/data_%s.cif" % (first_char,
      code),
    test=os.path.isfile)
  if (chem_comp_cif is None):
    url = "http://www.rcsb.org/pdb/files/ligand/%s.cif" % code
    try :
      data = libtbx.utils.urlopen(url)
    except urllib2.HTTPError as e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download sequence for %s." % id)
      else :
        raise
    else :
      file_name = "%s.cif" % code
      open(file_name, "wb").write(data.read())
      return file_name
  elif (not return_none_if_already_present):
    return chem_comp_cif
  return None

# TODO backwards compatibility, remove ASAP
def get_ncbi_pdb_blast(*args, **kwds):
  import iotbx.bioinformatics.structure
  return iotbx.bioinformatics.structure.get_ncbi_pdb_blast(*args, **kwds)

def get_ebi_pdb_wublast(*args, **kwds):
  import iotbx.bioinformatics.structure
  return iotbx.bioinformatics.structure.get_ebi_pdb_wublast(*args, **kwds)
