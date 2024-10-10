# TODO other PDB sites?
#
# See RCSB documentation at:
# https://www.rcsb.org/pages/download/http
#
# File format  Compression   Example URL
# PDB  uncompressed https://files.rcsb.org/download/4hhb.pdb
# CIF  uncompressed https://files.rcsb.org/download/4hhb.cif
# XML  uncompressed https://files.rcsb.org/download/4hhb.xml
# Data uncompressed https://files.rcsb.org/download/1btn-sf.cif
# CIF  uncompressed https://files.rcsb.org/ligands/download/HEM.cif
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
from collections import defaultdict
from libtbx.utils import Sorry, null_out
from libtbx import smart_open
from libtbx import Auto
import libtbx.utils
import libtbx.load_env
from six.moves.urllib.error import HTTPError
import gzip
import re
import os


def get_link(mirror, file_type, pdb_id=None, emdb_number=None):

  all_links_dict = {
      'rcsb': {
          'model_pdb': 'https://files.rcsb.org/pub/pdb/data/structures/divided/pdb/{mid_id}/pdb{pdb_id}.ent.gz',
          'model_cif': 'https://files.rcsb.org/pub/pdb/data/structures/divided/mmCIF/{mid_id}/{pdb_id}.cif.gz',
          'sequence': 'https://www.rcsb.org/fasta/entry/{pdb_id}',
          'sf': 'https://files.rcsb.org/download/{pdb_id}-sf.cif.gz',
          'map': 'https://files.rcsb.org/pub/emdb/structures/EMD-{emdb_number}/map/emd_{emdb_number}.map.gz',
          },
      'pdbe': {
          'model_pdb': 'https://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/pdb/{mid_id}/pdb{pdb_id}.ent.gz',
          'model_cif': 'https://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/mmCIF/{mid_id}/{pdb_id}.cif.gz',
          'sequence': 'https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta',
          'sf': 'https://www.ebi.ac.uk/pdbe/entry-files/download/r{pdb_id}sf.ent',
          'map': 'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{emdb_number}/map/emd_{emdb_number}.map.gz',
          },
      'pdbj': {
          'model_pdb': 'https://ftp.pdbj.org/pub/pdb/data/structures/divided/pdb/{mid_id}/pdb{pdb_id}.ent.gz',
          'model_cif': 'https://ftp.pdbj.org/pub/pdb/data/structures/divided/mmCIF/{mid_id}/{pdb_id}.cif.gz',
          'sequence': 'https://pdbj.org/rest/newweb/fetch/file?cat=pdb&type=fasta&id={pdb_id}',
          'sf': 'https://data.pdbjpw1.pdbj.org/pub/pdb/data/structures/divided/structure_factors/{mid_id}/r{pdb_id}sf.ent.gz',
          'map': 'https://ftp.pdbj.org/pub/emdb/structures/EMD-{emdb_number}/map/emd_{emdb_number}.map.gz',
          },
      # 'pdb-redo': {
      #     'model_pdb': 'https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.pdb',
      #     'model_cif': 'https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.cif',
      #     # these are from RCSB because PDB-redo does not have them
      #     'sequence': 'https://www.rcsb.org/fasta/entry/{pdb_id}',
      #     'sf': 'https://files.rcsb.org/download/{pdb_id}-sf.cif',
      #     'map': 'https://files.rcsb.org/pub/emdb/structures/EMD-{emdb_number}/map/emd_{emdb_number}.map.gz',
      #     },
  }


  assert mirror in ['rcsb', 'pdbe', 'pdbj']
  assert file_type in ['model_pdb', 'model_cif', 'sequence', 'sf', 'map']
  if file_type == 'map':
    assert emdb_number
  else:
    assert pdb_id
  mid_pdb_id = pdb_id[1:3]
  return all_links_dict[mirror][file_type].format(mid_id=mid_pdb_id, pdb_id=pdb_id, emdb_number=emdb_number)

def valid_pdb_id(id):
  return len(id) == 4 and re.match("[1-9]{1}[a-zA-Z0-9]{3}", id)

def get_file_type_from_params(data_type, format):
  """ Temporary, during refactoring. """
  conversion_dict = {
      "pdb":'model',
      "xray": 'sf',
      "fasta": 'sequence',
      "map": 'map',
      }
  result = conversion_dict[data_type]
  if result == 'model':
    if format == 'pdb':
      result = 'model_pdb'
    else:
      result = 'model_cif'
  return result

def fetch(id, data_type="pdb", format="pdb", mirror="rcsb", log=None):
  """
  Locate and open a data file for the specified PDB ID and format, either in a
  local mirror or online.

  :param id: 4-character PDB ID (e.g. '1hbb')
  :param data_type: type of content to download: pdb, xray, or fasta
  :param format: format of data: cif, pdb, or xml (or cif_or_pdb)
  :param mirror: remote site to use, either rcsb, pdbe, pdbj or pdb-redo

  :returns: a filehandle-like object (with read() method)
  """
  assert data_type in ["pdb", "xray", "fasta"]
  assert format in ["cif", "pdb", "xml", "cif_or_pdb"]
  assert mirror in ["rcsb", "pdbe", "pdbj", "pdb-redo"]
  id = id.lower()
  if not valid_pdb_id(id):
    raise Sorry("Invalid pdb id %s. Must be 4 characters, 1st is a number 1-9." % pdb_id)
  if (log is None) : log = null_out()

  file_type = get_file_type_from_params(data_type, format)
  url = get_link(mirror, file_type, pdb_id=id, emdb_number=None)
  need_to_decompress = url.split('.')[-1] == 'gz' and file_type != 'map'

  try :
    data = libtbx.utils.urlopen(url)
  except HTTPError as e :
    if e.getcode() == 404 :
      raise RuntimeError("Couldn't download %s for %s at %s." % (file_type, id, url))
    else :
      raise
  if need_to_decompress:
    return gzip.GzipFile(fileobj=data)
  return data

def write_data_to_disc(fname, data):
    with open(fname, "wb") as f:
      f.write(data.read())

def get_pdb(id, data_type, mirror, log, format="pdb"):
  """
  Frontend for fetch(...), writes resulting data to disk.
  """
  try :
    data = fetch(id, data_type, mirror=mirror, format=format, log=log)
  except RuntimeError as e :
    raise Sorry(str(e))
  default_value = (os.path.join(os.getcwd(), "{}.{}".format(id, format)), "Model")
  file_names_titles = defaultdict(lambda: default_value, {
      "xray":  (os.path.join(os.getcwd(), "{}-sf.cif".format(id)), "Structure factors"),
      "fasta": (os.path.join(os.getcwd(), "{}.fa".format(id)), "Sequence"),
  })
  file_name, title = file_names_titles[data_type]
  write_data_to_disc(file_name, data)
  print("%s saved to %s" % (title, file_name), file=log)
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
  chem_comp_cif = None
  if (chem_comp_cif is None):
    url = "https://files.rcsb.org/ligands/download/%s.cif" % code
    try :
      data = libtbx.utils.urlopen(url)
    except HTTPError as e :
      if e.getcode() == 404 :
        raise RuntimeError("Couldn't download sequence for %s." % id)
      else :
        raise
    else :
      file_name = "%s.cif" % code
      with open(file_name, "wb") as f:
        f.write(data.read())
      return file_name
  elif (not return_none_if_already_present):
    return chem_comp_cif
  return None

