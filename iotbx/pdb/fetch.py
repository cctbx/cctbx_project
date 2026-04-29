"""Fetch data from PDB"""
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
import libtbx.utils
import libtbx.load_env
from six.moves.urllib.error import HTTPError
import gzip
import re
import os


all_links_dict = {
    'rcsb': {
        'model_pdb': 'https://files.rcsb.org/pub/pdb/data/structures/divided/pdb/{mid_id}/pdb{pdb_id}.ent.gz',
        'model_cif': 'https://files.rcsb.org/pub/pdb/data/structures/divided/mmCIF/{mid_id}/{pdb_id}.cif.gz',
        'sequence': 'https://www.rcsb.org/fasta/entry/{pdb_id}',
        'sf': 'https://files.rcsb.org/download/{pdb_id}-sf.cif.gz',
        'em_map': 'https://files.rcsb.org/pub/emdb/structures/EMD-{emdb_number}/map/emd_{emdb_number}.map.gz',
        'em_half_map_1': 'https://files.rcsb.org/pub/emdb/structures/EMD-{emdb_number}/other/emd_{emdb_number}_half_map_1.map.gz',
        'em_half_map_2': 'https://files.rcsb.org/pub/emdb/structures/EMD-{emdb_number}/other/emd_{emdb_number}_half_map_2.map.gz',
        },
    'pdbe': {
        'model_pdb': 'https://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/pdb/{mid_id}/pdb{pdb_id}.ent.gz',
        'model_cif': 'https://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/mmCIF/{mid_id}/{pdb_id}.cif.gz',
        'sequence': 'https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta',
        'sf': 'https://www.ebi.ac.uk/pdbe/entry-files/download/r{pdb_id}sf.ent',
        'em_map': 'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{emdb_number}/map/emd_{emdb_number}.map.gz',
        'em_half_map_1': 'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{emdb_number}/other/emd_{emdb_number}_half_map_1.map.gz',
        'em_half_map_2': 'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{emdb_number}/other/emd_{emdb_number}_half_map_2.map.gz',
        },
    'pdbj': {
        'model_pdb': 'https://ftp.pdbj.org/pub/pdb/data/structures/divided/pdb/{mid_id}/pdb{pdb_id}.ent.gz',
        'model_cif': 'https://ftp.pdbj.org/pub/pdb/data/structures/divided/mmCIF/{mid_id}/{pdb_id}.cif.gz',
        'sequence': 'https://pdbj.org/rest/newweb/fetch/file?cat=pdb&type=fasta&id={pdb_id}',
        'sf': 'https://data.pdbjpw1.pdbj.org/pub/pdb/data/structures/divided/structure_factors/{mid_id}/r{pdb_id}sf.ent.gz',
        'em_map': 'https://ftp.pdbj.org/pub/emdb/structures/EMD-{emdb_number}/map/emd_{emdb_number}.map.gz',
        'em_half_map_1': 'https://ftp.pdbj.org/pub/databases/emdb/structures/EMD-{emdb_number}/other/emd_{emdb_number}_half_map_1.map.gz',
        'em_half_map_2': 'https://ftp.pdbj.org/pub/databases/emdb/structures/EMD-{emdb_number}/other/emd_{emdb_number}_half_map_2.map.gz',
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

def get_link(mirror, entity, pdb_id=None, emdb_number=None, link_templates=all_links_dict):
  assert mirror in link_templates.keys()
  if entity not in link_templates[mirror].keys():
    return None
  if entity.find('map') > 0:
    assert emdb_number
  else:
    assert pdb_id
  mid_pdb_id = pdb_id[1:3]
  return link_templates[mirror][entity].format(mid_id=mid_pdb_id, pdb_id=pdb_id, emdb_number=emdb_number)

def valid_pdb_id(id):
  return len(id) == 4 and re.match("[1-9]{1}[a-zA-Z0-9]{3}", id)

def fetch(id, entity='model_pdb', mirror="rcsb", emdb_number=None, link_templates=all_links_dict):
  """
  Locate and open a data file for the specified PDB ID and format, either in a
  local mirror or online.

  :param id: 4-character PDB ID (e.g. '1hbb')
  :param entity - one of 'model_pdb', 'model_cif', 'sequence', 'sf', 'em_map'
  :param mirror: remote site to use, either rcsb, pdbe, pdbj or pdb-redo

  :returns: a filehandle-like object (with read() method)
  """
  assert entity in ['model_pdb', 'model_cif', 'sequence', 'sf', 'em_map', 'em_half_map_1', 'em_half_map_2']
  assert mirror in ["rcsb", "pdbe", "pdbj"]
  id = id.lower()
  if not valid_pdb_id(id):
    raise Sorry("Invalid pdb id %s. Must be 4 characters, 1st is a number 1-9." % id)

  url = get_link(mirror, entity, pdb_id=id, emdb_number=emdb_number, link_templates=link_templates)
  need_to_decompress = url.split('.')[-1] == 'gz' and entity.find('map') < 0

  try :
    data = libtbx.utils.urlopen(url)
  except HTTPError as e :
    if e.getcode() == 404 or e.getcode() == 403 :
      raise RuntimeError("Couldn't download %s for %s at %s." % (entity, id, url))
    else :
      raise
  if need_to_decompress:
    return gzip.GzipFile(fileobj=data)
  return data

def write_data_to_disc(fname, data):
    with open(fname, "wb") as f:
      f.write(data.read())

def fetch_and_write(id, entity='model_pdb', mirror='rcsb', emdb_number=None, link_templates=all_links_dict, log=None):
  """
  Frontend for fetch(...), writes resulting data to disk.
  """
  try :
    data = fetch(id, entity, mirror=mirror, emdb_number=emdb_number, link_templates=link_templates)
  except RuntimeError as e :
    print(str(e),file=log)
    return None
  if (log is None) : log = null_out()
  default_value = (os.path.join(os.getcwd(), "{}.{}".format(id, format)), "Model")
  file_names_titles = defaultdict(lambda: default_value, {
      "model_pdb":  (os.path.join(os.getcwd(), "{}.pdb".format(id)), "Model in PDB format"),
      "model_cif":  (os.path.join(os.getcwd(), "{}.cif".format(id)), "Model in mmCIF format"),
      "sf":  (os.path.join(os.getcwd(), "{}-sf.cif".format(id)), "Structure factors"),
      "sequence": (os.path.join(os.getcwd(), "{}.fa".format(id)), "Sequence"),
      "em_map": (os.path.join(os.getcwd(), "emd_{}.map.gz".format(emdb_number)), "Cryo-EM map"),
      "em_half_map_1": (os.path.join(os.getcwd(), "emd_{}_half_map_1.map.gz".format(emdb_number)), "Cryo-EM half map 1"),
      "em_half_map_2": (os.path.join(os.getcwd(), "emd_{}_half_map_2.map.gz".format(emdb_number)), "Cryo-EM half map 2"),
  })
  file_name, title = file_names_titles[entity]
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

