
# return a list of all PDB entries with at least one alternate conformation
# in a protein chain.

from __future__ import division
from __future__ import print_function
from libtbx.easy_mp import pool_map
import os

def run():
  from mmtbx.wwpdb import rcsb_web_services
  assert ("PDB_MIRROR_PDB" in os.environ)
  f = open("pdb_with_alt_confs.txt", "w")
  f.write("# d_max <= 2.5, protein, xray, exp. data\n")
  protein_xray_structures_at_high_resolution = rcsb_web_services.post_query(
    query_xml=None,
    xray_only=True,
    d_max=2.51,
    protein_only=True,
    data_only=True) # XXX why trust anything else?
  all_results = pool_map(
    processes=8,
    chunksize=16,
    args=protein_xray_structures_at_high_resolution,
    func=get_nconfs)
  n_alts = 0
  for (pdb_id, n_confs) in all_results :
    if (n_confs > 1):
      n_alts += 1
      f.write(pdb_id + "\n")
  print("n_alts:", n_alts)
  f.close()

def get_nconfs(pdb_id):
  import iotbx.pdb.fetch
  assert ("PDB_MIRROR_PDB" in os.environ)
  n_confs = 0
  hierarchy, xray_structure = iotbx.pdb.fetch.load_pdb_structure(pdb_id,
    allow_unknowns=True)
  if (len(hierarchy.models()) > 1):
    n_confs = 0
  else :
    for chain in hierarchy.only_model().chains():
      if (chain.is_protein()):
        confs = chain.conformers()
        if (len(confs) > n_confs):
          n_confs = len(confs)
  return (pdb_id, n_confs)

if (__name__ == "__main__"):
  run()
