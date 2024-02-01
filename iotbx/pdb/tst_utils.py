from __future__ import absolute_import, division, print_function
import iotbx.pdb.utils

pdb_str_to_be_cif="""
CRYST1   40.339   36.116   46.266  90.00  90.00  90.00 P 1
ATOM      1  CA  ASP AXYB2      34.633  18.762  20.254  1.00 22.59           C
ATOM      2  CA  LYS AXYB3      36.047  17.704  23.610  1.00 19.79           C
ATOM      3  CA  ILE AXYB4      35.551  19.482  26.886  1.00 19.33           C
ATOM      4  CA AHIS AXYB5      38.649  21.223  28.218  0.50 19.79           C
ATOM      5  CA BHIS AXYB6      38.583  21.270  28.209  0.50 20.43           C
ATOM      6  CA  GLY A 138      38.261  15.285  27.690  1.00  6.80           C
ATOM      7  CA  ALA A 139      34.607  14.241  27.428  1.00  4.76           C
ATOM      8  CA ALEU A 140      33.091  14.490  23.937  0.50  5.08           C
ATOM      9  CA BLEU A 140      33.072  14.565  23.972  0.50  5.41           C
ATOM     10  CA  ASN A 141      30.271  17.061  23.474  1.00  5.65           C
"""

as_pdb = pdb_str_to_be_cif
# Convert pdb_str_hybrid_residues to mmcif:
from libtbx.test_utils import convert_pdb_to_cif_for_pdb_str
convert_pdb_to_cif_for_pdb_str(locals(), chain_addition="ZXLONG")
as_cif = pdb_str_to_be_cif # now it is cif

def exercise_all_chain_ids():
  ids = iotbx.pdb.utils.all_chain_ids()
  assert len(ids)==3906
  assert len(set(ids))==3906

def exercise_get_pdb_info():
  from iotbx.pdb.utils import get_pdb_info
  pdb_info_from_pdb = get_pdb_info(as_pdb)
  pdb_info_from_cif = get_pdb_info(as_cif)
  assert not pdb_info_from_pdb.hierarchy.is_similar_hierarchy(
    pdb_info_from_cif.hierarchy)
  for model in pdb_info_from_cif.hierarchy.models():
    for chain in model.chains():
      chain.id = chain.id.replace("ZXLONG","") # make it short again
  assert pdb_info_from_pdb.hierarchy.is_similar_hierarchy(
    pdb_info_from_cif.hierarchy)
  assert pdb_info_from_pdb.crystal_symmetry.is_similar_symmetry(
    pdb_info_from_cif.crystal_symmetry)

def run():
  exercise_get_pdb_info()
  exercise_all_chain_ids()
  print("OK")

if __name__ == '__main__':
  run()
