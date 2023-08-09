from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb

# ------------------------------------------------------------------------------

# from https://github.com/wwPDB/extended-wwPDB-identifier-examples
# https://github.com/wwPDB/extended-wwPDB-identifier-examples/blob/main/Models/7fgz-extended_CCD_code-model.cif
mmcif_str= '''data_phenix
loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.pdbx_PDB_model_num
  ATOM     1  N    .  LYS  A  278  ?    0.39900  -10.17100  39.80200  1.000  40.11000  N   ?  A  ?  1  1
  ATOM     2  CA   .  LYS  A  278  ?   -0.16900   -9.98800  41.17300  1.000  43.86000  C   ?  A  ?  1  1
  ATOM     3  C    .  LYS  A  278  ?    0.68700   -9.01100  41.99100  1.000  41.94000  C   ?  A  ?  1  1
  ATOM     4  O    .  LYS  A  278  ?    1.04400   -7.92000  41.55600  1.000  39.32000  O   ?  A  ?  1  1
  ATOM     5  CB   .  LYS  A  278  ?   -0.26000  -11.33600  41.90200  1.000  46.47000  C   ?  A  ?  1  1
  ATOM     6  CG   .  LYS  A  278  ?   -1.58300  -12.07400  41.71300  1.000  49.13000  C   ?  A  ?  1  1
  ATOM     7  CD   .  LYS  A  278  ?   -1.61100  -13.46800  42.31500  1.000  51.03000  C   ?  A  ?  1  1
  ATOM     8  CE   .  LYS  A  278  ?   -2.92300  -13.79900  42.99300  1.000  52.86000  C   ?  A  ?  1  1
  ATOM     9  N    .  LYS  A  279  ?    0.39900  -10.17100  39.80200  1.000  40.11000  N   ?  B  ?  1  1
  ATOM    10  CA   .  LYS  A  279  ?   -0.16900   -9.98800  41.17300  1.000  43.86000  C   ?  B  ?  1  1
  ATOM    11  C    .  LYS  A  279  ?    0.68700   -9.01100  41.99100  1.000  41.94000  C   ?  B  ?  1  1
  ATOM    12  O    .  LYS  A  279  ?    1.04400   -7.92000  41.55600  1.000  39.32000  O   ?  B  ?  1  1
  ATOM    13  CB   .  LYS  A  279  ?   -0.26000  -11.33600  41.90200  1.000  46.47000  C   ?  B  ?  1  1
  ATOM    14  CG   .  LYS  A  279  ?   -1.58300  -12.07400  41.71300  1.000  49.13000  C   ?  B  ?  1  1
  ATOM    15  CD   .  LYS  A  279  ?   -1.61100  -13.46800  42.31500  1.000  51.03000  C   ?  B  ?  1  1
  ATOM    16  CE   .  LYS  A  279  ?   -2.92300  -13.79900  42.99300  1.000  52.86000  C   ?  B  ?  1  1
  ATOM    17  NZ   .  LYS  A  279  ?   -3.20900  -12.85600  44.10000  1.000  54.19000  N   ?  B  ?  1  1
  HETATM  18  CA   .  CA   A  301  ?  -17.36200  -22.38500  28.04700  1.000  15.20000  CA  ?  C  ?  .  1
  HETATM  19  C10  .  7ZT  A  302  ?   -7.64600   -6.96500   5.79600  1.000  22.62000  C   ?  D  ?  .  1
  HETATM  20  C2   .  7ZT  A  302  ?   -8.46200   -5.53400   9.26500  1.000  16.68000  C   ?  D  ?  .  1
  HETATM  21  C3   .  7ZT  A  302  ?   -8.09200   -5.86500  10.71100  1.000  17.35000  C   ?  D  ?  .  1
  HETATM  22  C5   .  7ZT  A  302  ?   -6.78100   -7.68900  10.02700  1.000  17.11000  C   ?  D  ?  .  1
  HETATM  23  C6   .  7ZT  A  302  ?   -7.36300   -7.73000   8.63300  1.000  16.97000  C   ?  D  ?  .  1
  HETATM  24  C7   .  7ZT  A  302  ?   -8.07100   -7.95100  12.06400  1.000  16.77000  C   ?  D  ?  .  1
  HETATM  25  C8   .  7ZT  A  302  ?   -8.49000   -9.44000  12.02500  1.000  16.60000  C   ?  D  ?  .  1
  HETATM  26  C9   .  7ZT  A  302  ?   -8.43800   -6.31100   6.93500  1.000  20.13000  C   ?  D  ?  .  1
  HETATM  27  N1   .  7ZT  A  302  ?   -7.74300   -6.35500   8.24300  1.000  18.72000  N   ?  D  ?  .  1
  HETATM  28  N4   .  7ZT  A  302  ?   -7.97500   -7.33400  10.76700  1.000  17.04000  N   ?  D  ?  .  1
  HETATM  29  O1S  .  7ZT  A  302  ?   -7.57900   -7.37100   3.22400  1.000  25.47000  O   ?  D  ?  .  1
  HETATM  30  O2S  .  7ZT  A  302  ?   -8.53300   -5.21100   4.07800  1.000  25.01000  O   ?  D  ?  .  1
  HETATM  31  O3S  .  7ZT  A  302  ?   -9.69100   -7.32600   4.28200  1.000  26.42000  O   ?  D  ?  .  1
  HETATM  32  O8   .  7ZT  A  302  ?   -8.39100  -10.14000  10.78100  1.000  14.14000  O   ?  D  ?  .  1
  HETATM  33  S    .  7ZT  A  302  ?   -8.37200   -6.70800   4.27900  1.000  24.90000  S   ?  D  ?  .  1
  HETATM  34  C10  .  7ZT  A  303  ?  -10.12100   -0.28900  11.39600  1.000  75.92000  C   ?  E  ?  .  1
  HETATM  35  C2   .  7ZT  A  303  ?  -11.68200    2.49500   9.26400  1.000  77.03000  C   ?  E  ?  .  1
  HETATM  36  C3   .  7ZT  A  303  ?  -11.90700    3.97900   9.61700  1.000  76.82000  C   ?  E  ?  .  1
  HETATM  37  C5   .  7ZT  A  303  ?  -12.65400    3.12300  11.90800  1.000  69.90000  C   ?  E  ?  .  1
  HETATM  38  C6   .  7ZT  A  303  ?  -12.33300    1.71900  11.40900  1.000  69.07000  C   ?  E  ?  .  1
  HETATM  39  C7   .  7ZT  A  303  ?  -11.41500    5.23200  11.72500  1.000  71.46000  C   ?  E  ?  .  1
  HETATM  40  C8   .  7ZT  A  303  ?  -12.12600    5.97600  12.87100  1.000  70.74000  C   ?  E  ?  .  1
  HETATM  41  C9   .  7ZT  A  303  ?  -10.75600    0.35700  10.14200  1.000  74.40000  C   ?  E  ?  .  1
  HETATM  42  N1   .  7ZT  A  303  ?  -11.23300    1.73200  10.44600  1.000  74.61000  N   ?  E  ?  .  1
  HETATM  43  N4   .  7ZT  A  303  ?  -12.33500    4.28200  11.02200  1.000  73.77000  N   ?  E  ?  .  1
  HETATM  44  O1S  .  7ZT  A  303  ?   -7.73500   -0.30500  10.52500  1.000  74.01000  O   ?  E  ?  .  1
  HETATM  45  O2S  .  7ZT  A  303  ?   -9.04200   -2.31300  10.14100  1.000  83.78000  O   ?  E  ?  .  1
  HETATM  46  O3S  .  7ZT  A  303  ?   -8.26400   -1.86400  12.32100  1.000  70.44000  O   ?  E  ?  .  1
  HETATM  47  O8   .  7ZT  A  303  ?  -11.36000    5.92100  14.09600  1.000  61.89000  O   ?  E  ?  .  1
  HETATM  48  S    .  7ZT  A  303  ?   -8.75700   -1.21400  11.09100  1.000  77.53000  S   ?  E  ?  .  1
  HETATM  49  O    .  HOH  A  401  ?   10.56700   -6.53900  28.06400  1.000  30.11000  O   ?  F  ?  .  1
  HETATM  50  O    .  HOH  A  402  ?    2.85600   -8.84800  40.69200  1.000  32.14000  O   ?  F  ?  .  1
  HETATM  51  O    .  HOH  A  403  ?  -10.57700  -29.76200  22.83300  1.000  22.07000  O   ?  F  ?  .  1
  HETATM  52  O    .  HOH  A  404  ?   14.70500   -5.43700  13.98800  1.000  28.06000  O   ?  F  ?  .  1

loop_
  _chem_comp.id
  7ZT
  CA
  HOH
  LYS

loop_
  _struct_asym.id
  A
  B
  C
  D
  E
  F

'''

def test1():
  """
  Test correct reading of mmCIF file with multiple label_asym_id (multiple chains with
    same chain_id) into hierarchy and writing back and forth to mmCIF.
  """
  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  h = inp.construct_hierarchy()
  text = h.as_mmcif_string()
  for a,b in zip(mmcif_str.splitlines(),text.splitlines()):
    if a != b:
      print("from mmCIF text: %s" %(a))
      print("From hierarchy : %s" %(b))
      assert a==b
  assert text == mmcif_str
  chains = list(h.overall_counts().chain_ids.keys())
  chains.sort()
  answer = ['A']
  assert (chains == answer), '%s %s' % (chains, answer)
  resnames = sorted(h.overall_counts().resnames.keys())
  assert resnames == [' CA', '7ZT', 'HOH', 'LYS'], resnames

  # write as mmCIF string
  mmcif_text = h.as_mmcif_string()

  new_h = iotbx.pdb.input(lines=mmcif_text.split("\n"), source_info=None).construct_hierarchy()
  # as PDB:
  for a,b in zip(h.as_pdb_string().splitlines(),new_h.as_pdb_string().splitlines()):
    if a != b:
      print("Original as PDB: %s" %(a))
      print("From CIF as PDB: %s" %(b))
      assert a==b
  # As mmCIF
  for a,b in zip(h.as_mmcif_string().splitlines(),new_h.as_mmcif_string().splitlines()):
    if a != b:
      print("Original as CIF: %s" %(a))
      print("From CIF as CIF: %s" %(b))
      assert a==b

  assert h.is_similar_hierarchy(new_h)
  for a,b in zip(new_h.as_mmcif_string().splitlines(),mmcif_str.splitlines()):
    if a != b:
      print("From hierarchy: %s" %(a))
      print("Original CIF  : %s" %(b))
      assert a == b
  assert new_h.as_mmcif_string() == mmcif_str


if (__name__ == "__main__"):
  t0 = time.time()
  test1()
  print("OK. Time: %8.3f"%(time.time()-t0))
