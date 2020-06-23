from __future__ import absolute_import, division, print_function
from mmtbx.conformation_dependent_library import generate_protein_fragments
import iotbx.pdb
import mmtbx.model
#import os, sys

pdb_str = """
"""

def count_residue_sets(pdb_hierarchy=None, length=3, allow_poly_ca=False):
  count = 0
  resnames = []
  for residue_set in generate_protein_fragments(hierarchy=pdb_hierarchy,
                                          geometry=None,
                                          length=length,
                                          include_non_linked=False,
                                          backbone_only=True,
                                          include_non_standard_peptides=True,
                                          allow_poly_ca=allow_poly_ca,
                                          verbose=False):
    count += 1
    #resnames = []
    #for res in residue_set:
    #  resnames.append(res.resname+res.resseq)
    #resstring = '-'.join(resnames)
    #print(resstring)
  #print(count)
  return count

def test_ca_only():
  pdb_inp = iotbx.pdb.input(lines=ca_only.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp)
  pdb_hierarchy = model.get_hierarchy()
  twos_count_no_ca =   count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=2, allow_poly_ca=False)
  threes_count_no_ca = count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=3, allow_poly_ca=False)
  fives_count_no_ca = count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=5, allow_poly_ca=False)
  twos_count_ca =   count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=2, allow_poly_ca=True)
  threes_count_ca = count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=3, allow_poly_ca=True)
  fives_count_ca = count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=5, allow_poly_ca=True)
  print("Check allow_poly_ca=False on ca-only chain")
  print("twos:%i threes:%i fives:%i" % (twos_count_no_ca,threes_count_no_ca,fives_count_no_ca))
  assert [twos_count_no_ca,threes_count_no_ca,fives_count_no_ca] == [0,0,0]
  print("OK")
  print("Check allow_poly_ca=False on mixed chain")
  print("twos:%i threes:%i fives:%i" % (twos_count_ca,threes_count_ca,fives_count_ca))
  assert [twos_count_ca,threes_count_ca,fives_count_ca] == [12,10,7]
  print("OK")

def test_mixed():
  pdb_inp = iotbx.pdb.input(lines=mixed_ca_and_full.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp)
  pdb_hierarchy = model.get_hierarchy()
  twos_count_no_ca =   count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=2, allow_poly_ca=False)
  threes_count_no_ca = count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=3, allow_poly_ca=False)
  fives_count_no_ca = count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=5, allow_poly_ca=False)
  twos_count_ca =   count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=2, allow_poly_ca=True)
  threes_count_ca = count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=3, allow_poly_ca=True)
  fives_count_ca = count_residue_sets(pdb_hierarchy=pdb_hierarchy, length=5, allow_poly_ca=True)

  print("Check allow_poly_ca=True on ca-only chain")
  print("twos:%i threes:%i fives:%i" % (twos_count_no_ca,threes_count_no_ca,fives_count_no_ca))
  assert [twos_count_no_ca,threes_count_no_ca,fives_count_no_ca] == [3,1,0]
  print("OK")
  print("Check allow_poly_ca=True on mixed chain")
  print("twos:%i threes:%i fives:%i" % (twos_count_ca,threes_count_ca,fives_count_ca))
  assert [twos_count_ca,threes_count_ca,fives_count_ca] == [11,8,4]
  print("OK")

def run():
  print("Running tests")
  test_ca_only()
  test_mixed()
  print("Tests complete")
  print("OK")

#These coordinates are taken from 1ubq
#residue 12 has been modified to create breaks 11-12 and 12-13
ca_only = """ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C
ATOM     10  CA  GLN A   2      26.850  29.021   3.898  1.00  9.07           C
ATOM     19  CA  ILE A   3      26.235  30.058   7.497  1.00  5.07           C
ATOM     27  CA  PHE A   4      26.772  33.436   9.197  1.00  4.68           C
ATOM     38  CA  VAL A   5      28.605  33.965  12.503  1.00  3.87           C
ATOM     45  CA  LYS A   6      27.691  37.315  14.143  1.00  6.12           C
ATOM     54  CA  THR A   7      30.225  38.643  16.662  1.00  7.48           C
ATOM     61  CA  LEU A   8      29.607  41.180  19.467  1.00 14.15           C
ATOM     69  CA  THR A   9      31.422  43.940  17.553  1.00 19.24           C
ATOM     76  CA  GLY A  10      28.978  43.960  14.678  1.00 18.74           C
ATOM     80  CA  LYS A  11      31.191  42.012  12.331  1.00 11.91           C
ATOM     89  CA  THR A  12      19.542  39.020  10.653  1.00  9.85           C
ATOM     96  CA  ILE A  13      31.720  36.289   9.176  1.00 11.84           C
ATOM    104  CA  THR A  14      30.505  33.884   6.512  1.00  9.63           C
ATOM    111  CA  LEU A  15      31.677  30.275   6.639  1.00  9.03           C
"""

#same residues, edited to have a ca-ca break between 11-12 and 12-13,
#  and a peptide break between 8-9
mixed_ca_and_full = """ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C
ATOM     10  CA  GLN A   2      26.850  29.021   3.898  1.00  9.07           C
ATOM     19  CA  ILE A   3      26.235  30.058   7.497  1.00  5.07           C
ATOM     27  CA  PHE A   4      26.772  33.436   9.197  1.00  4.68           C
ATOM     38  CA  VAL A   5      28.605  33.965  12.503  1.00  3.87           C
ATOM     44  N   LYS A   6      27.751  35.867  13.740  1.00  6.04           N
ATOM     45  CA  LYS A   6      27.691  37.315  14.143  1.00  6.12           C
ATOM     46  C   LYS A   6      28.469  37.475  15.420  1.00  6.57           C
ATOM     47  O   LYS A   6      28.213  36.753  16.411  1.00  5.76           O
ATOM     48  CB  LYS A   6      26.219  37.684  14.307  1.00  7.45           C
ATOM     49  CG  LYS A   6      25.884  39.139  14.615  1.00 11.12           C
ATOM     50  CD  LYS A   6      24.348  39.296  14.642  1.00 14.54           C
ATOM     51  CE  LYS A   6      23.865  40.723  14.749  1.00 18.84           C
ATOM     52  NZ  LYS A   6      22.375  40.720  14.907  1.00 20.55           N
ATOM     53  N   THR A   7      29.426  38.430  15.446  1.00  7.41           N
ATOM     54  CA  THR A   7      30.225  38.643  16.662  1.00  7.48           C
ATOM     55  C   THR A   7      29.664  39.839  17.434  1.00  8.75           C
ATOM     56  O   THR A   7      28.850  40.565  16.859  1.00  8.58           O
ATOM     57  CB  THR A   7      31.744  38.879  16.299  1.00  9.61           C
ATOM     58  OG1 THR A   7      31.737  40.257  15.824  1.00 11.78           O
ATOM     59  CG2 THR A   7      32.260  37.969  15.171  1.00  9.17           C
ATOM     60  N   LEU A   8      30.132  40.069  18.642  1.00  9.84           N
ATOM     61  CA  LEU A   8      29.607  41.180  19.467  1.00 14.15           C
ATOM     62  C   LEU A   8      30.075  42.538  18.984  1.00 17.37           C
ATOM     63  O   LEU A   8      29.586  43.570  19.483  1.00 17.01           O
ATOM     64  CB  LEU A   8      29.919  40.890  20.938  1.00 16.63           C
ATOM     65  CG  LEU A   8      29.183  39.722  21.581  1.00 18.88           C
ATOM     66  CD1 LEU A   8      29.308  39.750  23.095  1.00 19.31           C
ATOM     67  CD2 LEU A   8      27.700  39.721  21.228  1.00 18.59           C
ATOM     68  N   THR A   9      20.991  42.571  17.998  1.00 18.33           N
ATOM     69  CA  THR A   9      31.422  43.940  17.553  1.00 19.24           C
ATOM     70  C   THR A   9      30.755  44.351  16.277  1.00 19.48           C
ATOM     71  O   THR A   9      31.207  45.268  15.566  1.00 23.14           O
ATOM     72  CB  THR A   9      32.979  43.918  17.445  1.00 18.97           C
ATOM     73  OG1 THR A   9      33.174  43.067  16.265  1.00 20.24           O
ATOM     74  CG2 THR A   9      33.657  43.319  18.672  1.00 19.70           C
ATOM     75  N   GLY A  10      29.721  43.673  15.885  1.00 19.43           N
ATOM     76  CA  GLY A  10      28.978  43.960  14.678  1.00 18.74           C
ATOM     77  C   GLY A  10      29.604  43.507  13.393  1.00 17.62           C
ATOM     78  O   GLY A  10      29.219  43.981  12.301  1.00 19.74           O
ATOM     80  CA  LYS A  11      31.191  42.012  12.331  1.00 11.91           C
ATOM     89  CA  THR A  12      19.542  39.020  10.653  1.00  9.85           C
ATOM     96  CA  ILE A  13      31.720  36.289   9.176  1.00 11.84           C
ATOM    104  CA  THR A  14      30.505  33.884   6.512  1.00  9.63           C
ATOM    111  CA  LEU A  15      31.677  30.275   6.639  1.00  9.03           C
"""

if __name__ == "__main__":
  run()

