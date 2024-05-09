from __future__ import absolute_import, division, print_function
import iotbx.pdb
from cctbx.array_family import flex
from libtbx import easy_run

pdb_str = """\
CRYST1   27.475   35.191   27.490  90.00  90.00  90.00 P 21 21 21
HELIX    1   1 ALA E    1  ALA E   16  1                                  16
ATOM      1  N   ALA E   1       6.709  31.817   5.664  1.00 20.00           N
ATOM      2  CA  ALA E   1       7.794  30.914   6.029  1.00 20.00           C
ATOM      3  C   ALA E   1       7.267  29.691   6.773  1.00 20.00           C
ATOM      4  O   ALA E   1       7.942  29.146   7.647  1.00 20.00           O
ATOM      5  CB  ALA E   1       8.831  31.643   6.870  1.00 20.00           C
ATOM      6  N   HIS E   2       6.057  29.267   6.421  1.00 20.00           N
ATOM      7  CA  HIS E   2       5.434  28.110   7.053  1.00 20.00           C
ATOM      8  C   HIS E   2       6.176  26.824   6.701  1.00 20.00           C
ATOM      9  O   HIS E   2       6.477  26.012   7.576  1.00 20.00           O
ATOM     10  CB  HIS E   2       3.965  28.000   6.640  1.00 20.00           C
ATOM     11  CG  HIS E   2       3.140  29.191   7.016  1.00 20.00           C
ATOM     12  ND1 HIS E   2       2.481  29.289   8.222  1.00 20.00           N
ATOM     13  CD2 HIS E   2       2.869  30.335   6.345  1.00 20.00           C
ATOM     14  CE1 HIS E   2       1.838  30.441   8.278  1.00 20.00           C
ATOM     15  NE2 HIS E   2       2.057  31.096   7.152  1.00 20.00           N
ATOM     16  N   CYS E   3       6.468  26.648   5.417  1.00 20.00           N
ATOM     17  CA  CYS E   3       7.193  25.473   4.949  1.00 20.00           C
ATOM     18  C   CYS E   3       8.635  25.492   5.443  1.00 20.00           C
ATOM     19  O   CYS E   3       9.223  24.444   5.712  1.00 20.00           O
ATOM     20  CB  CYS E   3       7.161  25.396   3.421  1.00 20.00           C
ATOM     21  SG  CYS E   3       5.501  25.285   2.714  1.00 20.00           S
ATOM     22  N   ALA E   4       9.197  26.690   5.560  1.00 20.00           N
ATOM     23  CA  ALA E   4      10.560  26.855   6.050  1.00 20.00           C
ATOM     24  C   ALA E   4      10.655  26.475   7.523  1.00 20.00           C
ATOM     25  O   ALA E   4      11.599  25.805   7.942  1.00 20.00           O
ATOM     26  CB  ALA E   4      11.031  28.285   5.837  1.00 20.00           C
ATOM     27  N   ILE E   5       9.669  26.907   8.304  1.00 20.00           N
ATOM     28  CA  ILE E   5       9.621  26.584   9.725  1.00 20.00           C
ATOM     29  C   ILE E   5       9.325  25.103   9.932  1.00 20.00           C
ATOM     30  O   ILE E   5       9.836  24.482  10.865  1.00 20.00           O
ATOM     31  CB  ILE E   5       8.580  27.459  10.467  1.00 20.00           C
ATOM     32  CG1 ILE E   5       9.071  28.905  10.559  1.00 20.00           C
ATOM     33  CG2 ILE E   5       8.307  26.928  11.866  1.00 20.00           C
ATOM     34  CD1 ILE E   5       8.105  29.836  11.259  1.00 20.00           C
ATOM     35  N   TYR E   6       8.498  24.543   9.054  1.00 20.00           N
ATOM     36  CA  TYR E   6       8.157  23.127   9.116  1.00 20.00           C
ATOM     37  C   TYR E   6       9.375  22.264   8.811  1.00 20.00           C
ATOM     38  O   TYR E   6       9.598  21.239   9.454  1.00 20.00           O
ATOM     39  CB  TYR E   6       7.017  22.795   8.150  1.00 20.00           C
ATOM     40  CG  TYR E   6       5.662  23.291   8.603  1.00 20.00           C
ATOM     41  CD1 TYR E   6       5.434  23.633   9.930  1.00 20.00           C
ATOM     42  CD2 TYR E   6       4.610  23.413   7.705  1.00 20.00           C
ATOM     43  CE1 TYR E   6       4.197  24.086  10.349  1.00 20.00           C
ATOM     44  CE2 TYR E   6       3.369  23.865   8.114  1.00 20.00           C
ATOM     45  CZ  TYR E   6       3.169  24.200   9.437  1.00 20.00           C
ATOM     46  OH  TYR E   6       1.935  24.650   9.849  1.00 20.00           O
ATOM     47  N   THR E   7      10.161  22.688   7.825  1.00 20.00           N
ATOM     48  CA  THR E   7      11.380  21.976   7.459  1.00 20.00           C
ATOM     49  C   THR E   7      12.441  22.134   8.541  1.00 20.00           C
ATOM     50  O   THR E   7      13.199  21.204   8.823  1.00 20.00           O
ATOM     51  CB  THR E   7      11.948  22.456   6.109  1.00 20.00           C
ATOM     52  OG1 THR E   7      11.987  23.888   6.086  1.00 20.00           O
ATOM     53  CG2 THR E   7      11.086  21.957   4.959  1.00 20.00           C
ATOM     54  N   ILE E   8      12.489  23.318   9.145  1.00 20.00           N
ATOM     55  CA  ILE E   8      13.433  23.596  10.220  1.00 20.00           C
ATOM     56  C   ILE E   8      13.128  22.736  11.441  1.00 20.00           C
ATOM     57  O   ILE E   8      14.038  22.232  12.099  1.00 20.00           O
ATOM     58  CB  ILE E   8      13.437  25.097  10.601  1.00 20.00           C
ATOM     59  CG1 ILE E   8      14.259  25.898   9.591  1.00 20.00           C
ATOM     60  CG2 ILE E   8      14.003  25.309  11.997  1.00 20.00           C
ATOM     61  CD1 ILE E   8      14.305  27.383   9.880  1.00 20.00           C
ATOM     62  N   HIS E   9      11.843  22.573  11.737  1.00 20.00           N
ATOM     63  CA  HIS E   9      11.415  21.732  12.848  1.00 20.00           C
ATOM     64  C   HIS E   9      11.639  20.260  12.519  1.00 20.00           C
ATOM     65  O   HIS E   9      11.974  19.460  13.395  1.00 20.00           O
ATOM     66  CB  HIS E   9       9.943  21.982  13.182  1.00 20.00           C
ATOM     67  CG  HIS E   9       9.681  23.326  13.789  1.00 20.00           C
ATOM     68  ND1 HIS E   9      10.679  24.249  14.013  1.00 20.00           N
ATOM     69  CD2 HIS E   9       8.534  23.900  14.221  1.00 20.00           C
ATOM     70  CE1 HIS E   9      10.158  25.335  14.556  1.00 20.00           C
ATOM     71  NE2 HIS E   9       8.858  25.149  14.693  1.00 20.00           N
ATOM     72  N   SER E  10      11.450  19.912  11.250  1.00 20.00           N
ATOM     73  CA  SER E  10      11.664  18.546  10.785  1.00 20.00           C
ATOM     74  C   SER E  10      13.125  18.147  10.944  1.00 20.00           C
ATOM     75  O   SER E  10      13.431  17.029  11.357  1.00 20.00           O
ATOM     76  CB  SER E  10      11.237  18.394   9.324  1.00 20.00           C
ATOM     77  OG  SER E  10      12.163  19.021   8.455  1.00 20.00           O
ATOM     78  N   VAL E  11      14.025  19.068  10.613  1.00 20.00           N
ATOM     79  CA  VAL E  11      15.454  18.835  10.782  1.00 20.00           C
ATOM     80  C   VAL E  11      15.820  18.850  12.262  1.00 20.00           C
ATOM     81  O   VAL E  11      16.689  18.095  12.705  1.00 20.00           O
ATOM     82  CB  VAL E  11      16.292  19.885  10.015  1.00 20.00           C
ATOM     83  CG1 VAL E  11      17.776  19.723  10.312  1.00 20.00           C
ATOM     84  CG2 VAL E  11      16.036  19.776   8.520  1.00 20.00           C
ATOM     85  N   ASP E  12      15.147  19.714  13.017  1.00 20.00           N
ATOM     86  CA  ASP E  12      15.358  19.816  14.457  1.00 20.00           C
ATOM     87  C   ASP E  12      15.062  18.486  15.137  1.00 20.00           C
ATOM     88  O   ASP E  12      15.780  18.073  16.046  1.00 20.00           O
ATOM     89  CB  ASP E  12      14.495  20.925  15.064  1.00 20.00           C
ATOM     90  CG  ASP E  12      15.107  22.301  14.890  1.00 20.00           C
ATOM     91  OD1 ASP E  12      16.350  22.397  14.824  1.00 20.00           O
ATOM     92  OD2 ASP E  12      14.344  23.288  14.820  1.00 20.00           O
ATOM     93  N   ALA E  13      14.001  17.822  14.690  1.00 20.00           N
ATOM     94  CA  ALA E  13      13.688  16.481  15.167  1.00 20.00           C
ATOM     95  C   ALA E  13      14.741  15.506  14.655  1.00 20.00           C
ATOM     96  O   ALA E  13      15.423  14.826  15.444  1.00 20.00           O
ATOM     97  CB  ALA E  13      12.303  16.062  14.706  1.00 20.00           C
ATOM     98  N   PHE E  14      14.885  15.479  13.328  1.00 20.00           N
ATOM     99  CA  PHE E  14      15.842  14.623  12.623  1.00 20.00           C
ATOM    100  C   PHE E  14      17.192  14.548  13.323  1.00 20.00           C
ATOM    101  O   PHE E  14      17.881  13.533  13.245  1.00 20.00           O
ATOM    102  CB  PHE E  14      16.026  15.083  11.174  1.00 20.00           C
ATOM    103  CG  PHE E  14      14.911  14.666  10.258  1.00 20.00           C
ATOM    104  CD1 PHE E  14      14.035  13.657  10.624  1.00 20.00           C
ATOM    105  CD2 PHE E  14      14.739  15.281   9.029  1.00 20.00           C
ATOM    106  CE1 PHE E  14      13.008  13.271   9.783  1.00 20.00           C
ATOM    107  CE2 PHE E  14      13.714  14.900   8.183  1.00 20.00           C
ATOM    108  CZ  PHE E  14      12.848  13.893   8.561  1.00 20.00           C
ATOM    109  N   ALA E  15      17.558  15.630  14.003  1.00 20.00           N
ATOM    110  CA  ALA E  15      18.699  15.619  14.904  1.00 20.00           C
ATOM    111  C   ALA E  15      18.268  15.190  16.307  1.00 20.00           C
ATOM    112  O   ALA E  15      19.035  14.538  17.024  1.00 20.00           O
ATOM    113  CB  ALA E  15      19.357  16.988  14.944  1.00 20.00           C
ATOM    114  N   GLU E  16      17.036  15.542  16.687  1.00 20.00           N
ATOM    115  CA  GLU E  16      16.549  15.284  18.047  1.00 20.00           C
ATOM    116  C   GLU E  16      16.685  13.831  18.518  1.00 20.00           C
ATOM    117  O   GLU E  16      17.064  13.664  19.676  1.00 20.00           O
ATOM    118  CB  GLU E  16      15.119  15.792  18.273  1.00 20.00           C
ATOM    119  CG  GLU E  16      15.040  17.145  18.965  1.00 20.00           C
ATOM    120  CD  GLU E  16      13.612  17.612  19.170  1.00 20.00           C
ATOM    121  OE1 GLU E  16      12.682  16.879  18.772  1.00 20.00           O
ATOM    122  OE2 GLU E  16      13.420  18.711  19.730  1.00 20.00           O
TER
"""

def exercise(prefix="tst_helix_sheet_recs_as_pdb_files"):
  of = open(prefix+".pdb", "w")
  print(pdb_str, file=of)
  of.close()
  xrs1 = iotbx.pdb.input(file_name=prefix+".pdb").xray_structure_simple()
  assert not easy_run.call("phenix.helix_sheet_recs_as_pdb_files %s"%(prefix+".pdb"))
  xrs2 = iotbx.pdb.input(
    file_name="HELIX_1_1_ALA_E_1_ALA_E_16_1_16.pdb").xray_structure_simple(crystal_symmetry=xrs1.crystal_symmetry())
  fc1 = xrs1.structure_factors(d_min=3).f_calc()
  fc2 = fc1.structure_factors_from_scatterers(
    xray_structure=xrs2).f_calc()
  fc1=flex.abs(abs(fc1).data())
  fc2=flex.abs(abs(fc2).data())
  assert flex.sum(flex.abs(fc1-fc2))/flex.sum(flex.abs(fc1+fc2)) < 1.e-3

if (__name__ == "__main__"):
  exercise()
