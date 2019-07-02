from __future__ import absolute_import, division, print_function

from mmtbx.command_line.fix_rotamer_outliers import run
import iotbx.pdb
import mmtbx.model

pdb_str_1 = """\
CRYST1   18.486   29.189   16.919  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A  58       6.848  19.603  11.185  1.00  7.73           N
ATOM      2  CA  TYR A  58       7.439  18.674  10.233  1.00  8.65           C
ATOM      3  C   TYR A  58       7.228  17.251  10.706  1.00  9.84           C
ATOM      4  O   TYR A  58       7.212  16.994  11.919  1.00  8.58           O
ATOM      5  CB  TYR A  58       8.929  18.959  10.045  1.00 20.00           C
ATOM      6  CG  TYR A  58       9.344  20.355  10.457  1.00 20.00           C
ATOM      7  CD1 TYR A  58       8.393  21.339  10.700  1.00 20.00           C
ATOM      8  CD2 TYR A  58      10.683  20.685  10.603  1.00 20.00           C
ATOM      9  CE1 TYR A  58       8.768  22.614  11.076  1.00 20.00           C
ATOM     10  CE2 TYR A  58      11.068  21.959  10.980  1.00 20.00           C
ATOM     11  CZ  TYR A  58      10.107  22.918  11.215  1.00 20.00           C
ATOM     12  OH  TYR A  58      10.483  24.189  11.588  1.00 20.00           O
ATOM     13  N   SER A  59       7.076  16.327   9.756  1.00  9.65           N
ATOM     14  CA  SER A  59       6.909  14.912  10.080  1.00  5.80           C
ATOM     15  C   SER A  59       7.795  14.069   9.169  1.00 10.35           C
ATOM     16  O   SER A  59       8.161  14.502   8.075  1.00 10.56           O
ATOM     17  CB  SER A  59       5.446  14.487   9.943  1.00 20.00           C
ATOM     18  OG  SER A  59       5.000  14.619   8.604  1.00 20.00           O
"""

pdb_str_2 = """\
CRYST1  130.848  130.848  149.639  90.00  90.00 120.00 P 64 2 2     12
ATOM     95  N   ARG A  30     -18.833  17.433 -14.879  1.00191.24           N
ATOM     96  CA  ARG A  30     -18.615  17.619 -16.290  1.00191.04           C
ATOM     97  C   ARG A  30     -18.032  19.002 -16.408  1.00190.51           C
ATOM     98  O   ARG A  30     -18.720  19.978 -16.161  1.00190.07           O
ATOM     99  CB  ARG A  30     -19.941  17.558 -17.076  1.00191.84           C
ATOM    100  CG  ARG A  30     -20.909  16.521 -16.584  1.00192.91           C
ATOM    101  CD  ARG A  30     -22.176  16.530 -17.408  1.00193.90           C
ATOM    102  NE  ARG A  30     -23.271  15.875 -16.703  1.00195.21           N
ATOM    103  CZ  ARG A  30     -24.498  15.744 -17.196  1.00195.78           C
ATOM    104  NH1 ARG A  30     -24.776  16.224 -18.399  1.00196.03           N
ATOM    105  NH2 ARG A  30     -25.449  15.149 -16.488  1.00195.90           N
ATOM    106  N   LYS A  31     -16.776  19.084 -16.815  1.00190.39           N
ATOM    107  CA  LYS A  31     -16.164  20.381 -17.016  1.00190.65           C
ATOM    108  C   LYS A  31     -17.139  21.024 -18.016  1.00191.02           C
ATOM    109  O   LYS A  31     -17.240  22.238 -18.137  1.00190.75           O
ATOM    110  CB  LYS A  31     -14.782  20.200 -17.673  1.00190.39           C
ATOM    111  CG  LYS A  31     -13.794  21.338 -17.444  1.00189.92           C
ATOM    112  CD  LYS A  31     -12.448  20.983 -18.063  1.00189.44           C
ATOM    113  CE  LYS A  31     -11.369  21.951 -17.632  1.00189.51           C
ATOM    114  NZ  LYS A  31     -10.046  21.567 -18.180  1.00189.58           N
ATOM    115  N   ASP A  32     -17.885  20.148 -18.681  1.00192.02           N
ATOM    116  CA  ASP A  32     -18.842  20.462 -19.742  1.00192.58           C
ATOM    117  C   ASP A  32     -19.950  21.497 -19.544  1.00191.35           C
ATOM    118  O   ASP A  32     -20.097  22.367 -20.384  1.00191.52           O
ATOM    119  CB  ASP A  32     -19.488  19.164 -20.226  1.00195.14           C
ATOM    120  CG  ASP A  32     -18.506  17.994 -20.259  1.00197.19           C
ATOM    121  OD1 ASP A  32     -18.792  16.987 -20.944  1.00197.86           O
ATOM    122  OD2 ASP A  32     -17.446  18.072 -19.596  1.00198.95           O
ATOM    123  N   ILE A  33     -20.744  21.422 -18.481  1.00189.52           N
ATOM    124  CA  ILE A  33     -21.803  22.408 -18.347  1.00187.36           C
ATOM    125  C   ILE A  33     -21.351  23.724 -17.684  1.00186.08           C
ATOM    126  O   ILE A  33     -22.173  24.531 -17.271  1.00185.90           O
ATOM    127  CB  ILE A  33     -23.073  21.768 -17.669  1.00186.97           C
ATOM    128  CG1 ILE A  33     -22.681  20.853 -16.498  1.00186.29           C
ATOM    129  CG2 ILE A  33     -23.832  20.928 -18.700  1.00186.44           C
ATOM    130  CD1 ILE A  33     -23.867  20.142 -15.833  1.00184.51           C
ATOM    131  N   SER A  34     -20.030  23.933 -17.618  1.00184.40           N
ATOM    132  CA  SER A  34     -19.424  25.157 -17.069  1.00182.95           C
ATOM    133  C   SER A  34     -19.038  25.962 -18.306  1.00182.29           C
ATOM    134  O   SER A  34     -18.590  25.403 -19.300  1.00182.43           O
ATOM    135  CB  SER A  34     -18.199  24.839 -16.225  1.00182.36           C
ATOM    136  OG  SER A  34     -16.990  25.070 -16.946  1.00182.11           O
ATOM    137  N   GLU A  35     -19.219  27.269 -18.236  1.00181.35           N
ATOM    138  CA  GLU A  35     -19.013  28.138 -19.394  1.00180.75           C
ATOM    139  C   GLU A  35     -17.702  28.735 -19.927  1.00179.75           C
ATOM    140  O   GLU A  35     -16.588  28.226 -19.772  1.00178.29           O
ATOM    141  CB  GLU A  35     -19.990  29.296 -19.259  1.00181.98           C
ATOM    142  CG  GLU A  35     -21.338  29.024 -19.904  1.00183.50           C
ATOM    143  CD  GLU A  35     -21.656  27.558 -20.136  1.00184.28           C
ATOM    144  OE1 GLU A  35     -21.062  26.961 -21.059  1.00183.95           O
ATOM    145  OE2 GLU A  35     -22.521  27.021 -19.411  1.00185.10           O
ATOM    146  N   ASN A  36     -17.914  29.823 -20.663  1.00179.92           N
ATOM    147  CA  ASN A  36     -16.844  30.610 -21.229  1.00180.62           C
ATOM    148  C   ASN A  36     -16.281  31.058 -19.926  1.00181.23           C
ATOM    149  O   ASN A  36     -15.103  31.341 -19.819  1.00181.62           O
ATOM    150  CB  ASN A  36     -17.388  31.832 -21.941  1.00180.82           C
ATOM    151  CG  ASN A  36     -17.645  31.611 -23.392  1.00182.08           C
ATOM    152  OD1 ASN A  36     -16.729  31.377 -24.171  1.00182.92           O
ATOM    153  ND2 ASN A  36     -18.907  31.705 -23.778  1.00183.27           N
ATOM    154  N   ALA A  37     -17.158  31.116 -18.928  1.00181.85           N
ATOM    155  CA  ALA A  37     -16.773  31.516 -17.583  1.00182.29           C
ATOM    156  C   ALA A  37     -15.531  30.754 -17.159  1.00182.76           C
ATOM    157  O   ALA A  37     -14.583  31.328 -16.615  1.00182.32           O
ATOM    158  CB  ALA A  37     -17.923  31.251 -16.606  1.00181.85           C
ATOM    159  N   LEU A  38     -15.547  29.456 -17.435  1.00183.67           N
ATOM    160  CA  LEU A  38     -14.427  28.579 -17.111  1.00184.12           C
ATOM    161  C   LEU A  38     -13.148  29.073 -17.763  1.00184.05           C
ATOM    162  O   LEU A  38     -12.040  28.809 -17.283  1.00183.42           O
ATOM    163  CB  LEU A  38     -14.738  27.147 -17.572  1.00184.55           C
ATOM    164  CG  LEU A  38     -13.625  26.095 -17.588  1.00184.94           C
ATOM    165  CD1 LEU A  38     -14.247  24.709 -17.608  1.00185.45           C
ATOM    166  CD2 LEU A  38     -12.721  26.297 -18.797  1.00184.46           C
ATOM    167  N   LYS A  39     -13.327  29.801 -18.859  1.00184.14           N
ATOM    168  CA  LYS A  39     -12.226  30.363 -19.628  1.00183.59           C
ATOM    169  C   LYS A  39     -11.302  31.159 -18.711  1.00182.98           C
ATOM    170  O   LYS A  39     -10.093  30.912 -18.655  1.00183.18           O
ATOM    171  CB  LYS A  39     -12.790  31.263 -20.749  1.00183.64           C
ATOM    172  CG  LYS A  39     -11.851  31.587 -21.912  1.00184.17           C
ATOM    173  CD  LYS A  39     -12.620  32.210 -23.080  1.00184.66           C
ATOM    174  CE  LYS A  39     -11.678  32.675 -24.183  1.00185.01           C
ATOM    175  NZ  LYS A  39     -12.408  33.271 -25.337  1.00185.25           N
"""

def get_necessary_inputs(pdb_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
  return mmtbx.model.manager(
      model_input = pdb_inp,
      build_grm= True)



def exercise_1():
  """ 58 is outlier """
  model = get_necessary_inputs(pdb_str_1)
  # pdb_h.write_pdb_file("fix_rot_out_ex1_start.pdb")
  model = run(
      args=[],
      model=model)
  rotamers = []
  # pdb_h.write_pdb_file("fix_rot_out_ex1_end.pdb")
  for res in model.get_hierarchy().only_chain().only_conformer().residues():
    rotamers.append(model.get_rotamer_manager().evaluate_residue(res))
  # print rotamers
  assert rotamers == ['m-80', 'p'], rotamers

def exercise_2():
  model = get_necessary_inputs(pdb_str_2)
  model = run(
      args=[],
      model = model)
  rotamers = []
  for res in model.get_hierarchy().only_chain().only_conformer().residues():
    rotamers.append(model.get_rotamer_manager().evaluate_residue(res))
  assert rotamers == ['mtt180', 'tttt', 'm-30', 'pt', 'p', 'mm-30', 'm-40',
    'EXCEPTION', 'tt', 'tttt'], rotamers

if (__name__ == "__main__"):
  exercise_1()
  exercise_2()
  print("OK")
