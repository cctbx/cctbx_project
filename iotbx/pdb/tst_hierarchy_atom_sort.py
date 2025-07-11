"""Test sorting of hierarchy"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from six.moves import zip


def validate_result(h1, h2):
  assert h1.is_similar_hierarchy(h2)
  atoms1 = h1.atoms()
  atoms2 = h2.atoms()
  for a1, a2 in zip(atoms1, atoms2):
    # print dir(a1)
    # print a1.id_str()
    # print a1.format_anisou_record()
    assert a1.id_str() == a2.id_str(), "%s != %s" % (a1.id_str(), a2.id_str())
    # cutting ANISOU because we do atoms_reset_serial()
    assert a1.format_anisou_record()[13:] == a2.format_anisou_record()[13:],"%s != %s"%(
        a1.format_anisou_record(), a2.format_anisou_record())

def exercise_1():
  """
  simple one residue test
  """
  in_pdb_str = """\
ATOM      1  N   THR C   1      12.121   9.329  18.086  1.00 10.00           N
ATOM      3  C   THR C   1      13.707  10.622  16.718  1.00 10.00           C
ATOM      6  OG1 THR C   1      10.087  11.287  17.482  1.00 10.00           O
ATOM      4  O   THR C   1      14.493  10.814  17.645  1.00 10.00           O
ATOM      2  CA  THR C   1      12.245  10.284  16.991  1.00 10.00           C
ATOM      5  CB  THR C   1      11.474  11.584  17.284  1.00 10.00           C
ATOM      7  CG2 THR C   1      11.619  12.563  16.129  1.00 10.00           C
  """
  out_pdb_str = """\
ATOM      1  N   THR C   1      12.121   9.329  18.086  1.00 10.00           N
ATOM      2  CA  THR C   1      12.245  10.284  16.991  1.00 10.00           C
ATOM      3  C   THR C   1      13.707  10.622  16.718  1.00 10.00           C
ATOM      4  O   THR C   1      14.493  10.814  17.645  1.00 10.00           O
ATOM      5  CB  THR C   1      11.474  11.584  17.284  1.00 10.00           C
ATOM      6  OG1 THR C   1      10.087  11.287  17.482  1.00 10.00           O
ATOM      7  CG2 THR C   1      11.619  12.563  16.129  1.00 10.00           C
  """
  in_pdb_str1 = """\
ATOM    233  N   ARG A  16      18.484  -1.972   5.177  1.00 57.75           N
ANISOU  233  N   ARG A  16     7670   7185   7086    504    345   1028       N
ATOM    234  CA  ARG A  16      19.270  -3.198   5.125  1.00 40.85           C
ANISOU  234  CA  ARG A  16     5601   5022   4897    530    438   1131       C
ATOM    235  CB  ARG A  16      19.983  -3.314   3.779  1.00 71.44           C
ANISOU  235  CB  ARG A  16     9452   8902   8790    573    399   1208       C
ATOM    236  CG  ARG A  16      21.082  -4.355   3.747  1.00 75.49           C
ANISOU  236  CG  ARG A  16    10032   9401   9249    612    465   1342       C
ATOM    237  CD  ARG A  16      21.602  -4.564   2.336  1.00 98.55           C
ANISOU  237  CD  ARG A  16    12919  12326  12199    634    431   1371       C
ATOM    238  NE  ARG A  16      20.667  -5.352   1.531  1.00130.20           N
ANISOU  238  NE  ARG A  16    16937  16300  16233    621    526   1310       N
ATOM    239  CZ  ARG A  16      19.766  -4.853   0.687  1.00155.16           C
ANISOU  239  CZ  ARG A  16    20034  19465  19455    600    486   1204       C
ATOM    240  NH1 ARG A  16      19.658  -3.545   0.494  1.00194.42           N
ANISOU  240  NH1 ARG A  16    24939  24468  24465    602    359   1173       N
ATOM    241  NH2 ARG A  16      18.975  -5.678   0.012  1.00136.93           N
ANISOU  241  NH2 ARG A  16    17732  17127  17169    577    574   1128       N
ATOM    242  C   ARG A  16      20.243  -3.255   6.301  1.00 31.70           C
ANISOU  242  C   ARG A  16     4497   3873   3675    553    440   1238       C
ATOM    243  O   ARG A  16      20.431  -4.315   6.908  1.00 56.41           O
ANISOU  243  O   ARG A  16     7702   6975   6754    552    549   1275       O
ATOM    244  H   ARG A  16      18.543  -1.461   4.487  1.00 69.30           H
ATOM    245  HA  ARG A  16      18.669  -3.955   5.205  1.00 49.02           H
ATOM    246  HB2 ARG A  16      19.331  -3.552   3.101  1.00 85.73           H
ATOM    247  HB3 ARG A  16      20.382  -2.457   3.563  1.00 85.73           H
ATOM    248  HG2 ARG A  16      21.820  -4.059   4.303  1.00 90.58           H
ATOM    249  HG3 ARG A  16      20.733  -5.200   4.072  1.00 90.58           H
ATOM    250  HD2 ARG A  16      21.720  -3.702   1.907  1.00118.26           H
ATOM    251  HD3 ARG A  16      22.446  -5.040   2.373  1.00118.26           H
ATOM    252  HE  ARG A  16      20.703  -6.208   1.609  1.00156.24           H
ATOM    253 HH11 ARG A  16      20.168  -3.003   0.926  1.00233.31           H
ATOM    254 HH12 ARG A  16      19.073  -3.238  -0.057  1.00233.31           H
ATOM    255 HH21 ARG A  16      19.040  -6.527   0.127  1.00164.32           H
ATOM    256 HH22 ARG A  16      18.393  -5.362  -0.538  1.00164.32           H
  """
  out_pdb_str1 = """\
ATOM      1  N   ARG A  16      18.484  -1.972   5.177  1.00 57.75           N
ANISOU    1  N   ARG A  16     7670   7185   7086    504    345   1028       N
ATOM      2  CA  ARG A  16      19.270  -3.198   5.125  1.00 40.85           C
ANISOU    2  CA  ARG A  16     5601   5022   4897    530    438   1131       C
ATOM      3  C   ARG A  16      20.243  -3.255   6.301  1.00 31.70           C
ANISOU    3  C   ARG A  16     4497   3873   3675    553    440   1238       C
ATOM      4  O   ARG A  16      20.431  -4.315   6.908  1.00 56.41           O
ANISOU    4  O   ARG A  16     7702   6975   6754    552    549   1275       O
ATOM      5  CB  ARG A  16      19.983  -3.314   3.779  1.00 71.44           C
ANISOU    5  CB  ARG A  16     9452   8902   8790    573    399   1208       C
ATOM      6  CG  ARG A  16      21.082  -4.355   3.747  1.00 75.49           C
ANISOU    6  CG  ARG A  16    10032   9401   9249    612    465   1342       C
ATOM      7  CD  ARG A  16      21.602  -4.564   2.336  1.00 98.55           C
ANISOU    7  CD  ARG A  16    12919  12326  12199    634    431   1371       C
ATOM      8  NE  ARG A  16      20.667  -5.352   1.531  1.00130.20           N
ANISOU    8  NE  ARG A  16    16937  16300  16233    621    526   1310       N
ATOM      9  CZ  ARG A  16      19.766  -4.853   0.687  1.00155.16           C
ANISOU    9  CZ  ARG A  16    20034  19465  19455    600    486   1204       C
ATOM     10  NH1 ARG A  16      19.658  -3.545   0.494  1.00194.42           N
ANISOU   10  NH1 ARG A  16    24939  24468  24465    602    359   1173       N
ATOM     11  NH2 ARG A  16      18.975  -5.678   0.012  1.00136.93           N
ANISOU   11  NH2 ARG A  16    17732  17127  17169    577    574   1128       N
ATOM     12  H   ARG A  16      18.543  -1.461   4.487  1.00 69.30           H
ATOM     13  HA  ARG A  16      18.669  -3.955   5.205  1.00 49.02           H
ATOM     14  HB2 ARG A  16      19.331  -3.552   3.101  1.00 85.73           H
ATOM     15  HB3 ARG A  16      20.382  -2.457   3.563  1.00 85.73           H
ATOM     16  HG2 ARG A  16      21.820  -4.059   4.303  1.00 90.58           H
ATOM     17  HG3 ARG A  16      20.733  -5.200   4.072  1.00 90.58           H
ATOM     18  HD2 ARG A  16      21.720  -3.702   1.907  1.00118.26           H
ATOM     19  HD3 ARG A  16      22.446  -5.040   2.373  1.00118.26           H
ATOM     20  HE  ARG A  16      20.703  -6.208   1.609  1.00156.24           H
ATOM     21 HH11 ARG A  16      20.168  -3.003   0.926  1.00233.31           H
ATOM     22 HH12 ARG A  16      19.073  -3.238  -0.057  1.00233.31           H
ATOM     23 HH21 ARG A  16      19.040  -6.527   0.127  1.00164.32           H
ATOM     24 HH22 ARG A  16      18.393  -5.362  -0.538  1.00164.32           H
TER
  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=out_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)

  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str1).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=out_pdb_str1).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)


def exercise_2():
  """
  multi chain test
  """
  in_pdb_str = """\
ATOM    674  N   PHE A  89     -13.997 -12.169 -16.080  1.00 30.39           N
ATOM    675  CA  PHE A  89     -14.205 -13.374 -15.271  1.00 31.46           C
ATOM    676  C   PHE A  89     -15.660 -13.679 -15.012  1.00 32.05           C
ATOM    677  O   PHE A  89     -16.062 -14.845 -15.199  1.00 33.50           O
ATOM    678  CB  PHE A  89     -13.457 -13.295 -13.964  1.00 30.50           C
ATOM    679  CG  PHE A  89     -12.122 -13.907 -14.019  1.00 29.32           C
ATOM    680  CD1 PHE A  89     -10.985 -13.176 -13.696  1.00 29.37           C
ATOM    681  CD2 PHE A  89     -11.991 -15.236 -14.376  1.00 31.82           C
ATOM    682  CE1 PHE A  89      -9.715 -13.749 -13.742  1.00 31.30           C
ATOM    683  CE2 PHE A  89     -10.753 -15.839 -14.430  1.00 29.53           C
ATOM    684  CZ  PHE A  89      -9.588 -15.093 -14.112  1.00 33.18           C
TER     685      PHE A  89
ATOM    686  N   GLN B   3       4.525  19.816   1.629  1.00 52.67           N
ATOM    687  CA  GLN B   3       5.708  20.537   1.020  1.00 53.13           C
ATOM    688  C   GLN B   3       5.683  20.664  -0.493  1.00 52.62           C
ATOM    689  O   GLN B   3       6.154  21.665  -1.066  1.00 51.62           O
ATOM    690  CB  GLN B   3       7.034  19.911   1.432  1.00 53.47           C
ATOM    691  CG  GLN B   3       7.831  20.744   2.462  1.00 55.56           C
ATOM    692  CD  GLN B   3       7.158  20.850   3.831  1.00 55.64           C
ATOM    693  OE1 GLN B   3       7.843  21.020   4.834  1.00 56.71           O
ATOM    694  NE2 GLN B   3       5.812  20.762   3.876  1.00 54.75           N
  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)


def exercise_3():
  """
  multi model test
  """
  in_pdb_str = """\
MODEL        1
ATOM    686  N   GLN B   3       4.525  19.816   1.629  1.00 52.67           N
ATOM    687  CA  GLN B   3       5.708  20.537   1.020  1.00 53.13           C
ATOM    688  C   GLN B   3       5.683  20.664  -0.493  1.00 52.62           C
ATOM    689  O   GLN B   3       6.154  21.665  -1.066  1.00 51.62           O
ATOM    690  CB  GLN B   3       7.034  19.911   1.432  1.00 53.47           C
ATOM    691  CG  GLN B   3       7.831  20.744   2.462  1.00 55.56           C
ATOM    692  CD  GLN B   3       7.158  20.850   3.831  1.00 55.64           C
ATOM    693  OE1 GLN B   3       7.843  21.020   4.834  1.00 56.71           O
ATOM    694  NE2 GLN B   3       5.812  20.762   3.876  1.00 54.75           N
ENDMDL
MODEL        2
ATOM    695  N   LEU B   4       5.127  19.654  -1.142  1.00 51.70           N
ATOM    696  CA  LEU B   4       4.838  19.817  -2.545  1.00 50.44           C
ATOM    697  C   LEU B   4       3.872  20.969  -2.808  1.00 50.31           C
ATOM    698  O   LEU B   4       4.335  21.988  -3.298  1.00 51.34           O
ATOM    699  CB  LEU B   4       4.446  18.508  -3.254  1.00 49.68           C
ATOM    700  CG  LEU B   4       4.223  18.577  -4.753  1.00 45.03           C
ATOM    701  CD1 LEU B   4       5.232  19.472  -5.480  1.00 41.90           C
ATOM    702  CD2 LEU B   4       4.267  17.162  -5.271  1.00 42.51           C
ENDMDL
  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)

def exercise_4():
  """
  alternative conformations test
  """
  in_pdb_str = """\
ATOM    146  N   TRP A   9      20.043  18.668  31.573  1.00  5.76           N
ATOM    147  CA  TRP A   9      19.105  18.066  32.480  1.00  5.43           C
ATOM    148  C   TRP A   9      19.223  16.553  32.483  1.00  6.09           C
ATOM    149  O   TRP A   9      19.315  15.939  31.413  1.00  8.92           O
ATOM    150  CB ATRP A   9      17.698  18.302  31.969  0.75  5.70           C
ATOM    151  CB BTRP A   9      17.649  18.637  32.191  0.25  6.08           C
ATOM    152  CG ATRP A   9      17.381  19.736  31.708  0.75  6.14           C
ATOM    153  CG BTRP A   9      17.332  20.104  32.297  0.25  6.49           C
ATOM    154  CD1ATRP A   9      17.385  20.336  30.500  0.75  7.28           C
ATOM    155  CD1BTRP A   9      17.086  20.842  33.421  0.25  7.86           C
ATOM    156  CD2ATRP A   9      16.969  20.687  32.646  0.75  7.31           C
ATOM    157  CD2BTRP A   9      17.137  21.008  31.214  0.25  6.86           C
ATOM    158  NE1ATRP A   9      17.005  21.638  30.598  0.75  9.09           N
ATOM    159  NE1BTRP A   9      16.817  22.137  33.085  0.25  9.75           N
ATOM    160  CE2ATRP A   9      16.727  21.889  31.902  0.75  8.75           C
ATOM    161  CE2BTRP A   9      16.827  22.288  31.730  0.25  8.56           C
ATOM    162  CE3ATRP A   9      16.752  20.668  34.034  0.75  8.75           C
ATOM    163  CE3BTRP A   9      17.208  20.883  29.825  0.25  7.27           C
ATOM    164  CZ2ATRP A   9      16.306  23.039  32.535  0.75 11.64           C
ATOM    165  CZ2BTRP A   9      16.597  23.398  30.929  0.25  9.81           C
ATOM    166  CZ3ATRP A   9      16.317  21.845  34.590  0.75 11.46           C
ATOM    167  CZ3BTRP A   9      16.994  21.984  29.037  0.25  8.89           C
ATOM    168  CH2ATRP A   9      16.103  22.991  33.865  0.75 12.15           C
ATOM    169  CH2BTRP A   9      16.699  23.238  29.593  0.25 10.35           C
ATOM    170  H   TRP A   9      19.784  18.842  30.771  1.00  6.91           H
ATOM    171  HA ATRP A   9      19.213  18.427  33.385  0.50  6.51           H
ATOM    172  HD1 TRP A   9      17.101  20.510  34.289  0.25  9.43           H
ATOM    173  HZ2ATRP A   9      16.167  23.824  32.057  0.75 13.97           H
ATOM    174  HZ2BTRP A   9      16.380  24.223  31.300  0.25 11.77           H
ATOM    175  HZ3ATRP A   9      16.158  21.868  35.506  0.75 13.75           H
ATOM    176  HH2ATRP A   9      16.950  22.201  29.951  0.75 10.90           H
  """
  out_pdb_str = """\
ATOM      1  N   TRP A   9      20.043  18.668  31.573  1.00  5.76           N
ATOM      2  CA  TRP A   9      19.105  18.066  32.480  1.00  5.43           C
ATOM      3  C   TRP A   9      19.223  16.553  32.483  1.00  6.09           C
ATOM      4  O   TRP A   9      19.315  15.939  31.413  1.00  8.92           O
ATOM      5  H   TRP A   9      19.784  18.842  30.771  1.00  6.91           H
ATOM      6  HD1 TRP A   9      17.101  20.510  34.289  0.25  9.43           H
ATOM      7  CB ATRP A   9      17.698  18.302  31.969  0.75  5.70           C
ATOM      8  CG ATRP A   9      17.381  19.736  31.708  0.75  6.14           C
ATOM      9  CD1ATRP A   9      17.385  20.336  30.500  0.75  7.28           C
ATOM     10  CD2ATRP A   9      16.969  20.687  32.646  0.75  7.31           C
ATOM     11  NE1ATRP A   9      17.005  21.638  30.598  0.75  9.09           N
ATOM     12  CE2ATRP A   9      16.727  21.889  31.902  0.75  8.75           C
ATOM     13  CE3ATRP A   9      16.752  20.668  34.034  0.75  8.75           C
ATOM     14  CZ2ATRP A   9      16.306  23.039  32.535  0.75 11.64           C
ATOM     15  CZ3ATRP A   9      16.317  21.845  34.590  0.75 11.46           C
ATOM     16  CH2ATRP A   9      16.103  22.991  33.865  0.75 12.15           C
ATOM     17  HA ATRP A   9      19.213  18.427  33.385  0.50  6.51           H
ATOM     18  HZ2ATRP A   9      16.167  23.824  32.057  0.75 13.97           H
ATOM     19  HZ3ATRP A   9      16.158  21.868  35.506  0.75 13.75           H
ATOM     20  HH2ATRP A   9      16.950  22.201  29.951  0.75 10.90           H
ATOM     21  CB BTRP A   9      17.649  18.637  32.191  0.25  6.08           C
ATOM     22  CG BTRP A   9      17.332  20.104  32.297  0.25  6.49           C
ATOM     23  CD1BTRP A   9      17.086  20.842  33.421  0.25  7.86           C
ATOM     24  CD2BTRP A   9      17.137  21.008  31.214  0.25  6.86           C
ATOM     25  NE1BTRP A   9      16.817  22.137  33.085  0.25  9.75           N
ATOM     26  CE2BTRP A   9      16.827  22.288  31.730  0.25  8.56           C
ATOM     27  CE3BTRP A   9      17.208  20.883  29.825  0.25  7.27           C
ATOM     28  CZ2BTRP A   9      16.597  23.398  30.929  0.25  9.81           C
ATOM     29  CZ3BTRP A   9      16.994  21.984  29.037  0.25  8.89           C
ATOM     30  CH2BTRP A   9      16.699  23.238  29.593  0.25 10.35           C
ATOM     31  HZ2BTRP A   9      16.380  24.223  31.300  0.25 11.77           H
  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=out_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)


def exercise_5():
  """
  insertion codes test + small atom names + oxt
  """
  in_pdb_str = """\
ATOM   2527  N   THR B  97A     -2.797 -21.183  13.289  1.00166.13           N
ATOM   2528  CA  THR B  97A     -3.395 -20.517  14.441  1.00156.77           C
ATOM   2529  cb  THR B  97A     -3.659 -19.020  14.175  1.00154.53           C
ATOM   2530  Og1 THR B  97A     -4.691 -18.877  13.191  1.00160.68           O
ATOM   2531  Cg2 THR B  97A     -2.392 -18.317  13.706  1.00141.48           C
ATOM   2532  C   THR B  97A     -4.707 -21.194  14.819  1.00135.28           C
ATOM   2533  O   THR B  97A     -5.201 -22.091  14.130  1.00138.82           O
ATOM   2534  OxT THR B  97A     -5.302 -20.848  15.835  1.00125.49           O
ATOM   2535  H   THR B  97A     -3.367 -21.412  12.686  1.00199.36           H
ATOM   2536  HA  THR B  97A     -2.790 -20.587  15.197  1.00188.13           H
ATOM   2537  HB  THR B  97A     -3.948 -18.599  15.000  1.00185.44           H
ATOM   2538  HG1 THR B  97A     -4.837 -18.063  13.045  1.00192.81           H
ATOM   2539 HG21 THR B  97A     -2.574 -17.378  13.544  1.00169.78           H
ATOM   2540 HG22 THR B  97A     -1.701 -18.391  14.383  1.00169.78           H
ATOM   2541 HG23 THR B  97A     -2.073 -18.723  12.885  1.00169.78           H
  """
  out_pdb_str = """\
ATOM      1  N   THR B  97A     -2.797 -21.183  13.289  1.00166.13           N
ATOM      2  CA  THR B  97A     -3.395 -20.517  14.441  1.00156.77           C
ATOM      3  C   THR B  97A     -4.707 -21.194  14.819  1.00135.28           C
ATOM      4  O   THR B  97A     -5.201 -22.091  14.130  1.00138.82           O
ATOM      5  cb  THR B  97A     -3.659 -19.020  14.175  1.00154.53           C
ATOM      6  Og1 THR B  97A     -4.691 -18.877  13.191  1.00160.68           O
ATOM      7  Cg2 THR B  97A     -2.392 -18.317  13.706  1.00141.48           C
ATOM      8  OxT THR B  97A     -5.302 -20.848  15.835  1.00125.49           O
ATOM      9  H   THR B  97A     -3.367 -21.412  12.686  1.00199.36           H
ATOM     10  HA  THR B  97A     -2.790 -20.587  15.197  1.00188.13           H
ATOM     11  HB  THR B  97A     -3.948 -18.599  15.000  1.00185.44           H
ATOM     12  HG1 THR B  97A     -4.837 -18.063  13.045  1.00192.81           H
ATOM     13 HG21 THR B  97A     -2.574 -17.378  13.544  1.00169.78           H
ATOM     14 HG22 THR B  97A     -1.701 -18.391  14.383  1.00169.78           H
ATOM     15 HG23 THR B  97A     -2.073 -18.723  12.885  1.00169.78           H
"""
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=out_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)

def exercise_6():
  """
  segids test
  """
  in_pdb_str = """\
ATOM   1270  N   ASP B  78      15.391  19.770  -5.786  1.00 23.52      seg1 N
ATOM   1271  CA  ASP B  78      16.810  19.893  -6.020  1.00 23.79      seg1 C
ATOM   1272  C   ASP B  78      17.210  21.335  -6.226  1.00 23.90      seg1 C
ATOM   1273  O   ASP B  78      18.173  21.750  -5.659  1.00 22.55      seg1 O
ATOM   1274  CB  ASP B  78      17.224  19.124  -7.214  1.00 24.08      seg1 C
ATOM   1275  CG  ASP B  78      17.134  17.633  -6.986  1.00 27.72      seg1 C
ATOM   1276  OD1 ASP B  78      17.487  16.952  -7.962  1.00 29.41      seg1 O
ATOM   1277  OD2 ASP B  78      16.728  17.173  -5.867  1.00 26.25      seg1 O
  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)


def exercise_7():
  """
  nucleic acids: DNA
  """
  in_pdb_str = """\
HEADER    DNA                                     23-JUL-13   4LTG
CRYST1   44.300   44.300   27.550  90.00  90.00  90.00 P 43 21 2     8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.022573  0.000000  0.000000        0.00000
SCALE2      0.000000  0.022573  0.000000        0.00000
SCALE3      0.000000  0.000000  0.036298        0.00000
ATOM      1  O5'  DT A   1     -32.333  29.688  19.797  1.00 93.95           O
ATOM      2  C5'  DT A   1     -31.362  28.983  18.995  1.00 72.88           C
ATOM      3  C4'  DT A   1     -30.925  29.834  17.826  1.00 71.64           C
ATOM      4  O4'  DT A   1     -30.444  31.085  18.365  1.00 82.53           O
ATOM      5  C3'  DT A   1     -29.767  29.260  17.021  1.00 54.91           C
ATOM      6  O3'  DT A   1     -30.292  28.508  15.934  1.00 73.97           O
ATOM      7  C2'  DT A   1     -28.986  30.493  16.584  1.00 88.67           C
ATOM      8  C1'  DT A   1     -29.345  31.571  17.616  1.00 77.37           C
ATOM      9  N1   DT A   1     -28.284  31.961  18.577  1.00 61.15           N
ATOM     10  C2   DT A   1     -27.981  33.304  18.661  1.00 41.64           C
ATOM     11  O2   DT A   1     -28.513  34.157  17.947  1.00 42.99           O
ATOM     12  N3   DT A   1     -26.971  33.606  19.543  1.00 27.38           N
ATOM     13  C4   DT A   1     -26.320  32.730  20.389  1.00 26.18           C
ATOM     14  O4   DT A   1     -25.463  33.133  21.131  1.00 31.78           O
ATOM     15  C5   DT A   1     -26.697  31.350  20.271  1.00 39.53           C
ATOM     16  C7   DT A   1     -26.007  30.342  21.135  1.00 48.69           C
ATOM     17  C6   DT A   1     -27.647  31.030  19.377  1.00 69.53           C
ATOM     18  H5'  DT A   1     -31.780  28.040  18.638  1.00 70.41           H
ATOM     19 H5''  DT A   1     -30.480  28.765  19.599  1.00 67.70           H
ATOM     20  H4'  DT A   1     -31.767  30.356  17.220  1.00 72.07           H
ATOM     21  H3'  DT A   1     -29.137  28.640  17.674  1.00 69.73           H
ATOM     22  H2'  DT A   1     -27.913  30.293  16.600  1.00 84.80           H
ATOM     23 H2''  DT A   1     -29.292  30.806  15.584  1.00 81.56           H
ATOM     24  H1'  DT A   1     -29.682  32.451  17.051  1.00 78.27           H
ATOM     25  H3   DT A   1     -26.719  34.607  19.617  1.00 28.74           H
ATOM     26  H71  DT A   1     -25.780  30.778  22.086  1.00 50.18           H
ATOM     27  H72  DT A   1     -26.652  29.499  21.273  1.00 48.02           H
ATOM     28  H73  DT A   1     -25.102  30.028  20.658  1.00 45.93           H
ATOM     29  H6   DT A   1     -27.940  29.991  19.294  1.00 61.95           H
ATOM     30 HO5'  DT A   1     -32.955  30.331  19.367  1.00 93.75           H
ATOM     31  P    DC A   2     -29.642  27.048  15.520  1.00 47.77           P
ATOM     32  OP1  DC A   2     -29.822  26.059  16.634  1.00 59.41           O
ATOM     33  OP2  DC A   2     -28.302  27.352  14.910  1.00 45.59           O
ATOM     34  O5'  DC A   2     -30.586  26.654  14.265  1.00 27.88           O
ATOM     35  C5'  DC A   2     -30.832  27.623  13.249  1.00 20.14           C
ATOM     36  C4'  DC A   2     -32.073  27.255  12.487  1.00 13.85           C
ATOM     37  O4'  DC A   2     -33.193  27.169  13.377  1.00 15.88           O
ATOM     38  C3'  DC A   2     -32.430  28.259  11.406  1.00 13.35           C
ATOM     39  O3'  DC A   2     -31.925  27.742  10.159  1.00 14.20           O
ATOM     40  C2'  DC A   2     -33.937  28.349  11.445  1.00 14.27           C
ATOM     41  C1'  DC A   2     -34.339  27.805  12.827  1.00 14.48           C
ATOM     42  N1   DC A   2     -34.783  28.809  13.819  1.00 12.22           N
ATOM     43  C2   DC A   2     -35.998  28.635  14.444  1.00 11.84           C
ATOM     44  O2   DC A   2     -36.679  27.636  14.165  1.00 13.62           O
ATOM     45  N3   DC A   2     -36.422  29.564  15.346  1.00 12.86           N
ATOM     46  C4   DC A   2     -35.666  30.624  15.614  1.00 13.36           C
ATOM     47  N4   DC A   2     -36.094  31.476  16.526  1.00 14.20           N
ATOM     48  C5   DC A   2     -34.419  30.834  14.961  1.00 13.19           C
ATOM     49  C6   DC A   2     -34.038  29.926  14.064  1.00 13.40           C
ATOM     50  H5'  DC A   2     -29.981  27.660  12.566  1.00 21.31           H
ATOM     51 H5''  DC A   2     -30.958  28.615  13.686  1.00 20.99           H
ATOM     52  H4'  DC A   2     -31.902  26.279  12.014  1.00 14.55           H
ATOM     53  H3'  DC A   2     -31.994  29.238  11.648  1.00 14.07           H
ATOM     54  H2'  DC A   2     -34.247  29.388  11.343  1.00 14.31           H
ATOM     55 H2''  DC A   2     -34.378  27.741  10.656  1.00 14.39           H
ATOM     56  H1'  DC A   2     -35.131  27.067  12.668  1.00 14.11           H
ATOM     57  H41  DC A   2     -36.985  31.333  16.978  1.00 14.22           H
ATOM     58  H42  DC A   2     -35.550  32.300  16.741  1.00 13.91           H
ATOM     59  H5   DC A   2     -33.793  31.694  15.176  1.00 13.09           H
ATOM     60  H6   DC A   2     -33.083  30.064  13.567  1.00 13.08           H
ATOM     61  P    DG A   3     -31.936  28.584   8.801  1.00 14.73           P
ATOM     62  OP1  DG A   3     -30.959  27.886   7.930  1.00 18.11           O
ATOM     63  OP2  DG A   3     -31.811  30.056   9.107  1.00 17.37           O
ATOM     64  O5'  DG A   3     -33.403  28.383   8.256  1.00 13.71           O
ATOM     65  C5'  DG A   3     -33.806  27.149   7.721  1.00 14.25           C
ATOM     66  C4'  DG A   3     -35.301  27.101   7.518  1.00 14.02           C
ATOM     67  O4'  DG A   3     -35.915  27.169   8.821  1.00 13.42           O
ATOM     68  C3'  DG A   3     -35.912  28.247   6.739  1.00 11.88           C
ATOM     69  O3'  DG A   3     -35.873  27.902   5.343  1.00 12.75           O
ATOM     70  C2'  DG A   3     -37.316  28.343   7.307  1.00 13.13           C
ATOM     71  C1'  DG A   3     -37.129  27.931   8.759  1.00 12.95           C
ATOM     72  N9   DG A   3     -36.979  29.070   9.665  1.00 11.96           N
ATOM     73  C8   DG A   3     -36.070  30.076   9.550  1.00 12.32           C
ATOM     74  N7   DG A   3     -36.133  30.941  10.523  1.00 11.47           N
ATOM     75  C5   DG A   3     -37.163  30.475  11.339  1.00 10.74           C
ATOM     76  C6   DG A   3     -37.705  30.997  12.534  1.00 10.48           C
ATOM     77  O6   DG A   3     -37.359  31.998  13.163  1.00 12.00           O
ATOM     78  N1   DG A   3     -38.779  30.242  12.984  1.00 11.07           N
ATOM     79  C2   DG A   3     -39.269  29.134  12.373  1.00 11.34           C
ATOM     80  N2   DG A   3     -40.300  28.553  12.949  1.00 12.77           N
ATOM     81  N3   DG A   3     -38.726  28.587  11.294  1.00 11.91           N
ATOM     82  C4   DG A   3     -37.695  29.317  10.822  1.00 11.19           C
ATOM     83  H5'  DG A   3     -33.506  26.345   8.396  1.00 14.44           H
ATOM     84 H5''  DG A   3     -33.306  26.988   6.764  1.00 12.69           H
ATOM     85  H4'  DG A   3     -35.563  26.153   7.027  1.00 13.88           H
ATOM     86  H3'  DG A   3     -35.366  29.177   6.938  1.00 12.05           H
ATOM     87  H2'  DG A   3     -37.696  29.364   7.241  1.00 13.37           H
ATOM     88 H2''  DG A   3     -37.990  27.654   6.793  1.00 12.65           H
ATOM     89  H1'  DG A   3     -37.982  27.308   9.062  1.00 12.43           H
ATOM     90  H8   DG A   3     -35.318  30.101   8.777  1.00 11.96           H
ATOM     91  H1   DG A   3     -39.244  30.570  13.853  1.00 10.57           H
ATOM     92  H21  DG A   3     -40.671  27.696  12.564  1.00 12.52           H
ATOM     93  H22  DG A   3     -40.682  28.919  13.808  1.00 12.23           H
ATOM    283  P    DA A  10     -25.244  44.397  -9.360  1.00 20.54           P
ATOM    284  OP1  DA A  10     -24.552  45.551  -9.942  1.00 29.48           O
ATOM    285  OP2  DA A  10     -26.370  43.671 -10.101  1.00 23.35           O
ATOM    286  O5'  DA A  10     -24.175  43.292  -8.979  1.00 20.55           O
ATOM    287  C5'  DA A  10     -23.000  43.584  -8.216  1.00 17.24           C
ATOM    288  C4'  DA A  10     -22.255  42.296  -7.935  1.00 18.54           C
ATOM    289  O4'  DA A  10     -23.137  41.343  -7.300  1.00 18.14           O
ATOM    290  C3'  DA A  10     -21.729  41.564  -9.158  1.00 16.43           C
ATOM    291  O3'  DA A  10     -20.487  42.076  -9.563  1.00 19.10           O
ATOM    292  C2'  DA A  10     -21.592  40.152  -8.649  1.00 17.51           C
ATOM    293  C1'  DA A  10     -22.786  40.023  -7.750  1.00 17.31           C
ATOM    294  N9   DA A  10     -23.917  39.473  -8.490  1.00 18.85           N
ATOM    295  C8   DA A  10     -24.650  40.048  -9.504  1.00 18.32           C
ATOM    296  N7   DA A  10     -25.646  39.314  -9.905  1.00 18.71           N
ATOM    297  C5   DA A  10     -25.485  38.125  -9.209  1.00 17.10           C
ATOM    298  C6   DA A  10     -26.202  36.927  -9.222  1.00 17.26           C
ATOM    299  N6   DA A  10     -27.256  36.711 -10.018  1.00 19.50           N
ATOM    300  N1   DA A  10     -25.863  35.988  -8.317  1.00 17.68           N
ATOM    301  C2   DA A  10     -24.804  36.206  -7.532  1.00 16.95           C
ATOM    302  N3   DA A  10     -24.039  37.281  -7.445  1.00 16.23           N
ATOM    303  C4   DA A  10     -24.456  38.228  -8.300  1.00 16.52           C
ATOM    304  H5'  DA A  10     -23.280  44.063  -7.275  1.00 17.63           H
ATOM    305 H5''  DA A  10     -22.358  44.270  -8.772  1.00 16.90           H
ATOM    306  H4'  DA A  10     -21.411  42.519  -7.268  1.00 18.75           H
ATOM    307  H3'  DA A  10     -22.450  41.620  -9.983  1.00 17.81           H
ATOM    308 HO3'  DA A  10     -19.996  41.652 -10.293  0.00 18.85           H
ATOM    309  H2'  DA A  10     -21.633  39.427  -9.461  1.00 17.10           H
ATOM    310 H2''  DA A  10     -20.666  40.033  -8.083  1.00 17.47           H
ATOM    311  H1'  DA A  10     -22.526  39.387  -6.894  1.00 17.04           H
ATOM    312  H8   DA A  10     -24.478  41.051  -9.866  1.00 18.12           H
ATOM    313  H61  DA A  10     -27.753  35.835  -9.976  1.00 19.52           H
ATOM    314  H62  DA A  10     -27.549  37.425 -10.670  1.00 20.35           H
ATOM    315  H2   DA A  10     -24.549  35.393  -6.866  1.00 18.44           H
TER     316       DA A  10
"""
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)

def exercise_8():
  """
  non-standard ligands
  """
  in_pdb_str = """\
HETATM  317 RU   RKL A 101     -45.719  30.968   8.877  1.00 11.71          RU
HETATM  318  C1  RKL A 101     -48.034  30.452  10.475  1.00 12.70           C
HETATM  319  N1  RKL A 101     -45.981  29.249   9.992  1.00 12.74           N
HETATM  320  C2  RKL A 101     -48.390  32.534   9.306  1.00 13.81           C
HETATM  321  N2  RKL A 101     -47.647  31.370   9.494  1.00 12.55           N
HETATM  322  C3  RKL A 101     -49.604  32.642   9.975  1.00 15.19           C
HETATM  323  N3  RKL A 101     -50.738  29.726  12.856  1.00 17.18           N
HETATM  324  C4  RKL A 101     -50.018  31.694  10.923  1.00 14.13           C
HETATM  325  N4  RKL A 101     -49.115  27.508  13.253  1.00 15.38           N
HETATM  326  C5  RKL A 101     -49.247  30.585  11.175  1.00 13.25           C
HETATM  327  N5  RKL A 101     -44.888  31.813  10.555  1.00 11.82           N
HETATM  328  C6  RKL A 101     -49.612  29.548  12.140  1.00 13.75           C
HETATM  329  N6  RKL A 101     -43.259  32.778  12.684  1.00 12.22           N
HETATM  330  C7  RKL A 101     -48.821  28.435  12.352  1.00 14.81           C
HETATM  331  N7  RKL A 101     -40.923  30.190   8.359  1.00 13.01           N
HETATM  332  C8  RKL A 101     -47.602  28.226  11.544  1.00 12.98           C
HETATM  333  N8  RKL A 101     -43.748  30.550   8.401  1.00 11.93           N
HETATM  334  C9  RKL A 101     -46.771  27.130  11.684  1.00 15.35           C
HETATM  335  N9  RKL A 101     -46.411  30.168   7.131  1.00 12.05           N
HETATM  336  C10 RKL A 101     -47.205  29.296  10.658  1.00 12.10           C
HETATM  337  N10 RKL A 101     -47.092  29.427   4.483  1.00 15.48           N
HETATM  338  C11 RKL A 101     -45.601  27.081  10.913  1.00 15.50           C
HETATM  339  N11 RKL A 101     -45.706  34.628   5.712  1.00 14.17           N
HETATM  340  C12 RKL A 101     -45.267  28.146  10.066  1.00 13.34           C
HETATM  341  N12 RKL A 101     -45.589  32.653   7.699  1.00 11.42           N
HETATM  342  C13 RKL A 101     -50.267  27.698  13.961  1.00 16.97           C
HETATM  343  C14 RKL A 101     -50.609  26.747  14.957  1.00 19.35           C
HETATM  344  C15 RKL A 101     -51.059  28.835  13.799  1.00 16.91           C
HETATM  345  C16 RKL A 101     -52.245  29.002  14.481  1.00 20.28           C
HETATM  346  C17 RKL A 101     -52.569  28.079  15.430  1.00 23.24           C
HETATM  347  C18 RKL A 101     -51.794  26.934  15.688  1.00 21.56           C
HETATM  348  C19 RKL A 101     -43.507  31.662  10.584  1.00 10.71           C
HETATM  349  C20 RKL A 101     -45.376  32.439  11.617  1.00 11.87           C
HETATM  350  C21 RKL A 101     -44.567  32.925  12.651  1.00 12.51           C
HETATM  351  C22 RKL A 101     -42.684  32.140  11.626  1.00 11.32           C
HETATM  352  C23 RKL A 101     -41.303  31.988  11.513  1.00 12.15           C
HETATM  353  C24 RKL A 101     -40.725  31.347  10.488  1.00 13.05           C
HETATM  354  C25 RKL A 101     -41.510  30.844   9.435  1.00 11.47           C
HETATM  355  C26 RKL A 101     -42.921  30.994   9.454  1.00 11.14           C
HETATM  356  C27 RKL A 101     -41.729  29.778   7.405  1.00 13.49           C
HETATM  357  C28 RKL A 101     -43.083  29.942   7.445  1.00 13.01           C
HETATM  358  C29 RKL A 101     -46.319  31.052   6.034  1.00 12.42           C
HETATM  359  C30 RKL A 101     -46.800  28.962   6.755  1.00 13.67           C
HETATM  360  C31 RKL A 101     -47.105  28.610   5.449  1.00 13.99           C
HETATM  361  C32 RKL A 101     -46.768  30.756   4.742  1.00 12.91           C
HETATM  362  C33 RKL A 101     -46.730  31.799   3.746  1.00 14.57           C
HETATM  363  C34 RKL A 101     -46.337  33.061   4.020  1.00 14.44           C
HETATM  364  C35 RKL A 101     -45.965  33.375   5.400  1.00 12.98           C
HETATM  365  C36 RKL A 101     -45.993  32.378   6.378  1.00 12.30           C
HETATM  366  C37 RKL A 101     -45.377  34.905   6.995  1.00 13.47           C
HETATM  367  C38 RKL A 101     -45.343  33.911   7.972  1.00 12.15           C
HETATM  368  H2  RKL A 101     -48.084  33.276   8.582  1.00 13.86           H
HETATM  369  H3  RKL A 101     -50.188  33.545   9.851  1.00 15.11           H
HETATM  370  H4  RKL A 101     -50.905  31.887  11.513  1.00 13.58           H
HETATM  371  H9  RKL A 101     -46.944  26.390  12.455  1.00 15.01           H
HETATM  372  H11 RKL A 101     -44.895  26.277  11.066  1.00 15.11           H
HETATM  373  H12 RKL A 101     -44.365  28.076   9.472  1.00 13.69           H
HETATM  374  H14 RKL A 101     -50.016  25.860  15.060  1.00 21.31           H
HETATM  375  H16 RKL A 101     -52.865  29.872  14.310  1.00 20.48           H
HETATM  376  H17 RKL A 101     -53.510  28.186  15.956  1.00 21.59           H
HETATM  377  H18 RKL A 101     -52.149  26.171  16.370  1.00 20.40           H
HETATM  378  H20 RKL A 101     -46.445  32.574  11.671  1.00 11.96           H
HETATM  379  H21 RKL A 101     -45.048  33.419  13.484  1.00 12.20           H
HETATM  380  H23 RKL A 101     -40.679  32.343  12.322  1.00 11.86           H
HETATM  381  H24 RKL A 101     -39.653  31.221  10.464  1.00 12.19           H
HETATM  382  H27 RKL A 101     -41.296  29.269   6.554  1.00 12.47           H
HETATM  383  H28 RKL A 101     -43.662  29.545   6.621  1.00 13.39           H
HETATM  384  H30 RKL A 101     -46.863  28.209   7.519  1.00 14.54           H
HETATM  385  H31 RKL A 101     -47.380  27.588   5.252  1.00 13.70           H
HETATM  386  H33 RKL A 101     -47.023  31.556   2.739  1.00 13.71           H
HETATM  387  H34 RKL A 101     -46.334  33.821   3.256  1.00 13.85           H
HETATM  388  H37 RKL A 101     -45.110  35.918   7.266  1.00 13.33           H
HETATM  389  H38 RKL A 101     -45.080  34.177   8.985  1.00 12.40           H
  """
  out_pdb_str = """\
HETATM    1  C1  RKL A 101     -48.034  30.452  10.475  1.00 12.70           C
HETATM    2  C10 RKL A 101     -47.205  29.296  10.658  1.00 12.10           C
HETATM    3  C11 RKL A 101     -45.601  27.081  10.913  1.00 15.50           C
HETATM    4  C12 RKL A 101     -45.267  28.146  10.066  1.00 13.34           C
HETATM    5  C13 RKL A 101     -50.267  27.698  13.961  1.00 16.97           C
HETATM    6  C14 RKL A 101     -50.609  26.747  14.957  1.00 19.35           C
HETATM    7  C15 RKL A 101     -51.059  28.835  13.799  1.00 16.91           C
HETATM    8  C16 RKL A 101     -52.245  29.002  14.481  1.00 20.28           C
HETATM    9  C17 RKL A 101     -52.569  28.079  15.430  1.00 23.24           C
HETATM   10  C18 RKL A 101     -51.794  26.934  15.688  1.00 21.56           C
HETATM   11  C19 RKL A 101     -43.507  31.662  10.584  1.00 10.71           C
HETATM   12  C2  RKL A 101     -48.390  32.534   9.306  1.00 13.81           C
HETATM   13  C20 RKL A 101     -45.376  32.439  11.617  1.00 11.87           C
HETATM   14  C21 RKL A 101     -44.567  32.925  12.651  1.00 12.51           C
HETATM   15  C22 RKL A 101     -42.684  32.140  11.626  1.00 11.32           C
HETATM   16  C23 RKL A 101     -41.303  31.988  11.513  1.00 12.15           C
HETATM   17  C24 RKL A 101     -40.725  31.347  10.488  1.00 13.05           C
HETATM   18  C25 RKL A 101     -41.510  30.844   9.435  1.00 11.47           C
HETATM   19  C26 RKL A 101     -42.921  30.994   9.454  1.00 11.14           C
HETATM   20  C27 RKL A 101     -41.729  29.778   7.405  1.00 13.49           C
HETATM   21  C28 RKL A 101     -43.083  29.942   7.445  1.00 13.01           C
HETATM   22  C29 RKL A 101     -46.319  31.052   6.034  1.00 12.42           C
HETATM   23  C3  RKL A 101     -49.604  32.642   9.975  1.00 15.19           C
HETATM   24  C30 RKL A 101     -46.800  28.962   6.755  1.00 13.67           C
HETATM   25  C31 RKL A 101     -47.105  28.610   5.449  1.00 13.99           C
HETATM   26  C32 RKL A 101     -46.768  30.756   4.742  1.00 12.91           C
HETATM   27  C33 RKL A 101     -46.730  31.799   3.746  1.00 14.57           C
HETATM   28  C34 RKL A 101     -46.337  33.061   4.020  1.00 14.44           C
HETATM   29  C35 RKL A 101     -45.965  33.375   5.400  1.00 12.98           C
HETATM   30  C36 RKL A 101     -45.993  32.378   6.378  1.00 12.30           C
HETATM   31  C37 RKL A 101     -45.377  34.905   6.995  1.00 13.47           C
HETATM   32  C38 RKL A 101     -45.343  33.911   7.972  1.00 12.15           C
HETATM   33  C4  RKL A 101     -50.018  31.694  10.923  1.00 14.13           C
HETATM   34  C5  RKL A 101     -49.247  30.585  11.175  1.00 13.25           C
HETATM   35  C6  RKL A 101     -49.612  29.548  12.140  1.00 13.75           C
HETATM   36  C7  RKL A 101     -48.821  28.435  12.352  1.00 14.81           C
HETATM   37  C8  RKL A 101     -47.602  28.226  11.544  1.00 12.98           C
HETATM   38  C9  RKL A 101     -46.771  27.130  11.684  1.00 15.35           C
HETATM   39  N1  RKL A 101     -45.981  29.249   9.992  1.00 12.74           N
HETATM   40  N10 RKL A 101     -47.092  29.427   4.483  1.00 15.48           N
HETATM   41  N11 RKL A 101     -45.706  34.628   5.712  1.00 14.17           N
HETATM   42  N12 RKL A 101     -45.589  32.653   7.699  1.00 11.42           N
HETATM   43  N2  RKL A 101     -47.647  31.370   9.494  1.00 12.55           N
HETATM   44  N3  RKL A 101     -50.738  29.726  12.856  1.00 17.18           N
HETATM   45  N4  RKL A 101     -49.115  27.508  13.253  1.00 15.38           N
HETATM   46  N5  RKL A 101     -44.888  31.813  10.555  1.00 11.82           N
HETATM   47  N6  RKL A 101     -43.259  32.778  12.684  1.00 12.22           N
HETATM   48  N7  RKL A 101     -40.923  30.190   8.359  1.00 13.01           N
HETATM   49  N8  RKL A 101     -43.748  30.550   8.401  1.00 11.93           N
HETATM   50  N9  RKL A 101     -46.411  30.168   7.131  1.00 12.05           N
HETATM   51 RU   RKL A 101     -45.719  30.968   8.877  1.00 11.71          RU
HETATM   52  H2  RKL A 101     -48.084  33.276   8.582  1.00 13.86           H
HETATM   53  H3  RKL A 101     -50.188  33.545   9.851  1.00 15.11           H
HETATM   54  H11 RKL A 101     -44.895  26.277  11.066  1.00 15.11           H
HETATM   55  H12 RKL A 101     -44.365  28.076   9.472  1.00 13.69           H
HETATM   56  H14 RKL A 101     -50.016  25.860  15.060  1.00 21.31           H
HETATM   57  H16 RKL A 101     -52.865  29.872  14.310  1.00 20.48           H
HETATM   58  H17 RKL A 101     -53.510  28.186  15.956  1.00 21.59           H
HETATM   59  H18 RKL A 101     -52.149  26.171  16.370  1.00 20.40           H
HETATM   60  H20 RKL A 101     -46.445  32.574  11.671  1.00 11.96           H
HETATM   61  H21 RKL A 101     -45.048  33.419  13.484  1.00 12.20           H
HETATM   62  H23 RKL A 101     -40.679  32.343  12.322  1.00 11.86           H
HETATM   63  H24 RKL A 101     -39.653  31.221  10.464  1.00 12.19           H
HETATM   64  H27 RKL A 101     -41.296  29.269   6.554  1.00 12.47           H
HETATM   65  H28 RKL A 101     -43.662  29.545   6.621  1.00 13.39           H
HETATM   66  H30 RKL A 101     -46.863  28.209   7.519  1.00 14.54           H
HETATM   67  H31 RKL A 101     -47.380  27.588   5.252  1.00 13.70           H
HETATM   68  H33 RKL A 101     -47.023  31.556   2.739  1.00 13.71           H
HETATM   69  H34 RKL A 101     -46.334  33.821   3.256  1.00 13.85           H
HETATM   70  H37 RKL A 101     -45.110  35.918   7.266  1.00 13.33           H
HETATM   71  H38 RKL A 101     -45.080  34.177   8.985  1.00 12.40           H
HETATM   72  H4  RKL A 101     -50.905  31.887  11.513  1.00 13.58           H
HETATM   73  H9  RKL A 101     -46.944  26.390  12.455  1.00 15.01           H
TER
  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=out_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)

def exercise_9():
  """
  HETATM+strange residue+altloc
  """
  in_pdb_str = """\
HETATM  177  N   DLE A  10      19.161  15.952  33.645  1.00  6.08           N
HETATM  178  CA  DLE A  10      19.055  14.520  33.727  1.00  8.34           C
HETATM  179  CB ADLE A  10      17.599  14.234  34.534  0.40  7.56           C
HETATM  180  CB BDLE A  10      17.810  13.938  33.307  0.30  5.95           C
HETATM  181  CB CDLE A  10      17.574  14.136  33.980  0.30 11.05           C
HETATM  182  CG ADLE A  10      16.464  14.576  33.593  0.40  8.13           C
HETATM  183  CG BDLE A  10      16.719  14.134  34.392  0.30  8.43           C
HETATM  184  CG CDLE A  10      16.476  14.206  32.929  0.30 12.34           C
HETATM  185  CD1ADLE A  10      16.341  13.585  32.446  0.40 13.59           C
HETATM  186  CD1BDLE A  10      16.718  13.157  35.542  0.30  9.56           C
HETATM  187  CD1CDLE A  10      16.348  12.835  32.290  0.30 15.22           C
HETATM  188  CD2ADLE A  10      15.102  14.630  34.276  0.40 11.11           C
HETATM  189  CD2BDLE A  10      16.647  15.478  35.059  0.30 13.57           C
HETATM  190  CD2CDLE A  10      15.254  14.620  33.685  0.30 12.74           C
HETATM  191  C   DLE A  10      19.895  13.952  34.844  1.00  6.92           C
HETATM  192  O   DLE A  10      19.946  14.484  35.950  1.00  9.33           O
HETATM  193  H   DLE A  10      19.184  16.418  34.367  1.00  7.29           H
HETATM  194  HA  DLE A  10      19.146  14.053  32.870  1.00 10.00           H
HETATM  195  HB1ADLE A  10      17.543  13.303  34.800  0.40  9.07           H
HETATM  196  HB1BDLE A  10      17.582  13.221  34.303  0.30 13.26           H
HETATM  197  HB1CDLE A  10      17.934  12.991  33.141  0.30  7.14           H
HETATM  198  HB2ADLE A  10      17.547  14.784  35.332  0.40  9.07           H
HETATM  199  HG ADLE A  10      16.645  15.459  33.210  0.40  9.75           H
HETATM  200  HG BDLE A  10      17.274  14.684  34.722  0.30 13.26           H
HETATM  201  HG CDLE A  10      17.521  14.354  32.479  0.30  7.14           H
HETATM  202 HD11ADLE A  10      15.596  13.829  31.891  0.40 20.38           H
HETATM  203 HD11BDLE A  10      16.681  12.261  35.201  0.30 14.35           H
HETATM  204 HD11CDLE A  10      15.891  12.243  32.892  0.30 22.83           H
HETATM  205 HD12ADLE A  10      16.202  12.703  32.798  0.40 20.38           H
HETATM  206 HD12BDLE A  10      15.953  13.319  36.099  0.30 14.35           H
HETATM  207 HD12CDLE A  10      15.850  12.908  31.473  0.30 22.83           H
HETATM  208 HD13ADLE A  10      17.146  13.599  31.924  0.40 20.38           H
HETATM  209 HD13BDLE A  10      17.519  13.270  36.058  0.30 14.35           H
HETATM  210 HD13CDLE A  10      17.223  12.487  32.102  0.30 22.83           H
HETATM  211 HD21ADLE A  10      14.426  14.832  33.625  0.40 16.66           H
HETATM  212 HD21BDLE A  10      17.384  15.574  35.666  0.30 20.36           H
HETATM  213 HD21CDLE A  10      15.860  14.014  33.935  0.30 10.12           H
HETATM  214 HD22ADLE A  10      15.109  15.312  34.951  0.40 16.66           H
HETATM  215 HD22BDLE A  10      16.689  16.168  34.393  0.30 20.36           H
HETATM  216 HD22CDLE A  10      14.534  14.774  33.069  0.30 19.10           H
HETATM  217 HD23ADLE A  10      14.912  13.781  34.681  0.40 16.66           H
HETATM  218 HD23BDLE A  10      15.821  15.551  35.543  0.30 20.36           H
HETATM  219 HD23CDLE A  10      15.007  13.925  34.299  0.30 19.10           H
  """
  out_pdb_str = """\
HETATM    1  N   DLE A  10      19.161  15.952  33.645  1.00  6.08           N
HETATM    2  CA  DLE A  10      19.055  14.520  33.727  1.00  8.34           C
HETATM    3  C   DLE A  10      19.895  13.952  34.844  1.00  6.92           C
HETATM    4  O   DLE A  10      19.946  14.484  35.950  1.00  9.33           O
HETATM    5  H   DLE A  10      19.184  16.418  34.367  1.00  7.29           H
HETATM    6  HA  DLE A  10      19.146  14.053  32.870  1.00 10.00           H
HETATM    7  CB ADLE A  10      17.599  14.234  34.534  0.40  7.56           C
HETATM    8  CG ADLE A  10      16.464  14.576  33.593  0.40  8.13           C
HETATM    9  CD1ADLE A  10      16.341  13.585  32.446  0.40 13.59           C
HETATM   10  CD2ADLE A  10      15.102  14.630  34.276  0.40 11.11           C
HETATM   11  HB1ADLE A  10      17.543  13.303  34.800  0.40  9.07           H
HETATM   12  HB2ADLE A  10      17.547  14.784  35.332  0.40  9.07           H
HETATM   13  HG ADLE A  10      16.645  15.459  33.210  0.40  9.75           H
HETATM   14 HD11ADLE A  10      15.596  13.829  31.891  0.40 20.38           H
HETATM   15 HD12ADLE A  10      16.202  12.703  32.798  0.40 20.38           H
HETATM   16 HD13ADLE A  10      17.146  13.599  31.924  0.40 20.38           H
HETATM   17 HD21ADLE A  10      14.426  14.832  33.625  0.40 16.66           H
HETATM   18 HD22ADLE A  10      15.109  15.312  34.951  0.40 16.66           H
HETATM   19 HD23ADLE A  10      14.912  13.781  34.681  0.40 16.66           H
HETATM   20  CB BDLE A  10      17.810  13.938  33.307  0.30  5.95           C
HETATM   21  CG BDLE A  10      16.719  14.134  34.392  0.30  8.43           C
HETATM   22  CD1BDLE A  10      16.718  13.157  35.542  0.30  9.56           C
HETATM   23  CD2BDLE A  10      16.647  15.478  35.059  0.30 13.57           C
HETATM   24  HB1BDLE A  10      17.582  13.221  34.303  0.30 13.26           H
HETATM   25  HG BDLE A  10      17.274  14.684  34.722  0.30 13.26           H
HETATM   26 HD11BDLE A  10      16.681  12.261  35.201  0.30 14.35           H
HETATM   27 HD12BDLE A  10      15.953  13.319  36.099  0.30 14.35           H
HETATM   28 HD13BDLE A  10      17.519  13.270  36.058  0.30 14.35           H
HETATM   29 HD21BDLE A  10      17.384  15.574  35.666  0.30 20.36           H
HETATM   30 HD22BDLE A  10      16.689  16.168  34.393  0.30 20.36           H
HETATM   31 HD23BDLE A  10      15.821  15.551  35.543  0.30 20.36           H
HETATM   32  CB CDLE A  10      17.574  14.136  33.980  0.30 11.05           C
HETATM   33  CG CDLE A  10      16.476  14.206  32.929  0.30 12.34           C
HETATM   34  CD1CDLE A  10      16.348  12.835  32.290  0.30 15.22           C
HETATM   35  CD2CDLE A  10      15.254  14.620  33.685  0.30 12.74           C
HETATM   36  HB1CDLE A  10      17.934  12.991  33.141  0.30  7.14           H
HETATM   37  HG CDLE A  10      17.521  14.354  32.479  0.30  7.14           H
HETATM   38 HD11CDLE A  10      15.891  12.243  32.892  0.30 22.83           H
HETATM   39 HD12CDLE A  10      15.850  12.908  31.473  0.30 22.83           H
HETATM   40 HD13CDLE A  10      17.223  12.487  32.102  0.30 22.83           H
HETATM   41 HD21CDLE A  10      15.860  14.014  33.935  0.30 10.12           H
HETATM   42 HD22CDLE A  10      14.534  14.774  33.069  0.30 19.10           H
HETATM   43 HD23CDLE A  10      15.007  13.925  34.299  0.30 19.10           H  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=out_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)

def exercise_10():
  """
  nucleic acids: RNA
  """
  in_pdb_str = """\
HEADER    RNA                                     26-MAR-97   1MIS
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000
ATOM      1  O5'   G A   1     -39.305 107.866 -51.789  1.00  0.00           O
ATOM      2  C5'   G A   1     -38.172 107.974 -50.957  1.00  0.00           C
ATOM      3  C4'   G A   1     -37.351 109.206 -51.345  1.00  0.00           C
ATOM      4  O4'   G A   1     -36.833 109.097 -52.674  1.00  0.00           O
ATOM      5  C3'   G A   1     -38.169 110.487 -51.328  1.00  0.00           C
ATOM      6  O3'   G A   1     -38.370 110.983 -50.015  1.00  0.00           O
ATOM      7  C2'   G A   1     -37.260 111.365 -52.146  1.00  0.00           C
ATOM      8  O2'   G A   1     -36.123 111.712 -51.383  1.00  0.00           O
ATOM      9  C1'   G A   1     -36.827 110.411 -53.252  1.00  0.00           C
ATOM     10  N9    G A   1     -37.778 110.451 -54.399  1.00  0.00           N
ATOM     11  C8    G A   1     -38.788 109.567 -54.744  1.00  0.00           C
ATOM     12  N7    G A   1     -39.381 109.860 -55.867  1.00  0.00           N
ATOM     13  C5    G A   1     -38.736 111.011 -56.303  1.00  0.00           C
ATOM     14  C6    G A   1     -38.952 111.785 -57.481  1.00  0.00           C
ATOM     15  O6    G A   1     -39.762 111.589 -58.384  1.00  0.00           O
ATOM     16  N1    G A   1     -38.098 112.882 -57.536  1.00  0.00           N
ATOM     17  C2    G A   1     -37.146 113.193 -56.587  1.00  0.00           C
ATOM     18  N2    G A   1     -36.416 114.294 -56.811  1.00  0.00           N
ATOM     19  N3    G A   1     -36.930 112.460 -55.487  1.00  0.00           N
ATOM     20  C4    G A   1     -37.761 111.387 -55.408  1.00  0.00           C
ATOM     21  H5'   G A   1     -37.557 107.080 -51.058  1.00  0.00           H
ATOM     22 H5''   G A   1     -38.499 108.071 -49.921  1.00  0.00           H
ATOM     23  H4'   G A   1     -36.512 109.317 -50.658  1.00  0.00           H
ATOM     24  H3'   G A   1     -39.114 110.342 -51.848  1.00  0.00           H
ATOM     25  H2'   G A   1     -37.743 112.247 -52.531  1.00  0.00           H
ATOM     26 HO2'   G A   1     -36.410 112.213 -50.618  1.00  0.00           H
ATOM     27  H1'   G A   1     -35.830 110.704 -53.596  1.00  0.00           H
ATOM     28  H8    G A   1     -39.075 108.696 -54.159  1.00  0.00           H
ATOM     29  H1    G A   1     -38.189 113.490 -58.338  1.00  0.00           H
ATOM     30  H21   G A   1     -36.579 114.849 -57.640  1.00  0.00           H
ATOM     31  H22   G A   1     -35.704 114.571 -56.150  1.00  0.00           H
ATOM     32 HO5'   G A   1     -39.804 107.089 -51.525  1.00  0.00           H
ATOM     33  P     C A   2     -39.575 112.019 -49.689  1.00  0.00           P
ATOM     34  OP1   C A   2     -39.502 112.375 -48.255  1.00  0.00           O
ATOM     35  OP2   C A   2     -40.825 111.449 -50.239  1.00  0.00           O
ATOM     36  O5'   C A   2     -39.205 113.330 -50.550  1.00  0.00           O
ATOM     37  C5'   C A   2     -38.206 114.229 -50.120  1.00  0.00           C
ATOM     38  C4'   C A   2     -38.024 115.328 -51.169  1.00  0.00           C
ATOM     39  O4'   C A   2     -37.635 114.785 -52.430  1.00  0.00           O
ATOM     40  C3'   C A   2     -39.302 116.088 -51.460  1.00  0.00           C
ATOM     41  O3'   C A   2     -39.648 117.004 -50.437  1.00  0.00           O
ATOM     42  C2'   C A   2     -38.900 116.766 -52.741  1.00  0.00           C
ATOM     43  O2'   C A   2     -37.943 117.772 -52.478  1.00  0.00           O
ATOM     44  C1'   C A   2     -38.214 115.608 -53.457  1.00  0.00           C
ATOM     45  N1    C A   2     -39.213 114.844 -54.247  1.00  0.00           N
ATOM     46  C2    C A   2     -39.629 115.388 -55.459  1.00  0.00           C
ATOM     47  O2    C A   2     -39.204 116.482 -55.828  1.00  0.00           O
ATOM     48  N3    C A   2     -40.513 114.676 -56.213  1.00  0.00           N
ATOM     49  C4    C A   2     -40.964 113.480 -55.805  1.00  0.00           C
ATOM     50  N4    C A   2     -41.837 112.816 -56.576  1.00  0.00           N
ATOM     51  C5    C A   2     -40.527 112.907 -54.565  1.00  0.00           C
ATOM     52  C6    C A   2     -39.665 113.627 -53.828  1.00  0.00           C
ATOM     53  H5'   C A   2     -37.265 113.700 -49.985  1.00  0.00           H
ATOM     54 H5''   C A   2     -38.505 114.679 -49.173  1.00  0.00           H
ATOM     55  H4'   C A   2     -37.260 116.032 -50.837  1.00  0.00           H
ATOM     56  H3'   C A   2     -40.104 115.385 -51.660  1.00  0.00           H
ATOM     57  H2'   C A   2     -39.744 117.164 -53.296  1.00  0.00           H
ATOM     58 HO2'   C A   2     -38.340 118.429 -51.904  1.00  0.00           H
ATOM     59  H1'   C A   2     -37.448 115.962 -54.119  1.00  0.00           H
ATOM     60  H41   C A   2     -42.140 113.214 -57.454  1.00  0.00           H
ATOM     61  H42   C A   2     -42.186 111.915 -56.282  1.00  0.00           H
ATOM     62  H5    C A   2     -40.825 111.939 -54.208  1.00  0.00           H
ATOM     63  H6    C A   2     -39.333 113.231 -52.894  1.00  0.00           H
ATOM    132  P     A A   5     -49.440 124.595 -55.324  1.00  0.00           P
ATOM    133  OP1   A A   5     -49.855 126.011 -55.204  1.00  0.00           O
ATOM    134  OP2   A A   5     -49.735 123.644 -54.229  1.00  0.00           O
ATOM    135  O5'   A A   5     -50.048 124.002 -56.693  1.00  0.00           O
ATOM    136  C5'   A A   5     -51.166 123.141 -56.672  1.00  0.00           C
ATOM    137  C4'   A A   5     -51.363 122.530 -58.060  1.00  0.00           C
ATOM    138  O4'   A A   5     -50.271 121.673 -58.404  1.00  0.00           O
ATOM    139  C3'   A A   5     -52.607 121.669 -58.101  1.00  0.00           C
ATOM    140  O3'   A A   5     -53.761 122.442 -58.379  1.00  0.00           O
ATOM    141  C2'   A A   5     -52.264 120.728 -59.232  1.00  0.00           C
ATOM    142  O2'   A A   5     -52.381 121.398 -60.470  1.00  0.00           O
ATOM    143  C1'   A A   5     -50.787 120.462 -58.947  1.00  0.00           C
ATOM    144  N9    A A   5     -50.642 119.381 -57.952  1.00  0.00           N
ATOM    145  C8    A A   5     -50.039 119.417 -56.717  1.00  0.00           C
ATOM    146  N7    A A   5     -50.053 118.272 -56.091  1.00  0.00           N
ATOM    147  C5    A A   5     -50.688 117.412 -56.982  1.00  0.00           C
ATOM    148  C6    A A   5     -50.989 116.037 -56.932  1.00  0.00           C
ATOM    149  N6    A A   5     -50.698 115.257 -55.883  1.00  0.00           N
ATOM    150  N1    A A   5     -51.588 115.490 -58.008  1.00  0.00           N
ATOM    151  C2    A A   5     -51.867 116.259 -59.057  1.00  0.00           C
ATOM    152  N3    A A   5     -51.638 117.555 -59.224  1.00  0.00           N
ATOM    153  C4    A A   5     -51.036 118.077 -58.126  1.00  0.00           C
ATOM    154  H5'   A A   5     -52.054 123.706 -56.387  1.00  0.00           H
ATOM    155 H5''   A A   5     -51.000 122.338 -55.954  1.00  0.00           H
ATOM    156  H4'   A A   5     -51.445 123.311 -58.809  1.00  0.00           H
ATOM    157  H3'   A A   5     -52.700 121.122 -57.163  1.00  0.00           H
ATOM    158  H2'   A A   5     -52.868 119.824 -59.228  1.00  0.00           H
ATOM    159 HO2'   A A   5     -51.824 122.180 -60.448  1.00  0.00           H
ATOM    160  H1'   A A   5     -50.247 120.181 -59.845  1.00  0.00           H
ATOM    161  H8    A A   5     -49.574 120.307 -56.306  1.00  0.00           H
ATOM    162  H61   A A   5     -50.933 114.274 -55.905  1.00  0.00           H
ATOM    163  H62   A A   5     -50.244 115.650 -55.071  1.00  0.00           H
ATOM    164  H2    A A   5     -52.341 115.756 -59.891  1.00  0.00           H
TER     262        C A   8
"""
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)

def exercise_11():
  """
  neutron structure
  """
  in_pdb_str = """\
ATOM    257  N   TRP A  16      36.612  25.839  63.233  1.00 11.33           N
ATOM    258  CA  TRP A  16      37.749  26.091  62.353  1.00 12.58           C
ATOM    259  C   TRP A  16      38.792  26.998  62.988  1.00 15.16           C
ATOM    260  O   TRP A  16      39.948  27.034  62.500  1.00 16.86           O
ATOM    261  CB  TRP A  16      37.252  26.708  61.028  1.00 12.69           C
ATOM    262  CG  TRP A  16      36.696  28.080  61.195  1.00 18.76           C
ATOM    263  CD1 TRP A  16      37.392  29.288  61.087  1.00 19.25           C
ATOM    264  CD2 TRP A  16      35.350  28.447  61.577  1.00 17.43           C
ATOM    265  NE1 TRP A  16      36.537  30.385  61.341  1.00 19.14           N
ATOM    266  CE2 TRP A  16      35.280  29.860  61.634  1.00 18.98           C
ATOM    267  CE3 TRP A  16      34.222  27.670  61.850  1.00 12.03           C
ATOM    268  CZ2 TRP A  16      34.127  30.548  61.929  1.00 14.01           C
ATOM    269  CZ3 TRP A  16      33.065  28.394  62.130  1.00 11.86           C
ATOM    270  CH2 TRP A  16      33.007  29.767  62.217  1.00 10.81           C
ATOM    271  H  ATRP A  16      35.888  26.290  63.061  0.52 13.60           H
ATOM    272  D  BTRP A  16      35.888  26.290  63.061  0.48 13.60           D
ATOM    273  HA  TRP A  16      38.175  25.232  62.150  1.00 15.10           H
ATOM    274  HB2 TRP A  16      37.990  26.741  60.399  1.00 15.23           H
ATOM    275  HB3 TRP A  16      36.567  26.134  60.651  1.00 15.23           H
ATOM    276  HD1 TRP A  16      38.294  29.359  60.875  1.00 23.10           H
ATOM    277  HE1ATRP A  16      36.753  31.217  61.320  0.77 22.97           H
ATOM    278  DE1BTRP A  16      36.753  31.217  61.320  0.23 22.97           D
ATOM    279  HE3 TRP A  16      34.242  26.740  61.846  1.00 14.43           H
ATOM    280  HZ2 TRP A  16      34.096  31.477  61.936  1.00 16.81           H
ATOM    281  HZ3 TRP A  16      32.278  27.917  62.268  1.00 14.23           H
ATOM    282  HH2 TRP A  16      32.212  30.179  62.470  1.00 12.97           H
  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=False)
  # print in_h.as_pdb_string()
  # in_h.sort_atoms_in_place()
  # print "="*50
  # print in_h.as_pdb_string()
  validate_result(in_h, out_h)

def exercise_12():
  """
  Beginning of a chain, CA with 3 hydrogens: H1, H2, H3
  """
  in_pdb_str = """\
ATOM      1  N   PHE A   1       2.305   1.496  17.812  1.00  1.22           N
ATOM      2  CA  PHE A   1       1.374   2.648  17.647  1.00  1.07           C
ATOM      3  C   PHE A   1       2.067   3.806  16.938  1.00  1.09           C
ATOM      4  O   PHE A   1       3.294   3.909  16.951  1.00  1.16           O
ATOM      5  CB  PHE A   1       0.850   3.108  19.009  1.00  1.41           C
ATOM      6  CG  PHE A   1       1.930   3.539  19.961  1.00  1.13           C
ATOM      7  CD1 PHE A   1       2.323   4.866  20.029  1.00  1.42           C
ATOM      8  CD2 PHE A   1       2.552   2.618  20.787  1.00  1.48           C
ATOM      9  CE1 PHE A   1       3.316   5.265  20.903  1.00  1.70           C
ATOM     10  CE2 PHE A   1       3.546   3.011  21.663  1.00  1.44           C
ATOM     11  CZ  PHE A   1       3.928   4.336  21.721  1.00  1.76           C
ATOM     12  H1  PHE A   1       2.118   1.070  18.571  1.00  1.40           H
ATOM     13  H2  PHE A   1       2.210   0.937  17.126  1.00  1.40           H
ATOM     14  H3  PHE A   1       3.143   1.794  17.838  1.00  1.40           H
ATOM     15  HA  PHE A   1       0.616   2.372  17.108  1.00  1.25           H
  """
  in_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=True)
  out_h = iotbx.pdb.input(source_info=None, lines=in_pdb_str).construct_hierarchy(sort_atoms=False)
  validate_result(in_h, out_h)

if (__name__ == "__main__"):
  exercise_1()
  exercise_2()
  exercise_3()
  exercise_4()
  exercise_5()
  exercise_6()
  exercise_7()
  exercise_8()
  exercise_9()
  exercise_10()
  exercise_11()
  print("OK")
