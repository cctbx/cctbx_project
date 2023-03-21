
from __future__ import absolute_import, division, print_function
from iotbx.pdb import hierarchy
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.monomer_library import server
from mmtbx.validation import model_properties
import iotbx.pdb
from libtbx.test_utils import show_diff, approx_equal
from libtbx.easy_pickle import loads, dumps
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO

def exercise_1():
  pdb_raw = """\
ATOM   1134  N   LYS A  82       5.933  36.285  21.572  1.00 70.94           N
ATOM   1135  CA  LYS A  82       6.564  37.423  20.931  1.00 76.69           C
ATOM   1136  C   LYS A  82       5.553  38.547  20.756  1.00 78.75           C
ATOM   1137  O   LYS A  82       5.325  39.038  19.654  1.00 86.47           O
ATOM   1138  CB  LYS A  82       7.179  37.024  19.583  1.00 82.32           C
ATOM   1139  CG  LYS A  82       8.190  38.035  19.048  0.00 70.34           C
ATOM   1140  CD  LYS A  82       9.429  38.129  19.944  0.00 67.69           C
ATOM   1141  CE  LYS A  82       9.983  39.545  20.014  0.00 64.44           C
ATOM   1142  NZ  LYS A  82      10.933  39.832  18.908  0.00 61.45           N
ATOM   1143  H   LYS A  82       5.139  36.115  21.291  1.00 85.12           H
ATOM   1144  HA  LYS A  82       7.279  37.749  21.501  1.00 92.03           H
ATOM   1145  HB2 LYS A  82       6.469  36.939  18.928  1.00 98.78           H
ATOM   1146  HB3 LYS A  82       7.636  36.175  19.687  1.00 98.78           H
ATOM   1147  HG2 LYS A  82       8.476  37.762  18.163  0.00 84.41           H
ATOM   1148  HG3 LYS A  82       7.775  38.912  19.011  0.00 84.41           H
ATOM   1149  HD2 LYS A  82       9.193  37.853  20.843  0.00 81.23           H
ATOM   1150  HD3 LYS A  82      10.122  37.551  19.589  0.00 81.23           H
ATOM   1151  HE2 LYS A  82       9.249  40.177  19.952  0.00 77.33           H
ATOM   1152  HE3 LYS A  82      10.453  39.662  20.854  0.00 77.33           H
ATOM   1153  HZ1 LYS A  82      11.237  40.666  18.977  0.00 73.75           H
ATOM   1154  HZ2 LYS A  82      10.523  39.738  18.123  0.00 73.75           H
ATOM   1155  HZ3 LYS A  82      11.621  39.269  18.944  0.00 73.75           H
ATOM   1156  N   LYS A  83       4.936  38.927  21.866  1.00 75.79           N
ATOM   1157  CA  LYS A  83       4.177  40.172  21.966  1.00 82.80           C
ATOM   1158  C   LYS A  83       4.081  40.508  23.460  1.00 86.23           C
ATOM   1159  O   LYS A  83       2.978  40.521  24.017  1.00 79.81           O
ATOM   1160  CB  LYS A  83       2.790  40.044  21.332  1.00 79.16           C
ATOM   1161  CG  LYS A  83       2.038  41.342  21.175  0.00 70.42           C
ATOM   1162  CD  LYS A  83       2.072  41.803  19.735  0.00 66.90           C
ATOM   1163  CE  LYS A  83       1.295  43.089  19.552  0.00 62.46           C
ATOM   1164  NZ  LYS A  83       1.004  43.350  18.118  0.00 60.73           N
ATOM   1165  H   LYS A  83       4.940  38.470  22.594  1.00 90.95           H
ATOM   1166  HA  LYS A  83       4.658  40.885  21.518  1.00 99.36           H
ATOM   1167  HB2 LYS A  83       2.251  39.459  21.887  1.00 95.00           H
ATOM   1168  HB3 LYS A  83       2.890  39.655  20.449  1.00 95.00           H
ATOM   1169  HG2 LYS A  83       1.113  41.213  21.435  0.00 84.51           H
ATOM   1170  HG3 LYS A  83       2.453  42.024  21.726  0.00 84.51           H
ATOM   1171  HD2 LYS A  83       2.992  41.962  19.471  0.00 80.28           H
ATOM   1172  HD3 LYS A  83       1.672  41.123  19.171  0.00 80.28           H
ATOM   1173  HE2 LYS A  83       0.452  43.024  20.027  0.00 74.95           H
ATOM   1174  HE3 LYS A  83       1.818  43.830  19.896  0.00 74.95           H
ATOM   1175  HZ1 LYS A  83       0.521  42.683  17.780  0.00 72.87           H
ATOM   1176  HZ2 LYS A  83       1.764  43.417  17.661  0.00 72.87           H
ATOM   1177  HZ3 LYS A  83       0.548  44.109  18.034  0.00 72.87           H
ATOM   3630  N   ASN A 242      -5.454  -3.027   1.145  0.00 67.69           N
ATOM   3631  CA  ASN A 242      -4.759  -2.535  -0.037  0.00 65.44           C
ATOM   3632  C   ASN A 242      -5.734  -2.397  -1.208  0.00 63.57           C
ATOM   3633  O   ASN A 242      -6.425  -3.357  -1.552  0.00 63.94           O
ATOM   3634  CB  ASN A 242      -3.626  -3.503  -0.392  0.00 63.13           C
ATOM   3635  CG  ASN A 242      -2.802  -3.044  -1.576  0.00 63.58           C
ATOM   3636  OD1 ASN A 242      -2.524  -1.862  -1.731  0.00 65.52           O
ATOM   3637  ND2 ASN A 242      -2.399  -3.988  -2.416  0.00 62.17           N
ATOM   3638  H   ASN A 242      -5.562  -3.880   1.129  0.00 81.22           H
ATOM   3639  HA  ASN A 242      -4.375  -1.665   0.151  0.00 78.53           H
ATOM   3640  HB2 ASN A 242      -3.032  -3.587   0.370  0.00 75.76           H
ATOM   3641  HB3 ASN A 242      -4.007  -4.368  -0.611  0.00 75.76           H
ATOM   3642 HD21 ASN A 242      -1.929  -3.779  -3.104  0.00 74.60           H
ATOM   3643 HD22 ASN A 242      -2.609  -4.810  -2.272  0.00 74.60           H
ATOM      2  CA ALYS A  32      10.574   8.177  11.768  0.40 71.49           C
ATOM      3  CB ALYS A  32       9.197   8.686  12.246  0.40 74.71           C
ATOM      2  CA BLYS A  32      10.574   8.177  11.768  0.40 71.49           C
ATOM      3  CB BLYS A  32       9.197   8.686  12.246  0.40 74.71           C
ATOM      5  CA AVAL A  33      11.708   5.617  14.332  0.50 71.42           C
ATOM      6  CB AVAL A  33      11.101   4.227  14.591  0.50 71.47           C
ATOM      5  CA BVAL A  33      11.708   5.617  14.332  0.40 71.42           C
ATOM      6  CB BVAL A  33      11.101   4.227  14.591  0.40 71.47           C
TER
ATOM      1  N   GLU X  18     -13.959  12.159  -6.598  1.00260.08           N
ATOM      2  CA  GLU X  18     -13.297  13.465  -6.628  1.00269.83           C
ATOM      3  C   GLU X  18     -11.946  13.282  -7.309  1.00269.18           C
ATOM      4  CB  GLU X  18     -13.128  14.035  -5.210  1.00261.96           C
ATOM      5  CG  GLU X  18     -14.455  14.401  -4.522  1.00263.56           C
ATOM      6  CD  GLU X  18     -14.291  15.239  -3.242  1.00264.89           C
ATOM      7  OE1 GLU X  18     -14.172  14.646  -2.143  1.00264.24           O
ATOM      8  OE2 GLU X  18     -14.309  16.498  -3.306  1.00264.37           O1-
HETATM  614  S   SO4 B 101      14.994  20.601  10.862  0.00  7.02           S
HETATM  615  O1  SO4 B 101      14.234  20.194  12.077  0.00  7.69           O
HETATM  616  O2  SO4 B 101      14.048  21.062   9.850  0.00  9.28           O
HETATM  617  O3  SO4 B 101      15.905  21.686  11.261  0.00  8.01           O
HETATM  618  O4  SO4 B 101      15.772  19.454  10.371  0.00  8.18           O
TER
HETATM  122  O   HOH S   1       5.334   8.357   8.032  1.00  0.00           O
HETATM  123  O   HOH S   2       5.396  15.243  10.734  1.00202.95           O
HETATM  124  O   HOH S   3     -25.334  18.357  18.032  0.00 20.00           O
"""
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_raw)
  xrs = pdb_in.xray_structure_simple()
  hierarchy = pdb_in.construct_hierarchy()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    raw_records=hierarchy.as_pdb_string(crystal_symmetry=xrs),
    crystal_symmetry=xrs,
    log=null_out())
  hierarchy.atoms().reset_i_seq()
  mstats = model_properties.model_statistics(
    pdb_hierarchy=hierarchy,
    xray_structure=xrs,
    all_chain_proxies=processed_pdb_file.all_chain_proxies,
    ignore_hd=True)
  out = StringIO()
  mstats.show(out=out)
  #print out.getvalue()
  assert not show_diff(out.getvalue(), """\
Overall:
  Number of atoms = 50  (anisotropic = 0)
  B_iso: mean =  96.0  max = 269.8  min =   0.0
  Occupancy: mean = 0.47  max = 1.00  min = 0.00
    warning: 22 atoms with zero occupancy
  69 total B-factor or occupancy problem(s) detected
  Atoms or residues with zero occupancy:
   LYS A  82   CG    occ=0.00
   LYS A  82   CD    occ=0.00
   LYS A  82   CE    occ=0.00
   LYS A  82   NZ    occ=0.00
   LYS A  83   CG    occ=0.00
   LYS A  83   CD    occ=0.00
   LYS A  83   CE    occ=0.00
   LYS A  83   NZ    occ=0.00
   ASN A 242  (all)  occ=0.00
   SO4 B 101  (all)  occ=0.00
   HOH S   3   O     occ=0.00
Macromolecules:
  Number of atoms = 42  (anisotropic = 0)
  B_iso: mean = 108.0  max = 269.8  min =  60.7
  Occupancy: mean = 0.51  max = 1.00  min = 0.00
    warning: 16 atoms with zero occupancy
  59 total B-factor or occupancy problem(s) detected
Ligands:
  Number of atoms = 5  (anisotropic = 0)
  B_iso: mean =   8.0  max =   9.3  min =   7.0
  Occupancy: mean = 0.00  max = 0.00  min = 0.00
    warning: 5 atoms with zero occupancy
  6 total B-factor or occupancy problem(s) detected
Waters:
  Number of atoms = 3  (anisotropic = 0)
  B_iso: mean =  74.3  max = 202.9  min =   0.0
  Occupancy: mean = 0.67  max = 1.00  min = 0.00
    warning: 1 atoms with zero occupancy
  4 total B-factor or occupancy problem(s) detected
(Hydrogen atoms not included in overall counts.)
""")
  assert (len(mstats.all.bad_adps) == 1)
  assert (mstats.all.n_zero_b == 1)
  mstats2 = loads(dumps(mstats))
  out1 = StringIO()
  out2 = StringIO()
  mstats.show(out=out1)
  mstats2.show(out=out2)
  assert (out1.getvalue() == out2.getvalue())
  # now with ignore_hd=False
  mstats3 = model_properties.model_statistics(
    pdb_hierarchy=hierarchy,
    xray_structure=xrs,
    all_chain_proxies=processed_pdb_file.all_chain_proxies,
    ignore_hd=False)
  out2 = StringIO()
  mstats3.show(out=out2)
  assert (out2.getvalue() != out.getvalue())
  assert ("""   LYS A  83   HZ3   occ=0.00""" in out2.getvalue())
  outliers = mstats3.all.as_gui_table_data(include_zoom=True)
  assert (len(outliers) == 86)
  # test with all_chain_proxies undefined
  mstats4 = model_properties.model_statistics(
    pdb_hierarchy=hierarchy,
    xray_structure=xrs,
    all_chain_proxies=None,
    ignore_hd=False)
  outliers = mstats4.all.as_gui_table_data(include_zoom=True)
  assert (len(outliers) == 86)

# corner case: deuterium as ligand (from 3qza)
def exercise_2():
  pdb_raw = """\
ATOM   6407  N   GLY A 388      -0.783   9.368 -16.436  1.00 51.96           N
ATOM   6408  CA  GLY A 388      -0.227   9.888 -15.197  1.00 54.04           C
ATOM   6409  C   GLY A 388      -0.637  11.320 -14.897  1.00 55.86           C
ATOM   6410  O   GLY A 388      -1.728  11.738 -15.347  1.00 56.70           O
ATOM   6411  OXT GLY A 388       0.129  12.024 -14.203  1.00 56.98           O
ATOM   6412  D   GLY A 388      -0.460   9.727 -17.309  1.00 51.44           D
ATOM   6413  HA2 GLY A 388      -0.561   9.258 -14.385  1.00 54.07           H
ATOM   6414  HA3 GLY A 388       0.843   9.835 -15.243  1.00 54.13           H
TER    6415      GLY A 388
HETATM 6416  D   D8U A 401     -12.236 -13.695 -42.992  1.00 15.23           D
HETATM 6417  O   DOD A1001      -4.151  -5.107 -38.592  1.00 13.40           O
HETATM 6418  D1  DOD A1001      -4.760  -5.026 -39.326  1.00 15.45           D
HETATM 6419  D2  DOD A1001      -4.625  -4.741 -37.845  1.00 14.81           D
"""
  mstats = get_mstats(pdb_raw)
  out = StringIO()
  mstats.show(out=out)
  assert ("Ligands:" in out.getvalue())
  assert ("B_iso: mean =  15.2  max =  15.2  min =  15.2" in out.getvalue())

# explicitly specified ligand selection
def exercise_3():
  pdb_raw = """\
ATOM      1  CA  GLY A   1      -0.227   9.888 -15.197  1.00 54.04           C
ATOM      2  CA  GLY A   2      -0.227   9.888 -15.197  1.00 54.04           C
ATOM      3  CA  GLY A   3      -0.227   9.888 -15.197  1.00 54.04           C
ATOM      4  CA  GLY A   4      -0.227   9.888 -15.197  1.00 54.04           C
ATOM      5  CA  GLY A   5      -0.227   9.888 -15.197  1.00 54.04           C
ATOM      6  CA  GLY A   6      -0.227   9.888 -15.197  1.00 54.04           C
ATOM      7  CA  GLY A   7      -0.227   9.888 -15.197  1.00 54.04           C
ATOM      8  CA  GLY A   8      -0.227   9.888 -15.197  1.00 54.04           C
ATOM      9  CA  GLY A   9      -0.227   9.888 -15.197  1.00 54.04           C
ATOM     10  CA  GLY A  10      -0.227   9.888 -15.197  1.00 54.04           C
HETATM   11  N   SEP A  11      -2.112   0.368  -0.991  1.00 20.00      A    N
HETATM   12  CA  SEP A  11      -0.692   0.284  -0.951  1.00 20.00      A    C
HETATM   13  CB  SEP A  11      -0.234   0.166   0.485  1.00 20.00      A    C
HETATM   14  OG  SEP A  11       1.130  -0.184   0.515  1.00 20.00      A    O
HETATM   15  C   SEP A  11      -0.237  -0.930  -1.727  1.00 20.00      A    C
HETATM   16  O   SEP A  11      -0.767  -2.051  -1.509  1.00 20.00      A    O
HETATM   18  P   SEP A  11       1.922  -0.008   1.871  1.00 20.00      A    P
HETATM   19  O1P SEP A  11       2.139   1.462   2.140  1.00 20.00      A    O
HETATM   20  O2P SEP A  11       3.259  -0.703   1.767  1.00 20.00      A    O-1
HETATM   21  O3P SEP A  11       1.127  -0.614   3.002  1.00 20.00      A    O-1
END"""
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_raw)
  xrs = pdb_in.xray_structure_simple()
  hierarchy = pdb_in.construct_hierarchy()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    raw_records=hierarchy.as_pdb_string(crystal_symmetry=xrs),
    crystal_symmetry=xrs,
    log=null_out())
  hierarchy.atoms().reset_i_seq()
  ligand_sel = hierarchy.atom_selection_cache().selection("resname SEP")
  mstats = model_properties.model_statistics(
    pdb_hierarchy=hierarchy,
    xray_structure=xrs,
    all_chain_proxies=processed_pdb_file.all_chain_proxies,
    ligand_selection=ligand_sel,
    ignore_hd=True)
  out = StringIO()
  mstats.show(out=out)
  assert (mstats.n_protein == 10)
  assert ("Ligands:" in out.getvalue())
  assert approx_equal(mstats.macromolecules.b_mean, 54.04)
  # now with just the raw selection string
  mstats = model_properties.model_statistics(
    pdb_hierarchy=hierarchy,
    xray_structure=xrs,
    all_chain_proxies=processed_pdb_file.all_chain_proxies,
    ligand_selection="resname SEP",
    ignore_hd=True)
  out = StringIO()
  mstats.show(out=out)
  assert (mstats.n_protein == 10)
  assert ("Ligands:" in out.getvalue())

def get_mstats(pdb_raw):
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_raw)
  xrs = pdb_in.xray_structure_simple()
  hierarchy = pdb_in.construct_hierarchy()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    raw_records=hierarchy.as_pdb_string(crystal_symmetry=xrs),
    crystal_symmetry=xrs,
    log=null_out())
  hierarchy.atoms().reset_i_seq()
  mstats = model_properties.model_statistics(
    pdb_hierarchy=hierarchy,
    xray_structure=xrs,
    all_chain_proxies=processed_pdb_file.all_chain_proxies,
    ignore_hd=True)
  return mstats

def test_zero_occupancy():
  pdb_raw = '''\
ATOM    100  N   LYS A  82       5.933  36.285  21.572  0.00 70.94           N
ATOM    101  CA  LYS A  82       6.564  37.423  20.931  1.00 76.69           C
ATOM    102  C   LYS A  82       5.553  38.547  20.756  1.00 78.75           C
ATOM    103  O   LYS A  82       5.325  39.038  19.654  1.00 86.47           O
ATOM    104  CB  LYS A  82       7.179  37.024  19.583  1.00 82.32           C
ATOM    105  CG  LYS A  82       8.190  38.035  19.048  1.00 70.34           C
ATOM    106  CD  LYS A  82       9.429  38.129  19.944  1.00 67.69           C
ATOM    107  CE  LYS A  82       9.983  39.545  20.014  1.00 64.44           C
ATOM    108  NZ  LYS A  82      10.933  39.832  18.908  1.00 61.45           N
ATOM    109  H   LYS A  82       5.139  36.115  21.291  1.00 85.12           H
ATOM    110  HA  LYS A  82       7.279  37.749  21.501  1.00 92.03           H
ATOM    111  HB2 LYS A  82       6.469  36.939  18.928  1.00 98.78           H
ATOM    112  HB3 LYS A  82       7.636  36.175  19.687  1.00 98.78           H
ATOM    113  HG2 LYS A  82       8.476  37.762  18.163  1.00 84.41           H
ATOM    114  HG3 LYS A  82       7.775  38.912  19.011  1.00 84.41           H
ATOM    115  HD2 LYS A  82       9.193  37.853  20.843  1.00 81.23           H
ATOM    116  HD3 LYS A  82      10.122  37.551  19.589  1.00 81.23           H
ATOM    117  HE2 LYS A  82       9.249  40.177  19.952  1.00 77.33           H
ATOM    118  HE3 LYS A  82      10.453  39.662  20.854  1.00 77.33           H
ATOM    119  HZ1 LYS A  82      11.237  40.666  18.977  1.00 73.75           H
ATOM    120  HZ2 LYS A  82      10.523  39.738  18.123  1.00 73.75           H
ATOM    121  HZ3 LYS A  82      11.621  39.269  18.944  1.00 73.75           H
'''
  mstats = get_mstats(pdb_raw)
  assert (mstats.all.n_outliers == 3)           # atom 100 has 0 occupancy
  assert (len(mstats.all.zero_occ) == 1)        # zero occupancy
  assert (len(mstats.all.partial_occ) == 1)     # occupancy < 1
  assert (len(mstats.all.different_occ) == 1)   # occupancies in residue differ
  return True

def test_partial_occupancy():
  pdb_raw = '''\
ATOM    100  N   LYS A  82       5.933  36.285  21.572  1.00 70.94           N
ATOM    101  CA  LYS A  82       6.564  37.423  20.931  1.00 76.69           C
ATOM    102  C   LYS A  82       5.553  38.547  20.756  1.00 78.75           C
ATOM    103  O   LYS A  82       5.325  39.038  19.654  1.00 86.47           O
ATOM    104  CB ALYS A  82       7.179  37.024  19.583  0.40 82.32           C
ATOM    105  CG ALYS A  82       8.190  38.035  19.048  0.40 70.34           C
ATOM    106  CD ALYS A  82       9.429  38.129  19.944  0.40 67.69           C
ATOM    107  CE ALYS A  82       9.983  39.545  20.014  0.40 64.44           C
ATOM    108  NZ ALYS A  82      10.933  39.832  18.908  0.40 61.45           N
ATOM    109  H  ALYS A  82       5.139  36.115  21.291  0.40 85.12           H
ATOM    110  HA ALYS A  82       7.279  37.749  21.501  0.40 92.03           H
ATOM    111  HB2ALYS A  82       6.469  36.939  18.928  0.40 98.78           H
ATOM    112  HB3ALYS A  82       7.636  36.175  19.687  0.40 98.78           H
ATOM    113  HG2ALYS A  82       8.476  37.762  18.163  0.40 84.41           H
ATOM    114  HG3ALYS A  82       7.775  38.912  19.011  0.40 84.41           H
ATOM    115  HD2ALYS A  82       9.193  37.853  20.843  0.40 81.23           H
ATOM    116  HD3ALYS A  82      10.122  37.551  19.589  0.40 81.23           H
ATOM    117  HE2ALYS A  82       9.249  40.177  19.952  0.40 77.33           H
ATOM    118  HE3ALYS A  82      10.453  39.662  20.854  0.40 77.33           H
ATOM    119  HZ1ALYS A  82      11.237  40.666  18.977  0.40 73.75           H
ATOM    120  HZ2ALYS A  82      10.523  39.738  18.123  0.40 73.75           H
ATOM    121  HZ3ALYS A  82      11.621  39.269  18.944  0.40 73.75           H
ATOM    122  CB BLYS A  82       7.179  37.024  19.583  0.60 82.32           C
ATOM    123  CG BLYS A  82       8.190  38.035  19.048  0.60 70.34           C
ATOM    124  CD BLYS A  82       9.429  38.129  19.944  0.60 67.69           C
ATOM    125  CE BLYS A  82       9.983  39.545  20.014  0.60 64.44           C
ATOM    126  NZ BLYS A  82      10.933  39.832  18.908  0.60 61.45           N
ATOM    127  H  BLYS A  82       5.139  36.115  21.291  0.60 85.12           H
ATOM    128  HA BLYS A  82       7.279  37.749  21.501  0.60 92.03           H
ATOM    129  HB2BLYS A  82       6.469  36.939  18.928  0.60 98.78           H
ATOM    130  HB3BLYS A  82       7.636  36.175  19.687  0.60 98.78           H
ATOM    131  HG2BLYS A  82       8.476  37.762  18.163  0.60 84.41           H
ATOM    132  HG3BLYS A  82       7.775  38.912  19.011  0.60 84.41           H
ATOM    133  HD2BLYS A  82       9.193  37.853  20.843  0.60 81.23           H
ATOM    134  HD3BLYS A  82      10.122  37.551  19.589  0.60 81.23           H
ATOM    135  HE2BLYS A  82       9.249  40.177  19.952  0.60 77.33           H
ATOM    136  HE3BLYS A  82      10.453  39.662  20.854  0.60 77.33           H
ATOM    137  HZ1BLYS A  82      11.237  40.666  18.977  0.60 73.75           H
ATOM    138  HZ2BLYS A  82      10.523  39.738  18.123  0.60 73.75           H
ATOM    139  HZ3BLYS A  82      11.621  39.269  18.944  0.50 73.75           H
'''
  mstats = get_mstats(pdb_raw)                 # atom 121, 139 sum < 1
  assert (mstats.all.n_outliers == 3)
  assert (len(mstats.all.partial_occ) == 2)    # occupancy < 1
  assert (len(mstats.all.different_occ) == 1)  # occupancies in residue differ
  return True

if (__name__ == "__main__"):
  exercise_1()
  exercise_2()
  exercise_3()
  test_zero_occupancy()
  test_partial_occupancy()
  print("OK")
