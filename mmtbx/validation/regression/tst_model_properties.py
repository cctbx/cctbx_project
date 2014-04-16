
# TODO test with waters, ligands, full PDB interpretation

from __future__ import division
from mmtbx.validation import model_properties
import iotbx.pdb.hierarchy
from libtbx.test_utils import show_diff
from cStringIO import StringIO

def exercise () :
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
"""
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string=pdb_raw)
  xrs = pdb_in.input.xray_structure_simple()
  mstats = model_properties.model_statistics(
    pdb_hierarchy=pdb_in.hierarchy,
    xray_structure=xrs,
    ignore_hd=True)
  out = StringIO()
  mstats.show(out=out)
  assert not show_diff(out.getvalue(), """\
Overall:
  Number of atoms = 26  (anisotropic = 0)
  B_iso: mean =  70.7  max =  86.5  min =  60.7
  Occupancy: mean = 0.38  max = 1.00  min = 0.00
    warning: 16 atoms with zero occupancy
  16 total B-factor or occupancy problems detected
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
   ASN A 242  (all)  occ=0.00
   ASN A 242  (all)  occ=0.00
   ASN A 242  (all)  occ=0.00
   ASN A 242  (all)  occ=0.00
   ASN A 242  (all)  occ=0.00
   ASN A 242  (all)  occ=0.00
   ASN A 242  (all)  occ=0.00
(Hydrogen atoms not included in overall counts.)
""")
  mstats2 = model_properties.model_statistics(
    pdb_hierarchy=pdb_in.hierarchy,
    xray_structure=xrs,
    ignore_hd=False)
  out2 = StringIO()
  mstats2.show(out=out2)
  assert (out2.getvalue() != out.getvalue())
  assert ("""   LYS A  83   HZ3   occ=0.00""" in out2.getvalue())

if (__name__ == "__main__") :
  exercise()
  print "OK"
