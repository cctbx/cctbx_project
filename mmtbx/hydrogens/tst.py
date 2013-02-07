from __future__ import division
import iotbx.pdb
from mmtbx import hydrogens
from mmtbx import monomer_library

m1_str = """\
CRYST1    9.756    9.585    9.568  90.00  90.00  90.00 P 1
ATOM      1  N   SER H   7       4.184   4.600   6.227  1.00 10.00           N
ATOM      2  CA  SER H   7       4.535   3.950   4.970  1.00 10.00           C
ATOM      3  C   SER H   7       6.049   3.863   4.811  1.00 10.00           C
ATOM      4  O   SER H   7       6.790   4.650   5.400  1.00 10.00           O
ATOM      5  CB  SER H   7       3.931   4.713   3.791  1.00 10.00           C
ATOM      6  OG  SER H   7       4.412   6.045   3.745  1.00 10.00           O
ATOM      8  HA  SER H   7       4.176   3.049   4.965  1.00 10.00           H
ATOM      9  HB2 SER H   7       4.173   4.263   2.966  1.00 10.00           H
ATOM     10  HB3 SER H   7       2.966   4.730   3.887  1.00 10.00           H
ATOM     11  HG  SER H   7       4.200   6.456   4.459  1.00 10.00           H
"""

m2_str = """\
CRYST1    9.756    9.585    9.568  90.00  90.00  90.00 P 1
ATOM      1  N   SER H   7       4.184   4.600   6.227  1.00 10.00           N
ATOM      2  CA  SER H   7       4.535   3.950   4.970  1.00 10.00           C
ATOM      3  C   SER H   7       6.049   3.863   4.811  1.00 10.00           C
ATOM      4  O   SER H   7       6.790   4.650   5.400  1.00 10.00           O
ATOM      5  CB  SER H   7       3.931   4.713   3.791  1.00 10.00           C
ATOM      6  OG  SER H   7       4.412   6.045   3.745  1.00 10.00           O
ATOM      8  HA  SER H   7       4.176   3.049   4.965  1.00 10.00           H
ATOM      9  HB2 SER H   7       4.173   4.263   2.966  1.00 10.00           H
ATOM     10  HB3 SER H   7       2.966   4.730   3.887  1.00 10.00           H
ATOM     11  HG ASER H   7       4.200   6.456   4.459  0.50 10.00           H
ATOM     11  DG BSER H   7       4.200   6.456   4.459  0.50 10.00           D
"""

m3_str = """
CRYST1    9.756    9.585    9.568  90.00  90.00  90.00 P 1
ATOM      1  N   SER H   7       4.184   4.600   6.227  1.00 10.00           N
ATOM      2  CA  SER H   7       4.535   3.950   4.970  1.00 10.00           C
ATOM      3  C   SER H   7       6.049   3.863   4.811  1.00 10.00           C
ATOM      4  O   SER H   7       6.790   4.650   5.400  1.00 10.00           O
ATOM      5  CB  SER H   7       3.931   4.713   3.791  1.00 10.00           C
ATOM      6  OG ASER H   7       4.412   6.045   3.745  0.50 10.00           O
ATOM      6  OG BSER H   7       2.515   4.664   3.823  0.50 10.00           O
ATOM      0  HA  SER H   7       4.173   3.050   4.984  1.00 10.00           H
ATOM      0  HG ASER H   7       4.483   6.348   4.536  0.50 10.00           H
ATOM      0  HG BSER H   7       2.253   4.558   4.625  0.50 10.00           H
"""

m4_str = """
CRYST1    8.942    9.538    9.274  90.00  90.00  90.00 P 1
ATOM      1  N   THR H  21       6.000   3.103   5.518  1.00 10.00           N
ATOM      2  CA  THR H  21       5.114   4.258   5.439  1.00 10.00           C
ATOM      3  C   THR H  21       3.895   4.074   6.337  1.00 10.00           C
ATOM      4  O   THR H  21       3.236   3.035   6.301  1.00 10.00           O
ATOM      5  CB  THR H  21       4.646   4.511   3.992  1.00 10.00           C
ATOM      6  OG1 THR H  21       5.786   4.660   3.137  1.00 10.00           O
ATOM      7  CG2 THR H  21       3.792   5.769   3.916  1.00 10.00           C
ATOM      8 DG21 THR H  21       3.708   6.055   2.993  1.00 10.00           D
ATOM      9 DG22 THR H  21       2.906   5.595   4.272  1.00 10.00           D
ATOM     10 DG23 THR H  21       4.205   6.482   4.429  1.00 10.00           D
"""

m5_str = """
CRYST1   19.756   19.585   19.568  90.00  90.00  90.00 P 1
ATOM      1  N   SER H   7       7.561   4.409   2.370  1.00 10.00           N
ATOM      2  CA  SER H   7       6.281   4.192   3.051  1.00 10.00           C
ATOM      3  C   SER H   7       6.112   3.781   4.520  1.00 10.00           C
ATOM      4  O   SER H   7       6.620   2.744   4.947  1.00 10.00           O
ATOM      5  CB  SER H   7       5.067   3.972   2.131  1.00 10.00           C
ATOM      6  HA  SER H   7       6.473   3.250   2.923  1.00 10.00           H
ATOM      7  OG ASER H   7       3.884   3.774   2.885  0.50 10.00           O
ATOM      8  HG ASER H   7       3.729   4.456   3.369  0.50 10.00           H
ATOM      9  OG BSER H   7       5.247   2.836   1.302  0.50 10.00           O
ATOM     10  HG BSER H   7       5.189   2.125   1.764  0.50 10.00           H
ATOM     11  N   THR H   8       5.397   4.602   5.281  1.00 10.00           N
ATOM     12  CA  THR H   8       5.161   4.328   6.694  1.00 10.00           C
ATOM     13  C   THR H   8       3.963   3.404   6.887  1.00 10.00           C
ATOM     14  O   THR H   8       4.113   2.257   7.310  1.00 10.00           O
ATOM     15  CB  THR H   8       4.938   5.627   7.495  1.00 10.00           C
ATOM     16  OG1 THR H   8       4.686   5.308   8.869  1.00 10.00           O
ATOM     17  CG2 THR H   8       3.757   6.410   6.935  1.00 10.00           C
ATOM     18 DG21 THR H   8       3.891   6.587   5.991  1.00 10.00           D
ATOM     19 DG22 THR H   8       2.937   5.904   7.048  1.00 10.00           D
ATOM     20 DG23 THR H   8       3.670   7.255   7.403  1.00 10.00           D
"""

loop = [(m1_str, [([4,5], [9])]),
        (m2_str, [([4,5], [9]), ([4,5], [10])]),
        (m3_str, [([4,6], [7]), ([4,8], [9])]),
        (m4_str, [([4,6], [7,8,9])]),
        (m5_str, [([4,6], [7]), ([4,8], [9]), ([14,16], [17,18,19])])]

def run():
  mon_lib_srv = monomer_library.server.server()
  for l in loop:
    pdb_inp = iotbx.pdb.input(source_info=None, lines=l[0])
    if 0: pdb_inp.write_pdb_file(file_name = "m1.pdb")
    ph = pdb_inp.construct_hierarchy()
    ph.atoms().reset_i_seq()
    xrs_answer = pdb_inp.xray_structure_simple()
    sel = hydrogens.rotatable(pdb_hierarchy=ph, mon_lib_srv=mon_lib_srv)
    assert sel == l[1]

if (__name__ == "__main__"):
  run()
  print "OK"
