from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.test_utils import approx_equal

pdb_6exy_ARG_168_D    = """\
CRYST1   37.314   58.352   63.867  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ARG A 168      22.726 -15.648  12.491  1.00 14.41           N
ATOM      2  CA  ARG A 168      22.811 -14.357  13.143  1.00 15.09           C
ATOM      3  C   ARG A 168      21.569 -13.549  12.811  1.00 13.95           C
ATOM      4  O   ARG A 168      20.857 -13.814  11.839  1.00 14.06           O
ATOM      5  CB  ARG A 168      24.057 -13.582  12.689  1.00 17.54           C
ATOM      6  CG  ARG A 168      25.343 -14.366  12.775  1.00 19.91           C
ATOM      7  CD  ARG A 168      25.613 -14.857  14.174  1.00 20.97           C
ATOM      8  NE  ARG A 168      26.737 -15.790  14.190  1.00 21.81           N
ATOM      9  CZ  ARG A 168      26.658 -17.059  14.578  1.00 21.76           C
ATOM     10  NH1 ARG A 168      27.742 -17.813  14.539  1.00 21.87           N
ATOM     11  NH2 ARG A 168      25.521 -17.587  15.017  1.00 21.76           N
ATOM     12  D   ARG A 168      22.668 -15.599  11.475  1.00 14.15           D
ATOM     13  DA  ARG A 168      22.844 -14.490  14.224  1.00 15.42           D
ATOM     14  DB2 ARG A 168      23.915 -13.295  11.650  1.00 17.96           D
ATOM     15  DB3 ARG A 168      24.174 -12.692  13.305  1.00 17.77           D
ATOM     16  DD2 ARG A 168      25.863 -14.006  14.804  1.00 21.25           D
ATOM     17  DD3 ARG A 168      24.731 -15.354  14.564  1.00 21.33           D
ATOM     18  DE  ARG A 168      27.476 -15.591  13.512  1.00 21.85           D
ATOM     19  DG2 ARG A 168      25.288 -15.225  12.102  1.00 20.44           D
ATOM     20  DG3 ARG A 168      26.173 -13.708  12.500  1.00 20.06           D
ATOM     21 DH11 ARG A 168      28.617 -17.413  14.207  1.00 22.11           D
ATOM     22 DH12 ARG A 168      27.696 -18.788  14.835  1.00 21.99           D
ATOM     23 DH21 ARG A 168      25.514 -18.567  15.308  1.00 21.74           D
ATOM     24 DH22 ARG A 168      24.658 -17.044  15.059  1.00 21.55           D
TER
END
"""

pdb_6exy_ARG_168_H    = """\
CRYST1   37.314   58.352   63.867  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ARG A 168      22.726 -15.648  12.491  1.00 14.41           N
ATOM      2  CA  ARG A 168      22.811 -14.357  13.143  1.00 15.09           C
ATOM      3  C   ARG A 168      21.569 -13.549  12.811  1.00 13.95           C
ATOM      4  O   ARG A 168      20.857 -13.814  11.839  1.00 14.06           O
ATOM      5  CB  ARG A 168      24.057 -13.582  12.689  1.00 17.54           C
ATOM      6  CG  ARG A 168      25.343 -14.366  12.775  1.00 19.91           C
ATOM      7  CD  ARG A 168      25.613 -14.857  14.174  1.00 20.97           C
ATOM      8  NE  ARG A 168      26.737 -15.790  14.190  1.00 21.81           N
ATOM      9  CZ  ARG A 168      26.658 -17.059  14.578  1.00 21.76           C
ATOM     10  NH1 ARG A 168      27.742 -17.813  14.539  1.00 21.87           N
ATOM     11  NH2 ARG A 168      25.521 -17.587  15.017  1.00 21.76           N
ATOM     12  H   ARG A 168      22.668 -15.599  11.475  1.00 14.15           H
ATOM     13  HA  ARG A 168      22.844 -14.490  14.224  1.00 15.42           H
ATOM     14  HB2 ARG A 168      23.915 -13.295  11.650  1.00 17.96           H
ATOM     15  HB3 ARG A 168      24.174 -12.692  13.305  1.00 17.77           H
ATOM     16  HD2 ARG A 168      25.863 -14.006  14.804  1.00 21.25           H
ATOM     17  HD3 ARG A 168      24.731 -15.354  14.564  1.00 21.33           H
ATOM     18  HE  ARG A 168      27.476 -15.591  13.512  1.00 21.85           H
ATOM     19  HG2 ARG A 168      25.288 -15.225  12.102  1.00 20.44           H
ATOM     20  HG3 ARG A 168      26.173 -13.708  12.500  1.00 20.06           H
ATOM     21 HH11 ARG A 168      28.617 -17.413  14.207  1.00 22.11           H
ATOM     22 HH12 ARG A 168      27.696 -18.788  14.835  1.00 21.99           H
ATOM     23 HH21 ARG A 168      25.514 -18.567  15.308  1.00 21.74           H
ATOM     24 HH22 ARG A 168      24.658 -17.044  15.059  1.00 21.55           H
TER
END
"""

pdb_6exy_ARG_168_HD   = """\
CRYST1   37.314   58.352   63.867  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ARG A 168      22.726 -15.648  12.491  1.00 14.41           N
ATOM      2  CA  ARG A 168      22.811 -14.357  13.143  1.00 15.09           C
ATOM      3  C   ARG A 168      21.569 -13.549  12.811  1.00 13.95           C
ATOM      4  O   ARG A 168      20.857 -13.814  11.839  1.00 14.06           O
ATOM      5  CB  ARG A 168      24.057 -13.582  12.689  1.00 17.54           C
ATOM      6  CG  ARG A 168      25.343 -14.366  12.775  1.00 19.91           C
ATOM      7  CD  ARG A 168      25.613 -14.857  14.174  1.00 20.97           C
ATOM      8  NE  ARG A 168      26.737 -15.790  14.190  1.00 21.81           N
ATOM      9  CZ  ARG A 168      26.658 -17.059  14.578  1.00 21.76           C
ATOM     10  NH1 ARG A 168      27.742 -17.813  14.539  1.00 21.87           N
ATOM     11  NH2 ARG A 168      25.521 -17.587  15.017  1.00 21.76           N

ATOM     12  H  AARG A 168      22.668 -15.599  11.475  0.50 14.15           H
ATOM     13  HA AARG A 168      22.844 -14.490  14.224  0.50 15.42           H
ATOM     14  HB2AARG A 168      23.915 -13.295  11.650  0.50 17.96           H
ATOM     15  HB3AARG A 168      24.174 -12.692  13.305  0.50 17.77           H
ATOM     16  HD2AARG A 168      25.863 -14.006  14.804  0.50 21.25           H
ATOM     17  HD3AARG A 168      24.731 -15.354  14.564  0.50 21.33           H
ATOM     18  HE AARG A 168      27.476 -15.591  13.512  0.50 21.85           H
ATOM     19  HG2AARG A 168      25.288 -15.225  12.102  0.50 20.44           H
ATOM     20  HG3AARG A 168      26.173 -13.708  12.500  0.50 20.06           H
ATOM     21 HH11AARG A 168      28.617 -17.413  14.207  0.50 22.11           H
ATOM     22 HH12AARG A 168      27.696 -18.788  14.835  0.50 21.99           H
ATOM     23 HH21AARG A 168      25.514 -18.567  15.308  0.50 21.74           H
ATOM     24 HH22AARG A 168      24.658 -17.044  15.059  0.50 21.55           H

ATOM     12  D  BARG A 168      22.668 -15.599  11.475  0.50 14.15           D
ATOM     13  DA BARG A 168      22.844 -14.490  14.224  0.50 15.42           D
ATOM     14  DB2BARG A 168      23.915 -13.295  11.650  0.50 17.96           D
ATOM     15  DB3BARG A 168      24.174 -12.692  13.305  0.50 17.77           D
ATOM     16  DD2BARG A 168      25.863 -14.006  14.804  0.50 21.25           D
ATOM     17  DD3BARG A 168      24.731 -15.354  14.564  0.50 21.33           D
ATOM     18  DE BARG A 168      27.476 -15.591  13.512  0.50 21.85           D
ATOM     19  DG2BARG A 168      25.288 -15.225  12.102  0.50 20.44           D
ATOM     20  DG3BARG A 168      26.173 -13.708  12.500  0.50 20.06           D
ATOM     21 DH11BARG A 168      28.617 -17.413  14.207  0.50 22.11           D
ATOM     22 DH12BARG A 168      27.696 -18.788  14.835  0.50 21.99           D
ATOM     23 DH21BARG A 168      25.514 -18.567  15.308  0.50 21.74           D
ATOM     24 DH22BARG A 168      24.658 -17.044  15.059  0.50 21.55           D
TER
END
"""

pdb_6exy_ARG_168_ac_D = """\
CRYST1   37.314   58.352   63.867  90.00  90.00  90.00 P 21 21 21
ATOM      1  N  AARG A 168      22.726 -15.648  12.491  0.46 14.41           N
ATOM      2  CA AARG A 168      22.811 -14.357  13.143  0.46 15.09           C
ATOM      3  C  AARG A 168      21.569 -13.549  12.811  0.46 13.95           C
ATOM      4  O  AARG A 168      20.857 -13.814  11.839  0.46 14.06           O
ATOM      5  CB AARG A 168      24.057 -13.582  12.689  0.46 17.54           C
ATOM      6  CG AARG A 168      25.343 -14.366  12.775  0.46 19.91           C
ATOM      7  CD AARG A 168      25.613 -14.857  14.174  0.46 20.97           C
ATOM      8  NE AARG A 168      26.737 -15.790  14.190  0.46 21.81           N
ATOM      9  CZ AARG A 168      26.658 -17.059  14.578  0.46 21.76           C
ATOM     10  NH1AARG A 168      27.742 -17.813  14.539  0.46 21.87           N
ATOM     11  NH2AARG A 168      25.521 -17.587  15.017  0.46 21.76           N
ATOM     12  D  AARG A 168      22.668 -15.599  11.475  0.46 14.15           D
ATOM     13  DA AARG A 168      22.844 -14.490  14.224  0.46 15.42           D
ATOM     14  DB2AARG A 168      23.915 -13.295  11.650  0.46 17.96           D
ATOM     15  DB3AARG A 168      24.174 -12.692  13.305  0.46 17.77           D
ATOM     16  DD2AARG A 168      25.863 -14.006  14.804  0.46 21.25           D
ATOM     17  DD3AARG A 168      24.731 -15.354  14.564  0.46 21.33           D
ATOM     18  DE AARG A 168      27.476 -15.591  13.512  0.46 21.85           D
ATOM     19  DG2AARG A 168      25.288 -15.225  12.102  0.46 20.44           D
ATOM     20  DG3AARG A 168      26.173 -13.708  12.500  0.46 20.06           D
ATOM     21 DH11AARG A 168      28.617 -17.413  14.207  0.46 22.11           D
ATOM     22 DH12AARG A 168      27.696 -18.788  14.835  0.46 21.99           D
ATOM     23 DH21AARG A 168      25.514 -18.567  15.308  0.46 21.74           D
ATOM     24 DH22AARG A 168      24.658 -17.044  15.059  0.46 21.55           D
ATOM     25  N  BARG A 168      22.711 -15.649  12.474  0.54 14.14           N
ATOM     26  CA BARG A 168      22.916 -14.337  13.062  0.54 14.70           C
ATOM     27  C  BARG A 168      21.754 -13.431  12.679  0.54 13.29           C
ATOM     28  O  BARG A 168      21.184 -13.538  11.588  0.54 13.90           O
ATOM     29  CB BARG A 168      24.218 -13.716  12.533  0.54 17.49           C
ATOM     30  CG BARG A 168      25.439 -14.599  12.700  0.54 20.99           C
ATOM     31  CD BARG A 168      26.144 -14.335  14.005  0.54 22.91           C
ATOM     32  NE BARG A 168      27.298 -15.226  14.175  0.54 24.26           N
ATOM     33  CZ BARG A 168      28.570 -14.887  13.964  0.54 24.87           C
ATOM     34  NH1BARG A 168      29.519 -15.802  14.139  0.54 24.92           N
ATOM     35  NH2BARG A 168      28.905 -13.657  13.596  0.54 24.77           N
ATOM     36  D  BARG A 168      22.563 -15.624  11.467  0.54 13.85           D
ATOM     37  DA BARG A 168      22.959 -14.410  14.146  0.54 15.03           D
ATOM     38  DB2BARG A 168      24.096 -13.533  11.468  0.54 17.99           D
ATOM     39  DB3BARG A 168      24.418 -12.777  13.046  0.54 17.29           D
ATOM     40  DD2BARG A 168      26.474 -13.296  14.007  0.54 23.34           D
ATOM     41  DD3BARG A 168      25.450 -14.509  14.828  0.54 22.84           D
ATOM     42  DE BARG A 168      27.141 -16.106  14.675  0.54 24.26           D
ATOM     43  DG2BARG A 168      25.143 -15.649  12.671  0.54 21.23           D
ATOM     44  DG3BARG A 168      26.149 -14.391  11.896  0.54 21.33           D
ATOM     45 DH11BARG A 168      29.265 -16.742  14.428  0.54 24.77           D
ATOM     46 DH12BARG A 168      30.500 -15.563  13.994  0.54 25.10           D
ATOM     47 DH21BARG A 168      29.887 -13.425  13.448  0.54 24.89           D
ATOM     48 DH22BARG A 168      28.193 -12.942  13.455  0.54 24.59           D
TER
END
"""

pdb_6exy_ARG_168_ac_H = """\
CRYST1   37.314   58.352   63.867  90.00  90.00  90.00 P 21 21 21
SCALE1      0.026800  0.000000  0.000000        0.00000
SCALE2      0.000000  0.017137  0.000000        0.00000
SCALE3      0.000000  0.000000  0.015658        0.00000
ATOM      1  N  AARG A 168      22.726 -15.648  12.491  0.46 14.41           N
ATOM      2  CA AARG A 168      22.811 -14.357  13.143  0.46 15.09           C
ATOM      3  C  AARG A 168      21.569 -13.549  12.811  0.46 13.95           C
ATOM      4  O  AARG A 168      20.857 -13.814  11.839  0.46 14.06           O
ATOM      5  CB AARG A 168      24.057 -13.582  12.689  0.46 17.54           C
ATOM      6  CG AARG A 168      25.343 -14.366  12.775  0.46 19.91           C
ATOM      7  CD AARG A 168      25.613 -14.857  14.174  0.46 20.97           C
ATOM      8  NE AARG A 168      26.737 -15.790  14.190  0.46 21.81           N
ATOM      9  CZ AARG A 168      26.658 -17.059  14.578  0.46 21.76           C
ATOM     10  NH1AARG A 168      27.742 -17.813  14.539  0.46 21.87           N
ATOM     11  NH2AARG A 168      25.521 -17.587  15.017  0.46 21.76           N
ATOM     12  H  AARG A 168      22.668 -15.599  11.475  0.46 14.15           H
ATOM     13  HA AARG A 168      22.844 -14.490  14.224  0.46 15.42           H
ATOM     14  HB2AARG A 168      23.915 -13.295  11.650  0.46 17.96           H
ATOM     15  HB3AARG A 168      24.174 -12.692  13.305  0.46 17.77           H
ATOM     16  HD2AARG A 168      25.863 -14.006  14.804  0.46 21.25           H
ATOM     17  HD3AARG A 168      24.731 -15.354  14.564  0.46 21.33           H
ATOM     18  HE AARG A 168      27.476 -15.591  13.512  0.46 21.85           H
ATOM     19  HG2AARG A 168      25.288 -15.225  12.102  0.46 20.44           H
ATOM     20  HG3AARG A 168      26.173 -13.708  12.500  0.46 20.06           H
ATOM     21 HH11AARG A 168      28.617 -17.413  14.207  0.46 22.11           H
ATOM     22 HH12AARG A 168      27.696 -18.788  14.835  0.46 21.99           H
ATOM     23 HH21AARG A 168      25.514 -18.567  15.308  0.46 21.74           H
ATOM     24 HH22AARG A 168      24.658 -17.044  15.059  0.46 21.55           H
ATOM     25  N  BARG A 168      22.711 -15.649  12.474  0.54 14.14           N
ATOM     26  CA BARG A 168      22.916 -14.337  13.062  0.54 14.70           C
ATOM     27  C  BARG A 168      21.754 -13.431  12.679  0.54 13.29           C
ATOM     28  O  BARG A 168      21.184 -13.538  11.588  0.54 13.90           O
ATOM     29  CB BARG A 168      24.218 -13.716  12.533  0.54 17.49           C
ATOM     30  CG BARG A 168      25.439 -14.599  12.700  0.54 20.99           C
ATOM     31  CD BARG A 168      26.144 -14.335  14.005  0.54 22.91           C
ATOM     32  NE BARG A 168      27.298 -15.226  14.175  0.54 24.26           N
ATOM     33  CZ BARG A 168      28.570 -14.887  13.964  0.54 24.87           C
ATOM     34  NH1BARG A 168      29.519 -15.802  14.139  0.54 24.92           N
ATOM     35  NH2BARG A 168      28.905 -13.657  13.596  0.54 24.77           N
ATOM     36  H  BARG A 168      22.563 -15.624  11.467  0.54 13.85           H
ATOM     37  HA BARG A 168      22.959 -14.410  14.146  0.54 15.03           H
ATOM     38  HB2BARG A 168      24.096 -13.533  11.468  0.54 17.99           H
ATOM     39  HB3BARG A 168      24.418 -12.777  13.046  0.54 17.29           H
ATOM     40  HD2BARG A 168      26.474 -13.296  14.007  0.54 23.34           H
ATOM     41  HD3BARG A 168      25.450 -14.509  14.828  0.54 22.84           H
ATOM     42  HE BARG A 168      27.141 -16.106  14.675  0.54 24.26           H
ATOM     43  HG2BARG A 168      25.143 -15.649  12.671  0.54 21.23           H
ATOM     44  HG3BARG A 168      26.149 -14.391  11.896  0.54 21.33           H
ATOM     45 HH11BARG A 168      29.265 -16.742  14.428  0.54 24.77           H
ATOM     46 HH12BARG A 168      30.500 -15.563  13.994  0.54 25.10           H
ATOM     47 HH21BARG A 168      29.887 -13.425  13.448  0.54 24.89           H
ATOM     48 HH22BARG A 168      28.193 -12.942  13.455  0.54 24.59           H
TER
END
"""

def _run_and_get_target(pdb_str):
  result = None
  file_name_prefix = "tst_pdb_interpretation_3"
  with open("%s.pdb"%file_name_prefix, "w") as fo:
    fo.write(pdb_str)
  assert not easy_run.call("phenix.pdb_interpretation %s.pdb > %s.log"%(
    file_name_prefix,file_name_prefix))
  with open("%s.log"%file_name_prefix, "r") as fo:
    for l in fo.readlines():
      if l.startswith("  target:"):
        l = float(l.strip().split()[1])
        result = l
        break
  assert result is not None
  return result

def run():
  r1 = _run_and_get_target(pdb_6exy_ARG_168_D)
  r2 = _run_and_get_target(pdb_6exy_ARG_168_H)
  r3 = _run_and_get_target(pdb_6exy_ARG_168_HD)
  r4 = _run_and_get_target(pdb_6exy_ARG_168_ac_D)
  r5 = _run_and_get_target(pdb_6exy_ARG_168_ac_H)
  assert approx_equal(r1,r2, 1.e-3)
  assert approx_equal(r1,675.542, 1.e-1)
  assert approx_equal(r4,r5, 1.e-3)
  assert approx_equal(r4,1332.54, 1.e-1)
  assert approx_equal(r3,1340.87, 1.e-1)

if (__name__ == "__main__"):
  run()
  print("OK")
