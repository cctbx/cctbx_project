from __future__ import absolute_import, division, print_function

import iotbx.pdb
import mmtbx.model
import libtbx.load_env
from libtbx.test_utils import show_diff, assert_lines_in_text
from six.moves import cStringIO as StringIO

pdb_str = """\
ATOM      1  N   GLY A   1      12.928   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      12.885   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      13.922   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      14.414   2.521   5.381  1.00 16.78           O
ATOM      5  H1  GLY A   1      11.977   4.525   6.514  0.00 16.77           H
ATOM      6  H2  GLY A   1      13.586   3.992   6.616  0.00 16.77           H
ATOM      7  H3  GLY A   1      13.250   5.598   6.177  0.00 16.77           H
ATOM      8  HA2 GLY A   1      11.900   3.815   4.398  0.00 16.57           H
ATOM      9  HA3 GLY A   1      13.100   5.065   4.014  0.00 16.57           H
ATOM     10  N   ASN A   2      14.281   2.923   3.155  1.00 15.02           N
ATOM     11  CA  ASN A   2      15.415   2.038   2.831  1.00 14.10           C
ATOM     12  C   ASN A   2      16.696   2.537   3.427  1.00 13.13           C
ATOM     13  O   ASN A   2      16.959   3.742   3.426  1.00 11.91           O
ATOM     14  CB  ASN A   2      15.591   1.881   1.341  1.00 15.38           C
ATOM     15  CG  ASN A   2      14.353   1.342   0.692  1.00 14.08           C
ATOM     16  OD1 ASN A   2      13.912   0.227   1.016  1.00 17.46           O
ATOM     17  ND2 ASN A   2      13.733   2.155  -0.169  1.00 11.72           N
ATOM     18  H   ASN A   2      13.820   3.334   2.343  0.00 15.02           H
ATOM     19  HA  ASN A   2      15.193   1.057   3.252  0.00 14.10           H
ATOM     20  HB2 ASN A   2      15.813   2.853   0.900  0.00 15.38           H
ATOM     21  HB3 ASN A   2      16.409   1.188   1.146  0.00 15.38           H
ATOM     22 HD21 ASN A   2      12.885   1.845  -0.643  0.00 11.72           H
ATOM     23 HD22 ASN A   2      14.108   3.086  -0.352  0.00 11.72           H
ATOM     24  N   ASN A   3      17.499   1.590   3.905  1.00 12.26           N
ATOM     25  CA  ASN A   3      18.744   1.904   4.589  1.00 11.74           C
ATOM     26  C   ASN A   3      19.982   1.332   3.895  1.00 11.10           C
ATOM     27  O   ASN A   3      20.065   0.119   3.648  1.00 10.42           O
ATOM     28  CB  ASN A   3      18.678   1.378   6.042  1.00 12.15           C
ATOM     29  CG  ASN A   3      19.931   1.739   6.861  1.00 12.82           C
ATOM     30  OD1 ASN A   3      20.235   2.925   7.072  1.00 15.05           O
ATOM     31  ND2 ASN A   3      20.666   0.715   7.306  1.00 13.48           N
ATOM     32  OXT ASN A   3      20.908   2.077   3.576  1.00 11.10           O
ATOM     33  H   ASN A   3      17.311   0.590   3.832  0.00 12.26           H
ATOM     34  HA  ASN A   3      18.863   2.987   4.586  0.00 11.74           H
ATOM     35  HB2 ASN A   3      17.812   1.815   6.539  0.00 12.15           H
ATOM     36  HB3 ASN A   3      18.588   0.292   6.023  0.00 12.15           H
ATOM     37 HD21 ASN A   3      21.508   0.892   7.854  0.00 13.48           H
ATOM     38 HD22 ASN A   3      20.385  -0.243   7.097  0.00 13.48           H
TER
ATOM     39  N   GLY B   1      12.928  -0.254   6.102  1.00 16.77           N
ATOM     40  CA  GLY B   1      12.885  -0.659   4.651  1.00 16.57           C
ATOM     41  C   GLY B   1      13.922  -1.726   4.419  1.00 16.16           C
ATOM     42  O   GLY B   1      14.414  -2.345   5.381  1.00 16.78           O
ATOM     43  H1  GLY B   1      11.977  -0.341   6.514  0.00 16.77           H
ATOM     44  H2  GLY B   1      13.586  -0.874   6.616  0.00 16.77           H
ATOM     45  H3  GLY B   1      13.250   0.732   6.177  0.00 16.77           H
ATOM     46  HA2 GLY B   1      11.900  -1.051   4.398  0.00 16.57           H
ATOM     47  HA3 GLY B   1      13.100   0.199   4.014  0.00 16.57           H
ATOM     48  N   ASN B   2      14.281  -1.943   3.155  1.00 15.02           N
ATOM     49  CA  ASN B   2      15.415  -2.828   2.831  1.00 14.10           C
ATOM     50  C   ASN B   2      16.696  -2.329   3.427  1.00 13.13           C
ATOM     51  O   ASN B   2      16.959  -1.124   3.426  1.00 11.91           O
ATOM     52  CB  ASN B   2      15.591  -2.985   1.341  1.00 15.38           C
ATOM     53  CG  ASN B   2      14.353  -3.524   0.692  1.00 14.08           C
ATOM     54  OD1 ASN B   2      13.912  -4.639   1.016  1.00 17.46           O
ATOM     55  ND2 ASN B   2      13.733  -2.711  -0.169  1.00 11.72           N
ATOM     56  H   ASN B   2      13.820  -1.532   2.343  0.00 15.02           H
ATOM     57  HA  ASN B   2      15.193  -3.809   3.252  0.00 14.10           H
ATOM     58  HB2 ASN B   2      15.813  -2.013   0.900  0.00 15.38           H
ATOM     59  HB3 ASN B   2      16.409  -3.678   1.146  0.00 15.38           H
ATOM     60 HD21 ASN B   2      12.885  -3.021  -0.643  0.00 11.72           H
ATOM     61 HD22 ASN B   2      14.108  -1.780  -0.352  0.00 11.72           H
ATOM     62  N   ASN B   3      17.499  -3.276   3.905  1.00 12.26           N
ATOM     63  CA  ASN B   3      18.744  -2.962   4.589  1.00 11.74           C
ATOM     64  C   ASN B   3      19.982  -3.534   3.895  1.00 11.10           C
ATOM     65  O   ASN B   3      20.065  -4.747   3.648  1.00 10.42           O
ATOM     66  CB  ASN B   3      18.678  -3.488   6.042  1.00 12.15           C
ATOM     67  CG  ASN B   3      19.931  -3.127   6.861  1.00 12.82           C
ATOM     68  OD1 ASN B   3      20.235  -1.941   7.072  1.00 15.05           O
ATOM     69  ND2 ASN B   3      20.666  -4.151   7.306  1.00 13.48           N
ATOM     70  OXT ASN B   3      20.908  -2.789   3.576  1.00 11.10           O
ATOM     71  H   ASN B   3      17.311  -4.276   3.832  0.00 12.26           H
ATOM     72  HA  ASN B   3      18.863  -1.879   4.586  0.00 11.74           H
ATOM     73  HB2 ASN B   3      17.812  -3.051   6.539  0.00 12.15           H
ATOM     74  HB3 ASN B   3      18.588  -4.574   6.023  0.00 12.15           H
ATOM     75 HD21 ASN B   3      21.508  -3.974   7.854  0.00 13.48           H
ATOM     76 HD22 ASN B   3      20.385  -5.109   7.097  0.00 13.48           H
TER
ATOM     77  N   GLY C   1      12.928   9.478   6.102  1.00 16.77           N
ATOM     78  CA  GLY C   1      12.885   9.073   4.651  1.00 16.57           C
ATOM     79  C   GLY C   1      13.922   8.006   4.419  1.00 16.16           C
ATOM     80  O   GLY C   1      14.414   7.387   5.381  1.00 16.78           O
ATOM     81  H1  GLY C   1      11.977   9.391   6.514  0.00 16.77           H
ATOM     82  H2  GLY C   1      13.586   8.858   6.616  0.00 16.77           H
ATOM     83  H3  GLY C   1      13.250  10.464   6.177  0.00 16.77           H
ATOM     84  HA2 GLY C   1      11.900   8.681   4.398  0.00 16.57           H
ATOM     85  HA3 GLY C   1      13.100   9.931   4.014  0.00 16.57           H
ATOM     86  N   ASN C   2      14.281   7.789   3.155  1.00 15.02           N
ATOM     87  CA  ASN C   2      15.415   6.904   2.831  1.00 14.10           C
ATOM     88  C   ASN C   2      16.696   7.403   3.427  1.00 13.13           C
ATOM     89  O   ASN C   2      16.959   8.608   3.426  1.00 11.91           O
ATOM     90  CB  ASN C   2      15.591   6.747   1.341  1.00 15.38           C
ATOM     91  CG  ASN C   2      14.353   6.208   0.692  1.00 14.08           C
ATOM     92  OD1 ASN C   2      13.912   5.093   1.016  1.00 17.46           O
ATOM     93  ND2 ASN C   2      13.733   7.021  -0.169  1.00 11.72           N
ATOM     94  H   ASN C   2      13.820   8.200   2.343  0.00 15.02           H
ATOM     95  HA  ASN C   2      15.193   5.923   3.252  0.00 14.10           H
ATOM     96  HB2 ASN C   2      15.813   7.719   0.900  0.00 15.38           H
ATOM     97  HB3 ASN C   2      16.409   6.054   1.146  0.00 15.38           H
ATOM     98 HD21 ASN C   2      12.885   6.711  -0.643  0.00 11.72           H
ATOM     99 HD22 ASN C   2      14.108   7.952  -0.352  0.00 11.72           H
ATOM    100  N   ASN C   3      17.499   6.456   3.905  1.00 12.26           N
ATOM    101  CA  ASN C   3      18.744   6.770   4.589  1.00 11.74           C
ATOM    102  C   ASN C   3      19.982   6.198   3.895  1.00 11.10           C
ATOM    103  O   ASN C   3      20.065   4.985   3.648  1.00 10.42           O
ATOM    104  CB  ASN C   3      18.678   6.244   6.042  1.00 12.15           C
ATOM    105  CG  ASN C   3      19.931   6.605   6.861  1.00 12.82           C
ATOM    106  OD1 ASN C   3      20.235   7.791   7.072  1.00 15.05           O
ATOM    107  ND2 ASN C   3      20.666   5.581   7.306  1.00 13.48           N
ATOM    108  OXT ASN C   3      20.908   6.943   3.576  1.00 11.10           O
ATOM    109  H   ASN C   3      17.311   5.456   3.832  0.00 12.26           H
ATOM    110  HA  ASN C   3      18.863   7.853   4.586  0.00 11.74           H
ATOM    111  HB2 ASN C   3      17.812   6.681   6.539  0.00 12.15           H
ATOM    112  HB3 ASN C   3      18.588   5.158   6.023  0.00 12.15           H
ATOM    113 HD21 ASN C   3      21.508   5.758   7.854  0.00 13.48           H
ATOM    114 HD22 ASN C   3      20.385   4.623   7.097  0.00 13.48           H
TER
ATOM    115  N   GLY D   1       9.009   2.179  -6.102  1.00 16.77           N
ATOM    116  CA  GLY D   1       9.052   1.774  -4.651  1.00 16.57           C
ATOM    117  C   GLY D   1       8.015   0.707  -4.419  1.00 16.16           C
ATOM    118  O   GLY D   1       7.523   0.088  -5.381  1.00 16.78           O
ATOM    119  H1  GLY D   1       9.960   2.092  -6.514  0.00 16.77           H
ATOM    120  H2  GLY D   1       8.351   1.559  -6.616  0.00 16.77           H
ATOM    121  H3  GLY D   1       8.687   3.165  -6.177  0.00 16.77           H
ATOM    122  HA2 GLY D   1      10.037   1.382  -4.398  0.00 16.57           H
ATOM    123  HA3 GLY D   1       8.837   2.632  -4.014  0.00 16.57           H
ATOM    124  N   ASN D   2       7.656   0.490  -3.155  1.00 15.02           N
ATOM    125  CA  ASN D   2       6.522  -0.395  -2.831  1.00 14.10           C
ATOM    126  C   ASN D   2       5.241   0.104  -3.427  1.00 13.13           C
ATOM    127  O   ASN D   2       4.978   1.309  -3.426  1.00 11.91           O
ATOM    128  CB  ASN D   2       6.346  -0.552  -1.341  1.00 15.38           C
ATOM    129  CG  ASN D   2       7.584  -1.091  -0.692  1.00 14.08           C
ATOM    130  OD1 ASN D   2       8.025  -2.206  -1.016  1.00 17.46           O
ATOM    131  ND2 ASN D   2       8.204  -0.278   0.169  1.00 11.72           N
ATOM    132  H   ASN D   2       8.117   0.901  -2.343  0.00 15.02           H
ATOM    133  HA  ASN D   2       6.744  -1.376  -3.252  0.00 14.10           H
ATOM    134  HB2 ASN D   2       6.124   0.420  -0.900  0.00 15.38           H
ATOM    135  HB3 ASN D   2       5.528  -1.245  -1.146  0.00 15.38           H
ATOM    136 HD21 ASN D   2       9.052  -0.588   0.643  0.00 11.72           H
ATOM    137 HD22 ASN D   2       7.829   0.653   0.352  0.00 11.72           H
TER
ATOM    138  N   ASN E   2       7.656   5.356  -3.155  1.00 15.02           N
ATOM    139  CA  ASN E   2       6.522   4.471  -2.831  1.00 14.10           C
ATOM    140  C   ASN E   2       5.241   4.970  -3.427  1.00 13.13           C
ATOM    141  O   ASN E   2       4.978   6.175  -3.426  1.00 11.91           O
ATOM    142  CB  ASN E   2       6.346   4.314  -1.341  1.00 15.38           C
ATOM    143  CG  ASN E   2       7.584   3.775  -0.692  1.00 14.08           C
ATOM    144  OD1 ASN E   2       8.025   2.660  -1.016  1.00 17.46           O
ATOM    145  ND2 ASN E   2       8.204   4.588   0.169  1.00 11.72           N
ATOM    146  H   ASN E   2       8.117   5.767  -2.343  0.00 15.02           H
ATOM    147  HA  ASN E   2       6.744   3.490  -3.252  0.00 14.10           H
ATOM    148  HB2 ASN E   2       6.124   5.286  -0.900  0.00 15.38           H
ATOM    149  HB3 ASN E   2       5.528   3.621  -1.146  0.00 15.38           H
ATOM    150 HD21 ASN E   2       9.052   4.278   0.643  0.00 11.72           H
ATOM    151 HD22 ASN E   2       7.829   5.519   0.352  0.00 11.72           H
TER
"""

def get_geometry_stats(lines, use_ncs=True):
  log_str=StringIO()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=lines)
  m = mmtbx.model.manager(model_input = pdb_inp, log = log_str)
  p = m.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_ncs_to_build_restraints = use_ncs
  m.process(make_restraints=True, pdb_interpretation_params=p)
  geom=StringIO()
  g = m.geometry_statistics()
  g.show(log=geom, exclude_protein_only_stats=True)
  return geom.getvalue(), log_str.getvalue()

def exercise_01():
  geom_ncs, log_ncs=get_geometry_stats(pdb_str, True)
  geom_no_ncs, log_no_ncs=get_geometry_stats(pdb_str, False)
  assert not show_diff(geom_ncs, geom_no_ncs)
  # assert not show_diff(log_ncs, log_no_ncs)
  assert_lines_in_text(log_ncs, """ Restraints were copied for chains:
    B, C
""")
  # print(log_ncs)

if(__name__ == "__main__"):
  if libtbx.env.find_in_repositories(relative_path="chem_data") is None:
    print("Skipping exercise_01(): chem_data directory not available")
  else:
    exercise_01()
    print('OK')
