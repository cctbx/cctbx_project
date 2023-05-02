from __future__ import absolute_import, division, print_function
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import format_cpu_times
import sys

params = monomer_library.pdb_interpretation.master_params.extract()
params.flip_symmetric_amino_acids = False

'''
Once this test is fixed, please move it to tst_pdb_interpretation.
The test is in a separate file to avoid accumulating failures in tst_pdb_interpretation
'''

def exercise_HOP2(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1   19.545   17.217   16.395  90.00  90.00  90.00 P 1
ATOM      1  P    DA A   7      12.372   5.356  10.324  1.00 21.25           P
ATOM      2  OP1  DA A   7      13.776   5.000  10.699  1.00 24.88           O
ATOM      3  OP2  DA A   7      11.412   5.668  11.395  1.00 21.95           O
ATOM      4  O5'  DA A   7      12.418   6.531   9.251  1.00 18.66           O
ATOM      5  C5'  DA A   7      13.460   6.584   8.265  1.00 19.04           C
ATOM      6  C4'  DA A   7      13.188   7.694   7.275  1.00 17.77           C
ATOM      7  O4'  DA A   7      11.855   7.547   6.729  1.00 18.07           O
ATOM      8  C3'  DA A   7      13.197   9.107   7.852  1.00 19.00           C
ATOM      9  O3'  DA A   7      14.545   9.611   7.912  1.00 18.59           O
ATOM     10  C2'  DA A   7      12.361   9.876   6.841  1.00 17.51           C
ATOM     11  C1'  DA A   7      11.311   8.840   6.438  1.00 17.34           C
ATOM     12  N9   DA A   7      10.040   8.980   7.158  1.00 17.04           N
ATOM     13  C8   DA A   7       9.563   8.289   8.257  1.00 15.09           C
ATOM     14  N7   DA A   7       8.347   8.641   8.609  1.00 13.14           N
ATOM     15  C5   DA A   7       8.018   9.653   7.711  1.00 14.61           C
ATOM     16  C6   DA A   7       6.875  10.442   7.551  1.00 14.42           C
ATOM     17  N6   DA A   7       5.792  10.332   8.309  1.00 13.77           N
ATOM     18  N1   DA A   7       6.875  11.362   6.554  1.00 14.01           N
ATOM     19  C2   DA A   7       7.954  11.462   5.782  1.00 14.66           C
ATOM     20  N3   DA A   7       9.083  10.774   5.830  1.00 14.29           N
ATOM     21  C4   DA A   7       9.050   9.873   6.822  1.00 15.47           C
ATOM     22  H5'  DA A   7      13.512   5.632   7.731  1.00 17.63           H
ATOM     23 H5''  DA A   7      14.417   6.770   8.760  1.00 18.60           H
ATOM     24  H4'  DA A   7      13.973   7.663   6.511  1.00 18.30           H
ATOM     25  H3'  DA A   7      12.753   9.126   8.850  1.00 17.97           H
ATOM     26  H2'  DA A   7      11.917  10.765   7.288  1.00 18.28           H
ATOM     27 H2''  DA A   7      12.977  10.220   6.005  1.00 16.30           H
ATOM     28  H1'  DA A   7      11.113   8.874   5.362  1.00 17.98           H
ATOM     29  H8   DA A   7      10.146   7.552   8.798  1.00 16.21           H
ATOM     30  H2   DA A   7       7.905  12.217   5.000  1.00 11.67           H
ATOM     31  D61A DA A   7       5.787   9.643   9.039  0.67 12.16           D
ATOM     32  D62A DA A   7       5.000  10.943   8.131  0.89 13.67           D
ATOM     33 DOP2A DA A   7      10.894   6.398  11.013  0.79 22.75           D
ATOM     34  H61B DA A   7       5.787   9.643   9.039  0.33 12.16           H
ATOM     35  H62B DA A   7       5.000  10.943   8.131  0.11 13.67           H
ATOM     36 HOP2B DA A   7      10.894   6.398  11.013  0.21 22.75           H
END
""".splitlines()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    params=params,
    raw_records=raw_records,
    force_symmetry=True)

def run(args):
  assert len(args) == 0
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  if 0: exercise_HOP2(mon_lib_srv, ener_lib)
  else: print('skipping until complete conversion to v3.2 atom names')
  print(format_cpu_times())

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
