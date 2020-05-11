from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run
from six.moves import range

pdbs = {
  'linking_BDP.pdb' : '''
CRYST1   30.794   82.051   32.278  90.00 117.89  90.00 P 1 21 1
SCALE1      0.032474  0.000000  0.017187        0.00000
SCALE2      0.000000  0.012188  0.000000        0.00000
SCALE3      0.000000  0.000000  0.035052        0.00000
ATOM      1  NZ  LYS A  76     -10.288 -11.233  -9.001  1.00 11.47           N
ATOM      2  HD3 LYS A  76      -8.299 -12.103  -7.539  1.00 13.18           H
ATOM      3  HZ1 LYS A  76     -10.482 -12.102  -8.722  0.00 11.31           H
ATOM      4  HZ3 LYS A  76      -9.860 -11.203  -9.840  0.00 11.31           H
ATOM      5  OH  TYR A  83       2.764 -11.498  13.475  1.00 14.89           O
ATOM      6  HH  TYR A  83       3.087 -11.482  14.250  0.00 12.06           H
ATOM      7  N   ALA A 102       3.459 -14.630  10.397  1.00  9.86           N
ATOM      8  H   ALA A 102       3.898 -13.982  10.750  1.00  9.75           H
ATOM      9  H   ALA A 103       3.966 -16.501  11.919  1.00 11.43           H
ATOM     10  CD  GLU A 131      -8.631  -9.699 -12.199  1.00 11.23           C
ATOM     11  OE2 GLU A 131      -8.579 -10.186 -11.066  1.00 15.02           O
ATOM     12  HG3 GLU A 131      -7.187  -9.899 -13.559  1.00  8.96           H
TER
HETATM   13  C3  NAG A1175       3.340 -15.138  17.816  1.00 18.88           C
HETATM   14  C7  NAG A1175       5.232 -15.904  20.454  1.00 23.06           C
HETATM   15  O3  NAG A1175       4.706 -15.282  17.480  1.00 19.56           O
HETATM   16  O4  NAG A1175       2.467 -14.929  15.523  1.00 18.37           O
HETATM   17  O7  NAG A1175       5.969 -15.238  21.198  1.00 27.13           O
HETATM   18  H3  NAG A1175       3.124 -14.101  18.073  1.00 19.18           H
HETATM   19  H82 NAG A1175       6.014 -16.957  18.773  1.00 21.76           H
HETATM   20  HO4 NAG A1175       3.094 -14.171  15.694  0.00 20.13           H
HETATM   21  C1  BDP A1176       5.315 -14.207  16.731  1.00 18.20           C
HETATM   22  C2  BDP A1176       6.818 -14.233  17.016  1.00 18.52           C
HETATM   23  C3  BDP A1176       7.628 -13.300  16.107  1.00 17.62           C
HETATM   24  C4  BDP A1176       7.296 -13.523  14.622  1.00 15.67           C
HETATM   25  C5  BDP A1176       5.772 -13.482  14.483  1.00 15.86           C
HETATM   26  C6  BDP A1176       5.352 -13.910  13.096  1.00 14.73           C
HETATM   27  O2  BDP A1176       7.014 -13.855  18.368  1.00 21.01           O
HETATM   28  O3  BDP A1176       9.012 -13.490  16.376  1.00 18.59           O
HETATM   29  O4  BDP A1176       7.901 -12.567  13.765  1.00 13.68           O
HETATM   30  O5  BDP A1176       5.100 -14.402  15.349  1.00 17.16           O
HETATM   31  O6A BDP A1176       5.577 -15.096  12.735  1.00 15.20           O
HETATM   32  O6B BDP A1176       4.762 -13.054  12.384  1.00 15.71           O
HETATM   33  H1  BDP A1176       4.902 -13.247  17.038  1.00 18.05           H
HETATM   34  H2  BDP A1176       7.167 -15.256  16.878  1.00 18.76           H
HETATM   35  H3  BDP A1176       7.358 -12.285  16.370  1.00 17.48           H
HETATM   36  HB  BDP A1176       6.144 -13.679  18.812  0.00 18.97           H
HETATM   37  H4  BDP A1176       7.635 -14.523  14.351  1.00 15.66           H
HETATM   38  H5  BDP A1176       5.414 -12.471  14.677  1.00 15.52           H
HETATM   39  HC  BDP A1176       9.118 -14.232  17.068  0.00 17.73           H
HETATM   40  C1  NAG A1177       8.992 -12.956  12.918  1.00 13.58           C
HETATM   41  C2  NAG A1177       9.213 -11.887  11.849  1.00 12.17           C
HETATM   42  N2  NAG A1177       7.954 -11.763  11.092  1.00 12.56           N
HETATM   43  O5  NAG A1177      10.174 -13.136  13.662  1.00 14.08           O
HETATM   44  H1  NAG A1177       8.740 -13.896  12.428  1.00 13.21           H
HETATM   45  H2  NAG A1177       9.424 -10.949  12.357  1.00 12.00           H
HETATM   46  H82 NAG A1177       5.243 -11.105  10.751  1.00 10.52           H
HETATM   47  HN2 NAG A1177       7.399 -12.601  11.034  1.00 11.60           H
''',
        }

links = {
  'linking_BDP.pdb' : '''
  Links applied
    BETA1-3
      " NAG A1175 " - " BDP A1176 "
    BETA1-4
      " BDP A1176 " - " NAG A1177 "''',
  }

def run(only_i=None):
  try: only_i=int(only_i)
  except ValueError: only_i=None
  except TypeError: only_i=None
  cifs = ""
  for pdb in pdbs:
    f=open(pdb, "w")
    f.write(pdbs[pdb])
    f.close()
    if pdb.endswith(".cif"): cifs += " %s" % pdb
    links[pdb] = links[pdb].splitlines()
  j=0
  for pdb in sorted(pdbs):
    #break
    if pdb.endswith(".cif"): continue
    if pdb.endswith(".params"): continue
    print('pdb',pdb)
    j+=1
    if only_i is not None and only_i!=j: continue
    # log_filename = "%s.log" % (pdb)
    cmd = "phenix.pdb_interpretation %s write_geo_file=False" % pdb
    result = easy_run.fully_buffered(cmd).raise_if_errors()
    assert (result.return_code == 0)
    for line in result.stdout_lines:
      if (line in links[pdb]):
        links[pdb].remove(line)
    if links[pdb]:
      raise RuntimeError("Missing expected log output %s" % pdb)
    print("OK")

if __name__=="__main__":
  import sys
  run(*tuple(sys.argv[1:]))
