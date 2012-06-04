from libtbx import easy_run
import sys

pdb1exr_fragment = """\
CRYST1   25.015   29.415   52.761  89.54  86.10  82.39 P 1           1
ATOM     43  N   GLU A   7      53.533  13.702  12.116  1.00 20.13           N
ATOM     44  CA  GLU A   7      53.762  12.263  12.038  1.00 18.74           C
ATOM     45  C   GLU A   7      52.498  11.484  11.693  1.00 16.59           C
ATOM     46  O   GLU A   7      52.569  10.583  10.885  1.00 17.01           O
ATOM     47  CB  GLU A   7      54.372  11.761  13.341  1.00 20.40           C
ATOM     48  CG  GLU A   7      54.593  10.273  13.376  1.00 25.03           C
ATOM     49  CD  GLU A   7      55.520   9.719  12.310  1.00 28.18           C
ATOM     50  OE1 GLU A   7      56.482  10.353  11.807  1.00 41.67           O
ATOM     51  OE2 GLU A   7      55.305   8.552  11.930  1.00 33.78           O
ATOM     52  N   GLN A   8      51.345  11.827  12.274  1.00 16.10           N
ATOM     53  CA  GLN A   8      50.125  11.194  11.846  1.00 15.79           C
ATOM     54  C   GLN A   8      49.860  11.431  10.346  1.00 14.53           C
ATOM     55  O   GLN A   8      49.477  10.500   9.631  1.00 14.74           O
ATOM     56  CB  GLN A   8      48.901  11.729  12.636  1.00 17.49           C
ATOM     57  CG  GLN A   8      48.898  11.267  14.104  1.00 19.61           C
ATOM     58  CD  GLN A   8      47.721  11.984  14.800  1.00 22.58           C
ATOM     59  OE1 GLN A   8      47.751  13.197  14.980  1.00 23.84           O
ATOM     60  NE2 GLN A   8      46.715  11.193  15.121  1.00 25.44           N
ATOM     61  N  AILE A   9      50.042  12.619   9.863  0.60 15.47           N
ATOM     62  N  BILE A   9      50.064  12.669   9.955  0.40 15.79           N
ATOM     63  CA AILE A   9      49.760  12.880   8.443  0.60 16.16           C
ATOM     64  CA BILE A   9      49.894  13.068   8.557  0.40 14.97           C
ATOM     65  C  AILE A   9      50.740  12.096   7.568  0.60 16.68           C
ATOM     66  C  BILE A   9      50.739  12.182   7.629  0.40 15.01           C
ATOM     67  O  AILE A   9      50.344  11.527   6.549  0.60 13.30           O
ATOM     68  O  BILE A   9      50.215  11.661   6.634  0.40 13.83           O
ATOM     69  CB AILE A   9      49.712  14.388   8.233  0.60 17.94           C
ATOM     70  CB BILE A   9      50.201  14.546   8.306  0.40 15.90           C
ATOM     71  CG1AILE A   9      48.529  15.065   8.881  0.60 19.71           C
ATOM     72  CG1BILE A   9      49.111  15.466   8.868  0.40 14.83           C
ATOM     73  CG2AILE A   9      49.644  14.632   6.712  0.60 18.98           C
ATOM     74  CG2BILE A   9      50.339  14.857   6.814  0.40 13.92           C
ATOM     75  CD1AILE A   9      48.524  16.563   9.024  0.60 23.19           C
ATOM     76  CD1BILE A   9      49.568  16.907   8.941  0.40 21.29           C
ATOM     77  N   ALA A  10      52.006  12.040   7.976  1.00 16.22           N
ATOM     78  CA  ALA A  10      52.970  11.271   7.207  1.00 15.86           C
ATOM     79  C   ALA A  10      52.525   9.824   7.134  1.00 14.08           C
ATOM     80  O   ALA A  10      52.586   9.200   6.066  1.00 14.53           O
ATOM     81  CB  ALA A  10      54.346  11.404   7.820  1.00 18.40           C
HETATM 1515  O  AHOH  2041      44.207   4.663 -11.870  0.33  8.04           O
HETATM 1516  O  BHOH  2042      44.917   5.551 -12.815  0.33 11.39           O
HETATM 1517  O  CHOH  2043      43.481   3.315 -12.617  0.33 10.33           O
END
"""

def exercise(args):
  assert len(args) == 0
  open("tmp_altloc_chain_break.pdb", "w").write(pdb1exr_fragment)
  command = " ".join([
    "phenix.pdbtools",
    "tmp_altloc_chain_break.pdb",
    "--geometry_regularization",
    "pdb_interpretation.nonbonded_weight=16",
    "geometry_minimization.macro_cycles=2"])
  gm_out = easy_run.fully_buffered(command=command) \
    .raise_if_errors() \
    .stdout_lines
  target_values = []
  for line in gm_out:
    if (line.startswith("target: ")):
      target_values.append(float(line.split()[-1]))
  assert len(target_values) == 4
  assert target_values[0] > 100
  assert target_values[3] < 1
  print "OK"

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
