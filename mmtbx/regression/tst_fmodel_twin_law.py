
from __future__ import division
from libtbx import easy_run
import os.path

def exercise (file_base="tmp_fmodel_twin_law") :
  open("%s.pdb" % file_base, "w").write("""\
CRYST1   12.000    5.000   12.000  90.00  90.00  90.00 P 1           1
ATOM     39  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
""")
  params = """
pdb_file = %s.pdb
high_resolution = 1.5
add_sigmas = True
twin_law = l,-k,h
twin_fraction = 0.4
output {
  label = I
  type = *real complex
  file_name = %s.mtz
}
""" % (file_base, file_base)
  open("%s.eff" % file_base, "w").write(params)
  result = easy_run.fully_buffered(
    "phenix.fmodel %s.eff" % file_base).raise_if_errors()
  assert ("""Artifically twinning the data with fraction 0.40""" in
          result.stdout_lines)
  assert ("""using twin law (l,-k,h)""" in result.stdout_lines)
  assert os.path.isfile("%s.mtz" % file_base)
  return ("%s.pdb" % file_base, "%s.mtz" % file_base)

if (__name__ == "__main__") :
  exercise()
  print "OK"
