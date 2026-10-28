"""reduce2 with qmready=True stops with a Sorry when the finished model has a
non-bonded contact below 1.0 A, H included. Model: 1yjp Asn 3 and Gln 4 with
C3-N4 = 3.1 A, too long for a peptide bond, so qmready adds OXT to Asn 3 and
NH3 to Gln 4; OXT ASN A 3 - H3 GLN A 4 = 0.892 A.
"""
from __future__ import absolute_import, division, print_function
import os
from libtbx.utils import null_out, Sorry
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

pdb_str_gap = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM     13  N   ASN A   3      23.342  28.731  29.139  1.00 12.26           N
ATOM     14  CA  ASN A   3      24.587  29.045  29.823  1.00 11.74           C
ATOM     15  C   ASN A   3      25.825  28.473  29.129  1.00 11.10           C
ATOM     16  O   ASN A   3      25.908  27.260  28.882  1.00 10.42           O
ATOM     17  CB  ASN A   3      24.521  28.519  31.276  1.00 12.15           C
ATOM     18  CG  ASN A   3      25.774  28.880  32.095  1.00 12.82           C
ATOM     19  OD1 ASN A   3      26.078  30.066  32.306  1.00 15.05           O
ATOM     20  ND2 ASN A   3      26.509  27.856  32.540  1.00 13.48           N
ATOM     21  N   GLN A   4      28.025  30.547  28.442  1.00 10.29           N
ATOM     22  CA  GLN A   4      29.414  30.207  28.043  1.00 10.53           C
ATOM     23  C   GLN A   4      30.465  30.925  28.932  1.00 10.24           C
ATOM     24  O   GLN A   4      30.577  32.162  28.959  1.00  8.86           O
ATOM     25  CB  GLN A   4      29.686  30.467  26.555  1.00  9.80           C
ATOM     26  CG  GLN A   4      30.974  29.777  26.057  1.00 10.25           C
ATOM     27  CD  GLN A   4      31.534  30.363  24.755  1.00 12.43           C
ATOM     28  OE1 GLN A   4      31.774  31.587  24.654  1.00 14.62           O
ATOM     29  NE2 GLN A   4      31.780  29.480  23.753  1.00  9.05           N
END
"""

# HEAD's qmready add output for pdb_str_gap (before the check saw H)
pdb_str_gap_h = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A   3      23.342  28.731  29.139  1.00 12.26           N
ATOM      2  CA  ASN A   3      24.587  29.045  29.823  1.00 11.74           C
ATOM      3  C   ASN A   3      25.825  28.473  29.129  1.00 11.10           C
ATOM      4  O   ASN A   3      25.908  27.260  28.882  1.00 10.42           O
ATOM      5  CB  ASN A   3      24.521  28.519  31.276  1.00 12.15           C
ATOM      6  CG  ASN A   3      25.774  28.880  32.095  1.00 12.82           C
ATOM      7  OD1 ASN A   3      26.078  30.066  32.306  1.00 15.05           O
ATOM      8  ND2 ASN A   3      26.509  27.856  32.540  1.00 13.48           N
ATOM      9  OXT ASN A   3      26.748  29.224  28.813  1.00 10.42           O
ATOM     10  H1  ASN A   3      23.385  29.084  28.162  1.00 12.26           H
ATOM     11  H2  ASN A   3      23.202  27.701  29.129  1.00 12.26           H
ATOM     12  H3  ASN A   3      22.549  29.184  29.637  1.00 12.26           H
ATOM     13  HA  ASN A   3      24.706  30.128  29.820  1.00 11.74           H
ATOM     14  HB2 ASN A   3      23.655  28.956  31.773  1.00 12.15           H
ATOM     15  HB3 ASN A   3      24.431  27.433  31.257  1.00 12.15           H
ATOM     16 HD21 ASN A   3      27.351  28.033  33.088  1.00 13.48           H
ATOM     17 HD22 ASN A   3      26.228  26.898  32.331  1.00 13.48           H
ATOM     18  N   GLN A   4      28.025  30.547  28.442  1.00 10.29           N
ATOM     19  CA  GLN A   4      29.414  30.207  28.043  1.00 10.53           C
ATOM     20  C   GLN A   4      30.465  30.925  28.932  1.00 10.24           C
ATOM     21  O   GLN A   4      30.577  32.162  28.959  1.00  8.86           O
ATOM     22  CB  GLN A   4      29.686  30.467  26.555  1.00  9.80           C
ATOM     23  CG  GLN A   4      30.974  29.777  26.057  1.00 10.25           C
ATOM     24  CD  GLN A   4      31.534  30.363  24.755  1.00 12.43           C
ATOM     25  OE1 GLN A   4      31.774  31.587  24.654  1.00 14.62           O
ATOM     26  NE2 GLN A   4      31.780  29.480  23.753  1.00  9.05           N
ATOM     27  OXT GLN A   4      31.213  30.245  29.634  1.00  8.86           O
ATOM     28  H1  GLN A   4      28.045  31.150  29.289  1.00 10.29           H
ATOM     29  H2  GLN A   4      27.556  31.056  27.666  1.00 10.29           H
ATOM     30  H3  GLN A   4      27.501  29.673  28.651  1.00 10.29           H
ATOM     31  HA  GLN A   4      29.535  29.134  28.190  1.00 10.53           H
ATOM     32  HB2 GLN A   4      28.850  30.085  25.969  1.00  9.80           H
ATOM     33  HB3 GLN A   4      29.794  31.540  26.396  1.00  9.80           H
ATOM     34  HG2 GLN A   4      31.743  29.878  26.823  1.00 10.25           H
ATOM     35  HG3 GLN A   4      30.760  28.723  25.880  1.00 10.25           H
ATOM     36 HE21 GLN A   4      32.153  29.810  22.863  1.00  9.05           H
ATOM     37 HE22 GLN A   4      31.592  28.487  23.891  1.00  9.05           H
END
"""

clash = ('pdb=" OXT ASN A   3 " - pdb=" H3  GLN A   4 " 0.892 A')

def run(model_text, prefix, extra_args=()):
  fn_in = prefix + ".pdb"
  fn_out = prefix + "_out.pdb"
  with open(fn_in, "w") as f:
    f.write(model_text)
  if os.path.exists(fn_out): os.remove(fn_out)
  args = [fn_in, "output.filename=" + fn_out, "overwrite=True"]
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=args + list(extra_args))
  return fn_out

def check_sorry(model_text, prefix, extra_args):
  try:
    fn_out = run(model_text, prefix, extra_args)
  except Sorry as e:
    assert "Bad clash(es) prevents qmready=True" in str(e), str(e)
    assert clash in str(e), str(e)
    assert not os.path.exists(prefix + "_out.pdb")
  else:
    raise AssertionError("no Sorry, %s written" % fn_out)

def test_add():
  check_sorry(pdb_str_gap, "tst_reduce2_qmready_clash_1", ["qmready=True"])

def test_optimize():
  check_sorry(pdb_str_gap_h, "tst_reduce2_qmready_clash_2",
    ["qmready=True", "approach=optimize"])

def test_add_no_qmready():
  # no termini added, no clash, no check
  fn_out = run(pdb_str_gap, "tst_reduce2_qmready_clash_3")
  assert os.path.isfile(fn_out)

if __name__ == "__main__":
  test_add()
  test_optimize()
  test_add_no_qmready()
  print("OK")
