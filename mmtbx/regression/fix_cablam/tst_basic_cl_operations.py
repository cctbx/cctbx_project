from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run
import libtbx.load_env

pdb_str = """\
CRYST1   32.501   39.502   44.640  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   ASN A   1      12.388  -5.203  29.298  1.00 13.23           N
ATOM      2  CA  ASN A   1      12.398  -3.931  28.500  1.00  6.85           C
ATOM      3  C   ASN A   1      12.384  -4.463  27.091  1.00  6.81           C
ATOM      5  O   ASN A   1      12.645  -5.635  26.750  1.00  8.28           O
ATOM      7  CB  ASN A   1      13.493  -3.033  29.043  1.00  7.32           C
ATOM      8  CG  ASN A   1      13.414  -1.559  28.658  1.00  7.93           C
ATOM      9  OD1 ASN A   1      14.388  -0.782  28.842  1.00 10.99           O
ATOM     10  ND2 ASN A   1      12.315  -1.201  28.134  1.00 11.12           N
"""

def exercise_out(prefix="tst_fix_cablam_basic_cl_operations_pdb_out"):
  in_fname = '%s_input.pdb' % prefix
  with open(in_fname, 'w') as f:
    f.write(pdb_str)
  cmd = " ".join([
      "phenix.cablam_idealization",
      in_fname,
      ])
  print(cmd)
  assert not easy_run.call(cmd)
  outfname = prefix + "_input_cablam_fixed.pdb"
  print("checking %s" % (outfname))
  assert os.path.isfile(outfname)


if __name__ == '__main__':
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping: probe not configured")
  else:
    exercise_out()
