
from __future__ import absolute_import, division, print_function
from mmtbx.disorder import analyze_model
from mmtbx.validation import molprobity
import iotbx.pdb
from libtbx.utils import null_out
import libtbx.load_env
import mmtbx.model
from six.moves import cStringIO as StringIO

def analyze_fragment(pdb_str):
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(pdb_in)
  validation = molprobity.molprobity(model, outliers_only=False)
  result = analyze_model.process_pdb_hierarchy(
    pdb_hierarchy=model.get_hierarchy(),
    validation=validation,
    log=null_out())
  return result

def exercise():
  # single alternate, split at C-alpha (from 1ytt)
  pdb_str = """\
CRYST1   65.508   72.216   45.035  90.00  90.00  90.00 P 21 21 21    8
ATOM    844  N   VAL A 216      11.617  50.561  30.170  1.00  5.99           N
ATOM    845  CA  VAL A 216      12.023  51.936  29.899  1.00  6.80           C
ATOM    846  C   VAL A 216      11.349  52.370  28.604  1.00  8.26           C
ATOM    847  O   VAL A 216      11.479  51.707  27.569  1.00  6.98           O
ATOM    848  CB  VAL A 216      13.550  52.097  29.765  1.00  6.32           C
ATOM    849  CG1 VAL A 216      13.899  53.570  29.583  1.00  6.65           C
ATOM    850  CG2 VAL A 216      14.248  51.530  31.002  1.00  6.72           C
ATOM    851  N   CYS A 217      10.629  53.485  28.685  1.00 10.22           N
ATOM    852  CA ACYS A 217       9.916  54.032  27.534  0.54 10.72           C
ATOM    853  CA BCYS A 217       9.883  54.049  27.565  0.45 11.26           C
ATOM    854  C   CYS A 217      10.525  55.369  27.120  1.00 11.35           C
ATOM    855  O   CYS A 217      11.061  56.102  27.951  1.00 11.69           O
ATOM    856  CB ACYS A 217       8.430  54.204  27.849  0.54 11.93           C
ATOM    857  CB BCYS A 217       8.437  54.279  28.023  0.45 12.95           C
ATOM    858  SG ACYS A 217       7.621  52.711  28.510  0.54 12.38           S
ATOM    859  SG BCYS A 217       7.253  54.846  26.764  0.45 16.59           S
ATOM    860  N   GLU A 218      10.471  55.668  25.825  1.00 11.30           N
ATOM    861  CA  GLU A 218      11.046  56.908  25.302  1.00 12.00           C
ATOM    862  C   GLU A 218       9.980  57.809  24.705  1.00 12.95           C
ATOM    863  O   GLU A 218       8.952  57.329  24.209  1.00 12.30           O
ATOM    864  CB  GLU A 218      12.127  56.628  24.248  1.00 13.87           C
ATOM    865  CG  GLU A 218      11.638  55.952  22.978  1.00 15.94           C
ATOM    866  CD  GLU A 218      12.723  55.827  21.912  1.00 17.80           C
ATOM    867  OE1 GLU A 218      13.918  55.750  22.262  1.00 17.06           O
ATOM    868  OE2 GLU A 218      12.382  55.827  20.711  1.00 20.67           O
END
"""
  result = analyze_fragment(pdb_str)
  out = StringIO()
  result.show(out=out)
  assert ("atom_group=B CYS  occ=0.45 phi=-111.6 psi=146.7  rot=t" in
    out.getvalue())

if (__name__ == "__main__"):
  if (not libtbx.env.has_module("probe")):
    print("Probe not configured, skipping test")
  else :
    exercise()
    print("OK")
