from __future__ import absolute_import, division, print_function
import io, sys
from math import cos, sin, radians
from scitbx.matrix import col, dihedral_angle
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import iotbx.pdb
from mmtbx.validation import cbetadev

# a2.pdb chain C1 10-17, backbone and CB: 11-16 have phi/psi corrections
pdb_str = """\
ATOM      1  N   HISC1  10      -5.933  -7.199  40.232  1.00 66.79           N
ATOM      2  CA  HISC1  10      -6.171  -7.818  38.931  1.00 79.24           C
ATOM      3  CB  HISC1  10      -6.253  -9.349  39.054  1.00 73.16           C
ATOM      9  C   HISC1  10      -7.432  -7.276  38.260  1.00 52.81           C
ATOM     10  O   HISC1  10      -8.471  -7.107  38.891  1.00 32.51           O
ATOM     11  N   METC1  11      -7.315  -7.025  36.964  1.00 47.49           N
ATOM     12  CA  METC1  11      -8.416  -6.540  36.142  1.00 44.78           C
ATOM     13  CB  METC1  11      -8.062  -6.738  34.660  1.00 48.84           C
ATOM     17  C   METC1  11      -9.741  -7.242  36.452  1.00 41.27           C
ATOM     18  O   METC1  11     -10.771  -6.612  36.648  1.00 42.99           O
ATOM     19  N   SERC1  12      -9.700  -8.562  36.495  1.00 45.91           N
ATOM     20  CA  SERC1  12     -10.911  -9.360  36.590  1.00 37.87           C
ATOM     21  CB  SERC1  12     -10.548 -10.836  36.502  1.00 27.69           C
ATOM     23  C   SERC1  12     -11.711  -9.113  37.868  1.00 51.16           C
ATOM     24  O   SERC1  12     -12.913  -9.402  37.914  1.00 56.64           O
ATOM     25  N   ILEC1  13     -11.053  -8.567  38.890  1.00 41.17           N
ATOM     26  CA  ILEC1  13     -11.584  -8.601  40.261  1.00 42.49           C
ATOM     27  CB  ILEC1  13     -10.464  -8.418  41.316  1.00 37.03           C
ATOM     31  C   ILEC1  13     -12.735  -7.644  40.553  1.00 39.93           C
ATOM     32  O   ILEC1  13     -12.647  -6.457  40.278  1.00 43.38           O
ATOM     33  N   ILEC1  14     -13.808  -8.181  41.126  1.00 35.12           N
ATOM     34  CA  ILEC1  14     -14.962  -7.383  41.492  1.00 39.52           C
ATOM     35  CB  ILEC1  14     -16.124  -8.262  41.960  1.00 38.38           C
ATOM     39  C   ILEC1  14     -14.629  -6.427  42.614  1.00 50.82           C
ATOM     40  O   ILEC1  14     -14.223  -6.862  43.695  1.00 48.23           O
ATOM     41  N   THRC1  15     -14.828  -5.131  42.356  1.00 54.80           N
ATOM     42  CA  THRC1  15     -14.499  -4.070  43.315  1.00 50.18           C
ATOM     43  CB  THRC1  15     -13.450  -3.066  42.750  1.00 45.76           C
ATOM     46  C   THRC1  15     -15.724  -3.271  43.752  1.00 58.74           C
ATOM     47  O   THRC1  15     -15.610  -2.367  44.586  1.00 57.15           O
ATOM     48  N   ASPC1  16     -16.889  -3.596  43.192  1.00 52.46           N
ATOM     49  CA  ASPC1  16     -18.091  -2.813  43.481  1.00 47.59           C
ATOM     50  CB  ASPC1  16     -17.974  -1.412  42.838  1.00 54.77           C
ATOM     54  C   ASPC1  16     -19.427  -3.486  43.099  1.00 36.40           C
ATOM     55  O   ASPC1  16     -19.668  -3.846  41.943  1.00 27.65           O
ATOM     56  N   VALC1  17     -20.284  -3.648  44.101  1.00 33.55           N
ATOM     57  CA  VALC1  17     -21.686  -3.977  43.884  1.00 29.34           C
ATOM     58  CB  VALC1  17     -22.070  -5.277  44.590  1.00 39.01           C
ATOM     61  C   VALC1  17     -22.546  -2.850  44.461  1.00 37.97           C
ATOM     62  O   VALC1  17     -22.254  -2.336  45.552  1.00 40.83           O
"""

def get_hierarchy(n_models=1):
  lines = pdb_str.splitlines()
  if n_models > 1:
    models = []
    for i in range(n_models):
      models += ['MODEL     %4d' % (i+1)] + lines + ['ENDMDL']
    lines = models
  return iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()

def run_cbetadev(hierarchy, out=None, quiet=True, **kwds):
  if out is None: out = null_out()
  return cbetadev.cbetadev(hierarchy, out=out, quiet=quiet, **kwds)

def exercise_frame():
  # a point at angle t in the plane perpendicular to CA->ideal is at
  # dihedral N-CA-ideal-point t
  N, CA, ideal = col((1.2, 0.9, 0.3)), col((0, 0, 0)), col((0.1, -0.4, 1.5))
  e1, e2 = cbetadev.deviation_frame(N, CA, ideal)
  for t in [30, 120, -75, 179]:
    p = ideal + (e1*cos(radians(t)) + e2*sin(radians(t)))*0.2
    assert approx_equal(dihedral_angle(
      sites=[N.elems, CA.elems, ideal.elems, p.elems], deg=True), t, eps=1e-6)

def exercise_ideal_xyz():
  h = get_hierarchy()
  uncorrected = run_cbetadev(h, collect_ideal=True)
  corrected = run_cbetadev(h, apply_phi_psi_correction=True, collect_ideal=True)
  n_moved = 0
  for u, c in zip(uncorrected.results, corrected.results):
    assert u.id_str() == c.id_str()
    # the modelled CB is the reported deviation from the reported ideal
    assert approx_equal(abs(col(u.xyz) - col(u.ideal_xyz)), u.deviation, eps=1e-6)
    assert approx_equal(abs(col(c.xyz) - col(c.ideal_xyz)), c.deviation, eps=1e-6)
    assert corrected.get_beta_ideal()[c.id_str()] == c.ideal_xyz
    if c.deviation == u.deviation:
      assert c.ideal_xyz == u.ideal_xyz # no correction, no move
    else:
      n_moved += 1
      assert abs(col(c.ideal_xyz) - col(u.ideal_xyz)) > 1e-3
  assert n_moved == 6, n_moved

def exercise_summary():
  h = get_hierarchy(n_models=2)
  stdout = sys.stdout
  for quiet in [False, True]:
    sys.stdout = captured = io.StringIO()
    out = io.StringIO()
    try:
      run_cbetadev(h, out=out, quiet=quiet, apply_phi_psi_correction=True)
    finally:
      sys.stdout = stdout
    assert captured.getvalue() == '', captured.getvalue()
    n = out.getvalue().count('Outliers removed')
    assert n == (0 if quiet else 1), (quiet, n)
    if not quiet:
      assert 'Num. of residues :    16' in out.getvalue(), out.getvalue()

def run():
  exercise_frame()
  exercise_ideal_xyz()
  exercise_summary()
  print('OK')

if __name__ == '__main__':
  run()
