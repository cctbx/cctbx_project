from __future__ import absolute_import, division, print_function
from math import pi
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
from libtbx import group_args
import iotbx.pdb
from mmtbx.validation import cbetadev
from mmtbx.validation import cbd_utils

# 3a0f chain A 126-132: GLU 129 (phi 71, psi 100) is in an empty bin of the
# others table with populated neighbours
pdb_str = """\
ATOM    961  N   ASN A 126       2.784  75.120  30.098  1.00 26.44           N
ATOM    962  CA  ASN A 126       1.877  74.929  28.975  1.00 26.67           C
ATOM    963  C   ASN A 126       1.498  76.271  28.390  1.00 26.91           C
ATOM    964  O   ASN A 126       1.081  76.355  27.239  1.00 25.80           O
ATOM    965  CB  ASN A 126       0.630  74.193  29.446  1.00 26.95           C
ATOM    966  CG  ASN A 126       0.129  74.698  30.792  1.00 28.31           C
ATOM    967  OD1 ASN A 126       0.743  74.447  31.828  1.00 29.99           O
ATOM    968  ND2 ASN A 126      -0.975  75.429  30.778  1.00 28.06           N
ATOM    969  N   GLY A 127       1.650  77.307  29.231  1.00 27.61           N
ATOM    970  CA  GLY A 127       1.468  78.704  28.886  1.00 27.35           C
ATOM    971  C   GLY A 127       2.684  79.251  28.181  1.00 28.00           C
ATOM    972  O   GLY A 127       3.751  78.632  28.164  1.00 27.60           O
ATOM    973  N   GLY A 128       2.538  80.492  27.745  1.00 28.46           N
ATOM    974  CA  GLY A 128       2.634  80.783  26.353  1.00 29.73           C
ATOM    975  C   GLY A 128       3.383  81.937  25.808  1.00 30.49           C
ATOM    976  O   GLY A 128       2.804  82.885  25.339  1.00 30.62           O
ATOM    977  N   GLU A 129       4.691  81.723  25.811  1.00 31.64           N
ATOM    978  CA  GLU A 129       5.793  82.570  25.348  1.00 30.66           C
ATOM    979  C   GLU A 129       6.118  83.798  26.171  1.00 29.80           C
ATOM    980  O   GLU A 129       5.450  84.808  26.089  1.00 29.78           O
ATOM    981  CB  GLU A 129       5.875  82.734  23.838  1.00 31.48           C
ATOM    982  CG  GLU A 129       6.226  81.401  23.024  1.00 31.98           C
ATOM    983  CD  GLU A 129       7.082  80.344  23.764  1.00 32.69           C
ATOM    984  OE1 GLU A 129       8.140  80.683  24.358  1.00 32.33           O
ATOM    985  OE2 GLU A 129       6.690  79.148  23.728  1.00 31.23           O
ATOM    986  N   ARG A 130       7.160  83.650  26.985  1.00 28.30           N
ATOM    987  CA  ARG A 130       7.562  84.671  27.935  1.00 27.32           C
ATOM    988  C   ARG A 130       8.662  85.565  27.355  1.00 27.42           C
ATOM    989  O   ARG A 130       8.688  86.761  27.660  1.00 28.27           O
ATOM    990  CB  ARG A 130       8.039  84.022  29.253  1.00 27.63           C
ATOM    991  CG  ARG A 130       6.971  83.207  29.999  1.00 25.87           C
ATOM    992  CD  ARG A 130       7.609  82.323  31.071  1.00 25.27           C
ATOM    993  NE  ARG A 130       6.597  81.504  31.716  1.00 23.96           N
ATOM    994  CZ  ARG A 130       5.846  81.882  32.751  1.00 23.93           C
ATOM    995  NH1 ARG A 130       5.993  83.102  33.308  1.00 21.64           N
ATOM    996  NH2 ARG A 130       4.938  81.027  33.232  1.00 20.73           N
ATOM    997  N   LEU A 131       9.561  84.987  26.547  1.00 25.82           N
ATOM    998  CA  LEU A 131      10.603  85.741  25.873  1.00 25.37           C
ATOM    999  C   LEU A 131      10.113  86.241  24.535  1.00 25.47           C
ATOM   1000  O   LEU A 131       9.497  85.503  23.802  1.00 25.27           O
ATOM   1001  CB  LEU A 131      11.813  84.852  25.595  1.00 25.25           C
ATOM   1002  CG  LEU A 131      13.234  85.195  26.041  1.00 23.59           C
ATOM   1003  CD1 LEU A 131      14.200  84.458  25.169  1.00 22.57           C
ATOM   1004  CD2 LEU A 131      13.579  86.703  26.142  1.00 21.21           C
ATOM   1005  N   ALA A 132      10.408  87.494  24.203  1.00 25.99           N
ATOM   1006  CA  ALA A 132      10.001  88.067  22.911  1.00 25.74           C
ATOM   1007  C   ALA A 132      10.995  89.140  22.443  1.00 25.83           C
ATOM   1008  O   ALA A 132      11.483  89.939  23.253  1.00 26.63           O
ATOM   1009  CB  ALA A 132       8.592  88.624  22.992  1.00 24.97           C
TER
END
"""

def get_deviations(**kwds):
  hierarchy = iotbx.pdb.input(lines=pdb_str.splitlines(),
                              source_info=None).construct_hierarchy()
  validation = cbetadev.cbetadev(pdb_hierarchy=hierarchy,
                                 out=null_out(),
                                 **kwds)
  return dict(((r.resseq.strip(), r.resname), r.deviation)
              for r in validation.results)

def tst_01():
  """
  Neighbour estimate on a small table.
  """
  table = {(60, 90)  : (0.1, 0.),
           (60, 100) : (0.1, pi/2),
           (70, 90)  : (0.2, 0.),
           }
  # mean of the three vectors: (0.1, 0.) + (0., 0.1) + (0.2, 0.)
  r, t = cbd_utils.get_neighbour_correction(table, (70, 100))
  assert approx_equal(r, ((0.3/3)**2+(0.1/3)**2)**0.5)
  assert approx_equal(t, 0.32175, eps=1e-4)
  # only two populated neighbours
  assert cbd_utils.get_neighbour_correction(table, (60, 80)) is None
  assert cbd_utils.get_neighbour_correction(table, (70, 100),
                                            min_neighbours=4) is None
  # psi 180 and -180 are the same angle and either key can be in a table
  table = {(60, 180)  : (0.1, 0.),
           (70, -180) : (0.1, 0.),
           (80, 170)  : (0.1, 0.),
           }
  r, t = cbd_utils.get_neighbour_correction(table, (70, 170))
  assert approx_equal((r, t), (0.1, 0.))
  table = {(60, 180)  : (0.1, 0.),
           (70, 180)  : (0.1, 0.),
           (80, -160) : (0.1, 0.),
           }
  r, t = cbd_utils.get_neighbour_correction(table, (70, -170))
  assert approx_equal((r, t), (0.1, 0.))

def tst_02():
  """
  Only residues with 0 < phi < 150 are filled.
  """
  table = cbd_utils.others_phi_psi
  assert (70, 100) not in table
  assert cbd_utils.get_neighbour_correction(table, (70, 100)) is not None
  # an empty bin with negative phi and enough populated neighbours
  negative = None
  for phi in range(-170, 0, 10):
    for psi in range(-170, 180, 10):
      if (phi, psi) in table: continue
      if cbd_utils.get_neighbour_correction(table, (phi, psi)) is not None:
        negative = (phi, psi)
        break
    if negative: break
  assert negative is not None
  result = group_args(deviation=0.3, dihedral=45.)
  residue = group_args(resname='GLU')
  for phi_psi, filled in [((71., 100.), True), (negative, False)]:
    assert cbd_utils.get_phi_psi_correction(result, residue, phi_psi) is None
    rc = cbd_utils.get_phi_psi_correction(result, residue, phi_psi,
                                          fill_empty_phi_psi_bins=True)
    assert (rc is not None) == filled, phi_psi

def tst_03():
  """
  Corrections of the fragment with and without fill_empty_phi_psi_bins.
  """
  default = get_deviations()
  assert approx_equal(default[('129', 'GLU')], 0.4690, eps=1e-4)
  for mode in ['apply_positive_phi_correction', 'apply_phi_psi_correction']:
    off = get_deviations(**{mode : True})
    on = get_deviations(fill_empty_phi_psi_bins=True, **{mode : True})
    assert approx_equal(off[('129', 'GLU')], 0.4690, eps=1e-4), mode
    assert approx_equal(on[('129', 'GLU')], 0.3807, eps=1e-4), mode
    # the other residues are unchanged by the fill
    for key in off:
      if key == ('129', 'GLU'): continue
      assert approx_equal(on[key], off[key]), (mode, key)

def tst_04():
  """
  fill_empty_phi_psi_bins needs one of the corrections.
  """
  from iotbx.cli_parser import run_program
  from mmtbx.programs import cbetadev as cbetadev_program
  file_name = 'tst_cbetadev_04.pdb'
  with open(file_name, 'w') as f:
    f.write(pdb_str)
  for args in [['fill_empty_phi_psi_bins=True'],
               ['fill_empty_phi_psi_bins=True', 'mahalanobis_z_score=True']]:
    try:
      run_program(program_class=cbetadev_program.Program,
                  args=[file_name] + args,
                  logger=null_out())
    except Sorry as e:
      assert str(e).find('fill_empty_phi_psi_bins requires') > -1
    else:
      raise AssertionError('Sorry not raised')
  run_program(program_class=cbetadev_program.Program,
              args=[file_name,
                    'apply_positive_phi_correction=True',
                    'fill_empty_phi_psi_bins=True'],
              logger=null_out())

if (__name__ == "__main__"):
  tst_01()
  tst_02()
  tst_03()
  tst_04()
  print("OK")
