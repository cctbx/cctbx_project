from __future__ import absolute_import, division, print_function

import pytest

@pytest.mark.parametrize("mosaic", (True, False))
def test_crystal(mosaic):
  from dxtbx.model import CrystalFactory
  from scitbx import matrix

  real_space_a = matrix.col((35.2402102454, -7.60002142787, 22.080026774))
  real_space_b = matrix.col((22.659572494, 1.47163505925, -35.6586361881))
  real_space_c = matrix.col((5.29417246554, 38.9981792999, 4.97368666613))

  if mosaic:
    from dxtbx.model import MosaicCrystalKabsch2010
    c1 = MosaicCrystalKabsch2010(
        real_space_a=real_space_a,
        real_space_b=real_space_b,
        real_space_c=real_space_c,
        space_group_symbol="P 1 2/m 1")
    c1.set_mosaicity(0.1)
  else:
    from dxtbx.model import Crystal
    c1 = Crystal(
        real_space_a=real_space_a,
        real_space_b=real_space_b,
        real_space_c=real_space_c,
        space_group_symbol="P 1 2/m 1")

  d = c1.to_dict()
  c2 = CrystalFactory.from_dict(d)
  eps = 1e-7
  assert abs(matrix.col(d['real_space_a']) - real_space_a) <= eps
  assert abs(matrix.col(d['real_space_b']) - real_space_b) <= eps
  assert abs(matrix.col(d['real_space_c']) - real_space_c) <= eps
  assert d['space_group_hall_symbol'] == "-P 2y"
  if mosaic:
    assert d['mosaicity'] == 0.1
  assert c1 == c2

def test_crystal_with_scan_points():
  from dxtbx.model import Crystal, CrystalFactory
  from scitbx import matrix

  real_space_a = matrix.col((35.2402102454, -7.60002142787, 22.080026774))
  real_space_b = matrix.col((22.659572494, 1.47163505925, -35.6586361881))
  real_space_c = matrix.col((5.29417246554, 38.9981792999, 4.97368666613))

  c1 = Crystal(
      real_space_a=real_space_a,
      real_space_b=real_space_b,
      real_space_c=real_space_c,
      space_group_symbol="P 1 2/m 1")

  A = c1.get_A()
  c1.set_A_at_scan_points([A for i in range(5)])

  # Set the B covariance. The values are nonsense, just ensure they are
  # all different
  from scitbx.array_family import flex
  cov_B = flex.double(range((9*9))) * 1e-5
  c1.set_B_covariance(cov_B)
  cov_B.reshape(flex.grid(1, 9, 9))
  cov_B_array = flex.double(flex.grid(5, 9, 9))
  for i in range(5):
    cov_B_array[i:(i+1), :, :] = cov_B
  c1.set_B_covariance_at_scan_points(cov_B_array)
  cov_B = c1.get_B_covariance()

  d = c1.to_dict()
  c2 = CrystalFactory.from_dict(d)
  eps = 1e-9
  for Acomp in d['A_at_scan_points']:
    for e1, e2 in zip(A, Acomp):
      assert abs(e1 - e2) <= eps
  for covBcomp in d['B_covariance_at_scan_points']:
    for e1, e2 in zip(cov_B, covBcomp):
      assert abs(e1 - e2) <= eps

  assert c1 == c2
