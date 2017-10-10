from __future__ import absolute_import, division
from __future__ import print_function

from builtins import zip
from builtins import range
from builtins import object
class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_crystal()
    self.tst_crystal(mosaic=True)
    self.tst_crystal_with_scan_points()

  def tst_crystal(self, mosaic=False):
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
    assert(abs(matrix.col(d['real_space_a']) - real_space_a) <= eps)
    assert(abs(matrix.col(d['real_space_b']) - real_space_b) <= eps)
    assert(abs(matrix.col(d['real_space_c']) - real_space_c) <= eps)
    assert(d['space_group_hall_symbol'] == "-P 2y")
    if mosaic:
      assert(d['mosaicity'] == 0.1)
    assert(c1 == c2)
    print('OK')

  def tst_crystal_with_scan_points(self):
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

    d = c1.to_dict()
    c2 = CrystalFactory.from_dict(d)
    eps = 1e-9
    for Acomp in (d['A_at_scan_points']):
      for e1, e2 in zip(A, Acomp):
        assert(abs(e1 - e2) <= eps)
    assert(c1 == c2)
    print('OK')

if __name__ == '__main__':
  test = Test()
  test.run()
