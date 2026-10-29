from __future__ import absolute_import, division, print_function
from math import pi, cos, sin
from libtbx.test_utils import approx_equal
from mmtbx.validation.cbd_utils import radial_deviation

def as_xy(v):
  return (v.r*cos(v.t), v.r*sin(v.t))

def exercise_angle_wrapping():
  # 0.3 at 175 deg minus 0.1 at 90 deg lands at ~193.9 deg if not wrapped
  rc = radial_deviation(0.3, 175) - radial_deviation(0.1, 90)
  assert -pi < rc.t <= pi, rc.t/pi*180
  assert approx_equal(rc.t/pi*180, 193.88072251867072-360, eps=1e-6)
  # wrapping must not change the vector itself
  a = radial_deviation(0.3, 175)
  b = radial_deviation(0.1, 90)
  assert approx_equal(as_xy(rc), (as_xy(a)[0]-as_xy(b)[0],
                                  as_xy(a)[1]-as_xy(b)[1]), eps=1e-9)
  # sweep: results always in [-180, 180] (+/-180 equivalent, fp noise)
  for t1 in range(-180, 181, 15):
    for t2 in range(-180, 181, 15):
      for op in ['add', 'sub']:
        x = radial_deviation(0.3, t1)
        y = radial_deviation(0.2, t2)
        z = x+y if op=='add' else x-y
        assert -pi-1e-12 <= z.t <= pi+1e-12, (t1, t2, op, z.t/pi*180)

def run():
  exercise_angle_wrapping()
  print('OK')

if __name__ == '__main__':
  run()
