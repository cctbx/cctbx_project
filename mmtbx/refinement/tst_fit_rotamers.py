def exercise_rotate_point_around_axis():
  from scitbx.matrix import col, rotate_point_around_axis
  cb = col([7.767, 5.853, 7.671])
  cg = col([6.935, 5.032, 8.622])
  ca = col([7.000, 7.000, 7.000])
  count_unchanged = 0
  for angle_i in range(-360,361,15):
    cg_r = col(rotate_point_around_axis(
      axis_point_1=ca, axis_point_2=cb, point=cg, angle=angle_i, deg=True))
    assert abs(abs(cb-cg)-abs(cb-cg_r)) < 1.e-6
    assert abs(abs(ca-cg)-abs(ca-cg_r)) < 1.e-6
    if (angle_i in [-360,0,360]):
      assert abs(cg_r-cg) < 1.e-6
      count_unchanged += 1
    cg_r_alt = ca.rt_for_rotation_around_axis_through(
      point=cb, angle=angle_i, deg=True) * cg
    assert abs(cg_r_alt - cg_r) < 1.e-6
  assert count_unchanged == 3

def run(args):
  assert len(args) == 0
  exercise_rotate_point_around_axis()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
