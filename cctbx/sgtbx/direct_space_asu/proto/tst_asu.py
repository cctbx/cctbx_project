from cctbx.sgtbx.direct_space_asu.proto import direct_space_asu
from cctbx.sgtbx.direct_space_asu.proto import cut
from boost.rational import int as rint


def run():
  r1 = rint(1)
  r2 = rint(1,2)
  c = cut( (1,2,3), r1 )
  point = (r1,r2,r1)
  b = c.evaluate(point)
  print b
  s = direct_space_asu("19")
  print s.hall_symbol
  print s.n_faces()
  b = s.is_inside(point)
  print b
  s.show_comprehensive_summary()


if (__name__ == "__main__"):
  run()

