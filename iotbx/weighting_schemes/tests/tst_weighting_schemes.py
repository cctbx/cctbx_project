from iotbx import weighting_schemes
from libtbx.tst_utils import show_diff

def exercise():
  unit_weighting = weighting_schemes.unit_weighting()
  assert unit_weighting.type() == "unit"
  assert str(unit_weighting) == "w=1"
  shelx_weighting = weighting_schemes.mainstream_shelx_weighting(0.1234,
                                                                 0.5678)
  assert shelx_weighting.type() == "calc"
  assert not show_diff(
    str(shelx_weighting),
    "w=1/[\s^2^(Fo^2^)+(0.1234P)^2^+0.5678P] where P=(Fo^2^+2Fc^2^)/3")

def run():
  exercise()
  print "OK"

if __name__ == '__main__':
  run()
