from __future__ import absolute_import, division, print_function
def exercise():
  from cctbx.web import change_basis
  from libtbx.utils import Sorry
  from libtbx.test_utils import Exception_expected
  def check(string, expected):
    rt = change_basis.p_from_string(string)
    ex = change_basis.p_from_string(expected)
    assert rt.r == ex.r
    assert rt.t == ex.t
  check("", "a,b,c")
  check("a,b,c", "a,b,c")
  check("{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}", "a,b,c")
  check("2, 0, 0, 0, 2, 0, 0, 0, 2", "2a,2b,2c")
  check("{{1/2, 0, 0}, {0, 1/2, 0}, {0, 0, 1/2}}", "1/2a,1/2b,1/2c")
  check("{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 1/2, 1/3, 1/4",
    "a+1/2,b+1/3,c+1/4")
  check("(1 0 0 | 1/2) (0 1 0 | 1/3) (0 0 1 | 1/4)",
    "a+1/2,b+1/3,c+1/4")
  check("1/2,1/3,1/4", "a+1/2,b+1/3,c+1/4")
  try:
    change_basis.p_from_string("0")
  except Sorry as e:
    assert str(e) == \
      'Uninterpretable expression for change-of-basis matrix'
  else: raise Exception_expected
  #
  rt = change_basis.w_from_string("1/2,1/3,1/4")
  assert rt.r.elems == (0,)*9
  assert str(rt.t.elems) == "(1/2, 1/3, 1/4)"
  try:
    change_basis.w_from_string("0")
  except Sorry as e:
    assert str(e) == \
      'Uninterpretable expression for symmetry matrix'
  else: raise Exception_expected
  #
  xyz = change_basis.xyz_from_string("0.5,1/3,4")
  assert str(xyz) == "(0.5, 1/3, 4)"
  try:
    change_basis.xyz_from_string("0")
  except Sorry as e:
    assert str(e) == \
      'Uninterpretable expression for coordinates'
  else: raise Exception_expected

def run(args):
  assert len(args) == 0
  exercise()
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
