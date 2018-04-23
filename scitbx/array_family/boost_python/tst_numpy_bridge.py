from __future__ import division
from scitbx.array_family import flex
from libtbx.test_utils import Exception_expected, show_diff

try:
  import numpy as np
except ImportError:
  numpy = None

def exercise_type_conversions(flex_type, verbose):
  """Test converting various flex array types"""
  if (flex_type is flex.bool):
    z = False
  else:
    z = 0
  for n in xrange(10):
    fa = flex_type([z]*n)
    na = fa.as_numpy_array()
    assert na is not None
    if (n == 0 and verbose):
      print "flex.%s -> numpy %s" % (flex_type.__name__, na.dtype)
    assert na.shape == (n,)
    fna = flex_type(na)
    assert fna.all_eq(fa)

def test_int_conversions():
  """More exhaustive tests for converting int-type arrays"""
  fa = flex.int_range(1,7)
  na = fa.as_numpy_array()
  assert na.tolist() == list(fa)
  fna = flex.int(na)
  assert fna.all() == (6,)
  assert fna.origin() == (0,)
  assert fna.focus() == (6,)
  assert fna.all_eq(fa)
  fa[0] = 99
  assert na[0] == 1
  fa[0] = 1

  # Test that a 2D gridded versa gets properly converted
  fa.reshape(flex.grid(2,3))
  na = fa.as_numpy_array()
  assert na.tolist() == [[1, 2, 3], [4, 5, 6]]
  fna = flex.int(na)
  assert fna.all() == (2,3)
  assert fna.origin() == (0,0)
  assert fna.focus() == (2,3)
  assert fna.all_eq(fa)

  # Test that a gridded versa gets converted properly
  fa = flex.int_range(4*2*3) + 1
  fa.reshape(flex.grid(4,2,3))
  na = fa.as_numpy_array()
  assert na.tolist() == [
    [[1, 2, 3], [4, 5, 6]],
    [[7, 8, 9], [10, 11, 12]],
    [[13, 14, 15], [16, 17, 18]],
    [[19, 20, 21], [22, 23, 24]]]
  fna = flex.int(na)
  assert fna.all() == (4,2,3)
  assert fna.origin() == (0,0,0)
  assert fna.focus() == (4,2,3)
  assert fna.all_eq(fa)

  # Test converting a numpy array of a different type
  fa = flex.int(np.array([1,2], dtype=np.int64))
  assert fa[0] == 1
  assert fa[1] == 2

  # Test converting a very large array (memory leaks)
  bigarray = np.ones(1000000)
  fa = flex.int(bigarray)
  assert (fa == 1).count(True)
  # And, back to numpy
  bigarray = fa.as_numpy_array()
  assert np.sum(bigarray == 1) == len(bigarray)


def run(args):
  assert args in [[], ["--forever"]]
  verbose = True
  while True:
    if (flex.int().as_numpy_array(optional=True) is None):
      try:
        flex.int().as_numpy_array()
      except RuntimeError, e:
        assert not show_diff(str(e), "numpy API not available")
      else:
        raise Exception_expected
    else:
      for flex_type in [
            flex.bool,
            flex.int,
            flex.long,
            flex.float,
            flex.double,
            flex.complex_double,
            flex.size_t]:
        exercise_type_conversions(flex_type, verbose)
      test_int_conversions()
    if (len(args) == 0):
      break
    verbose = False
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
