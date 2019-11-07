from __future__ import absolute_import, division, print_function
import scitbx.math
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
import sys

class test_function_5_1:

  def __init__(self, beta=2):
    self.beta = beta

  def functional(self, x):
    assert x.size() == 1
    alpha = x[0]
    return -alpha / (alpha**2 + self.beta)

  def gradients(self, x):
    assert x.size() == 1
    alpha = x[0]
    return flex.double([(alpha**2 - self.beta) / (alpha**2 + self.beta)**2])

def exercise(verbose):
  line_search = scitbx.math.line_search_more_thuente_1994()
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  line_search.show_status(f=out)
  print(file=out)
  assert line_search.maxfev == 20
  line_search.gtol = 0.01
  assert line_search.gtol == 0.01
  x = flex.double([0])
  search_direction = flex.double([1])
  tf = test_function_5_1()
  line_search.start(
    x=x,
    functional=tf.functional(x),
    gradients=tf.gradients(x),
    search_direction=search_direction,
    initial_estimate_of_satisfactory_step_length=1)
  assert line_search.info_code == -1
  line_search.show_status(f=out)
  print("x:", list(x), file=out)
  print(file=out)
  while (line_search.info_code == -1):
    line_search.next(
      x=x,
      functional=tf.functional(x),
      gradients=tf.gradients(x))
    if (verbose):
      line_search.show_status(f=out)
      print("x:", list(x), file=out)
      print("functional:", tf.functional(x), file=out)
      print(file=out)
  assert line_search.info_code == 1
  assert line_search.nfev == 3
  assert approx_equal(line_search.stp, 1.41997705997)
  print("OK")

if (__name__ == "__main__"):
  exercise(verbose="--verbose" in sys.argv[1:])
