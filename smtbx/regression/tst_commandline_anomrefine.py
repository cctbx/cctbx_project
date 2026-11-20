from __future__ import division
from libtbx.easy_run import fully_buffered
from libtbx.test_utils import approx_equal
from smtbx.regression.test_data import fnames

ciffile=fnames.thpp_cif
insfile=fnames.thpp_ins
hklfile=fnames.thpp_hkl

commands = [
  "smtbx.anom_refine {} {} F -e13000 -t -O -d0 -s0 -c25".format(
      ciffile, hklfile),
  "smtbx.anom_refine {} {} F -e13000 -T -O -d0 -s0 -c25".format(
      insfile, hklfile)
  ]

expected_results = [
  [13000, .152, -.053, .743, 1.303],
  [13000, -.032, .001, .005, .008, .015, .016, 3.674, 4.279]
  ]

def numbers_in(out_lines):
  """ The numeric output, ignoring whatever else reached stdout.

  smtbx.anom_refine writes nothing but rows of numbers, but it does not have
  stdout to itself: importing fast_linalg initialises OpenBLAS and announces
  where it found it, and anything else on the way in is free to do the same.
  Taking a line only when every token on it is a number keeps those out
  without silently dropping a token from a line that is a result -- a result
  line with something unparsable on it fails here, as it should.
  """
  result = []
  for line in out_lines:
    values = []
    for val in line.split():
      try:
        values.append(float(val))
      except ValueError:
        values = None
        break
    if values:
      result.extend(values)
  return result


def run():
  for command, expected_result in zip(commands, expected_results):
    error_text = "Wrong result for '{}'".format(command)
    print(command)
    run_buf = fully_buffered(command)
    run_buf.raise_if_errors()
    result = numbers_in(run_buf.stdout_lines)
    assert len(result) == len(expected_result), error_text
    for x, y in zip(result, expected_result):
      assert approx_equal(x, y, eps=.0015), error_text

run()
