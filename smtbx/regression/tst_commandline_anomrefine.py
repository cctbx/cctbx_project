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

def run():
  for command, expected_result in zip(commands, expected_results):
    error_text = "Wrong result for '{}'".format(command)
    print(command)
    run_buf = fully_buffered(command)
    run_buf.raise_if_errors()
    out_lines = run_buf.stdout_lines
    result = []
    for line in out_lines:
      for val in line.split():
        result.append(float(val))
    assert len(result) == len(expected_result), error_text
    for x, y in zip(result, expected_result):
      assert approx_equal(x, y, eps=.0015), error_text

run()
