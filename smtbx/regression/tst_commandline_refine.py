from __future__ import division
import os
from smtbx.command_line import refine
from smtbx.refinement import model
from libtbx.test_utils import approx_equal
from smtbx.regression.test_data import fnames

insfile = fnames.thpp_ins
hklfile = fnames.thpp_hkl
ciffile = fnames.thpp_cif
outfile = fnames.thpp_out

class run_options:
  def __init__(self):
    self.overwrite = True
    self.profile = None
    self.max_cycles = 8
    self.stop_if_max_derivative_below = 1e-7
    self.stop_if_shift_norm_below = 1e-7

def check_result(value):
  xm = model.from_cif(model=outfile, reflections=hklfile+'=hklf4')
  assert approx_equal(xm.xray_structure.scatterers()[0].site[0], value)

def run_tests():
  run_args = [insfile, hklfile, outfile]
  refine.run(run_args, run_options())
  check_result(0.167193)
  run_args = [ciffile, hklfile+'=hklf4', outfile]
  refine.run(run_args, run_options())
  check_result(0.166661)
  os.remove(outfile)

run_tests()
