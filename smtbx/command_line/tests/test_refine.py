from __future__ import division
import os
from smtbx.command_line import refine
from smtbx.refinement import model
from libtbx.test_utils import approx_equal

working_dir = os.path.dirname(__file__)
insfile = os.path.join(working_dir, 'thpp.ins')
hklfile = os.path.join(working_dir, 'thpp.hkl')
ciffile = os.path.join(working_dir, 'thpp.cif')
outfile = os.path.join(working_dir, 'thpp_out.cif')

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
  check_result(0.166963)
  run_args = [ciffile, hklfile+'=hklf4', outfile]
  refine.run(run_args, run_options())
  check_result(0.166644)
  os.remove(outfile)

run_tests()
