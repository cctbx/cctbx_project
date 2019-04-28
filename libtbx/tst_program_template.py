from __future__ import absolute_import, division, print_function

import sys

import libtbx.phil
from libtbx.program_template import ProgramTemplate
from libtbx.utils import multi_out

# =============================================================================
class TestProgram(ProgramTemplate):

  master_phil = """
parameter_a = None
  .type = str
parameter_b = 0
  .type = int
parameter_c = None
  .type = float
"""

working_phil = libtbx.phil.parse("parameter_a = not None\nparameter_b = 5")

# -----------------------------------------------------------------------------
def test_phil():
  master_phil = libtbx.phil.parse(TestProgram.master_phil)
  params = master_phil.fetch(working_phil).extract()
  logger = multi_out()
  logger.register('stdout', sys.stdout)
  test_program = TestProgram(None, params, master_phil, logger)

  full_phil = libtbx.phil.parse(test_program.get_program_phil_str())
  full = master_phil.fetch(full_phil).extract()
  assert(full.parameter_a == 'not None')
  assert(full.parameter_b == 5)
  assert(full.parameter_c is None)
  assert('parameter_c' not in test_program.get_program_phil_str(True))

# =============================================================================
if __name__ == '__main__':
  test_phil()
