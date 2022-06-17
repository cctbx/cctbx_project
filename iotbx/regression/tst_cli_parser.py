from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry

from libtbx.test_utils import Exception_expected, Exception_not_expected

# =============================================================================
def test_dry_run():
  class testProgram(ProgramTemplate):
    def validate(self):
      raise Sorry('This is a test')

  try:
    run_program(program_class=testProgram, args=['--dry-run', '--quiet'])
  except Sorry:
    pass
  else:
    raise Exception_expected

  class testProgram(ProgramTemplate):
    def validate(self):
      pass

  try:
    run_program(program_class=testProgram, args=['--dry-run', '--quiet'])
  except Exception:
    raise Exception_not_expected

# =============================================================================
if __name__ == '__main__':
  test_dry_run()

  print("OK")
