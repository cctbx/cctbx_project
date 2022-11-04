from __future__ import absolute_import, division, print_function

import os

from iotbx.cli_parser import run_program
from iotbx.data_manager import DataManager
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

# -----------------------------------------------------------------------------
def test_label_parsing():
  class testProgram(ProgramTemplate):
    program_name = 'tst_cli_parser'
    datatypes = ['miller_array', 'model', 'phil', 'restraint']
    def validate(self):
      pass
    def run(self):
      pass

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_mtz = os.path.join(data_dir, 'data',
                          'insulin_unmerged_cutted_from_ccp4.mtz')
  labels = ['M_ISYM', 'BATCH', 'I,SIGI,merged', 'IPR,SIGIPR,merged',
            'FRACTIONCALC', 'XDET', 'YDET', 'ROT', 'WIDTH', 'LP', 'MPART',
            'FLAG', 'BGPKRATIOS']
  phil_filename = testProgram.program_name + '_all.eff'

  # check label matching and parsing
  for phil_args, label_args in [
    (['labels.name', 'labels.name', 'labels.name'], ['XD', 'SIGI', 'WIDTH']),
    (['user_selected_labels', 'user_selected_labels', 'user_selected_labels'], ['YD', 'M_ISYM', 'IPR']),
    (['labels.name', 'user_selected_labels', 'user_selected_labels'], ['BATCH', 'FRACTION', 'BGP']),
    (['user_selected_labels', 'labels.name', 'user_selected_labels', 'labels.name'], ['LP', 'ROT', 'FLAG', 'SIGI'])
  ]:

    combined_args = ['{}={}'.format(phil_arg, label_arg) for phil_arg, label_arg in zip(phil_args, label_args)]

    run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', '--write-all', data_mtz] + combined_args)

    dm = DataManager()
    dm.process_phil_file(phil_filename)
    working_phil = dm.master_phil.fetch(dm.get_phil())
    extract = working_phil.extract()
    for label in extract.data_manager.miller_array[0].labels:
      assert label.name in labels
    for label in label_args:
      assert label in extract.data_manager.miller_array[0].user_selected_labels

  os.remove(phil_filename)

  # check ambiguous label
  try:
    run_program(program_class=testProgram, args=['--quiet', data_mtz, 'labels.name=SIG'])
  except Sorry as s:
    assert 'The label, SIG' in str(s)
  try:
    run_program(program_class=testProgram, args=['--quiet', data_mtz, 'user_selected_label=SIG'])
  except Sorry as s:
    assert 'The label, SIG' in str(s)

# =============================================================================
if __name__ == '__main__':
  test_dry_run()
  test_label_parsing()

  print("OK")
