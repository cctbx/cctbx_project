from __future__ import absolute_import, division, print_function

import os
import sys

from six.moves import cStringIO as StringIO

from iotbx.cli_parser import run_program
from iotbx.data_manager import DataManager
from libtbx.program_template import ProgramTemplate
from libtbx.utils import multi_out, Sorry
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
    (['labels.name'], ['xd']),
    (['labels.name', 'labels.name', 'labels.name'], ['XD', 'SIGI', 'WIDTH']),
    (['user_selected_labels', 'user_selected_labels', 'user_selected_labels'], ['YD', 'M_ISYM', 'IPR']),
    (['labels.name', 'user_selected_labels', 'user_selected_labels'], ['BATCH', 'FRACTION', 'BGP']),
    (['user_selected_labels', 'labels.name', 'user_selected_labels', 'labels.name'], ['LP', 'ROT', 'FLAG', 'SIGI'])
  ]:

    combined_args = ['{}={}'.format(phil_arg, label_arg) for phil_arg, label_arg in zip(phil_args, label_args)]

    parser_log = StringIO()
    logger = multi_out()
    logger.register('parser_log', parser_log)

    run_program(
      program_class=testProgram,
      args=['--overwrite', '--write-all', data_mtz] + combined_args,
      logger=logger)

    parser_log.flush()
    text = parser_log.getvalue()
    assert 'Combined labels PHIL' in text

    logger.close()

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

# -----------------------------------------------------------------------------
def test_user_selected_labels():
  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_mtz = os.path.join(data_dir, 'data',
                          'insulin_unmerged_cutted_from_ccp4.mtz')

  class testProgram(ProgramTemplate):
    program_name = 'tst_cli_parser'
    datatypes = ['miller_array', 'phil']
    master_phil_str = '''
other_file = None
  .type = path
'''
    def validate(self):
      pass
    def run(self):
      self.result = self.data_manager.get_miller_array_user_selected_labels(data_mtz)
      params = self.data_manager.export_phil_scope(as_extract=True)
      assert params.data_manager.miller_array[0].user_selected_labels == self.result
    def get_results(self):
      return self.result

  phil_name = 'user_selected.eff'

  # phil file with just the DataManager scope
  phil_one = '''
data_manager {
  miller_array {
    file = %s
    user_selected_labels = IPR,SIGIPR,merged
    user_selected_labels = FRACTIONCALC
  }
}
''' % data_mtz

  # phil file with the same file as a separate parameter
  phil_two = '''
data_manager {
  miller_array {
    file = %s
    user_selected_labels = IPR,SIGIPR,merged
    user_selected_labels = FRACTIONCALC
  }
}
other_file = %s
''' % (data_mtz, data_mtz)

  for phil in [phil_one, phil_two]:
    with open(phil_name, 'w') as f:
      f.write(phil)

    result = run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', phil_name])

    assert len(result) == 2
    for label in result:
      assert label in ['FRACTIONCALC' ,'IPR,SIGIPR,merged']

  os.remove(phil_name)

# =============================================================================
if __name__ == '__main__':
  test_dry_run()
  test_label_parsing()
  test_user_selected_labels()

  print("OK")
