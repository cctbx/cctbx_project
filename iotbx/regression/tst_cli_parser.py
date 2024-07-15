from __future__ import absolute_import, division, print_function

import os
import json

from multiprocessing import Process
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
def test_model_type_parsing():
  data_dir = os.path.dirname(os.path.abspath(__file__))
  model_1yjp = os.path.join(data_dir, 'data', '1yjp.pdb')
  model_2erl = os.path.join(data_dir, 'data', '2ERL.pdb')

  class testProgram(ProgramTemplate):

    datatypes = ['map_coefficients', 'miller_array', 'model', 'phil', 'real_map']

    def validate(self):
      pass

    def run(self):
      if self.data_manager.get_model_type(filename=model_1yjp) != ['neutron']:
        raise Sorry('1jyp should be neutron')
      if self.data_manager.get_model_type(filename=model_2erl) != ['electron']:
        raise Sorry('2erl should be electron')

  # model.type order matches input model order
  run_program(program_class=testProgram, args=[model_1yjp, model_2erl,
    'model.type=neutron', 'model.type=electron', '--quiet'])

  # other order will fail
  try:
    run_program(program_class=testProgram, args=[model_1yjp, model_2erl,
      'model.type=electron', 'model.type=neutron', '--quiet'])
  except Sorry as e:
    assert '1jyp' in str(e)
    assert 'neutron' in str(e)

  # and there needs to be a model.type for each model
  try:
    run_program(program_class=testProgram, args=[model_1yjp, model_2erl,
      'model.type=neutron', '--quiet'])
  except Sorry as e:
    assert 'Please specify exactly one "model.type" for each model' in str(e)

# -----------------------------------------------------------------------------
def test_user_selected_labels():
  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_mtz = os.path.join(data_dir, 'data',
                          'insulin_unmerged_cutted_from_ccp4.mtz')

  class testProgram(ProgramTemplate):
    program_name = 'tst_cli_parser'
    datatypes = ['map_coefficients', 'miller_array', 'phil', 'real_map']
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

  # duplicate user labels
  try:
    result = run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'miller_array.user_selected=WID',
            'miller_array.user_selected=WIDT'])
  except Sorry as s:
    assert 'duplicate user_selected_labels' in str(s)

  # map_coefficients.user_selected will not work on files without map coefficients
  try:
    result = run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'miller_array.user_selected=WID',
            'map_coefficients.user_selected=WIDT'])
  except Sorry as s:
    assert 'does not seem to have map_coefficients data' in str(s)

  # map_coefficients.user_selected_labels should be mirrored in
  # miller_array.user_selected_labels
  data_mtz = os.path.join(data_dir, 'data',
                          'phaser_1.mtz')
  class testProgram(ProgramTemplate):

    datatypes = ['map_coefficients', 'miller_array', 'phil', 'real_map']

    def validate(self):
      pass

    def run(self):
      ma_labels = self.data_manager.get_miller_array_user_selected_labels()
      mc_labels = self.data_manager.get_map_coefficients_user_selected_labels()
      print(ma_labels, file=self.logger)
      print(mc_labels, file=self.logger)
      assert ma_labels == mc_labels

  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite', data_mtz,
          'map_coefficients.user_selected=FW']
  )

  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite', data_mtz,
          'miller_array.user_selected=FW']
  )

  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite', data_mtz,
          'map_coefficients.labels.name=FW']
  )

  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite', data_mtz,
          'miller_array.labels.name=FW']
  )

  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite', data_mtz,
          'map_coefficients.user_selected=FWT,PHIFWT']
  )

  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite', data_mtz,
          'miller_array.user_selected=FWT,PHIFWT']
  )

  # non map_coefficients should not be copied
  try:
    run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'miller_array.user_selected=FP']
    )
  except AssertionError:
    pass

  # or used to select non map_coefficients data
  try:
    run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'map_coefficients.user_selected=FP']
    )
  except Sorry as s:
    assert 'is not recognized to be map_coefficients data' in str(s)

  try:
    run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'map_coefficients.user_selected=WT',
            'miller_array.user_selected=PHIF']
    )
  except Sorry as s:
    assert 'duplicate user_selected_labels' in str(s)

  try:
    run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'miller_array.user_selected=WT',
            'miller_array.user_selected=PHIF']
    )
  except Sorry as s:
    assert 'duplicate user_selected_labels' in str(s)

  try:
    run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'map_coefficients.user_selected=WT',
            'map_coefficients.user_selected=PHIF']
    )
  except Sorry as s:
    assert 'duplicate user_selected_labels' in str(s)

  try:
    run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'map_coefficients.labels.name=FWT,PHIFWT',
            'map_coefficients.user_selected=PHIF']
    )
  except Sorry as s:
    assert 'duplicate user_selected_labels' in str(s)

  try:
    run_program(
      program_class=testProgram,
      args=['--quiet', '--overwrite', data_mtz,
            'miller_array.labels.name=FWT,PHIFWT',
            'map_coefficients.user_selected=PHIF']
    )
  except Sorry as s:
    assert 'duplicate user_selected_labels' in str(s)

# -----------------------------------------------------------------------------
def test_json():
  class testProgram(ProgramTemplate):
    program_name = 'tst_cli_parser'
    datatypes = ['model', 'phil']
    def validate(self):
      pass
    def run(self):
      pass

  data_dir = os.path.dirname(os.path.abspath(__file__))
  model_1yjp = os.path.join(data_dir, 'data', '1yjp.pdb')

  # check for get_results_as_JSON function
  parser_log = StringIO()
  logger = multi_out()
  logger.register('parser_log', parser_log)

  run_program(
    program_class=testProgram,
    args=['--overwrite', '--json', model_1yjp],
    logger=parser_log
  )

  parser_log.flush()
  text = parser_log.getvalue()
  assert 'WARNING: The get_results_as_JSON function has not been defined for this program' in text

  logger.close()

  # check that json file is output
  expected_result = {'key': 'value'}
  class testProgram(ProgramTemplate):
    program_name = 'tst_cli_parser'
    datatypes = ['model', 'phil']
    def validate(self):
      pass
    def run(self):
      pass
    def get_results_as_JSON(self):
      dummy_results = expected_result
      return json.dumps(dummy_results, indent=2)

  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite', '--json', model_1yjp]
  )

  expected_filename = 'tst_cli_parser_result.json'
  assert os.path.exists(expected_filename)
  with open(expected_filename, 'r') as f:
    result = json.loads(f.read())
  assert result == expected_result
  os.remove(expected_filename)

  # check non-default filename
  expected_filename = 'non-default.json'
  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite', '--json',
          '--json-filename', expected_filename, model_1yjp]
  )
  assert os.path.exists(expected_filename)
  with open(expected_filename, 'r') as f:
    result = json.loads(f.read())
  assert result == expected_result
  os.remove(expected_filename)

  # try only with --json_filename and without .json extension
  filename = 'non-default'
  run_program(
    program_class=testProgram,
    args=['--quiet', '--overwrite',
          '--json_filename', filename, model_1yjp]
  )
  assert os.path.exists(expected_filename)
  with open(expected_filename, 'r') as f:
    result = json.loads(f.read())
  assert result == expected_result
  os.remove(expected_filename)

# -----------------------------------------------------------------------------
# since --diff-params calls sys.exit, run in a separate process
class testProgram(ProgramTemplate):
  program_name = 'test_diff_params'
  master_phil_str = '''
diff_test_parameter = None
.type = str
'''
  def run():
    pass
  def validate(self):
    pass

def run_diff_program(args):
  return run_program(program_class=testProgram, args=args)

def run_function_in_process(args):
  p = Process(target=run_diff_program, args=[args])
  p.start()
  p.join()

def test_diff_params():

  expected_filename = 'test_diff_params_modified.eff'
  if os.path.exists(expected_filename):
    os.remove(expected_filename)

  # no diff
  args = ['--quiet', '--diff-params']
  run_function_in_process(args)
  assert not os.path.exists(expected_filename)

  # program diff
  args = ['--quiet', '--diff-params', 'diff_test=abc']
  run_function_in_process(args)
  with open(expected_filename, 'r') as f:
    text = f.read()
    assert 'diff_test_parameter = abc' in text.strip(), text

  # DataManager diff
  data_dir = os.path.dirname(os.path.abspath(__file__))
  model_1yjp = os.path.join(data_dir, 'data', '1yjp.pdb')
  args = [model_1yjp, '--quiet', '--diff-params']
  run_function_in_process(args)
  with open(expected_filename, 'r') as f:
    text = f.read()
    assert text.count('1yjp.pdb') == 2, text
    assert 'diff_test_parameter' not in text.strip(), text

  # both diff
  args = ['--quiet', '--diff-params', model_1yjp, 'diff_test=abc']
  run_function_in_process(args)
  with open(expected_filename, 'r') as f:
    text = f.read()
    assert text.count('1yjp.pdb') == 2, text
    assert 'diff_test_parameter' in text.strip(), text

  if os.path.exists(expected_filename):
    os.remove(expected_filename)

# =============================================================================
if __name__ == '__main__':
  test_dry_run()
  test_label_parsing()
  test_model_type_parsing()
  test_user_selected_labels()
  test_json()
  test_diff_params()

  print("OK")
