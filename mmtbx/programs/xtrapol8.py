"""Execute Xtrapol8 calculations"""
from __future__ import absolute_import, division, print_function
import sys
from phenix.program_template import ProgramTemplate
from libtbx.utils import Sorry
import mmtbx.maps.xtrapol8 as xtrapol8

master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
n_alpha = 10
  .type = int
  .help = number of alpha values to test
'''

#  =============================================================================

class Program(ProgramTemplate):

  description = """
Run xtrapol8.

  phenix.xtrapol8 reference_model.pdb reference.mtz triggered.mtz
"""

  datatypes = ['model', 'phil', 'miller_array', 'restraint']

  master_phil_str = master_phil_str

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validate inputs:', file = self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    fs  = self.data_manager.has_miller_arrays(
      raise_sorry = False,
      expected_n  = 2,
      exact_count = True)
    if (len(self.data_manager.get_miller_array_names()) != 2):
      raise Sorry('Exactly 2 reflection files are expected.')

    print('  ...all good.', file = self.logger)

  # ---------------------------------------------------------------------------

  def run(self):

    print('Using model file:', self.data_manager.get_default_model_name(),
      file=self.logger)
    print('Using reflection file(s):', self.data_manager.get_miller_array_names(),
      file=self.logger)

    model_reference = self.data_manager.get_model()
    hkls = self.data_manager.get_miller_array_names()
    # hack until data_manager.fmodel supports this
    fn_reference = hkls[0]
    fn_triggered = hkls[1]
    f_obs_reference = self.data_manager.get_miller_arrays(['FOBS'], fn_reference)[0]
    f_obs_triggered = self.data_manager.get_miller_arrays(['FOBS'], fn_triggered)[0]

    xtr = xtrapol8.manager(
      model_reference = model_reference,
      f_obs_reference = f_obs_reference,
      f_obs_triggered = f_obs_triggered,
      log             = sys.stdout)
    xtr.run()

