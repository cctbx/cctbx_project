from __future__ import division, print_function
from libtbx.program_template import ProgramTemplate
import mmtbx.nci.hbond

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.hbond: tool to find H bonds in an atomic model

Usage example:
  phenix.hbond model.pdb
  '''

  datatypes = ['model', 'phil']

  master_phil_str = mmtbx.nci.hbond.master_phil_str

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(raise_sorry=True)

  # ---------------------------------------------------------------------------
  def run(self):
    print('Using model: %s' % self.data_manager.get_default_model_name(),
      file=self.logger)
    model = self.data_manager.get_model()
    self.results = mmtbx.nci.hbond.find(model = model)
    self.results.show(log = self.logger)
    print("-"*79, file=self.logger)
    sk = self.results.get_theta_2_skew_and_kurtosis()
    if(sk.skew_all is not None):
      print("all: skew: %7.3f kurtosis: %7.3f"%(sk.skew_all, sk.kurtosis_all), file=self.logger)
    if(sk.skew_filtered is not None):
      print("filtered: skew: %7.3f kurtosis: %7.3f"%(sk.skew_filtered, sk.kurtosis_filtered), file=self.logger)
    self.results.as_pymol()

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results

# =============================================================================
# end
