from __future__ import absolute_import, division, print_function
from six.moves import range
from dials.array_family import flex
from xfel.merging.application.worker import worker

class error_modifier_ha14(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(error_modifier_ha14, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Adjust intensity errors -- ha14 method'

  def run(self, experiments, reflections):
    '''Modify intensity errors according to the ha14 error model'''
    assert self.params.merging.error.model == "ha14"

    new_reflections = flex.reflection_table()
    self.logger.log_step_time("ERROR_MODIFIER_HA14")
    if len(reflections) > 0:
      self.logger.log("Modifying intensity errors -- ha14 method...")
      for expt_id in range(len(experiments)):
        refls = reflections.select(reflections['id'] == expt_id)
        refls = self.modify_errors(refls)
        new_reflections.extend(refls)
    self.logger.log_step_time("ERROR_MODIFIER_HA14", True)

    return experiments, new_reflections

  def modify_errors(self, reflections):
    '''Formerly sdfac_auto, ha14 method applies sdfac to each-image data assuming negative intensities are normally distributed noise'''
    assert 0 == (reflections['intensity.sum.variance'] <= 0.0).count(True)
    I_over_sig = reflections['intensity.sum.value'] / flex.sqrt(reflections['intensity.sum.variance'])
    negative_I_over_sig = I_over_sig.select(I_over_sig < 0.)
    #assert that at least a few I/sigmas are less than zero
    negative_I_over_sig_count = I_over_sig.select(I_over_sig < 0.).size()
    if negative_I_over_sig_count > 2:
      # get a rough estimate for the SDFAC, assuming that negative measurements
      # represent false predictions and therefore normally distributed noise.
      no_signal = negative_I_over_sig
      for xns in range(len(no_signal)):
        no_signal.append(-no_signal[xns])
      Stats = flex.mean_and_variance(no_signal)
      SDFAC = Stats.unweighted_sample_standard_deviation()
    else:
      SDFAC = 1.
    reflections['intensity.sum.variance'] *= (SDFAC**2)
    if self.params.output.log_level == 0:
      self.logger.log("The applied SDFAC is %7.4f"%SDFAC)
    return reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(error_modifier)
