from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from xfel.merging.application.worker import worker
from xfel.merging.application.reflection_table_utils import reflection_table_utils

class error_modifier(worker):

  def __repr__(self):
    return 'Adjust intensity errors'

  def run(self, experiments, reflections):
    '''Modify intensity errors according to an error model'''
    self.logger.log_step_time("ERROR_MODIFIER")
    if len(reflections) > 0:
      self.logger.log("Modifying intensity errors...")
      if self.params.merging.error.model == "errors_from_sample_residuals":
        reflections = self.modify_errors_from_residuals(reflections)
    self.logger.log_step_time("ERROR_MODIFIER", True)

    return experiments, reflections

  def modify_errors_from_residuals(self, reflections):
    '''For each asu hkl replace the variances of individual measured intensities with the variance of the intensity distribution'''
    new_reflections = flex.reflection_table()
    number_of_hkls = 0
    number_of_multiply_measured_hkls = 0
    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections):
      number_of_hkls += 1
      number_of_measurements = len(refls)
      assert number_of_measurements > 0
      if number_of_measurements > 1:
        stats = flex.mean_and_variance(refls['intensity.sum.value'])
        refls['intensity.sum.variance'] = flex.double(number_of_measurements, stats.unweighted_sample_variance())
        number_of_multiply_measured_hkls += 1
      new_reflections.extend(refls)

    self.logger.log("Modified errors for (multiply-measured) %d symmetry-reduced HKLs out of %d symmetry-reduced HKLs"%(number_of_multiply_measured_hkls, number_of_hkls))

    return new_reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(error_modifier)
