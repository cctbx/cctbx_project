from __future__ import division
from six.moves import range
from xfel.merging.application.worker import worker
from dials.array_family import flex

class error_model(worker):

  def run(self, experiments, reflections):

    self.logger.log_step_time("ERROR_MODEL")
    if len(reflections) > 0:
      self.logger.log("Modifying intensity errors...")
      if self.params.merging.error.model == "errors_from_sample_residuals":
        reflections = self.modify_errors_from_residuals(reflections)
    self.logger.log_step_time("ERROR_MODEL", True)

    return experiments, reflections

  def modify_errors_from_residuals(self, reflections):
    '''For each asu hkl replace the variances of individual measured intensities with the variance of the intensity distribution'''
    new_reflections = flex.reflection_table()
    number_of_hkls = 0
    number_of_multiply_measured_hkls = 0
    for refls in self.get_next_hkl_reflection_table(reflections):
      number_of_hkls += 1
      number_of_measurements = len(refls)
      assert number_of_measurements > 0
      if number_of_measurements > 1:
        stats = flex.mean_and_variance(refls['intensity.sum.value'])
        refls['intensity.sum.variance'] = flex.double(number_of_measurements, stats.unweighted_sample_variance())
        number_of_multiply_measured_hkls += 1
      new_reflections.extend(refls)

    self.logger.log("Modified errors for (multiply-measured) %d symmetry-reduced hkl's out of %d symmetry-reduced hkl's"%(number_of_multiply_measured_hkls, number_of_hkls))

    return new_reflections

  def get_next_hkl_reflection_table(self, reflections):
    '''Generate asu hkl slices from an asu hkl-sorted reflection table'''
    i_begin = 0
    hkl_ref = reflections[0].get('miller_index_asymmetric')
    for i in range(len(reflections)):
      hkl = reflections[i].get('miller_index_asymmetric')
      if hkl == hkl_ref:
        continue
      else:
        yield reflections[i_begin:i]
        i_begin = i
        hkl_ref = hkl
    yield reflections[i_begin:i+1]
