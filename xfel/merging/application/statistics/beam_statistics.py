from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dials.array_family import flex

class beam_statistics(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(beam_statistics, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Beam statistics'

  def run(self, experiments, reflections):
    self.logger.log_step_time("BEAM_STATISTICS")
    f_wavelengths = flex.double([b.get_wavelength() for b in experiments.beams()])

    flex_all_wavelengths = self.mpi_helper.aggregate_flex(f_wavelengths, flex.double, root=None)
    average_wavelength = flex.mean(flex_all_wavelengths)
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Wavelength: %f"%average_wavelength)

    # save the average wavelength to the phil parameters
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Average wavelength (%f A) is saved to phil parameters"%average_wavelength)
    if not 'average_wavelength' in (self.params.statistics).__dict__:
      self.params.statistics.__inject__('average_wavelength', average_wavelength)
    else:
      self.params.statistics.__setattr__('average_wavelength', average_wavelength)

    self.logger.log_step_time("BEAM_STATISTICS", True)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(beam_statistics)
