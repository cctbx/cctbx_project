from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dxtbx.imageset import ImageSetFactory
from dxtbx.model.experiment_list import ExperimentList
from dials.array_family import flex
import os

from dials.command_line.stills_process import Processor
class integrate_only_processor(Processor):
  def __init__(self, params):
    self.params = params

class integrate(worker):
  """
  Calls the stills process version of dials.integrate
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(integrate, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Integrate reflections'

  def run(self, experiments, reflections):
    from dials.util import log
    self.logger.log_step_time("INTEGRATE")

    logfile = os.path.splitext(self.logger.rank_log_file_path)[0] + "_integrate.log"
    log.config(logfile=logfile)
    processor = integrate_only_processor(self.params)

    # Re-generate the image sets using their format classes so we can read the raw data
    # Integrate the experiments one at a time to not use up memory
    all_integrated_expts = ExperimentList()
    all_integrated_refls = None
    current_imageset = None
    current_imageset_path = None
    for expt_id, expt in enumerate(experiments):
      assert len(expt.imageset.paths()) == 1 and len(expt.imageset) == 1
      self.logger.log("Starting integration experiment %d"%expt_id)
      refls = reflections.select(reflections['id'] == expt_id)
      if expt.imageset.paths()[0] != current_imageset_path:
        current_imageset_path = expt.imageset.paths()[0]
        current_imageset = ImageSetFactory.make_imageset(expt.imageset.paths())
      idx = expt.imageset.indices()[0]
      expt.imageset = current_imageset[idx:idx+1]
      idents = refls.experiment_identifiers()
      del idents[expt_id]
      idents[0] = expt.identifier
      refls['id'] = flex.int(len(refls), 0)

      try:
        integrated = processor.integrate(experiments[expt_id:expt_id+1], refls)
      except RuntimeError:
        self.logger.log("Error integrating expt %d"%expt_id)
        continue

      all_integrated_expts.append(expt)
      if all_integrated_refls:
        all_integrated_refls = flex.reflection_table.concat([all_integrated_refls, integrated])
      else:
        all_integrated_refls = integrated

    if all_integrated_refls is None:
      all_integrated_refls = flex.reflection_table()
    self.logger.log("Integration done, %d experiments, %d reflections" %
                    (len(all_integrated_expts), len(all_integrated_refls)))
    return all_integrated_expts, all_integrated_refls


if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(integrate)
