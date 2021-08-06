from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dxtbx.imageset import ImageSetFactory
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

    # make sure that all ranks make an imageset
    first_imageset_path = experiments[0].imageset.paths()[0] if experiments else None
    all_first_imageset_paths = self.mpi_helper.comm.gather(first_imageset_path, root=0)
    if self.mpi_helper.rank==0:
      all_first_imageset_paths = [x for x in all_first_imageset_paths if x is not None]
      dummy_imageset_path = all_first_imageset_paths[0]
    else:
      dummy_imageset_path = None
    current_imageset_path = self.mpi_helper.comm.bcast(dummy_imageset_path, root=0)
    current_imageset = ImageSetFactory.make_imageset([current_imageset_path])



    # Re-generate the image sets using their format classes so we can read the raw data
    # Integrate the experiments one at a time to not use up memory
    all_integrated = flex.reflection_table()
#    current_imageset = None
#    current_imageset_path = None
    for expt_id, expt in enumerate(experiments):
      assert len(expt.imageset.paths()) == 1 and len(expt.imageset) == 1
      self.logger.log("Starting integration experiment %d"%expt_id)
      refls = reflections.select(reflections['exp_id'] == expt.identifier)
      if expt.imageset.paths()[0] != current_imageset_path:
        current_imageset_path = expt.imageset.paths()[0]
        current_imageset = ImageSetFactory.make_imageset(expt.imageset.paths())
      idx = expt.imageset.indices()[0]
      expt.imageset = current_imageset[idx:idx+1]
      idents = refls.experiment_identifiers()
      del idents[expt_id]
      idents[0] = expt.identifier
      refls['id'] = flex.int(len(refls), 0)

      integrated = processor.integrate(experiments[expt_id:expt_id+1], refls)

      idents = integrated.experiment_identifiers()
      del idents[0]
      idents[expt_id] = expt.identifier
      integrated['id'] = flex.int(len(integrated), expt_id)
      all_integrated.extend(integrated)

    self.logger.log("Integration done, %d experiments, %d reflections"%(len(experiments), len(all_integrated)))
    return experiments, all_integrated

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(integrate)
