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

    if self.params.mp.psana2_mode:
      self.check_psana2(split_comm=False)
      paths_ = list(set([p for iset in experiments.imagesets() for p in iset.paths()]))
      paths = list(set([p for plist in self.mpi_helper.comm.allgather(paths_) for p in plist]))
      assert len(paths) == 1

      current_imageset_path = paths[0]
      current_imageset = ImageSetFactory.make_imageset([current_imageset_path])

    for expt_id, expt in enumerate(experiments):
      assert len(expt.imageset.paths()) == 1 and len(expt.imageset) == 1
      self.logger.log("Starting integration experiment %d"%expt_id)
      if self.params.integration.recruitment.expand_nave_parameters: # NKS request predictions on larger envelope
        eta = expt.crystal.get_half_mosaicity_deg()
        deff = expt.crystal.get_domain_size_ang()
        factor = self.params.integration.recruitment.expansion_factor
        expt.crystal.set_domain_size_ang(deff/factor)
        expt.crystal.set_half_mosaicity_deg(eta*factor)
        self.logger.log("Expand focus experiment, half_mosaicity_deg=%8f-->%8f domain_size_ang=%.0f-->%.0f"%(
          eta,expt.crystal.get_half_mosaicity_deg(),deff,expt.crystal.get_domain_size_ang()))
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

      if self.params.integration.recruitment.expand_nave_parameters: # NKS Nave params on I/sigma spots, reapply envelope
        import math
        d_spacing = expt.crystal.get_unit_cell().d(integrated["miller_index"])
        isigma3 = (integrated["intensity.sum.value"]/flex.sqrt(integrated["intensity.sum.variance"])
                   )>self.params.integration.recruitment.significance_cutoff
        from serialtbx.mono_simulation.max_like import minimizer
        try:
          M = minimizer(
                d_i=d_spacing.select(isigma3),
                psi_i=integrated['delpsical.rad'].select(isigma3),
                eta_rad=(math.pi/180.) * (2.*expt.crystal.get_half_mosaicity_deg()),
                Deff=expt.crystal.get_domain_size_ang(),
          )
          new_half_mosaicity_deg = M.x[1] * 180.0 / (2.0 * math.pi)
          new_domain_size_ang = 2.0 / M.x[0]
          self.logger.log("  Refit Nave on 3 sigma, half_mosaicity_deg=%8f domain_size_ang=%.0f"%(
          new_half_mosaicity_deg,new_domain_size_ang))
          expt.crystal.set_half_mosaicity_deg(new_half_mosaicity_deg)
          expt.crystal.set_domain_size_ang(new_domain_size_ang)
          permitted_psi_envelope_rad = ((d_spacing/new_domain_size_ang) + (new_half_mosaicity_deg * math.pi/180.))
          accepted = flex.abs(integrated['delpsical.rad']) <= permitted_psi_envelope_rad
          integrated = integrated.select(accepted)

        except Exception as e:
          self.logger.log("Recovering "+str(e))
          expt.crystal.set_domain_size_ang(deff)
          expt.crystal.set_half_mosaicity_deg(eta)
          integrated = processor.integrate(experiments[expt_id:expt_id+1], refls)
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
