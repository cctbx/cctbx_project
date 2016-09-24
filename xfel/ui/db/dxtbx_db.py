from __future__ import division
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.experiment import Experiment, Event, Cell_Bin
from scitbx.array_family import flex

def log_frame(experiments, reflections, params, run, n_strong, timestamp = None, two_theta_low = None, two_theta_high = None):
  app = dxtbx_xfel_db_application(params)
  db_run = app.get_run(run_number=run)
  if params.input.trial is None:
    db_trial = app.get_trial(trial_id = params.input.trial_id)
    params.input.trial = db_trial.trial
  else:
    db_trial = app.get_trial(trial_number = params.input.trial)

  if params.input.rungroup is None:
    db_event = app.create_event(timestamp = timestamp, run_id = db_run.id, trial_id = db_trial.id, n_strong = n_strong,
      two_theta_low = two_theta_low, two_theta_high = two_theta_high)
  else:
    db_event = app.create_event(timestamp = timestamp, run_id = db_run.id, trial_id = db_trial.id, rungroup_id = params.input.rungroup, n_strong = n_strong,
      two_theta_low = two_theta_low, two_theta_high = two_theta_high)

  if experiments is not None:
    assert len(experiments) == 1
    db_experiment = app.create_experiment(experiments[0])
    app.link_imageset_frame(db_experiment.imageset, db_event)

    d = experiments[0].crystal.get_unit_cell().d(reflections['miller_index'])
    for db_bin in db_experiment.crystal.cell.bins: # will be [] if there are no isoforms
      sel = (d <= float(db_bin.d_max)) & (d > float(db_bin.d_min))
      sel &= reflections['intensity.sum.value'] > 0
      refls = reflections.select(sel)
      n_refls = len(refls)
      Cell_Bin(app,
               count = n_refls,
               bin_id = db_bin.id,
               crystal_id = db_experiment.crystal.id,
               avg_intensity = flex.mean(refls['intensity.sum.value']) if n_refls > 0 else None,
               avg_sigma = flex.mean(flex.sqrt(refls['intensity.sum.variance'])) if n_refls > 0 else None,
               avg_i_sigi = flex.mean(refls['intensity.sum.value'] /
                                      flex.sqrt(refls['intensity.sum.variance'])) if n_refls > 0 else None)

class dxtbx_xfel_db_application(xfel_db_application):
  def create_experiment(self, experiment):
    return Experiment(self, experiment=experiment)

  def create_event(self, **kwargs):
    return Event(self, **kwargs)

  def get_event(self, event_id):
    return Event(self, event_id=event_id)

  def link_imageset_frame(self, imageset, event):
    query = "INSERT INTO `%s_imageset_event` (imageset_id, event_id, event_run_id) VALUES (%d, %d, %d)" % (
      self.params.experiment_tag, imageset.id, event.id, event.run_id)
    self.execute_query(query, commit=True)
