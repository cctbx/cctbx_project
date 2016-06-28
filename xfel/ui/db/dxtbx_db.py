from __future__ import division
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.experiment import Experiment, Event, Cell_Bin

def log_frame(experiments, reflections, params, run, timestamp = None):
  assert len(experiments) == 1

  app = dxtbx_xfel_db_application(params)
  db_experiment = app.create_experiment(experiments[0])
  db_run = app.get_run(run_number=run)
  db_trial = app.get_trial(trial_number=params.input.trial)
  db_event = app.create_event(timestamp = timestamp, run_id = db_run.id, trial_id = db_trial.id)
  app.link_imageset_frame(db_experiment.imageset, db_event)

  d = experiments[0].crystal.get_unit_cell().d(reflections['miller_index'])
  for db_bin in db_experiment.crystal.cell.bins: # will be [] if there are no isoforms
    sel = (d <= float(db_bin.d_max)) & (d > float(db_bin.d_min))
    Cell_Bin(app,
             count = len(reflections.select(sel)),
             bin_id = db_bin.id,
             crystal_id = db_experiment.crystal.id)

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
