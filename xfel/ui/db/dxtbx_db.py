from __future__ import division
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.experiment import Experiment

def log_frame(experiments, reflections, params, timestamp = None):
  assert len(experiments) == 1

  app = dxtbx_xfel_db_application(params)
  db_experiment = app.create_experiment(experiments[0])
  pass

class dxtbx_xfel_db_application(xfel_db_application):
  def create_experiment(self, experiment):
    return Experiment(self, experiment=experiment)
