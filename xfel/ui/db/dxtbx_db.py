from __future__ import absolute_import, division, print_function
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.experiment import Experiment, Event, Bin, Cell_Bin
from scitbx.array_family import flex

def log_frame(experiments, reflections, params, run, n_strong, timestamp = None,
              two_theta_low = None, two_theta_high = None, db_event = None):
  app = dxtbx_xfel_db_application(params)
  db_run = app.get_run(run_number=run)
  if params.input.trial is None:
    db_trial = app.get_trial(trial_id = params.input.trial_id)
    params.input.trial = db_trial.trial
  else:
    db_trial = app.get_trial(trial_number = params.input.trial)

  if db_event is None:
    if params.input.rungroup is None:
      db_event = app.create_event(timestamp = timestamp,
                                  run_id = db_run.id,
                                  trial_id = db_trial.id,
                                  n_strong = n_strong,
                                  two_theta_low = two_theta_low,
                                  two_theta_high = two_theta_high)
    else:
      db_event = app.create_event(timestamp = timestamp,
                                  run_id = db_run.id,
                                  trial_id = db_trial.id,
                                  rungroup_id = params.input.rungroup,
                                  n_strong = n_strong,
                                  two_theta_low = two_theta_low,
                                  two_theta_high = two_theta_high)

  for i, experiment in enumerate(experiments or []):
    reflections_i = reflections.select(reflections['id']==i)
    db_experiment = app.create_experiment(experiment)
    app.link_imageset_frame(db_experiment.imageset, db_event)

    d = experiment.crystal.get_unit_cell().d(reflections['miller_index']).select(reflections['id']==i)
    if len(db_experiment.crystal.cell.bins) == 0: # will be [] if there are no isoforms and no target cells
      from cctbx.crystal import symmetry
      cs = symmetry(unit_cell = db_experiment.crystal.cell.unit_cell, space_group_symbol = db_experiment.crystal.cell.lookup_symbol)
      mset = cs.build_miller_set(anomalous_flag=False, d_min=db_trial.d_min)
      n_bins = 10 # FIXME use n_bins as an attribute on the trial table
      binner = mset.setup_binner(n_bins=n_bins)
      for i in binner.range_used():
        d_max, d_min = binner.bin_d_range(i)
        Bin(app, number = i, d_min = d_min, d_max = d_max,
            total_hkl = binner.counts_complete()[i], cell_id = db_experiment.crystal.cell.id)
      db_experiment.crystal.cell.bins = app.get_cell_bins(db_experiment.crystal.cell.id)
      assert len(db_experiment.crystal.cell.bins) == n_bins

    for db_bin in db_experiment.crystal.cell.bins:
      sel = (d <= float(db_bin.d_max)) & (d > float(db_bin.d_min))
      sel &= reflections_i['intensity.sum.value'] > 0
      refls = reflections_i.select(sel)
      n_refls = len(refls)
      Cell_Bin(app,
               count = n_refls,
               bin_id = db_bin.id,
               crystal_id = db_experiment.crystal.id,
               avg_intensity = flex.mean(refls['intensity.sum.value']) if n_refls > 0 else None,
               avg_sigma = flex.mean(flex.sqrt(refls['intensity.sum.variance'])) if n_refls > 0 else None,
               avg_i_sigi = flex.mean(refls['intensity.sum.value'] /
                                      flex.sqrt(refls['intensity.sum.variance'])) if n_refls > 0 else None)
  return db_event

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
