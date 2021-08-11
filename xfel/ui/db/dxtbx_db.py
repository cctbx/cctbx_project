from __future__ import absolute_import, division, print_function
from xfel.ui.db.xfel_db import xfel_db_application
from xfel.ui.db.experiment import Imageset, Experiment, Event, Bin, Cell_Bin, Cell, Detector, Crystal, Beam
from scitbx.array_family import flex

def log_frame(experiments, reflections, params, run, n_strong, timestamp = None,
              two_theta_low = None, two_theta_high = None, db_event = None, app = None, trial = None):
  if app is None:
    app = dxtbx_xfel_db_application(params, mode = 'cache_commits')
  else:
    app.mode = 'cache_commits'

  if isinstance(run, int) or isinstance(run, str):
    db_run = app.get_run(run_number=run)
  else:
    db_run = run

  if trial is None:
    if params.input.trial is None:
      db_trial = app.get_trial(trial_id = params.input.trial_id)
      params.input.trial = db_trial.trial
    else:
      db_trial = app.get_trial(trial_number = params.input.trial)
  else:
    db_trial = trial

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

  inserts = ""

  if app.last_query is None:
    app.last_query = ""

  def save_last_id(name):
    nonlocal inserts
    inserts += app.last_query + ";\n"
    inserts += "SELECT LAST_INSERT_ID() INTO @%s_id;\n"%name

  if experiments:
    save_last_id('event')
  else:
    inserts += app.last_query + ";\n"

  for i, experiment in enumerate(experiments or []):
    reflections_i = reflections.select(reflections['id']==i)

    imageset = Imageset(app)
    save_last_id('imageset')

    beam = Beam(app, beam = experiment.beam)
    save_last_id('beam')

    detector = Detector(app, detector = experiment.detector)
    save_last_id('detector')

    cell = Cell(app, crystal=experiment.crystal, isoform_id = None)
    save_last_id('cell')

    crystal = Crystal(app, crystal = experiment.crystal, make_cell = False, cell_id = "@cell_id")
    save_last_id('crystal')

    inserts += ("INSERT INTO `%s_experiment` (imageset_id, beam_id, detector_id, crystal_id, crystal_cell_id) " + \
                "VALUES (@imageset_id, @beam_id, @detector_id, @crystal_id, @cell_id);\n") % (
      params.experiment_tag)

    inserts += "INSERT INTO `%s_imageset_event` (imageset_id, event_id, event_run_id) VALUES (@imageset_id, @event_id, %d);\n" % (
      params.experiment_tag, db_run.id)

    d = experiment.crystal.get_unit_cell().d(reflections['miller_index']).select(reflections['id']==i)
    from cctbx.crystal import symmetry
    cs = symmetry(unit_cell = experiment.crystal.get_unit_cell(), space_group = experiment.crystal.get_space_group())
    mset = cs.build_miller_set(anomalous_flag=False, d_min=db_trial.d_min)
    n_bins = 10 # FIXME use n_bins as an attribute on the trial table
    binner = mset.setup_binner(n_bins=n_bins)
    for i in binner.range_used():
      d_max, d_min = binner.bin_d_range(i)
      Bin(app, number = i, d_min = d_min, d_max = d_max,
          total_hkl = binner.counts_complete()[i], cell_id = '@cell_id')
      save_last_id('bin')

      sel = (d <= float(d_max)) & (d > float(d_min))
      sel &= reflections_i['intensity.sum.value'] > 0
      refls = reflections_i.select(sel)
      n_refls = len(refls)
      Cell_Bin(app,
               count = n_refls,
               bin_id = '@bin_id',
               crystal_id = '@crystal_id',
               avg_intensity = flex.mean(refls['intensity.sum.value']) if n_refls > 0 else None,
               avg_sigma = flex.mean(flex.sqrt(refls['intensity.sum.variance'])) if n_refls > 0 else None,
               avg_i_sigi = flex.mean(refls['intensity.sum.value'] /
                                      flex.sqrt(refls['intensity.sum.variance'])) if n_refls > 0 else None)
      inserts += app.last_query + ";\n"
  app.mode = 'execute'
  return inserts

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
