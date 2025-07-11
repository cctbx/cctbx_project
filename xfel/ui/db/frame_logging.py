from __future__ import absolute_import, division, print_function

from dials.command_line.stills_process import Processor
from xfel.ui.db.dxtbx_db import log_frame, dxtbx_xfel_db_application
from xfel.ui.db.run import Run
from xfel.ui.db.trial import Trial

class DialsProcessorWithLogging(Processor):
  '''Overrides for steps of dials processing of stills with XFEL GUI database logging.'''

  def __init__(self, params, composite_tag = None, rank = 0):
    self.debug_file_handle = None
    super(DialsProcessorWithLogging, self).__init__(params, composite_tag, rank)
    self.tt_low = None
    self.tt_high = None
    if self.params.experiment_tag is None:
      return

    assert params.db.logging_batch_size >= 1

    from libtbx.mpi4py import MPI
    comm = MPI.COMM_WORLD

    self.queries = []
    self.rank = rank

    if comm.size > 1:
      if rank == 0:
        db_app = dxtbx_xfel_db_application(params, cache_connection=False)
        run = db_app.get_run(run_number=self.params.input.run_num)
        run_id = run.id
        rund = run._db_dict
        trial = db_app.get_trial(trial_number = params.input.trial)
        trial_id = trial.id
        triald = trial._db_dict
      else:
        db_app = None
        run_id = None
        rund = None
        trial_id = None
        triald = None
      self.db_app, run_id, rund, trial_id, triald = comm.bcast((db_app, run_id, rund, trial_id, triald), root=0)
      self.run = Run(self.db_app, run_id = run_id, **rund)
      self.trial = Trial(self.db_app, trial_id = trial_id, **triald)
    else:
      self.db_app = dxtbx_xfel_db_application(params, cache_connection=True)
      self.run = self.db_app.get_run(run_number=self.params.input.run_num)
      self.trial = self.db_app.get_trial(trial_number = params.input.trial)
    self.db_app.mode = 'cache_commits'
    self.n_strong = None

  def finalize(self):
    super(DialsProcessorWithLogging, self).finalize()
    if self.params.experiment_tag is None:
      return
    self.log_batched_frames()
    self.trial = None

  def log_batched_frames(self):
    current_run = self.params.input.run_num
    current_dbrun = self.run
    inserts = "BEGIN;\n" # start a transaction
    for q in self.queries:
      experiments, reflections, run, n_strong, timestamp, two_theta_low, two_theta_high, db_event = q
      if run != current_run:
        self.db_app.mode = "execute"
        current_run = run
        current_dbrun = self.db_app.get_run(run_number=run)
        self.db_app.mode = "cache_commits"

      inserts += log_frame(experiments, reflections, self.params, current_dbrun, n_strong, timestamp = timestamp,
                           two_theta_low = two_theta_low, two_theta_high = two_theta_high,
                           db_event = db_event, app = self.db_app, trial = self.trial)
    inserts += "COMMIT;\n"

    # patch up query so for example '@row_id' becomes @row_id
    newinserts = []
    for line in inserts.split('\n'):
      if '@' in line:
        newline = []
        for word in line.split(' '):
          if '@' in word:
            word = word.replace("'", "")
          newline.append(word)
        line = ' '.join(newline)
      newinserts.append(line)
    inserts = '\n'.join(newinserts)

    self.db_app.execute_query(inserts, commit=False) # transaction, so don't commit twice
    self.queries = []

  def log_frame(self, experiments, reflections, run, n_strong, timestamp = None,
                two_theta_low = None, two_theta_high = None, db_event = None):
    # update an existing db_event if db_event is not None
    if self.params.experiment_tag is None:
      return
    self.queries.append((experiments, reflections, run, n_strong, timestamp,
                         two_theta_low, two_theta_high,
                         db_event))
    if len(self.queries) >= self.params.db.logging_batch_size:
      self.log_batched_frames()
    return db_event

  def get_run_and_timestamp(self, obj):
    sets = obj.imagesets()
    assert len(sets) == 1
    imageset = sets[0]
    assert len(imageset) == 1
    format_obj = imageset.data().reader().format_class._current_instance_ # XXX
    try: # XTC specific version
      import psana
      run = str(format_obj.get_run_from_index(imageset.indices()[0]).run())
      timestamp = format_obj.get_psana_timestamp(imageset.indices()[0])
      evt = format_obj._get_event(imageset.indices()[0])
      if evt:
        fid = evt.get(psana.EventId).fiducials()
        timestamp += ", fid:" + str(fid)
      return run, timestamp
    except (ImportError, AttributeError): # General version
      run = self.params.input.run_num
      timestamp = self.tag
      return run, timestamp

  def pre_process(self, experiments):
    super(DialsProcessorWithLogging, self).pre_process(experiments)
    if self.params.radial_average.enable:
      try:
        self.radial_average(experiments)
      except Exception as e:
        run, timestamp = self.get_run_and_timestamp(experiments)
        self.log_frame(experiments, None, run, 0, timestamp = timestamp)
        raise e

  def radial_average(self, experiments):
      from dxtbx.command_line.radial_average import run as radial_run
      from scitbx.array_family import flex

      if self.params.radial_average.verbose:
        run, timestamp = self.get_run_and_timestamp(experiments)
        if timestamp == self.tag:
          print("Radial average of run %s, timestamp %s"%(str(run), self.tag))
        else:
          print("Radial average of run %s, tag %s, timestamp %s"%(str(run), self.tag, timestamp))

      imageset = experiments.imagesets()[0]
      two_thetas, radial_average_values = radial_run(self.params.radial_average, imageset = imageset)

      def get_closest_idx(data, val):
        deltas = flex.abs(data - val)
        return flex.first_index(deltas, flex.min(deltas))

      if self.params.radial_average.two_theta_low is not None:
        self.tt_low = radial_average_values[get_closest_idx(two_thetas, self.params.radial_average.two_theta_low)]

      if self.params.radial_average.two_theta_high is not None:
        self.tt_high = radial_average_values[get_closest_idx(two_thetas, self.params.radial_average.two_theta_high)]

      if self.params.radial_average.verbose:
        print("Two theta low and high for run %s, timestamp %s: %f, %f"%(str(run), timestamp, self.tt_low, self.tt_high))

  def find_spots(self, experiments):
    observed = super(DialsProcessorWithLogging, self).find_spots(experiments)
    self.n_strong = len(observed)
    if (
        not self.params.dispatch.index or
        len(observed) < self.params.dispatch.hit_finder.minimum_number_of_reflections
    ):
      run, timestamp = self.get_run_and_timestamp(experiments)
      self.log_frame(experiments, None, run, self.n_strong, timestamp = timestamp,
                     two_theta_low = self.tt_low, two_theta_high = self.tt_high)
    return observed

  def index(self, experiments, reflections):
    try:
      experiments, indexed = super(DialsProcessorWithLogging, self).index(experiments, reflections)
    except Exception as e:
      run, timestamp = self.get_run_and_timestamp(experiments)
      self.log_frame(experiments, None, run, self.n_strong, timestamp = timestamp,
                     two_theta_low = self.tt_low, two_theta_high = self.tt_high)
      raise e
    else:
      if not self.params.dispatch.integrate:
        run, timestamp = self.get_run_and_timestamp(experiments)
        self.log_frame(experiments, indexed, run, self.n_strong, timestamp = timestamp,
                       two_theta_low = self.tt_low, two_theta_high = self.tt_high)
    return experiments, indexed

  def integrate(self, experiments, indexed):
    # Results should always be logged after integration even if it fails.
    run, timestamp = self.get_run_and_timestamp(experiments)
    try:
      integrated = super(DialsProcessorWithLogging, self).integrate(experiments, indexed)
      self.log_frame(experiments, integrated, run, self.n_strong, timestamp = timestamp,
                     two_theta_low = self.tt_low, two_theta_high = self.tt_high)
      return integrated
    except Exception as e:
      self.log_frame(experiments, indexed, run, self.n_strong, timestamp = timestamp,
                     two_theta_low = self.tt_low, two_theta_high = self.tt_high)
      raise e
    else:
      self.log_frame(experiments, indexed, run, len(indexed), timestamp = timestamp,
                     two_theta_low = self.tt_low, two_theta_high = self.tt_high)
