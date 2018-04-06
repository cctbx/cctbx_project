from __future__ import division

from dials.command_line.stills_process import Processor

class DialsProcessorWithLogging(Processor):
  '''Overrides for steps of dials processing of stills with XFEL GUI database logging.
  A run object and a params instance attribute are required.'''

  def log_frame(self, experiments, reflections, run, n_strong, timestamp = None,
                two_theta_low = None, two_theta_high = None, db_event = None):
    # update an existing db_event if db_event is not None
    if self.params.experiment_tag is None:
      return
    try:
      from xfel.ui.db.dxtbx_db import log_frame
      log_frame(experiments, reflections, self.params, run, n_strong, timestamp = timestamp,
                two_theta_low = two_theta_low, two_theta_high = two_theta_high,
                db_event = db_event)
    except Exception, e:
      import traceback; traceback.print_exc()
      print str(e), "event", timestamp
      if reflections is None:
        self.debug_write("db_logging_failed", "fail")
      else:
        self.debug_write("db_logging_failed_%d" % len(reflections), "fail")

  def find_spots(self, datablock, run, timestamp = None,
                 two_theta_low = None, two_theta_high = None):
    observed = super(DialsProcessorWithLogging, self).find_spots(datablock)
    db_event = self.log_frame(None, None, run.run(), len(observed), timestamp = timestamp,
                              two_theta_low = two_theta_low, two_theta_high = two_theta_high)
    return (observed, db_event)

  def index(self, datablock, observed):
    # Results are never logged directly after indexing, only after the indexing solution
    # is refined.
    experiments, indexed = super(DialsProcessorWithLogging, self).index(datablock, observed)
    return (experiments, indexed)

  def refine(self, experiments, indexed, run, timestamp = None,
             two_theta_low = None, two_theta_high = None,
             db_event = None, log = False):
    # Default is log==False so that the integration step can log the result, *whether or not
    # it is successful.* This prevents duplicate indexing and integration entries that would
    # be misunderstood to be multiple lattices.
    try:
      experiments, indexed = super(DialsProcessorWithLogging, self).refine(experiments, indexed)
      if log:
        db_event = self.log_frame(experiments, indexed, run.run(), len(observed),
                                  timestamp = timestamp,
                                  two_theta_low = two_theta_low, two_theta_high = two_theta_high,
                                  db_event = db_event)
      return (experiments, indexed, db_event)
    except Exception as e:
      if log:
        db_event = self.log_frame(experiments, indexed, run.run(), len(observed),
                                  timestamp = timestamp,
                                  two_theta_low = two_theta_low, two_theta_high = two_theta_high,
                                  db_event = db_event)
      raise e

  def integrate(self, experiments, indexed, run, timestamp = None,
                two_theta_low = None, two_theta_high = None, db_event = None):
    # Results should always be logged after integration even if it fails.
    try:
      integrated = super(DialsProcessorWithLogging, self).integrate(experiments, indexed)
      db_event = self.log_frame(experiments, integrated, run.run(), len(), timestamp = timestamp,
                                two_theta_low = two_theta_low, two_theta_high = two_theta_high,
                                db_event = db_event)
      return (integrated, db_event)
    except Exception as e:
      db_event = self.log_frame(experiments, indexed, run.run(), len(indexed), timestamp = timestamp,
                                two_theta_low = two_theta_low, two_theta_high = two_theta_high,
                                db_event = db_event)
      raise e
