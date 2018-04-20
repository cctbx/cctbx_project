from __future__ import division

from dials.command_line.stills_process import Processor

class DialsProcessorWithLogging(Processor):
  '''Overrides for steps of dials processing of stills with XFEL GUI database logging.'''

  def __init__(self, params, composite_tag = None):
    super(DialsProcessorWithLogging, self).__init__(params, composite_tag)
    self.tt_low = None
    self.tt_high = None

  def log_frame(self, experiments, reflections, run, n_strong, timestamp = None,
                two_theta_low = None, two_theta_high = None, db_event = None):
    # update an existing db_event if db_event is not None
    if self.params.experiment_tag is None:
      return
    from xfel.ui.db.dxtbx_db import log_frame
    db_event = log_frame(experiments, reflections, self.params, run, n_strong, timestamp = timestamp,
                         two_theta_low = two_theta_low, two_theta_high = two_theta_high,
                         db_event = db_event)
    return db_event

  def get_run_and_timestamp(self, obj):
    try:
      sets = obj.extract_imagesets()
    except AttributeError:
      sets = obj.imagesets()
    assert len(sets) == 1
    imageset = sets[0]
    assert len(imageset) == 1
    format_obj = imageset.data().reader().format_class._current_instance_ # XXX
    run = format_obj.get_run(imageset.indices()[0])
    timestamp = format_obj.get_psana_timestamp(imageset.indices()[0])
    return run, timestamp

  def pre_process(self, datablock):
    super(DialsProcessorWithLogging, self).pre_process(datablock)

    if self.params.radial_average.enable:
      from dxtbx.command_line.radial_average import run as radial_run
      from scitbx.array_family import flex
      imageset = datablock.extract_imagesets()[0]
      two_thetas, radial_average_values = radial_run(self.params.radial_average, imageset = imageset)

      def get_closest_idx(data, val):
        deltas = flex.abs(data - val)
        return flex.first_index(deltas, flex.min(deltas))

      if self.params.radial_average.two_theta_low is not None:
        self.tt_low = radial_average_values[get_closest_idx(two_thetas, self.params.radial_average.two_theta_low)]

      if self.params.radial_average.two_theta_high is not None:
        self.tt_high = radial_average_values[get_closest_idx(two_thetas, self.params.radial_average.two_theta_high)]

  def find_spots(self, datablock):
    observed = super(DialsProcessorWithLogging, self).find_spots(datablock)
    run, timestamp = self.get_run_and_timestamp(datablock)
    self.db_event = self.log_frame(None, None, run, len(observed), timestamp = timestamp,
                                   two_theta_low = self.tt_low, two_theta_high = self.tt_high)
    return observed

  def integrate(self, experiments, indexed):
    # Results should always be logged after integration even if it fails.
    run, timestamp = self.get_run_and_timestamp(experiments)
    try:
      integrated = super(DialsProcessorWithLogging, self).integrate(experiments, indexed)
      self.log_frame(experiments, integrated, run, len(integrated), timestamp = timestamp,
                     two_theta_low = self.tt_low, two_theta_high = self.tt_high,
                     db_event = self.db_event)
      return integrated
    except Exception as e:
      self.log_frame(experiments, indexed, run, len(indexed), timestamp = timestamp,
                     two_theta_low = self.tt_low, two_theta_high = self.tt_high,
                     db_event = self.db_event)
      raise e
