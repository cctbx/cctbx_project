from __future__ import absolute_import, division, print_function

from libtbx.resource_monitor import ResourceMonitor
from xfel.merging.application.worker import worker


resource_monitor = None


class MonitorWorker(worker):
  """Monitor and plot CPU and GPU resources during `cctbx.xfel.merge`"""

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    kwargs = dict(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)
    super(MonitorWorker, self).__init__(**kwargs)

  def __repr__(self):
    return 'Monitor CPU and GPU resources'

  def run(self, experiments, reflections):
    global resource_monitor
    if resource_monitor is None:
      resource_monitor = ResourceMonitor(
        detail=self.params.monitor.detail,
        period=self.params.monitor.period,
        plot=self.params.monitor.plot,
        prefix=self.params.monitor.prefix,
        write=self.params.monitor.write,
      )
    if resource_monitor.active:
      self.logger.log('Stopping resource monitor')
      resource_monitor.stop()
    else:  # if resource_monitor.active == False
      self.logger.log('Starting resource monitor')
      resource_monitor.start()
    return experiments, reflections
