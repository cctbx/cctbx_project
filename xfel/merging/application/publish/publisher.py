from __future__ import division
from xfel.merging.application.worker import worker
from xfel.command_line.upload_mtz import run_with_preparsed as run_publish
import os

class publisher(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(publisher, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return "Publish output mtz and main log to a Google Drive folder"

  def run(self, experiments, reflections):

    mtz_fname = self.params.output.prefix + "_all.mtz"
    mtz_path = os.path.join(self.params.output.output_dir, mtz_fname)
    self.params.publish.input.mtz_file = mtz_path
    self.params.publish.input.log_file = None
    self.params.publish.input.version = None

    if self.mpi_helper.rank==0: run_publish(self.params.publish)

    return experiments, reflections
