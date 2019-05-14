from __future__ import absolute_import, division, print_function

"""
Base classes for the merging application
"""

from xfel.merging.application.mpi_helper import mpi_helper
from xfel.merging.application.mpi_logger import mpi_logger

class worker(object):
  """ Base class for the worker objects. Performs validation and does work using the run method """
  def __init__(self, params):
    self.params = params

    # create logger
    self.logger = mpi_logger(self.params)

    # create MPI helper
    self.mpi_helper = mpi_helper()

  def __repr__(self):
    return 'Unknown'

  def validate(self):
    """ Override to perform any validation of the input parameters """
    pass

  def run(self, experiments, reflections):
    """ Process the data """
    pass

class factory(object):
  """ Constructs worker objects """

  @staticmethod
  def from_parameters(param, additional_info=None):
    """ Construct a list of workers given the params object. The list contains all workers
        that comprise a single step, in the order that they will be executed """
    pass

def exercise_worker(worker_class):
  """ Boilerplate code for testing a worker class """
  from xfel.merging.application.phil.phil import phil_scope
  from dials.util.options import OptionParser
  from dxtbx.model.experiment_list import ExperimentListDumper
  # Create the parser
  parser = OptionParser(phil=phil_scope)

  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)

  # Load the data for the worker
  if 'simple_file_loader' in str(worker_class):
    experiments = reflections = None
  else:
    from xfel.merging.application.input.file_loader import simple_file_loader
    loader = simple_file_loader(params)
    loader.validate()
    experiments, reflections = loader.run(None, None)

  worker = worker_class(params)
  worker.validate()
  experiments, reflections = worker.run(experiments, reflections)

  prefix = worker_class.__name__
  reflections.as_msgpack_file(prefix + ".mpack")
  ExperimentListDumper(experiments).as_file(prefix + ".json")

