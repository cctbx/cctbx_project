from __future__ import division

"""
Base classes for the merging application
"""

class worker(object):
  """ Base class for the worker objects. Performs validation and does work using the run method """
  def __init__(self, params):
    self.params = params

  def validate(self):
    """ Override to perform any validation of the input parameters """
    pass

  def run(self, experiments, reflections):
    """ Process the data """
    pass

class factory(object):
  """ Constructs worker objects """

  @staticmethod
  def from_parameters(param):
    """ Construct a list of workers given the params object. The list contains all workers
        that comprise a single step, in the order that they will be executed """
    pass
