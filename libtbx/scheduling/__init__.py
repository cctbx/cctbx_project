"""
Generic job scheduling
"""

from __future__ import division


class SchedulingError(Exception):
  """
  Package exception
  """


class identifier(object):
  """
  Job identifier
  """

  def __init__(self):

    self.jobid = id( self )


  def __eq__(self, other):

    return self.jobid == other.jobid


  def __ne__(self, other):

    return not ( self == other )


  def __hash__(self):

    return hash( self.jobid )


  def __str__(self):

    return "identifier( jobid = %s )" % self.jobid


  def __repr__(self):

    return str( self )


class holder(object):
  """
  Context manager for a manager
  """

  def __init__(self, creator):

    self.creator = creator
    self.manager = None


  def __enter__(self):

    self.manager = self.creator.create()
    return self.manager


  def __exit__(self, exc_type, exc_val, exc_tb):

    self.creator.destroy( manager = self.manager )
    return False
