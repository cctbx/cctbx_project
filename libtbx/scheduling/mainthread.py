"""
Mainthread

Offers low-overhead execution on a single thread

Methods:
  has_results(): returns whether there are any results available
  results(): return an iterator yielding ( identifier, result ) tuples
  submit(target, args = (), kwargs = {}): submits job, return identifier
  is_empty(): returns whether there are no more jobs or results waiting
  is_full(): returns whether there are any jobs waiting
  shutdown(): no effect
  join(): finish processing all jobs
  terminate(): no effect
"""

from __future__ import absolute_import, division, print_function

from collections import deque

from libtbx.scheduling import identifier
from libtbx.scheduling import result


class manager(object):
  """
  Mainthread manager
  """

  def __init__(self):

    self.outqueue = deque()
    self.inqueue = deque()


  def job_count(self):

    return len( self.outqueue )


  def has_results(self):

    return self.inqueue


  def results(self):

    while not self.is_empty():
      while not self.has_results():
        self.poll()

      yield self.inqueue.popleft()


  def submit(self, target, args = (), kwargs = {}):

    jobid = identifier()
    self.outqueue.append( ( jobid, target, args, kwargs ) )
    return jobid


  def is_empty(self):

    return not self.outqueue and not self.inqueue


  def is_full(self):

    return self.outqueue


  def join(self):

    while self.outqueue:
      self.poll()


  def shutdown(self):

    pass


  def resume(self):

    pass


  def terminate(self):

    pass


  # Internal methods
  def poll(self):

    if self.outqueue:
      ( jobid, target, args, kwargs ) = self.outqueue.popleft()

      try:
        value = target( *args, **kwargs )

      except Exception as e:
        res = result.error( exception = e, traceback = result.get_traceback_info() )

      else:
        res = result.success( value = value )

      self.inqueue.append( ( jobid, res ) )


class creator(object):
  """
  Creator for manager

  Note this is a singleton object
  """

  @staticmethod
  def create():

    return manager()


  @staticmethod
  def destroy(manager):

    pass
