from __future__ import division

from libtbx.queuing_system_utils import communication

import threading
import Queue
import time

class SchedulerActor(object):
  """
  Interface to a scheduling.Manager object that is polled by another thread
  """

  def __init__(self, manager, waittime):

    self.manager = manager
    self.waittime = waittime
    self.jobs_in = Queue.Queue()
    self.jobs_out = Queue.Queue()

    self.identifier_for = {}
    self.active = True
    self.daemon = threading.Thread( target = self.run )
    self.daemon.daemon = True
    self.daemon.start()


  def run(self):

    while self.active:
      self.manager.poll()

      while self.manager.completed_results:
        current = self.manager.completed_results.popleft()
        assert current.identifier in self.identifier_for

        try:
          self.jobs_out.put(
            ( self.identifier_for[ current.identifier ], current ),
            block = False,
            )

        except Queue.Full:
          self.manager.completed_results.appendleft( current )
          break

        del self.identifier_for[ current.identifier ]

      while True:
        try:
          (identifier, ( target, args, kwargs ) ) = self.jobs_in.get( block = False )
          i = self.manager.submit( target = target, args = args, kwargs = kwargs )
          assert i not in self.identifier_for
          self.identifier_for[ i ] = identifier

        except Queue.Empty:
          break

      if not self.manager.completed_results:
        time.sleep( self.waittime )

    self.manager.join()

    for result in self.manager.results:
      assert result.identifier in self.identifier_for
      self.jobs_out.put( ( self.identifier_for[ result.identifier ], result ) )
      del self.identifier_for[ result.identifier ]


  def shutdown(self):

    self.active = False
    self.daemon.join()


class SubmitJob(communication.Command):
  """
  Request for job submission
  """

  def __init__(self, identifier, target, args, kwargs):

    self.identifier = identifier
    self.callparms = ( target, args, kwargs )


  def process(self, environment):

    environment.jobs_in.put( ( self.identifier, self.callparms ) )


class GetJob(communication.Command):

  def process(self, environment):

    ( identifier, result ) = environment.jobs_out.get( block = False )
    return ( identifier, result.value )


class SchedulerServer(communication.Server):
  """
  Server-side component
  """

  def __init__(self, instream, outstream, manager, waittime = 0.1):

    super( SchedulerServer, self ).__init__(
      instream = instream,
      outstream = outstream,
      environment = SchedulerActor( manager = manager, waittime = waittime ),
      )


class SchedulerClient(communication.Client):
  """
  Client-side connection to a potentially remote SchedulerServer
  """

  def __init__(self, instream, outstream):

    super(SchedulerClient, self).__init__(
      instream = instream,
      outstream = outstream,
      )
    self.running = set()
    self.result_for = {}


  def is_known(self, identifier):

    return identifier in self.running or identifier in self.result_for


  def submit(self, identifier, target, args, kwargs):

    assert not self.is_known( identifier = identifier )

    response = self.send(
      command = SubmitJob(
        identifier = identifier,
        target = target,
        args = args,
        kwargs = kwargs,
        ),
      )

    try:
      response() # successful if no exception raised

    except Exception, e:
      raise RuntimeError, "Submission failure: %s" % e

    self.running.add( identifier )


  def is_alive(self, identifier):

    self.poll()
    return identifier not in self.result_for


  def get_result(self, identifier):

    return self.result_for[ identifier ]


  def remove(self, identifier):

    del self.result_for[ identifier ]


  def poll(self):

    while True:
      response = self.send( command = GetJob() )

      try:
        ( identifier, value ) = response()

      except Queue.Empty:
        break

      assert identifier in self.running
      self.running.remove( identifier )

      assert identifier not in self.result_for
      self.result_for[ identifier ] = value


  def Job(self, target, args = (), kwargs = {}):

    return Job(
      handler = self,
      target = target,
      args = args,
      kwargs = kwargs,
      )


class PreSubmissionStatus(object):
  """
  A job that has not been submitted yet
  """

  @staticmethod
  def start(job):

    job.handler.submit(
      identifier = job.name,
      target = job.target,
      args = job.args,
      kwargs = job.kwargs,
      )

    job.status = RunningStatus


  @staticmethod
  def is_alive(job):

    raise RuntimeError, "job has not been submitted yet"


class RunningStatus(object):
  """
  A job that has been submitted
  """

  @staticmethod
  def start(job):

    raise RuntimeError, "start called second time"


  @staticmethod
  def is_alive(job):

    alive = job.handler.is_alive( identifier = job.name )

    if not alive:
      job.result = job.handler.get_result( identifier = job.name )
      job.handler.remove( identifier = job.name )
      job.status = FinishedStatus

    return alive


class FinishedStatus(object):
  """
  A job that has been submitted
  """

  @staticmethod
  def start(job):

    raise RuntimeError, "start called second time"


  @staticmethod
  def is_alive(job):

    return False


class Job(object):
  """
  Job object to execute function calls on remote machines accessible via
  a network channel

  Restrictions: target, args and kwargs has to be pickleable
  """

  def __init__(self, handler, target, args = (), kwargs = {}):

    self.handler = handler

    self.target = target
    self.args = args
    self.kwargs = kwargs

    self.status = PreSubmissionStatus
    self.result = None
    self.exitcode = 0


  @property
  def name(self):

    return "remote_job_%d" % id( self )


  def start(self):

    self.status.start( job = self )


  def is_alive(self):

    return self.status.is_alive( job = self )


  def join(self):

    while self.is_alive():
      time.sleep( 0.1 )


  def __str__(self):

    return "%s(name = '%s')" % ( self.__class__.__name__, self.name )

