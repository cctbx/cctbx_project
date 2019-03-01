from __future__ import absolute_import, division, print_function

from libtbx.queuing_system_utils import communication
from six.moves import cPickle as pickle
import threading
from six.moves import queue
import time

class SchedulerActor(object):
  """
  Interface to a scheduling.Manager object that is polled by another thread
  """

  def __init__(self, manager, waittime):

    self.manager = manager
    self.waittime = waittime
    self.jobs_in = queue.Queue()
    self.jobs_out = queue.Queue()

    self.identifier_for = {}
    self.active = True
    self.daemon = threading.Thread( target = self.run )
    self.daemon.daemon = True
    self.daemon.start()


  def run(self):

    while self.active:
      self.manager.poll()

      while self.manager.completed_results:
        current = self.manager.completed_results[0]
        assert current.identifier in self.identifier_for

        try:
          self.jobs_out.put(
            ( self.identifier_for[ current.identifier ], current ),
            block = False,
            )

        except queue.Full:
          break

        self.manager.completed_results.popleft()
        del self.identifier_for[ current.identifier ]

      while True:
        try:
          ( identifier, ( target, args, kwargs ) ) = self.jobs_in.get( block = False )

        except queue.Empty:
          break

        i = self.manager.submit( target = target, args = args, kwargs = kwargs )
        assert i not in self.identifier_for
        self.identifier_for[ i ] = identifier

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


class SubmitJobs(communication.Command):
  """
  Request for job submission
  """

  def __init__(self, job_infos):

    self.job_infos = job_infos


  def process(self, environment):

    for ( identifier, callparms ) in self.job_infos:
      environment.jobs_in.put( ( identifier, callparms ) )


class GetJobs(communication.Command):

  def process(self, environment):

    jobs = []

    while True:
      try:
        jres = environment.jobs_out.get( block = False )

      except queue.Empty:
        break

      jobs.append( jres )

    return jobs


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

  def __init__(self, instream, outstream, waittime = 1, submittime = 0.5, poolsize = 4):

    super(SchedulerClient, self).__init__(
      instream = instream,
      outstream = outstream,
      )
    self.running = set()
    self.result_for = {}
    self.waiting = []
    self.waittime = waittime
    self.submittime = submittime
    self.polltime = time.time() # assume queue is empty
    self.enqueuetime = self.polltime
    self.poolsize = poolsize


  def is_known(self, identifier):

    return identifier in self.running or identifier in self.result_for


  def submit(self, identifier, target, args, kwargs):

    assert not self.is_known( identifier = identifier )

    if not self.waiting:
      self.enqueuetime = time.time()

    self.waiting.append( ( identifier, ( target, args, kwargs ) ) )
    self.running.add( identifier )
    self.poll()


  def is_alive(self, identifier):

    self.poll()
    return identifier not in self.result_for


  def get_result(self, identifier):

    return self.result_for[ identifier ]


  def remove(self, identifier):

    del self.result_for[ identifier ]


  def poll(self):

    now = time.time()

    if self.running and ( self.waittime < ( now - self.polltime ) or now < self.polltime ):
      response = self.send( command = GetJobs() )

      jobs = response()

      for ( identifier, value ) in jobs:
        assert identifier in self.running
        self.running.remove( identifier )

        assert identifier not in self.result_for
        self.result_for[ identifier ] = value

      now = time.time()
      self.polltime = now

    if (
      ( self.poolsize <= len( self.waiting ) )
      or ( self.waiting
        and ( self.submittime < ( now - self.enqueuetime ) or now < self.enqueuetime )
        )
      ):
      response = self.send( command = SubmitJobs( job_infos = self.waiting ) )

      try:
        response()

      except Exception as e:
        raise RuntimeError("Submission failure: %s" % e)

      self.waiting = []


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

    raise RuntimeError("job has not been submitted yet")


class RunningStatus(object):
  """
  A job that has been submitted
  """

  @staticmethod
  def start(job):

    raise RuntimeError("start called second time")


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

    raise RuntimeError("start called second time")


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


class RemoteFactory(object):
  """
  Remote instance method factory. There is no check that the instance has
  a method passed along with the constructor
  """

  def __init__(self, calculation, method):

    from libtbx.object_oriented_patterns import lazy_initialization
    self.instance = lazy_initialization( func = calculation )
    self.method = method


  def __call__(self, *args, **kwargs):

    return getattr( self.instance(), self.method )( *args, **kwargs )


def object_to_argument(obj):

  return pickle.dumps( obj, 0 )


def argument_to_object(arg):

  return pickle.loads( arg.decode( "string-escape" ) )


def command_merge(cmdline):

  return "%s %r %r %s" % cmdline


def command_unmerge(cmdline):

  return cmdline


def server_process_command_line(
  job_factory,
  queue_factory,
  executable = "libtbx.remote_processing",
  folder = ".",
  transformation = command_merge,
  ):

  return transformation(
    cmdline = (
      executable,
      object_to_argument( obj = job_factory ),
      object_to_argument( obj = queue_factory ),
      "--folder=%s" % folder,
      ),
    )
