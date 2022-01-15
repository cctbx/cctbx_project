from __future__ import absolute_import, division, print_function

from six.moves.queue import Empty, Full
from libtbx.scheduling import result

import time
from collections import deque
from six.moves import range
from six.moves import zip

class Identifier(object):
  """
  Job identifier
  """

  def __init__(self, target, args, kwargs):

    self.target = target
    self.args = args
    self.kwargs = kwargs


  def __repr__(self):

    return "Identifier(id = %s)" % id( self )


class Result(object):
  """
  Empty calculation result
  """

  def __init__(self, identifier, value):

    self.identifier = identifier
    self.value = value


  def __repr__(self):

    return "Result(id = %s)" % id( self.identifier )


  def __call__(self):

    return self.value()


class ProcessingException(Exception):
  """
  An exception signalling temporary error with post-processing, e.g. a timeout
  """


class NullProcessor(object):
  """
  A singleton class to not do anything
  """

  @staticmethod
  def prepare(target):

    return target


  @staticmethod
  def finalize(identifier):

    return Result( identifier = identifier, value = result.success( value = None ) )


  @staticmethod
  def reset():

    pass


class RetrieveTarget(object):
  """
  A simple mangled target
  """

  def __init__(self, queue, target):

    self.queue = queue
    self.target = target


  def __call__(self, *args, **kwargs):

    result = self.target( *args, **kwargs )
    self.queue.put( result )


class RetrieveProcessor(object):
  """
  Adds a retrieve step
  """

  def __init__(self, queue, block = True, timeout = None, grace = None):

    self.queue = queue
    self.block = block
    self.timeout = timeout
    self.grace = max( timeout, grace )
    self.request_timeout = None


  def prepare(self, target):

    return RetrieveTarget( queue = self.queue, target = target )


  def finalize(self, job):

    try:
      job_result = self.queue.get( block = self.block, timeout = self.timeout )

    except Empty:
      now = time.time()

      if self.request_timeout is None:
        self.request_timeout = now
        raise ProcessingException("Timeout on %s" % self.queue)

      if ( now - self.request_timeout ) <= self.grace:
        raise ProcessingException("Timeout on %s" % self.queue)

      self.reset()
      raise RuntimeError("Timeout on %s" % self.queue)

    self.request_timeout = None

    return job_result


  def reset(self):

    self.request_timeout = None


class FetchProcessor(object):
  """
  Adds a fetch step. Only affects postprocessing. Can only be combined with
  job classes that save the return value as property
  """

  def __init__(self, transformation):

    self.transformation = transformation


  def prepare(self, target):

    return target



  def finalize(self, job):

    return self.transformation( job.result )


  @classmethod
  def Unpack(cls):

    return cls( transformation = lambda r: r() )


  @classmethod
  def Identity(cls):

    return cls( transformation = lambda r: r )


class ExecutionUnit(object):
  """
  Wrapper to make factory functions countable
  """

  def __init__(self, factory, priority = 0, processor = NullProcessor):

    self.factory = factory
    self.priority = priority
    self.processor = processor


  def start(self, target, args, kwargs):

    job = self.factory(
      target = self.processor.prepare( target = target ),
      args = args,
      kwargs = kwargs,
      )
    job.start()
    return job


  def finalize(self, job):

    return self.processor.finalize( job = job )


  def reset(self):

    self.processor.reset()


  def __str__(self):

    return "%s(id = %s)" % ( self.__class__.__name__, id( self ) )


class WorkerProcess(object):
  """
  A class that runs a worker process
  """

  def __call__(self, inqueue, outqueue):

    # Check-in
    outqueue.put( None )

    while True:
      data = inqueue.get()

      # Sentinel
      if not data:
        outqueue.put( None )
        break

      assert len( data ) == 3
      ( target, args, kwargs ) = data
      result = target( *args, **kwargs )
      outqueue.put( result )


class TerminatedException(Exception):
  """
  Exception class signalling that process has died
  """


class WorkerMaintenance(object):
  """
  Methods to create/destroy workers
  """

  def __init__(self, job_factory, queue_create, queue_teardown, grace = 10):

    self.job_factory = job_factory
    self.queue_create = queue_create
    self.queue_teardown = queue_teardown
    self.grace = grace


  def launch(self, executor, inqueue, outqueue):

    process = self.job_factory(
      target = executor,
      kwargs = { "inqueue": inqueue, "outqueue": outqueue },
      )
    process.start()
    return process


  def create(self, executor):

    q1 = self.queue_create()
    q2 = self.queue_create()
    process = self.launch( executor = executor, inqueue = q1, outqueue = q2 )
    return ( process, q2, q1 )


  def connect(self, process, inqueue):

    if process.is_alive():
      result = inqueue.get( block = False )
      assert result is None


  def functional(self, process):

    return process.is_alive()


  def cleanup(self, process, inqueue, outqueue):

    process.join()
    self.queue_teardown( inqueue )
    self.queue_teardown( outqueue )


  def shutdown(self, process, inqueue, outqueue):

    if process.is_alive():
      outqueue.put( None )

      try:
        inqueue.get( timeout = self.grace )

      except Empty:
        process.terminate()

    self.cleanup( process = process, inqueue = inqueue, outqueue = outqueue )


class Worker(object):
  """
  An existing process that runs jobs
  """

  def __init__(self, maintenance, priority = 0):

    self.maintenance = maintenance
    self.priority = priority
    self.startup()


  def startup(self):

    ( self.process, self.inqueue, self.outqueue ) = self.maintenance.create(
      executor = WorkerProcess()
      )
    self.connected = False


  def connect(self):

    assert not self.connected
    self.maintenance.connect( process = self.process, inqueue = self.inqueue )
    self.connected = True


  def shutdown(self):

    self.maintenance.shutdown(
      process = self.process,
      inqueue = self.inqueue,
      outqueue = self.outqueue,
      )


  def restart(self):

    self.shutdown()
    self.startup()


  def functional(self):

    return self.maintenance.functional( process = self.process )


  def start(self, target, args, kwargs):

    self.outqueue.put( ( target, args, kwargs ) )
    return WorkerJob( worker = self )


  def get(self, block):

    if not self.connected:
      self.connect()

    try:
      return self.inqueue.get( block = block )

    except Exception:
      if not self.functional():
        self.restart()
        raise TerminatedException("Worker died")

      else:
        raise


  def finalize(self, job):

    job.join()
    assert job.exitcode == 0

    return job.result


  def reset(self):

    pass


  def __str__(self):

    return "%s(id = %s)" % ( self.__class__.__name__, id( self ) )


class WorkerJob(object):
  """
  A job that is dispatched to a worker
  """

  def __init__(self, worker):

    self.worker = worker
    self.result = None
    self.exitcode = None


  def poll(self, block):

    if self.exitcode is None:
      try:
        self.result = self.worker.get( block = block )
        self.exitcode = 0

      except Empty:
        return True

      except Exception as e:
        self.result = result.error( exception = e )
        self.exitcode = 1
        self.err = e

    return False


  def is_alive(self):

    return self.poll( block = False )


  def join(self):

    self.poll( block = True )


class RunningState(object):
  """
  A job that is currently running
  """

  def __init__(self, job):

    self.job = job


  def is_alive(self):

    return self.job.is_alive()


  def initiate_postprocessing(self, job):

    self.job.join()

    exit_code = getattr( self.job, "exitcode", 0 ) # Thread has no "exitcode" attribute

    if exit_code == 0:
      job.status = PostprocessingState( job = self.job )

    else:
      err = getattr( self.job, "err", RuntimeError( "exit code = %s" % exit_code ) )
      job.status = ValueState( value = result.error( exception = err ) )


  def perform_postprocessing(self, job):

    self.initiate_postprocessing( job = job )
    job.perform_postprocessing()


  def get(self, job):

    self.initiate_postprocessing( job = job )
    return job.get()


  def __str__(self):

    return "running"


class PostprocessingState(object):
  """
  A job that is trying to postprocess
  """

  def __init__(self, job):

    self.job = job


  def is_alive(self):

    return False


  def initiate_postprocessing(self, job):

    raise RuntimeError("initiate_postprocess called second time")


  def perform_postprocessing(self, job):

    try:
      value = result.success( value = job.unit.finalize( job = self.job ) )

    except ProcessingException as e:
      raise

    except Exception as e:
      value = result.error( exception = e )

    job.status = ValueState( value = value )


  def get(self, job):

    self.perform_postprocessing( job = job )
    return job.get()


  def __str__(self):

    return "postprocessing"


class ValueState(object):
  """
  A job that either finished, exited with an error or failed postprocessing
  """

  def __init__(self, value):

    self.value = value


  def is_alive(self):

    return False


  def initiate_postprocessing(self, job):

    raise RuntimeError("initiate_postprocess called second time")


  def perform_postprocessing(self, job):

    pass


  def get(self, job):

    return self.value


  def __str__(self):

    return "finished"


class ExecutingJob(object):
  """
  A job that is executing on a Unit
  """

  def __init__(self, unit, identifier):

    self.unit = unit
    self.identifier = identifier
    self.status = RunningState(
      job = self.unit.start(
        target = identifier.target,
        args = identifier.args,
        kwargs = identifier.kwargs,
        )
      )


  def is_alive(self):

    return self.status.is_alive()


  def initiate_postprocessing(self):

    self.status.initiate_postprocessing( job = self )


  def perform_postprocessing(self):

    self.status.perform_postprocessing( job = self )


  def get(self):

    return Result(
      identifier = self.identifier,
      value = self.status.get( job = self ),
      )


  def __str__(self):

    return "%s(\n  unit = %s,\n  identifier = %s,\n  status = %s\n)" % (
      self.__class__.__name__,
      self.unit,
      self.identifier,
      self.status,
      )


# Manager implementation
"""
Process queue Manager interface

Manager properties:
  waiting_jobs
  running_jobs
  completed_jobs
  known_jobs
  completed_results

Manager iterators:
  results

Manager methods:
  submit(target, args = (), kwargs = {}):
  is_empty():
  is_full():
  wait():
  wait_for(identifier):
  result_for(identifier):
  join():
  poll():
"""


class Allocated(object):
  """
  The number of units is know upfront
  """

  def __init__(self, units):

    self.units = set( units )


  def empty(self):

    return not bool( self.units )


  def put(self, unit):

    self.units.add( unit )


  def get(self):

    unit = max( self.units, key = lambda s: s.priority )
    self.units.remove( unit )
    return unit


class Unlimited(object):
  """
  No limit on the number of execution units
  """

  def __init__(self, factory):

    self.factory = factory


  def empty(self):

    return False


  def put(self, unit):

    pass


  def get(self):

    return self.factory()


class UnlimitedCached(object):
  """
  No limit on the number of execution units. Created units are cached for reuse
  """

  def __init__(self, factory):

    self.factory = factory
    self.cache = set()


  def empty(self):

    return False


  def put(self, unit):

    self.cache.add( unit )


  def get(self):

    try:
      unit = self.cache.pop()

    except KeyError:
      unit = self.factory()

    return unit


class Scheduler(object):
  """
  Pure scheduler
  Job startup and pulldown is delegated to a handler object
  """

  def __init__(self, handler, polling_interval = 0.01):

    self.handler = handler
    self.executing = set()
    self.waiting_jobs = deque()
    self.postprocessing = deque()
    self.completed_results = deque()
    self.polling_interval = polling_interval


  @property
  def executing_jobs(self):

    return [ ej.identifier for ej in self.executing ]


  @property
  def postprocessing_jobs(self):

    return [ ej.identifier for ej in self.postprocessing ]


  @property
  def running_jobs(self):

    return self.executing_jobs + self.postprocessing_jobs


  @property
  def completed_jobs(self):

    return [ result.identifier for result in self.completed_results ]


  @property
  def known_jobs(self):

    return list( self.waiting_jobs ) + self.running_jobs + self.completed_jobs


  @property
  def results(self):

    while self.known_jobs:
      self.wait()
      yield self.completed_results.popleft()


  def submit(self, target, args = (), kwargs = {}):

    identifier = Identifier( target = target, args = args, kwargs = kwargs )
    self.waiting_jobs.append( identifier )
    self.poll()
    return identifier


  def is_empty(self):

    return not ( self.executing or self.postprocessing )


  def is_full(self):

    return self.handler.empty()


  def wait(self):

    while not self.completed_results and not self.is_empty():
      self.poll()
      time.sleep( self.polling_interval )


  def wait_for(self, identifier):

    if identifier not in self.known_jobs:
      raise RuntimeError("Job identifier not known")

    while identifier not in self.completed_jobs and not self.is_empty():
      self.poll()
      time.sleep( self.polling_interval )


  def result_for(self, identifier):

    self.wait_for( identifier = identifier )
    index = self.completed_jobs.index( identifier )
    result = self.completed_results[ index ]
    self.completed_results.remove( result )

    return result


  def join(self):

    while not self.is_empty():
      self.poll()
      time.sleep( self.polling_interval )


  def poll(self):

    # Process finished jobs
    for job in list( self.executing ):
      if not job.is_alive():
        self.executing.remove( job )
        job.initiate_postprocessing()
        self.postprocessing.append( job )

    # Postprocess jobs
    for i in range( len( self.postprocessing ) ):
      ej = self.postprocessing.popleft()

      try:
        ej.perform_postprocessing()

      except ProcessingException:
        self.postprocessing.append( ej )
        continue

      self.completed_results.append( ej.get() )
      self.handler.put( ej.unit )

    # Submit new jobs
    while not self.is_full() and self.waiting_jobs:
      identifier = self.waiting_jobs.popleft()
      unit = self.handler.get()
      job = ExecutingJob( unit = unit, identifier = identifier )
      self.executing.add( job )


# Backward compatibility
def Manager(units, polling_interval = 0.01):

  return Scheduler(
    handler = Allocated( units = units ),
    polling_interval = polling_interval,
    )


class Adapter(object):
  """
  Behaves like a manager, but uses an external resource to run the jobs

  Unintuitive feature:
    adapter.is_empty() == True and adapter.is_full() == True is possible
    simultaneously, even if the number of cpus is not zero. This is because
    adapter.is_empty() report the status of the adapter queue, while
    adapter.is_full() reports that of the manager. This gives the behaviour
    one normally expects, i.e. no submission if there are no free cpus, and
    empty status when jobs submitted through the adaptor have all been
    processed.
  """

  def __init__(self, manager):

    self.manager = manager
    self.completed_results = deque()
    self.active_jobs = set()


  @property
  def waiting_jobs(self):

    return [ j for j in self.manager.waiting_jobs if j in self.active_jobs ]


  @property
  def running_jobs(self):

    return [ j for j in self.manager.running_jobs if j in self.active_jobs ]


  @property
  def completed_jobs(self):

    return [ result.identifier for result in self.completed_results ]


  @property
  def known_jobs(self):

    return self.waiting_jobs + self.running_jobs + self.completed_jobs


  @property
  def results(self):

    while self.known_jobs:
      self.wait()
      yield self.completed_results.popleft()


  def submit(self, target, args = (), kwargs = {}):

    identifier = self.manager.submit( target = target, args = args, kwargs = kwargs )
    self.active_jobs.add( identifier )
    return identifier


  def is_empty(self):

    return ( not self.running_jobs and not self.waiting_jobs )


  def is_full(self):

    return self.manager.is_full()


  def wait(self):

    self.transfer_finished_jobs()

    while not self.completed_jobs and not self.is_empty():
      self.poll()
      time.sleep( self.manager.polling_interval )


  def wait_for(self, identifier):

    if identifier not in self.completed_jobs:
      if identifier not in self.known_jobs:
        raise RuntimeError("Job identifier not known")

      assert identifier in self.manager.known_jobs
      self.manager.wait_for( identifier = identifier )
      self.transfer_finished_jobs()

    assert identifier in self.completed_jobs


  def result_for(self, identifier):

    self.wait_for( identifier = identifier )
    index = self.completed_jobs.index( identifier )
    result = self.completed_results[ index ]
    self.completed_results.remove( result )

    return result


  def join(self):

    while not self.is_empty():
      self.poll()
      time.sleep( self.manager.polling_interval )


  def poll(self):

    self.manager.poll()
    self.transfer_finished_jobs()


  def transfer_finished_jobs(self):

    for result in list( self.manager.completed_results ):
      if result.identifier in self.active_jobs:
        self.manager.completed_results.remove( result )
        self.active_jobs.remove( result.identifier )
        self.completed_results.append( result )


class MainthreadJob(object):
  """
  Adaptor object that behaves like a job factory, but executes the calculation
  on the main thread. This can be used for debugging, or when only one CPU is
  available, as it has the lowest startup overhead.
  It is possible to use more than one, but it does not necessarily makes sense
  """

  def __init__(self, target, args, kwargs):

    self.target = target
    self.args = args
    self.kwargs = kwargs
    self.exception = None


  def start(self):

    try:
      self.target( *self.args, **self.kwargs )

    except Exception as e:
      self.exitcode = 1
      self.err = e


  def join(self):

    pass


  def is_alive(self):

    return False


class MainthreadQueue(object):
  """
  A queue to be used with MainthreadJob
  """

  def __init__(self):

    self.deque = deque()


  def put(self, obj, block = True, timeout = None):

    # NB: block and timeout are ignored, because it is not safe to use this
    #     with multiple threads
    self.deque.append( obj )


  def get(self, block = True, timeout = None):

    # NB: block and timeout are ignored, because it is not safe to use this
    #     with multiple threads
    return self.deque.popleft()


class MainthreadManager(object):
  """
  Dummy process queue executing on the main thread. Implements the Manager
  interface
  """

  def __init__(self, processor):

    self.waiting_jobs = deque()
    self.completed_results = deque()

    self.unit = ExecutionUnit( factory = MainthreadJob, processor = processor )


  @property
  def executing_jobs(self):

    return []


  @property
  def postprocessing_jobs(self):

    return []


  @property
  def running_jobs(self):

    return self.executing_jobs + self.postprocessing_jobs


  @property
  def completed_jobs(self):

    return [ result.identifier for result in self.completed_results ]


  @property
  def known_jobs(self):

    return list( self.waiting_jobs ) + self.running_jobs + self.completed_jobs


  @property
  def results(self):

    while self.known_jobs:
      self.wait()
      yield self.completed_results.popleft()


  def submit(self, target, args = (), kwargs = {}):

    identifier = Identifier( target = target, args = args, kwargs = kwargs )
    self.waiting_jobs.append( identifier )
    return identifier


  def is_empty(self):

    return not self.is_full()


  def is_full(self):

    return bool( self.waiting_jobs )


  def wait(self):

    if not self.completed_results:
      self.poll()


  def wait_for(self, identifier):

    if identifier not in self.known_jobs:
      raise RuntimeError("Job identifier not known")

    while identifier not in self.completed_jobs and not self.is_empty():
      self.poll()


  def result_for(self, identifier):

    self.wait_for( identifier = identifier )
    index = self.completed_jobs.index( identifier )
    result = self.completed_results[ index ]
    self.completed_results.remove( result )

    return result


  def join(self):

    while not self.is_empty():
      self.poll()


  def poll(self):

    if not self.is_empty():
      identifier = self.waiting_jobs.popleft()
      job = ExecutingJob( unit = self.unit, identifier = identifier )
      assert not job.is_alive()
      job.initiate_postprocessing()

      while True:
        try:
          job.perform_postprocessing()

        except ProcessingException:
          time.sleep( 0.1 )
          continue

        break

      self.completed_results.append( job.get() )


# Pool implementation
"""
Process queue Pool interface

Pools offer lower control over execution, as it is not possible to know the
state of individual workers

Pool properties:
  known_jobs
  completed_jobs
  completed_results

  job_count

Pool iterators:
  results

Pool methods:
  submit(target, args = (), kwargs = {}):
  is_empty():
  is_full():
  wait():
  join():
  poll():
  shutdown():
"""

class Immortal(object):
  """
  A worker that is not supposed to die
  """

  def is_alive(self):

    return True


  def job_start(self):

    pass


  def job_done(self):

    pass


class JobCount(object):
  """
  A worker that completes N jobs, and then exits
  """

  def __init__(self, maxjobs):

    self.jobs_to_go = maxjobs


  def is_alive(self):

    return 0 < self.jobs_to_go


  def job_start(self):

    pass


  def job_done(self):

    self.jobs_to_go -= 1


class MaxTime(object):
  """
  A worker that exits when a time limit has been reached
  """

  def __init__(self, seconds):

    self.endtime = time.time() + seconds


  def is_alive(self):

    return time.time() < self.endtime


  def job_start(self):

    pass


  def job_done(self):

    pass


class MaxWallClock(object):
  """
  A worker that exits when a wall clock time limit has been reached
  """

  def __init__(self, seconds):
    from libtbx.development.timers import work_clock
    self.maxtime = work_clock() + seconds


  def is_alive(self):

    from libtbx.development.timers import work_clock
    return work_clock() < self.maxtime


  def job_start(self):

    pass


  def job_done(self):

    pass


def pool_process_cycle(pid, inqueue, outqueue, waittime, lifeman, sentinel):

  manager = lifeman()
  outqueue.put( ( worker_startup_event, pid ) )

  while manager.is_alive():
    try:
      data = inqueue.get( timeout = waittime )

    except Empty:
      continue

    if data == sentinel:
      outqueue.put( ( worker_shutdown_event, pid ) )
      break

    assert len( data ) == 4
    ( jobid, target, args, kwargs ) = data
    outqueue.put( ( job_started_event, ( jobid, pid ) ) )
    manager.job_start()
    result = target( *args, **kwargs )
    outqueue.put( ( job_finished_event, ( jobid, pid, result ) ) )
    manager.job_done()

  else:
    outqueue.put( ( worker_termination_event, pid ) )


class ProcessRegister(object):
  """
  A simple class that encapsulates job management
  """

  def __init__(self):

    self.running_on = {}
    self.terminateds = deque()
    self.requested_shutdowns = deque()
    self.results = deque()


  @property
  def process_count(self):

    return len( self.running_on )


  def record_process_startup(self, pid):

    if pid in self.running_on:
      raise RuntimeError("Existing worker with identical processID")

    self.running_on[ pid ] = None


  def record_job_start(self, jobid, pid):

    if pid not in self.running_on:
      raise RuntimeError("Unknown processID")

    if self.running_on[ pid ] is not None:
      raise RuntimeError("Attempt to start process on busy worker")

    self.running_on[ pid ] = jobid


  def record_job_finish(self, jobid, pid, value):

    if pid not in self.running_on:
      raise RuntimeError("Unknown processID")

    if self.running_on[ pid ] != jobid:
      raise RuntimeError("Inconsistent register information: jobid/pid mismatch")

    self.running_on[ pid ] = None
    self.results.append( ( jobid, result.success( value = value ) ) )


  def record_process_shutdown(self, pid):

    self.record_process_exit( pid = pid, container = self.requested_shutdowns )


  def record_process_termination(self, pid):

    self.record_process_exit( pid = pid, container = self.terminateds )


  def record_process_exit(self, pid, container):

    if pid not in self.running_on:
      raise RuntimeError("Unknown processID")

    if self.running_on[ pid ] is not None:
      raise RuntimeError("Shutdown of busy worker")

    container.append( pid )
    del self.running_on[ pid ]


  def record_process_crash(self, pid, exception):

    if pid not in self.running_on:
      raise RuntimeError("Unknown processID")

    jobid = self.running_on[ pid ]

    if jobid is not None:
      self.results.append( ( jobid, result.error( exception = exception ) ) )

    del self.running_on[ pid ]
    self.terminateds.append( pid )


# Stand-alone functions for pickling
def job_started_event(register, data):

  ( jobid, pid ) = data
  register.record_job_start( jobid = jobid, pid = pid )


def job_finished_event(register, data):

  ( jobid, pid, value ) = data
  register.record_job_finish( jobid = jobid, pid = pid, value = value )


def worker_startup_event(register, data):

  register.record_process_startup( pid = data )


def worker_shutdown_event(register, data):

  register.record_process_shutdown( pid = data )


def worker_termination_event(register, data):

  register.record_process_termination( pid = data )


def worker_crash_event(register, data):

  ( pid, exception ) = data
  register.record_process_crash( pid = pid, exception = exception )


class ConstantStrategy(object):
  """
  Keep pool size constant
  """

  def __init__(self, capacity):

    self.capacity = capacity


  def target(self, job_count, process_count):

    return self.capacity


def normal_state(strategy, target, process_count):

  if target < process_count:
    strategy.state = downsize_state()
    target = process_count

  return target


class downsize_state(object):

  def __init__(self):

    self.starttime = time.time()


  def __call__(self, strategy, target, process_count):

    if process_count <= target or strategy.grace < time.time() - self.starttime:
      strategy.state = normal_state

    else:
      target = process_count

    return target


class BreathingStrategy(object):
  """
  Varies pool size according to load
  """

  def __init__(self, minimum, maximum, buffering, grace = 1):

    self.minimum = minimum
    self.maximum = maximum
    self.buffering = buffering
    self.grace = grace
    self.state = normal_state


  @property
  def capacity(self):

    return self.maximum


  def target(self, job_count, process_count):

    target = max(
      min( self.maximum, job_count + self.buffering ),
      self.minimum,
      )

    return self.state(
      strategy = self,
      target = target,
      process_count = process_count,
      )


class ProcessPool(object):
  """
  Process pool
  """

  SENTINEL = None

  def __init__(self,
    inqueue,
    outqueue,
    job_factory,
    loadmanager,
    lifemanager,
    waittime = 0.01,
    grace = 2,
    ):

    self.job_factory = job_factory
    self.loadmanager = loadmanager
    self.lifemanager = lifemanager

    self.inqueue = inqueue
    self.outqueue = outqueue

    self.waittime = waittime
    self.grace = grace

    self.register = ProcessRegister()

    from itertools import count

    self.pid_assigner = count()
    self.process_numbered_as = {}
    self.recycleds = deque()
    self.terminatings = set()
    self.unreporteds = deque()
    self.outstanding_shutdown_requests = 0

    self.identifier_for = {}
    self.completed_results = deque()

    self.manage()


  @property
  def job_count(self):

    return len( self.identifier_for )


  @property
  def completed_jobs(self):

    return [ result.identifier for result in self.completed_results ]


  @property
  def known_jobs(self):

    return list(self.identifier_for.values()) + self.completed_jobs


  @property
  def process_count(self):

    return max(
      ( len( self.process_numbered_as ) + len( self.terminatings )
        - self.outstanding_shutdown_requests ),
      0,
      )


  @property
  def results(self):

    while self.known_jobs:
      self.wait()
      yield self.completed_results.popleft()


  def submit(self, target, args = (), kwargs = {}):

    identifier = Identifier( target = target, args = args, kwargs = kwargs )
    jobid = id( identifier )
    self.outqueue.put( ( jobid, target, args, kwargs ) )
    self.identifier_for[ jobid ] = identifier
    return identifier


  def is_empty(self):

    return not self.identifier_for


  def is_full(self):

    return self.loadmanager.capacity <= self.job_count


  def wait(self):

    while not self.completed_results and not self.is_empty():
      self.poll()
      time.sleep( self.waittime )


  def poll(self):

    for ( pid, process ) in self.process_numbered_as.items():
      if not process.is_alive():
        process.join()

        exit_code = getattr( process, "exitcode", 0 ) # Thread has no "exitcode" attribute

        if exit_code != 0:
          err = getattr(
            process,
            "err",
            RuntimeError( "exit code = %s" % exit_code ),
            )
          self.unreporteds.append( ( worker_crash_event, ( pid, err ) ) )

        self.terminatings.add( pid )
        del self.process_numbered_as[ pid ]

    while self.unreporteds:
      try:
        self.inqueue.put( self.unreporteds[0], timeout = self.grace )

      except Full:
        break

      self.unreporteds.popleft()

    while True:
      try:
        ( event, data ) = self.inqueue.get( timeout = self.waittime )

      except Empty:
        break

      event( register = self.register, data = data )

    while self.register.terminateds:
      pid = self.register.terminateds.popleft()

      try:
        self.terminatings.remove( pid )

      except KeyError:
        self.register.terminateds.appendleft( pid )
        break

      self.recycleds.append( pid )

    while self.register.requested_shutdowns:
      pid = self.register.requested_shutdowns.popleft()

      try:
        self.terminatings.remove( pid )

      except KeyError:
        self.register.requested_shutdowns.appendleft( pid )
        break

      self.recycleds.append( pid )
      assert 0 < self.outstanding_shutdown_requests
      self.outstanding_shutdown_requests -= 1

    while self.register.results:
      ( jobid, value ) = self.register.results.popleft()
      self.completed_results.append(
        Result( identifier = self.identifier_for[ jobid ], value = value ),
        )
      del self.identifier_for[ jobid ]

    self.manage()


  def join(self):

    self.shutdown()

    while self.process_numbered_as:
      self.poll()
      time.sleep( self.waittime )

    self.poll()


  def shutdown(self):

    self.loadmanager = ConstantStrategy( capacity = 0 )
    self.manage()


  def terminate(self):

    self.shutdown()
    time.sleep( self.waittime )
    self.poll()

    for ( pid,  process ) in self.process_numbered_as.items():
      if process.is_alive():
        if hasattr( process, "terminate" ): # Thread has no terminate
          try:
            process.terminate()
            process.join()

          except Exception:
            pass

        self.terminatings.add( pid )
        del self.process_numbered_as[ pid ]

    self.poll()


  # Internal methods
  def start_process(self):

    try:
      pid = self.recycleds.popleft()

    except IndexError:
      pid = next(self.pid_assigner)

    process = self.job_factory(
      target = pool_process_cycle,
      kwargs = {
        "pid": pid,
        "inqueue": self.outqueue,
        "outqueue": self.inqueue,
        "waittime": self.waittime,
        "lifeman": self.lifemanager,
        "sentinel": self.SENTINEL,
        },
      )

    process.start()
    self.process_numbered_as[ pid ] = process


  def stop_process(self):

    self.outqueue.put( self.SENTINEL )
    self.outstanding_shutdown_requests += 1


  def manage(self):

    count = self.process_count
    target = self.loadmanager.target(
      job_count = self.job_count,
      process_count = count,
      )

    if count < target:
      while self.process_count < target:
        self.start_process()

    elif target < count:
      while target < self.process_count:
        self.stop_process()


class MainthreadPool(object):
  """
  Single process pool that executes on the main thread

  This is used for all single-CPU jobs, as it is more efficient than spawning a
  single job
  """

  def __init__(self):

    self.waiting_jobs = deque()
    self.completed_results = deque()


  @property
  def job_count(self):

    return len( self.waiting_jobs )


  @property
  def completed_jobs(self):

    return [ result.identifier for result in self.completed_results ]


  @property
  def known_jobs(self):

    return list( self.waiting_jobs ) + self.completed_jobs


  @property
  def process_count(self):

    return 1


  @property
  def results(self):

    while self.known_jobs:
      self.wait()
      yield self.completed_results.popleft()


  def submit(self, target, args = (), kwargs = {}):

    identifier = Identifier( target = target, args = args, kwargs = kwargs )
    self.waiting_jobs.append( identifier )
    return identifier


  def is_empty(self):

    return not self.waiting_jobs


  def is_full(self):

    return 1 <= self.job_count


  def wait(self):

    while not self.completed_results and not self.is_empty():
      self.poll()


  def poll(self):

    if self.waiting_jobs:
      current = self.waiting_jobs.popleft()

      try:
        value = current.target( *current.args, **current.kwargs )

      except Exception as e:
        res = result.error( exception = e )

      else:
        res = result.success( value = value )

      self.completed_results.append( Result( identifier = current, value = res ) )


  def join(self):

    while self.waiting_jobs:
      self.poll()


  def shutdown(self):

    pass


  def terminate(self):

    pass


# Parallel execution iterator
class RechargeableIterator(object):
  """
  Allows new elements to be filled in
  """

  def __init__(self):

    self.intermittent = deque()


  def next(self):

    try:
      return self.intermittent.popleft()

    except IndexError:
      raise StopIteration


  def has_next(self):

    return bool( self.intermittent )


  def append(self, value):

    self.intermittent.append( value )


  def extend(self, iterable):

    self.intermittent.extend( iterable )


class NoPooler(object):
  """
  Does not group jobs
  """

  @classmethod
  def pack(cls, calcsiter, manager):

    ( target, args, kwargs ) = next(calcsiter) # raise StopIteration
    identifier = manager.submit(
      target = target,
      args = args,
      kwargs = kwargs,
      )
    return [ identifier ]


  @classmethod
  def unpack(self, result, resiter):

    resiter.append( result )


class PooledRun(object):
  """
  Runs multiple jobs
  """

  def __init__(self, identifiers):

    self.identifiers = identifiers


  def __call__(self):

    results = []

    for iden in self.identifiers:
      try:
        res = iden.target( *iden.args, **iden.kwargs )

      except Exception as e:
        results.append( ( True, e ) )

      else:
        results.append( ( False, res ) )

    return results


  def results(self, raw):

    results = []

    try:
      res = raw()

    except Exception as e:
      results.extend(
        [
          Result(
            identifier = i,
            value = result.error( exception = e ),
            )
          for i in self.identifiers
          ]
        )

    else:
      assert len( res ) == len( self.identifiers )

      for ( identifier, ( failure, data ) ) in zip( self.identifiers, res ):
        if failure:
          results.append(
            Result( identifier = identifier, value = result.error( exception = data ) )
            )

        else:
          results.append(
            Result( identifier = identifier, value = result.success( value = data ) )
            )

    return results


class Pooler(object):
  """
  Pools up a number of jobs and runs this way

  pool - number of jobs to pool
  """

  def __init__(self, size):

    assert 0 < size
    self.size = size


  def pack(self, calcsiter, manager):

    identifiers = [
      Identifier( target = target, args = args, kwargs = kwargs )
      for ( index, ( target, args, kwargs ) )
      in zip( range(self.size), calcsiter )
      ]

    if not identifiers:
      raise StopIteration

    manager.submit(
      target = PooledRun( identifiers = identifiers ),
      args = (),
      kwargs = {},
      )

    return identifiers


  def unpack(self, result, resiter):

    resiter.extend( result.identifier.target.results( raw = result ) )


def get_pooler(size):

  if size == 1:
    return NoPooler

  elif 1 < size:
    return Pooler( size = size )

  else:
    raise ValueError("Invalid pool size: %s" % size)



class ParallelForIterator(object):
  """
  Creates an iterator that executes calls on a Manager-like object

  calculations - an iterable of calculations yielding ( target, args, kwargs ) tuples
  manager - execution manager
  """

  def __init__(self, calculations, manager, pool = 1):

    self.manager = manager
    self.pooler = get_pooler( size = pool )
    self.resiter = RechargeableIterator()
    self.calcsiter = iter( calculations )
    self.resume()


  def _ongoing(self, orderer):

    while not self.manager.is_full():
      try:
        identifiers = self.pooler.pack(
          calcsiter = self.calcsiter,
          manager = self.manager,
          )

      except StopIteration:
        self.suspend()
        break

      for identifier in identifiers:
        orderer.job_submitted( identifier = identifier )


  def _terminating(self, orderer):

    pass


  def next(self, orderer):

    if not self.resiter.has_next():
      result = next(self.manager.results) # raise StopIteration
      self.pooler.unpack( result = result, resiter = self.resiter )

    self.process( orderer = orderer )
    assert self.resiter.has_next()
    r = next(self.resiter)
    return ( r.identifier, r )


  def suspend(self):

    self.process = self._terminating


  def resume(self):

    self.process = self._ongoing


class Ordering(object):
  """
  Base class for ordering classes
  """

  def __init__(self, parallel_for):

    self.parallel_for = parallel_for
    self.parallel_for.process( orderer = self )


  def __del__(self):

    self.parallel_for.suspend()

    # Fetch all processing values
    for elem in self:
      pass

    self.parallel_for.resume()


  def __iter__(self):

    return self


class FinishingOrder(Ordering):
  """
  Results returned as jobs finish
  """

  def __init__(self, parallel_for):

    super( FinishingOrder, self ).__init__( parallel_for = parallel_for )


  def job_submitted(self, identifier):

    pass


  def next(self):

    ( identifier, result ) = self.parallel_for.next( orderer = self )
    return ( ( identifier.target, identifier.args, identifier.kwargs ), result )


class SubmissionOrder(Ordering):
  """
  Results returned as they are submitted
  """

  def __init__(self, parallel_for):

    self.submitteds = deque()
    self.result_for = {}
    super( SubmissionOrder, self ).__init__( parallel_for = parallel_for )


  def job_submitted(self, identifier):

    self.submitteds.append( identifier )


  def next(self):

    while not self.submitteds or self.submitteds[0] not in self.result_for:
      ( identifier, result ) = self.parallel_for.next( orderer = self )
      self.result_for[ identifier ] = result

    first = self.submitteds.popleft()
    result = self.result_for[ first ]
    del self.result_for[ first ]
    return ( ( first.target, first.args, first.kwargs ), result )
