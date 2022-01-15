"""
Process Pool

Pools offer lower control over execution, as it is not possible to know the
state of individual workers

Common methods:
  has_results(): returns whether there are any results available
  results(): return an iterator yielding ( identifier, result ) tuples
  submit(target, args = (), kwargs = {}): submits job, return identifier
  is_empty(): returns whether there are no more jobs or results waiting
  is_full(): returns whether the number of currently processing jobs is at maximum
  shutdown(): sets the number of workers to zero
  resume(): starts workers as necessary
  join(): finish processing all jobs, and then shutdown
  terminate(): kills all processing

Pool methods:
  job_count(): return number of jobs submitted (waiting + running)
  worker_count(): return an approximate number of active workers
"""

from __future__ import absolute_import, division, print_function

import time
from collections import deque
from six.moves.queue import Empty, Full

from libtbx.scheduling import SchedulingError, identifier
from libtbx.scheduling import result

# Accounting
class process_register(object):
  """
  A simple class that encapsulates job accounting
  """

  def __init__(self):

    self.running_on = {}
    self.terminateds = deque()
    self.requested_shutdowns = deque()
    self.results = deque()


  def process_count(self):

    return len( self.running_on )


  def record_process_startup(self, pid):

    if pid in self.running_on:
      raise SchedulingError("Existing worker with identical processID")

    self.running_on[ pid ] = None


  def record_job_start(self, jobid, pid):

    if pid not in self.running_on:
      raise SchedulingError("Unknown processID")

    if self.running_on[ pid ] is not None:
      raise SchedulingError("Attempt to start process on busy worker")

    self.running_on[ pid ] = jobid


  def record_job_finish(self, jobid, pid, value):

    if pid not in self.running_on:
      raise SchedulingError("Unknown processID")

    if self.running_on[ pid ] != jobid:
      raise SchedulingError("Inconsistent register information: jobid/pid mismatch")

    self.running_on[ pid ] = None
    self.results.append( ( jobid, value ) )


  def record_process_shutdown(self, pid):

    self.record_process_exit( pid = pid, container = self.requested_shutdowns )


  def record_process_termination(self, pid):

    self.record_process_exit( pid = pid, container = self.terminateds )


  def record_process_exit(self, pid, container):

    if pid not in self.running_on:
      raise SchedulingError("Unknown processID")

    if self.running_on[ pid ] is not None:
      raise SchedulingError("Shutdown of busy worker")

    container.append( pid )
    del self.running_on[ pid ]


  def record_process_crash(self, pid, exception, traceback):

    if pid not in self.running_on:
      raise SchedulingError("Unknown processID")

    jobid = self.running_on[ pid ]

    if jobid is not None:
      self.results.append(
        ( jobid, result.error( exception = exception, traceback = traceback ) )
        )

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

  ( pid, exception, traceback ) = data
  register.record_process_crash( pid = pid, exception = exception, traceback = traceback )


# Lifecycle controllers
class unlimited(object):
  """
  A worker that is not supposed to die
  """

  def active(self):

    return True


  def record_job_start(self):

    pass


  def record_job_end(self):

    pass


class jobcount_limit(object):
  """
  A worker that completes N jobs, and then exits
  """

  def __init__(self, maxjobs):

    self.jobs_to_go = maxjobs


  def active(self):

    return 0 < self.jobs_to_go


  def record_job_start(self):

    pass


  def record_job_end(self):

    self.jobs_to_go -= 1


class runtime_limit(object):
  """
  A worker that exits when a time limit has been reached
  """

  def __init__(self, seconds):

    self.endtime = time.time() + seconds


  def active(self):

    return time.time() < self.endtime


  def record_job_start(self):

    pass


  def record_job_done(self):

    pass


class wallclock_limit(object):
  """
  A worker that exits when a wall clock time limit has been reached
  """

  def __init__(self, seconds):
    from libtbx.development.timers import work_clock
    self.maxtime = work_clock() + seconds


  def active(self):
    from libtbx.development.timers import work_clock
    return work_clock() < self.maxtime


  def record_job_start(self):

    pass


  def record_job_done(self):

    pass


def pool_process_cycle(
  pid,
  inqueue,
  outqueue,
  waittime,
  lifecycle,
  termination_signal,
  idle_timeout,
  ):

  controller = lifecycle()
  outqueue.put( ( worker_startup_event, pid ) )
  last_activity = time.time()

  while controller.active():
    if last_activity + idle_timeout < time.time():
      outqueue.put( ( worker_termination_event, pid ) )
      break

    try:
      data = inqueue.get( timeout = waittime )

    except Empty:
      continue

    if data == termination_signal:
      outqueue.put( ( worker_shutdown_event, pid ) )
      break

    assert len( data ) == 4
    ( jobid, target, args, kwargs ) = data
    outqueue.put( ( job_started_event, ( jobid, pid ) ) )
    controller.record_job_start()

    try:
      value = target( *args, **kwargs )

    except Exception as e:
      res = result.error( exception = e, traceback = result.get_traceback_info() )

    else:
      res = result.success( value = value )

    outqueue.put( ( job_finished_event, ( jobid, pid, res ) ) )
    controller.record_job_end()
    last_activity = time.time()

  else:
    outqueue.put( ( worker_termination_event, pid ) )


# Autoscaling
class constant_capacity(object):
  """
  Keep pool size constant
  """

  def __init__(self, capacity):

    self.capacity = capacity


  def __call__(self, manager):

    if manager.worker_count() != self.capacity:
      manager.set_worker_count( target = self.capacity )


class upscaling_capacity(object):
  """
  Increase pool size when needed, rely on timeout for downscaling
  """

  def __init__(self, minimum, maximum, buffering):

    self.minimum = minimum
    self.maximum = maximum
    self.buffering = buffering


  @property
  def capacity(self):

    return self.maximum


  def __call__(self, manager):

    ideal = max(
      min( self.maximum, manager.job_count() + self.buffering ),
      self.minimum,
      )

    worker_count = manager.worker_count()

    if worker_count < ideal:
      manager.set_worker_count( target = ideal )


class breathing_capacity(object):
  """
  Varies pool size according to load
  """

  def __init__(self, minimum, maximum, buffering, downsize_step = 2, downsize_delay = 5):

    self.minimum = minimum
    self.maximum = maximum
    self.buffering = buffering
    self.downsize_step = downsize_step
    self.downsize_delay = downsize_delay

    self.downsize_wait_start = None


  @property
  def capacity(self):

    return self.maximum


  def __call__(self, manager):

    ideal = max(
      min( self.maximum, manager.job_count() + self.buffering ),
      self.minimum,
      )

    worker_count = manager.worker_count()

    if worker_count < ideal:
      manager.set_worker_count( target = ideal )
      self.downsize_wait_start = None

    elif ideal < worker_count:
      if self.downsize_wait_start is None:
        self.downsize_wait_start = time.time()

      elif self.downsize_delay < time.time() - self.downsize_wait_start:
        target = max( worker_count - self.downsize_step, ideal )
        manager.set_worker_count( target = target )
        self.downsize_wait_start = time.time()

    else:
      self.downsize_wait_start = None


class manager(object):
  """
  Process pool
  """

  TERMINATION_SIGNAL = None

  def __init__(self,
    inqueue,
    outqueue,
    job_factory,
    autoscaling,
    lifecycle,
    waittime = 0.01,
    stalltime = 2,
    idle_timeout = 120,
    ):

    self.job_factory = job_factory
    self.autoscaling = autoscaling
    self.lifecycle = lifecycle

    self.inqueue = inqueue
    self.outqueue = outqueue

    self.waittime = waittime
    self.stalltime = stalltime
    self.idle_timeout = idle_timeout

    self.register = process_register()

    from itertools import count

    self.pid_assigner = count()
    self.process_numbered_as = {}
    self.recycleds = deque()
    self.terminatings = set()
    self.unreporteds = deque()
    self.outstanding_shutdown_requests = 0

    self.running_jobs = set()
    self.completed_results = deque()

    self.manage()


  def job_count(self):

    return len( self.running_jobs )


  def worker_count(self):

    return max(
      ( len( self.process_numbered_as ) + len( self.terminatings )
        - self.outstanding_shutdown_requests ),
      0,
      )


  def is_empty(self):

    return not self.running_jobs and not self.completed_results


  def is_full(self):

    return self.autoscaling.capacity <= self.job_count()


  def has_results(self):

    return self.completed_results


  def results(self):

    while not self.is_empty():
      while not self.has_results():
        self.wait()
        self.poll()
        self.manage()

      yield self.completed_results.popleft()


  def submit(self, target, args = (), kwargs = {}):

    jobid = identifier()
    self.outqueue.put( ( jobid, target, args, kwargs ) )
    self.running_jobs.add( jobid )
    return jobid


  def shutdown(self):

    self.set_worker_count( target = 0 )


  def resume(self):

    self.manage()


  def join(self):

    while self.running_jobs:
      self.poll()
      self.wait()

    self.poll()


  def terminate(self):

    self.shutdown()
    self.wait()
    self.poll()

    for process in self.process_numbered_as.values():
      if process.is_alive():
        if hasattr( process, "terminate" ): # Thread has no terminate
          try:
            process.terminate()

          except Exception:
            pass

    while self.process_numbered_as:
      self.poll()
      self.wait()


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
        "lifecycle": self.lifecycle,
        "termination_signal": self.TERMINATION_SIGNAL,
        "idle_timeout": self.idle_timeout,
        },
      )

    process.start()
    self.process_numbered_as[ pid ] = process


  def stop_process(self):

    self.outqueue.put( self.TERMINATION_SIGNAL )
    self.outstanding_shutdown_requests += 1


  def wait(self):

    time.sleep( self.waittime )


  def poll(self):

    for (pid, process) in list(self.process_numbered_as.items()):
      if not process.is_alive():
        process.join()

        exit_code = getattr( process, "exitcode", 0 ) # Thread has no "exitcode" attribute

        if exit_code != 0:
          data = (
            pid,
            result.get_exception( process = process, exit_code = exit_code ),
            result.get_crash_info( process = process ),
            )
          self.unreporteds.append( ( worker_crash_event, data ) )

        self.terminatings.add( pid )
        del self.process_numbered_as[ pid ]

    while self.unreporteds:
      try:
        self.inqueue.put( self.unreporteds[0], timeout = self.stalltime )

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
      pid = self.register.terminateds[0]

      try:
        self.terminatings.remove( pid )

      except KeyError:
        break

      self.register.terminateds.popleft()
      self.recycleds.append( pid )

    while self.register.requested_shutdowns:
      pid = self.register.requested_shutdowns[0]

      try:
        self.terminatings.remove( pid )

      except KeyError:
        break

      self.register.requested_shutdowns.popleft()
      self.recycleds.append( pid )
      assert 0 < self.outstanding_shutdown_requests
      self.outstanding_shutdown_requests -= 1

    while self.register.results:
      ( jobid, res ) = self.register.results.popleft()
      self.completed_results.append( ( jobid, res ) )
      self.running_jobs.remove( jobid )


  def set_worker_count(self, target):

    while self.worker_count() < target:
      self.start_process()

    while target < self.worker_count():
      self.stop_process()


  def manage(self):

    self.autoscaling( manager = self )


class creator(object):
  """
  Information to create and destroy a manager
  """

  def __init__(
    self,
    job_factory,
    inq_factory,
    outq_factory,
    autoscaling,
    lifecycle,
    waittime = 0.01,
    stalltime = 2,
    idle_timeout = 120,
    ):

    self.job_factory = job_factory
    self.inq_factory = inq_factory
    self.outq_factory = outq_factory
    self.autoscaling = autoscaling
    self.lifecycle = lifecycle
    self.waittime = waittime
    self.stalltime = stalltime
    self.idle_timeout = idle_timeout


  def create(self):

    return manager(
      inqueue = self.inq_factory.create(),
      outqueue = self.outq_factory.create(),
      job_factory = self.job_factory,
      autoscaling = self.autoscaling,
      lifecycle = self.lifecycle,
      waittime = self.waittime,
      stalltime = self.stalltime,
      idle_timeout = self.idle_timeout,
      )


  def destroy(self, manager):

    manager.terminate()
    manager.join()
    self.inq_factory.destroy( manager.inqueue )
    self.outq_factory.destroy( manager.outqueue )
