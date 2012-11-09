from __future__ import division
import time
from collections import deque

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

  def __init__(self, identifier):

    self.identifier = identifier


  def __repr__(self):

    return "Result(id = %s)" % id( self.identifier )


  def __call__(self):

    return None


class ErrorEnding(object):
  """
  Error raised by the calculation result
  """

  def __init__(self, identifier, exception):

    self.identifier = identifier
    self.exception = exception


  def __repr__(self):

    return "ErrorEnding(id = %s)" % id( self.identifier )


  def __call__(self):

    raise self.exception


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

    return Result( identifier = identifier )


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


class RetrieveResult(object):
  """
  Result bundled with job
  """

  def __init__(self, identifier, result):

    self.identifier = identifier
    self.result = result


  def __repr__(self):

    return "RetrieveResult(id = %s)" % id( self.identifier )


  def __call__(self):

    return self.result


class RetrieveProcessor(object):
  """
  Adds a retrieve step
  """

  def __init__(self, queue, block = True, timeout = None):

    self.queue = queue
    self.block = block
    self.timeout = timeout


  def prepare(self, target):

    return RetrieveTarget( queue = self.queue, target = target )


  def finalize(self, identifier):

    from Queue import Empty

    try:
      result = self.queue.get( block = self.block, timeout = self.timeout )

    except Empty:
      raise ProcessingException, "Timeout on %s" % self.queue

    return RetrieveResult( identifier = identifier, result = result )


  def reset(self):

    pass


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


  def finalize(self, identifier):

    return self.processor.finalize( identifier = identifier )


  def error(self, identifier, exception):

    self.processor.reset()
    return ErrorEnding( identifier = identifier, exception = exception )


  def __str__(self):

    return "%s(id = %s)" % ( self.__class__.__name__, id( self ) )


class ExecutingJob(object):
  """
  A job that is executing on a Unit
  """

  def __init__(self, unit, identifier):

    self.unit = unit
    self.identifier = identifier
    self.job = self.unit.start(
      target = identifier.target,
      args = identifier.args,
      kwargs = identifier.kwargs,
      )


  def is_alive(self):

    return self.job.is_alive()


  def postprocess(self):

    try:
      self.job.join()

    except Exception, e:
      self.method = self.unit.error
      self.args = ( self.identifier, e )

    else:
      self.method = self.unit.finalize
      self.args = ( self.identifier, )


  def get(self):

    assert hasattr( self, "method" ) and hasattr( self, "args" )
    return self.method( *self.args )


  def __str__(self):

    return "%s(\n  unit = %s,\n  indentifier = %s,\n  job = %s\n)" % (
      self.__class__.__name__,
      self.unit,
      self.identifier,
      self.job,
      )


class ResultIterator(object):
  """
  Iterate over completed results
  """

  def __init__(self, manager):

    self.manager = manager


  def next(self):

    self.manager.wait()

    try:
      result = self.manager.completed_results.popleft()

    except IndexError:
      raise StopIteration

    return result


  def __iter__(self):

    return self


class OrderedResultIterator(object):
  """
  Iterate over completed results in specified order
  """

  def __init__(self, manager, identifiers):

    self.manager = manager
    self.identifiers = iter( identifiers )


  def next(self):

    identifier = self.identifiers.next() # raise StopIteration
    return self.manager.result_for( identifier = identifier )


  def __iter__(self):

    return self


class Manager(object):
  """
  Process queue

  Manager properties:
    waiting_jobs
    running_jobs
    completed_jobs
    known_jobs

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

  def __init__(self, units, polling_interval = 0.01):

    self.available_units = set( units )
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


  def results_in_order_of(self, identifiers):

    return OrderedResultIterator( manager = self, identifiers = identifiers )


  def submit(self, target, args = (), kwargs = {}):

    identifier = Identifier( target = target, args = args, kwargs = kwargs )
    self.waiting_jobs.append( identifier )
    self.poll()
    return identifier


  def is_empty(self):

    return not ( self.executing or self.postprocessing )


  def is_full(self):

    return not self.available_units


  def wait(self):

    while not self.completed_results and not self.is_empty():
      self.poll()
      time.sleep( self.polling_interval )


  def wait_for(self, identifier):

    if identifier not in self.known_jobs:
      raise RuntimeError, "Job identifier not known"

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
        job.postprocess()
        self.postprocessing.append( job )

    # Postprocess jobs
    for i in range( len( self.postprocessing ) ):
      ej = self.postprocessing.popleft()

      try:
        result = ej.get()

      except ProcessingException, e:
        self.postprocessing.append( ej )
        continue

      self.available_units.add( ej.unit )
      self.completed_results.append( result )

    # Submit new jobs
    while not self.is_full() and self.waiting_jobs:
      identifier = self.waiting_jobs.popleft()
      unit = max( self.available_units, key = lambda s: s.priority )
      self.available_units.remove( unit )
      job = ExecutingJob( unit = unit, identifier = identifier )
      self.executing.add( job )


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

    return ResultIterator( manager = self )


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
        raise RuntimeError, "Job identifier not known"

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

    except Exception, e:
      self.exception = e


  def join(self):

    if self.exception:
      raise self.exception


  def is_alive(self):

    return False


class MainthreadQueue(object):
  """
  A queue to be used with MainthreadJob
  """

  def __init__(self):

    self.deque = deque()


  def put(self, obj):

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

    return ResultIterator( manager = self )


  def results_in_order_of(self, identifiers):

    return OrderedResultIterator( manager = self, identifiers = identifiers )


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
      raise RuntimeError, "Job identifier not known"

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
      job.postprocess()
      result = job.get()
      self.completed_results.append( result )


# Parallel execution iterator
class ParallelForIterator(object):
  """
  Creates an iterator that executes a function with variable input argument on a
  Manager-like object

  function - the function to execute
  calculations - an iterable of calculations yielding ( target, args, kwargs ) tuples
  manager - execution manager
  """

  def __init__(self, calculations, manager):

    self.calcsiter = iter( calculations )
    self.manager = manager

    self.process = self._ongoing


  def _ongoing(self, orderer):

    while not self.manager.is_full():
      try:
        ( target, args, kwargs ) = self.calcsiter.next()

      except StopIteration:
        self.process = self._terminating
        break

      identifier = self.manager.submit(
          target = target,
          args = args,
          kwargs = kwargs,
          )
      orderer.job_submitted( identifier = identifier )


  def _terminating(self, orderer):

    pass


  def suspend(self):

    self.process = self._terminating


  def restart(self):

    self.process = self._ongoing


  def next(self, orderer):

    result = self.manager.results.next() # raise StopIteration
    self.process( orderer = orderer )
    return ( result.identifier, result )


class FinishingOrder(object):
  """
  Results returned as jobs finish
  """

  def __init__(self, parallel_for):

    self.parallel_for = parallel_for
    self.parallel_for.process( orderer = self )


  def job_submitted(self, identifier):

    pass


  def next(self):

    ( identifier, result ) = self.parallel_for.next( orderer = self )
    return ( ( identifier.target, identifier.args, identifier.kwargs ), result )


  def __iter__(self):

    return self


class SubmissionOrder(object):
  """
  Results returned as they are submitted
  """

  def __init__(self, parallel_for):

    self.parallel_for = parallel_for
    self.submitteds = deque()
    self.result_for = {}
    self.parallel_for.process( orderer = self )


  def job_submitted(self, identifier):

    self.submitteds.append( identifier )


  def next(self):

    if not self.submitteds:
      raise StopIteration

    first = self.submitteds.popleft()

    while first not in self.result_for:
      ( identifier, result ) = self.parallel_for.next( orderer = self )
      self.result_for[ identifier ] = result

    result = self.result_for[ first ]
    del self.result_for[ first ]
    return ( ( first.target, first.args, first.kwargs ), result )


  def __iter__(self):

    return self

