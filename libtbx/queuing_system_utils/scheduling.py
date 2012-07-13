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

    return RetrieveResult(
      identifier = identifier,
      result = self.queue.get( block = self.block, timeout = self.timeout ),
      )


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
  """

  def __init__(self, units, polling_interval = 0.01):

    self.available_units = set( units )
    self.execinfo_for = {}
    self.waiting_jobs = deque()
    self.completed_results = deque()
    self.polling_interval = polling_interval


  @property
  def running_jobs(self):

    return [ info[1] for info in self.execinfo_for.values() ]


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
    self.poll()
    return identifier


  def is_empty(self):

    return not self.execinfo_for


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
    for job in self.execinfo_for.keys():
      if not job.is_alive():
        job.join()
        ( unit, identifier ) = self.execinfo_for[ job ]
        result = unit.finalize( identifier = identifier )
        del self.execinfo_for[ job ]
        self.available_units.add( unit )
        self.completed_results.append( result )

    # Submit new jobs
    while not self.is_full() and self.waiting_jobs:
      identifier = self.waiting_jobs.popleft()
      unit = max( self.available_units, key = lambda s: s.priority )
      self.available_units.remove( unit )
      job = unit.start(
        target = identifier.target,
        args = identifier.args,
        kwargs = identifier.kwargs
        )
      self.execinfo_for[ job ] = ( unit, identifier )


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

