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

  def __init__(self, queue):

    self.queue = queue


  def prepare(self, target):

    return RetrieveTarget( queue = self.queue, target = target )


  def finalize(self, identifier):

    return RetrieveResult( identifier = identifier, result = self.queue.get() )


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


class JobIterator(object):
  """
  Iterate over completed jobs
  """

  def __init__(self, manager):

    self.manager = manager


  def next(self):

    self.manager.wait()

    try:
      identifier = self.manager.completed_jobs.popleft()

    except IndexError:
      raise StopIteration

    return identifier


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
    self.completed_jobs = deque()
    self.polling_interval = polling_interval


  @property
  def running_jobs(self):

    return [ info[1] for info in self.execinfo_for.values() ]


  @property
  def finished_jobs(self):

    return JobIterator( manager = self )


  def submit(self, target, args = (), kwargs = {}):

    identifier = Identifier( target = target, args = args, kwargs = kwargs )
    self.waiting_jobs.append( identifier )
    self.poll()
    return identifier


  def empty(self):

    return not self.execinfo_for


  def full(self):

    return not self.available_units


  def wait(self):

    while not self.completed_jobs and not self.empty():
      self.poll()
      time.sleep( self.polling_interval )


  def join(self):

    while not self.empty():
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
        self.completed_jobs.append( result )

    # Submit new jobs
    while not self.full() and self.waiting_jobs:
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
    It is possible to have adapter.empty() == True and adapter.full() == True
    simultaneously, even if the number of cpus is not zero. This is because
    adapter.empty() report the status of the adapter queue, while adapter.full()
    reports that of the manager. This gives the behaviour one normally expects,
    i.e. no submission if there are no free cpus, and empty status when
    jobs submitted through the adaptor have all been processed.
  """

  def __init__(self, manager):

    self.manager = manager
    self.completed_jobs = deque()
    self.active_jobs = set()


  @property
  def waiting_jobs(self):

    return [ j for j in self.manager.waiting_jobs if j in self.active_jobs ]


  @property
  def running_jobs(self):

    return [ j for j in self.manager.running_jobs if j in self.active_jobs ]


  @property
  def finished_jobs(self):

    return JobIterator( manager = self )


  def submit(self, target, args = (), kwargs = {}):

    identifier = self.manager.submit( target = target, args = args, kwargs = kwargs )
    self.active_jobs.add( identifier )
    return identifier


  def empty(self):

    return ( not self.running_jobs and not self.waiting_jobs )


  def full(self):

    return self.manager.full()


  def wait(self):

    self.transfer_finished_jobs()

    while not self.completed_jobs and not self.empty():
      self.poll()
      time.sleep( self.manager.polling_interval )


  def join(self):

    while not self.empty():
      self.poll()
      time.sleep( self.manager.polling_interval )


  def poll(self):

    self.manager.poll()
    self.transfer_finished_jobs()


  def transfer_finished_jobs(self):

    for ( index, result ) in enumerate( list( self.manager.completed_jobs ) ):
      if result.identifier in self.active_jobs:
        self.manager.completed_jobs.remove( result )
        self.active_jobs.remove( result.identifier )
        self.completed_jobs.append( result )

