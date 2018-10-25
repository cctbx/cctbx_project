"""
Virtual Manager

Behaves as a normal manager instance, but performs processing through an internally
held manager, and only handles jobs submitted through it.

As this is a virtual resource, there is no creator associated with it.

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
"""

from __future__ import absolute_import, division, print_function

from collections import deque


class server(object):
  """
  Holds the underlying manager and sorts finished jobs according to submitter
  """

  def __init__(self, resource):

    self.resource = resource
    self.results_submitted_through = {}
    self.unfinisheds_submitted_through = {}
    self.submitter_for = {}


  def manager(self):

    identifier = object()
    assert identifier not in self.results_submitted_through
    assert identifier not in self.unfinisheds_submitted_through
    self.results_submitted_through[ identifier ] = deque()
    self.unfinisheds_submitted_through[ identifier ] = 0
    return manager( server = self, identifier = identifier )


  # Internal methods
  def submit(self, identifier, target, args = (), kwargs = {}):

    jobid = self.resource.submit( target = target, args = args, kwargs = kwargs )
    assert jobid not in self.submitter_for
    self.submitter_for[ jobid ] = identifier
    self.unfinisheds_submitted_through[ identifier ] += 1

    return jobid


  def results_for(self, identifier):

    return self.results_submitted_through[ identifier ]


  def unfinisheds_for(self, identifier):

    return self.unfinisheds_submitted_through[ identifier ]


  def update(self):

    while self.resource.has_results():
      self.fetch()


  def fetch(self):

    ( jobid, result ) = next(self.resource.results())
    identifier = self.submitter_for[ jobid ]
    del self.submitter_for[ jobid ]
    self.results_submitted_through[ identifier ].append( ( jobid, result ) )
    self.unfinisheds_submitted_through[ identifier ] -= 1


  # Calls forwarded to underlying manager
  def is_full(self):

    return self.resource.is_full()


  def shutdown(self):

    self.resource.shutdown()


  def resume(self):

    self.resource.resume()


  def join(self):

    self.resource.join()


  def terminate(self):

    self.resource.terminate()


class manager(object):
  """
  Behaves like a manager, but uses an external resource to run the jobs

  Unintuitive feature:
    manager.is_empty() == True and manager.is_full() == True is possible
    simultaneously, even if the number of cpus is not zero. This is because
    adapter.is_empty() report the status of the adapter queue, while
    adapter.is_full() reports that of the manager. This gives the behaviour
    one normally expects, i.e. no submission if there are no free cpus, and
    empty status when jobs submitted through the adaptor have all been
    processed.
  """

  def __init__(self, server, identifier):

    self.server = server
    self.identifier = identifier


  @property
  def completed_results(self):

    return self.server.results_for( identifier = self.identifier )


  def has_results(self):

    self.server.update()
    return self.completed_results


  def results(self):

    self.server.update()
    results = self.completed_results

    while not self.is_empty():
      while not results:
        self.server.fetch()

      yield results.popleft()


  def submit(self, target, args = (), kwargs = {}):

    return self.server.submit(
      identifier = self.identifier,
      target = target,
      args = args,
      kwargs = kwargs,
      )


  def is_empty(self):

    return ( self.server.unfinisheds_for( identifier = self.identifier ) == 0
      and not self.completed_results )


  def is_full(self):

    return self.server.is_full()


  def shutdown(self):

    self.server.shutdown()


  def resume(self):

    self.server.resume()


  def join(self):

    self.server.join()


  def terminate(self):

    self.server.terminate()
