from __future__ import absolute_import, division, print_function

from builtins import object
from collections import deque
from six.moves import range, zip

class single_pooler(object):
  """
  Does not group jobs
  """

  @staticmethod
  def submit_one_job(calcsiter, manager):

    ( target, args, kwargs ) = next(calcsiter)
    return (
      manager.submit( target = target, args = args, kwargs = kwargs ),
      ( target, args, kwargs ),
      )


  @staticmethod
  def process_one_job(calculation, result):

    return ( calculation, result )


  @staticmethod
  def insert_result_into(result, resiter):

    resiter.append( result )


class pooled_run(object):
  """
  Runs multiple jobs
  """

  def __init__(self, calculations):

    self.calculations = calculations


  def __call__(self):

    from libtbx.scheduling import result
    results = []

    for ( target, args, kwargs ) in self.calculations:
      try:
        value = target( *args, **kwargs )

      except Exception as e:
        results.append( result.error( exception = e ) )

      else:
        results.append( result.success( value = value ) )

    return results


class multi_pooler(object):
  """
  Pools up a number of jobs and runs this way

  size - number of jobs to pool
  """

  def __init__(self, size):

    assert 0 < size
    self.size = size


  def submit_one_job(self, calcsiter, manager):

    calculations = [ c for ( i, c ) in zip( range(self.size), calcsiter ) ]

    if not calculations:
      raise StopIteration

    return (
      manager.submit( target = pooled_run( calculations = calculations ) ),
      calculations,
      )


  @staticmethod
  def process_one_job(calculation, result):

    try:
      results = result()

    except Exception:
      return [ ( c, result ) for c in calculation ]

    else:
      assert len( results ) == len( calculation )
      return [ ( c, r ) for ( c, r ) in zip( calculation, results ) ]


  @staticmethod
  def insert_result_into(result, resiter):

    resiter.extend( result )


class finishing_order(object):
  """
  Orders output finishing order
  """

  @staticmethod
  def next_submitted_job(identifier):

    pass


  @staticmethod
  def next_returned_result(identifier, result, pooler, resiter):

    pooler.insert_result_into( result = result, resiter = resiter )


class submission_order(object):
  """
  Orders output in submission order
  """

  def __init__(self):

    self.identifiers = deque()
    self.result_for = {}


  def next_submitted_job(self, identifier):

    self.identifiers.append( identifier )


  def next_returned_result(self, identifier, result, pooler, resiter):

    self.result_for[ identifier ] = result

    while self.identifiers and self.identifiers[0] in self.result_for:
      top = self.identifiers.popleft()
      pooler.insert_result_into( result = self.result_for[ top ], resiter = resiter )
      del self.result_for[ top ]


class ongoing_state(object):
  """
  The iteration is ongoing, iterable has not been exhausted

  Note this is a singleton
  """

  @staticmethod
  def iterate(pfi):

    ( identifier, calcdata ) = pfi.pooler.submit_one_job(
      calcsiter = pfi.calcsiter,
      manager = pfi.manager,
      )

    pfi.calculation_data_for[ identifier ] = calcdata
    pfi.orderer.next_submitted_job( identifier = identifier )


  @classmethod
  def fillup(cls, pfi):

    if pfi.manager.is_empty():
      cls.iterate( pfi = pfi )

    while not pfi.manager.is_full():
      cls.iterate( pfi = pfi )


  @staticmethod
  def emptyhandler():

    pass


class exhausted_state(object):
  """
  The iteration is finishing, iterable has been exhausted

  Note this is a singleton
  """

  @staticmethod
  def fillup(pfi):

    pass


  @staticmethod
  def emptyhandler():

    raise StopIteration


class iterator(object):
  """
  Creates an iterator that executes calls on a Manager-like object

  calculations - an iterable of calculations yielding ( target, args, kwargs ) tuples
  manager - execution manager
  """

  def __init__(self, calculations, manager, poolsize = 1, keep_input_order = False):

    self.manager = manager
    self.calcsiter = iter( calculations )
    self.resiter = deque()

    assert 0 < poolsize

    if poolsize == 1:
      self.pooler = single_pooler

    else:
      self.pooler = multi_pooler( size = poolsize )

    self.calculation_data_for = {}

    if keep_input_order:
      self.orderer = submission_order()

    else:
      self.orderer = finishing_order

    self.resume()


  def __iter__(self):

    return self


  def __next__(self):

    while not self.resiter:
      self.process_next_one()

      try:
        self.state.fillup( pfi = self )

      except StopIteration:
        self.state = exhausted_state

    return self.resiter.popleft()


  def suspend(self):

    self.state = exhausted_state


  def resume(self):

    self.state = ongoing_state


  # Internal
  def process_next_one(self):

    try:
      ( identifier, result ) = next(self.manager.results()) # raise StopIteration
      processed = self.pooler.process_one_job(
        calculation = self.calculation_data_for[ identifier ],
        result = result,
        )
      self.orderer.next_returned_result(
        identifier = identifier,
        result = processed,
        pooler = self.pooler,
        resiter = self.resiter,
        )

    except StopIteration:
      self.state.emptyhandler()
