from __future__ import absolute_import, division, print_function
from collections import Counter
from contextlib import contextmanager
from libtbx.mpi4py import MPI
import numpy as np


import sys
def system_exception_handler(exception_type, value, traceback):
  try:
    rank = MPI.COMM_WORLD.Get_rank()
    from traceback import print_exception
    sys.stderr.write("\nTrying to abort all MPI processes because of exception in process %d:\n"%rank)
    print_exception(exception_type, value, traceback)
    sys.stderr.write("\n")
    sys.stderr.flush()
  finally:
    try:
      MPI.COMM_WORLD.Abort(1)
    except Exception as e:
      sys.stderr.write("\nFailed to execute: MPI.COMM_WORLD.Abort(1)\n")
      sys.stderr.flush()
      raise e
sys.excepthook = system_exception_handler


@contextmanager
def adaptive_collective(rooted_variant, non_rooted_variant):
  """
  Declare and switch from rooted to non-rooted collective if root is None.
  Statement `with adaptive_collective(gather, allgather) as gather:` will yield
  `allgather` if root is None, and gather (with an appropriate root) otherwise.
  """
  def adaptive_collective_dispatcher(*args, **kwargs):
    root = kwargs.pop('root', 0)
    if root is None:
      collective = non_rooted_variant
    else:  # if root is an integer
      collective = rooted_variant
      kwargs['root'] = root
    return collective(*args, **kwargs)
  yield adaptive_collective_dispatcher


class mpi_helper(object):
  def __init__(self):
    self.MPI = MPI
    self.comm = self.MPI.COMM_WORLD
    self.rank = self.comm.Get_rank()
    self.size = self.comm.Get_size()
    self.error = (None,None) # (rank,description)

  def time(self):
    return self.MPI.Wtime()

  def finalize(self, mpi_finalize = True):
    if mpi_finalize:
      self.MPI.Finalize()

  def cumulative_flex(self, flex_array, flex_type=None, root=0):
    """
    Build a cumulative sum flex array out of multiple same-size flex arrays.
    Example: (a1,a2,a3) + (b1, b2, b3) = (a1+b1, a2+b2, a3+b3)
    """
    flex_type = flex_type if flex_type is not None else type(flex_array)
    with adaptive_collective(self.comm.gather, self.comm.allgather) as gather:
      list_of_all_flex_arrays = gather(flex_array, root=root)
    if list_of_all_flex_arrays is None:
      return None
    cumulative = flex_type(flex_array.size(), 0)
    for flex_array in list_of_all_flex_arrays:
      if flex_array is not None:
        cumulative += flex_array
    return cumulative

  def aggregate_flex(self, flex_array, flex_type=None, root=0):
    """
    Build an aggregate flex array out of multiple flex arrays
    Example: (a1,a2,a3) + (b1, b2, b3) = (a1, a2, a3, b1, b2, b3)
    """
    flex_type = flex_type if flex_type is not None else type(flex_array)
    with adaptive_collective(self.comm.gather, self.comm.allgather) as gather:
      list_of_all_flex_arrays = gather(flex_array, root=root)
    if list_of_all_flex_arrays is None:
      return None
    aggregate = flex_type()
    for flex_array in list_of_all_flex_arrays:
      if flex_array is not None:
        aggregate.extend(flex_array)
    return aggregate

  def count(self, data, root=0):
    """
    Return total `Counter` of occurrences of each element in data across ranks.
    Example: (a1, a1, a2) + (a1, a2, a3) = {a1: 3, a2: 2, a1: 1}
    """
    with adaptive_collective(self.comm.gather, self.comm.allgather) as gather:
      counters = gather(Counter(data), root=root)
    return sum(counters, Counter()) if counters is not None else None

  def sum(self, data, root=0):
    """
    Sum values of data across all ranks.
    Example: a1 + a2 + a3 = a1+a2+a3
    """
    with adaptive_collective(self.comm.reduce, self.comm.allreduce) as reduce:
      return reduce(data, self.MPI.SUM, root=root)

  def set_error(self, description):
    self.error = (self.rank, description)

  def check_errors(self):
    all_errors = self.comm.allreduce([self.error], self.MPI.SUM)
    actual_errors = [error for error in all_errors if error != (None,None)]
    if len(actual_errors) > 0:
      sys.stderr.write("\nAborting MPI process %d because of the following error(s):"%self.rank)
      for error in actual_errors:
        sys.stderr.write("\nError reported by process %d: %s\n"%(error[0], error[1]))
      sys.stderr.flush()
      self.comm.Abort(1)

  def gather_variable_length_numpy_arrays(self, send_arrays, root=0, dtype=float):
    with adaptive_collective(self.comm.gather, self.comm.allgather) as gather:
      lengths = gather(send_arrays.size, root=root)
    gathered_array = np.empty(np.sum(lengths), dtype=dtype) if lengths else None
    with adaptive_collective(self.comm.Gatherv, self.comm.Allgatherv) as gather_v:
      gather_v(sendbuf=send_arrays, recvbuf=(gathered_array, lengths), root=root)
    return gathered_array
