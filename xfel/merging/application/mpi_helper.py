from __future__ import absolute_import, division, print_function
from six.moves import range
from libtbx.mpi4py import MPI
from dials.array_family import flex

import sys
def system_exception_handler(exception_type, value, traceback):
  abort_all_mpi_processes()

sys.excepthook = system_exception_handler

def abort_all_mpi_processes():
  try:
      import libtbx.mpi4py as mpi4py
      rank = mpi4py.MPI.COMM_WORLD.Get_rank()
      from traceback import print_exception
      sys.stderr.write("\nAborting all MPI processes because of exception in process %d:\n"%rank)
      print_exception(exception_type, value, traceback)
      sys.stderr.write("\n")
      sys.stderr.flush()
  finally:
      try:
          import libtbx.mpi4py as mpi4py
          mpi4py.MPI.COMM_WORLD.Abort(1)
      except Exception as e:
          sys.stderr.write("\nFailed to abort MPI process\n")
          sys.stderr.flush()
          raise e

class mpi_helper(object):
  def __init__(self):
    self.MPI = MPI
    self.comm = self.MPI.COMM_WORLD
    self.rank = self.comm.Get_rank()
    self.size = self.comm.Get_size()
    self.error = 0

  def time(self):
    return self.MPI.Wtime()

  def finalize(self):
    self.MPI.Finalize()

  def cumulative_flex(self, flex_array, flex_type):
    '''Build a cumulative sum flex array out of multiple same-size flex arrays.'''
    # Example: (a1,a2,a3) + (b1, b2, b3) = (a1+b1, a2+b2, a3+b3)
    if self.rank == 0:
      cumulative = flex_type(flex_array.size(), 0)
    else:
      cumulative = None

    list_of_all_flex_arrays = self.comm.gather(flex_array, 0)

    if self.rank == 0:
      for i in range(len(list_of_all_flex_arrays)):
        flex_array = list_of_all_flex_arrays[i]
        if flex_array is not None:
          cumulative += flex_array

    return cumulative

  def aggregate_flex(self, flex_array, flex_type):
    '''Build an aggregate flex array out of multiple flex arrays'''
    # Example: (a1,a2,a3) + (b1, b2, b3) = (a1, a2, a3, b1, b2, b3)
    if self.rank == 0:
      aggregate = flex_type()
    else:
      aggregate = None

    list_of_all_flex_arrays = self.comm.gather(flex_array, 0)

    if self.rank == 0:
      for i in range(len(list_of_all_flex_arrays)):
        flex_array = list_of_all_flex_arrays[i]
        if flex_array is not None:
          aggregate.extend(flex_array)

    return aggregate

  def sum(self, data, root=0):
    return self.comm.reduce(data, self.MPI.SUM, root=root)

  def set_error(self, error):
    self.error = error

  def abort_on_error(self):
    all_errors = self.comm.allreduce([self.error], self.MPI.SUM)
    if all_errors.count(1) != 0:
      sys.stderr.write("\nError detected... Aborting MPI process %d\n"%self.rank)
      sys.stderr.flush()
      self.comm.Abort(1)
