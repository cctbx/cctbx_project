from __future__ import absolute_import, division, print_function
from six.moves import range
from libtbx.mpi4py import MPI
from dials.array_family import flex

class mpi_helper(object):
  def __init__(self):
    self.MPI = MPI
    self.comm = self.MPI.COMM_WORLD
    self.rank = self.comm.Get_rank()
    self.size = self.comm.Get_size()

  def time(self):
    return self.MPI.Wtime()

  def finalize(self):
    self.MPI.Finalize()

  def cumulative_flex(self, flex_array, flex_type):
    '''Build a cumulative sum flex array out of multiple flex arrays. All arays must be the same size.'''
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
