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

  def extend_flex(self, data, flex_type):
    '''Build an extended flex array out of multiple flex arrays.'''
    extended = None

    all_data = self.comm.gather(data, 0)

    if self.rank == 0:
      if len(all_data) > 0:
        extended = all_data[0]
        for i in range(1, self.size):
          extended.extend(all_data[i])

    return extended

  def cumulative_flex(self, data, flex_type):
    '''Build a cumulative sum flex array out of multiple flex arrays. All arays must be the same size.'''
    if self.rank == 0: # only rank 0 will actually get the cumulative array
      cumulative = flex_type(data.size(), 0)
    else:
      cumulative = None

    all_data = self.comm.gather(data, 0)

    if self.rank == 0:
      for i in range(data.size()):
        for j in range(self.size):
          cumulative[i] += all_data[j][i]

    return cumulative

  def sum(self, data):
    return self.comm.reduce(data, self.MPI.SUM, 0)
