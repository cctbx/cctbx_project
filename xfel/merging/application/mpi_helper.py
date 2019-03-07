from __future__ import division
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

  def cumulative_flex(self, data, flex_type):
    '''Build a cumulative flex array out of flex arrays from individual ranks'''
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
