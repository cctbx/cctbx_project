from __future__ import division

from libtbx.mpi4py import MPI

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

