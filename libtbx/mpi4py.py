from __future__ import absolute_import, division, print_function
import time

''' mpi4py wrapper: emulating mpi4py behavior for a single rank when the real mpi4py is not installed '''

class mpiCommEmulator(object):
  def Get_rank(self):
    return 0
  def Get_size(self):
    return 1
  def barrier(self):
    pass
  def bcast(self, transmitted, root):
    return transmitted
  def reduce(self, data, operation, root):
    if operation == mpiEmulator.SUM or operation == mpiEmulator.MAX or operation == mpiEmulator.MIN:
      return data
    else:
      assert False, "Unsupported MPI reduce operation %s"%(operation)
  def allreduce(self, data, operation):
    return self.reduce(data, operation, 0)
  def alltoall(self, items):
    return items
  def scatter(self, items):
    return items
  def gather(self, item, root):
    items = []
    items.append(item)
    return items
  def Abort(self,error):
    pass

class mpiEmulator(object):
  COMM_WORLD = mpiCommEmulator()

  SUM = "SUM"
  MAX = "MAX"
  MIN = "MIN"
  # TODO: implement more operations as needed

  def Wtime(self):
    return time.time()
  def Finalize(self):
    pass

try:
  from mpi4py import MPI
except ImportError:
  print ("\nWarning: could not import mpi4py. Running as a single process.\n")
  MPI = mpiEmulator()
