from __future__ import absolute_import, division, print_function
import time

''' mpi4py wrapper: emulating mpi4py behavior for a single rank when the real mpi4py is not installed '''

class mpiEmulator(object):

  SUM = "SUM"
  MAX = "MAX"
  MIN = "MIN"
  # TODO: implement more operations as needed

  def Wtime(self):
    return time.time()
  def Finalize(self):
    pass

class mpiCommEmulator(object):
  def Get_rank(self):
    return 0
  def Get_size(self):
    return 1
  def barrier(self):
    pass
  def bcast(self, obj, root=0):
    return obj
  def reduce(self, sendobj, op=mpiEmulator.SUM, root=0):
    if op == mpiEmulator.SUM or op == mpiEmulator.MAX or op == mpiEmulator.MIN:
      return sendobj
    else:
      assert False, "Unsupported MPI reduce operation %s"%(op)
  def allreduce(self, sendobj, op=mpiEmulator.SUM):
    return self.reduce(sendobj, op, 0)
  def alltoall(self, sendobj):
    return sendobj
  def scatter(self, sendobj, root=0):
    assert root == 0 and len(sendobj) == 1
    return sendobj[0]
  def gather(self, sendobj, root=0):
    items = []
    items.append(sendobj)
    return items
  def Abort(self,errorcode=0):
    import sys
    sys.exit()
  @property
  def rank(self):
    return self.Get_rank()
  @property
  def size(self):
    return self.Get_size()

mpiEmulator.COMM_WORLD = mpiCommEmulator()

try:
  from mpi4py import MPI
except ImportError:
  print ("\nWarning: could not import mpi4py. Running as a single process.\n")
  MPI = mpiEmulator()
