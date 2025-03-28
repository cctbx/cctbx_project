from __future__ import absolute_import, division, print_function
import numpy as np
import time
import sys

''' mpi4py wrapper: emulating mpi4py behavior for a single rank when the real mpi4py is not installed '''

class mpiEmulator(object):

  SUM = "SUM"
  MAX = "MAX"
  MIN = "MIN"
  LAND = "LAND"
  LOR = "LOR"
  LXOR = "LXOR"
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
  def Barrier(self):
    pass
  def bcast(self, obj, root=0):
    return obj
  def Bcast(self, buf, root=0):
    pass
  def reduce(self, sendobj, op=mpiEmulator.SUM, root=0):
    if op in {mpiEmulator.SUM, mpiEmulator.MAX, mpiEmulator.MIN,
              mpiEmulator.LAND, mpiEmulator.LOR, mpiEmulator.LXOR}:
      return sendobj
    else:
      assert False, "Unsupported MPI reduce operation %s"%(op)
  def Reduce(self, sendbuf, recvbuf, op=mpiEmulator.SUM, root=0):
    if isinstance(sendbuf, np.ndarray):  # numpy array
      recvbuf[:] = sendbuf[:]
    else:
      # Handle scalar values or other types
      if isinstance(recvbuf, list):
        for i, value in enumerate(sendbuf):
          recvbuf[i] = value
      else:
        assert False, "Unsupported MPI emulator data type."
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
  def Gatherv(self, sendbuf, recvbuf, root=0):
    assert len(recvbuf) == 2, "Other ways of using Gatherv are not implemented"
    rbuff, counts = recvbuf
    if len(counts) == 1:
      sendbuf = (sendbuf,)
    counter = 0
    for item, count in zip(sendbuf, counts):
      rbuff[counter:counter+count] = item
      counter += count
  def allgather(self, sendobj):
    return [sendobj]
  def Allgatherv(self, sendbuf, recvbuf):
    return self.Gatherv(sendbuf, recvbuf)
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

class MpiDisabledError(Exception):
  pass

try:
  import libtbx
  if libtbx.mpi_import_guard.disable_mpi:
    raise MpiDisabledError
  from mpi4py import MPI
  using_mpi = True
except ImportError:
  print ("\nWarning: could not import mpi4py. Running as a single process.\n")
  MPI = mpiEmulator()
  using_mpi = False
except MpiDisabledError:
  MPI = mpiEmulator()
  using_mpi = False


def mpi_abort_on_exception(func):
  """
  A decorator for functions that will be called in an MPI context. This ensures
  the MPI job will abort if a single rank raises an exception (or exits) out of
  the decorated function.
  """
  def wrapped_func(*args, **kwargs):
    try:
      return func(*args, **kwargs)
    except Exception:
      sys.excepthook(*sys.exc_info())
      MPI.COMM_WORLD.Abort(1)
    except SystemExit as e:
      if e.code:
        sys.stderr.write(e.code + '\n')
      MPI.COMM_WORLD.Abort(1)
  if using_mpi and MPI.COMM_WORLD.size > 1:
    return wrapped_func
  else:
    return func
