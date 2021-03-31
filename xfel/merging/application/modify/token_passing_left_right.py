from __future__ import division, print_function
from mpi4py import MPI
mpi_rank = MPI.COMM_WORLD.Get_rank()
mpi_size = MPI.COMM_WORLD.Get_size()
comm = MPI.COMM_WORLD

def token_passing_left_right(value):
  comm.barrier()
  src = mpi_rank - 1 if mpi_rank != 0 else mpi_size - 1
  dst = mpi_rank + 1 if mpi_rank != mpi_size - 1 else 0

  if mpi_rank % 2 == 0:
    comm.send(value, dest=dst)
    m = comm.recv(source=src)
  else:
    m = comm.recv(source=src)
    comm.send(value, dest=dst)
  tokens = [m,value]
  comm.barrier()

  src = mpi_rank + 1 if mpi_rank != mpi_size - 1 else 0
  dst = mpi_rank - 1 if mpi_rank != 0 else mpi_size - 1

  if mpi_rank % 2 == 0:
    comm.send(value, dest=dst)
    m = comm.recv(source=src)
  else:
    m = comm.recv(source=src)
    comm.send(value, dest=dst)
  tokens.append(m)
  comm.barrier()
  return tokens

if __name__=="__main__":
  T = token_passing_left_right(100+mpi_rank)
  print("On Rank %d the tokens are: %s"%(mpi_rank, str(T)))
  print("OK")

