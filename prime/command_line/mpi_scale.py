from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.mpi_scale
"""
Find initial scaling factors for all integration results
"""
from mpi4py import MPI
import sys
from prime.postrefine.mod_input import process_input, read_pickles
import numpy as np

# setup mpi
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
assert size>1

def master(frame_indices):
  for frame_index in frame_indices:
    rankreq, isdata = comm.recv(source=MPI.ANY_SOURCE)
    if not isdata:
      comm.send(frame_index, dest=rankreq)
  flag_done = np.zeros(size)
  while np.sum(flag_done) < size-1:
    rankreq, _ = comm.recv(source=MPI.ANY_SOURCE)
    comm.send('endrun', dest=rankreq)
    flag_done[rankreq] = 1

def client(frame_files, iparams):
  result = {}
  while True:
    comm.send((rank, False), dest=0)
    frame_index = comm.recv(source=0)
    if str(frame_index) == 'endrun':
      break
    #start scaling
    from prime.postrefine import postref_handler
    prh = postref_handler()
    pres, _ = prh.scale_frame_by_mean_I(frame_index, frame_files[frame_index], iparams, 0, 'average')
    result[frame_index] = pres
  return result

def run(argv):
  #broadcast parameters
  if rank == 0:
    iparams, txt_out_input = process_input(argv)
    frame_files = read_pickles(iparams.data)
  else:
    iparams = None
    frame_files = None
  frame_files = comm.bcast(frame_files, root=0)
  iparams = comm.bcast(iparams, root=0)
  #assign scaling task
  if rank == 0:
    master(range(len(frame_files)))
    result = {}
  else:
    result = client(frame_files, iparams)
  result = comm.gather(result, root=0)
  if rank == 0:
    print "Scaling is done on %d cores"%len(result)
  MPI.Finalize()

if __name__ == "__main__":
  argv = sys.argv[1:] if len(sys.argv) > 1 else None
  run(argv)
