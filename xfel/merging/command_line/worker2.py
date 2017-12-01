from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME mpi.worker2

print "Deprecated on 12/01/17. Code will be removed at a later date"
exit()

from mpi4py import MPI
from xfel.command_line.mpi_merge import scaling_manager_mpi

comm = MPI.Comm.Get_parent()
size = comm.Get_size()
rank = comm.Get_rank()

received_info = {}
received_info = comm.bcast(received_info, root=0)
file_names = received_info["file_names"]
sm = scaling_manager_mpi(received_info["miller_set"],
                         received_info["model"],
                         received_info["params"])

assert sm.params.backend == 'MySQL' # only option that makes sense
from xfel.merging.database.merging_database import manager
db_mgr = manager(sm.params)

for ix in xrange(len(file_names)):
  if ix%size == rank:
    print "Rank %d processing %s"%(rank, file_names[ix])
    sm.tar_to_scale_frame_adapter(tar_list=[file_names[ix],], db_mgr=db_mgr)

comm.gather(sm, root=0)
comm.Disconnect()
