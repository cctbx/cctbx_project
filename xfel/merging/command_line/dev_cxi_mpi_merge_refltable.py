from __future__ import division
from __future__ import print_function
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME dev.cxi.mpi_merge_refltable
#
# $Id$

from xfel.merging.command_line.dev_mpi_cluster_two_merge import scaling_manager_mpi, Script
from xfel.merging.command_line.dev_cxi_merge_refltable import refltable_scaling_manager

class refltable_scaling_manager_mpi(scaling_manager_mpi, refltable_scaling_manager):
  pass

if (__name__ == "__main__"):
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  script = Script(refltable_scaling_manager_mpi)
  result = script.run(comm=comm,timing=False)
  if rank == 0:
    script.show_plot(result)
  print("DONE")
