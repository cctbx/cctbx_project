from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mpi_fxs_launch
from xfel.amo.pnccd_ana.mpi_fxs_launch import launch
import sys
launch(sys.argv[1:])
