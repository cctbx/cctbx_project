"""Density modification to enhance chain ends"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.grow_density

from mmtbx import grow_density
import sys

if(__name__ == "__main__"):
  grow_density.cmd_run(
    args         = sys.argv[1:],
    command_name = "phenix.grow_density")

