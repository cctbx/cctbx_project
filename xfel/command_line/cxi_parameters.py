from __future__ import absolute_import, division, print_function
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.parameters
#
# Utility for printing out default parameters for cctbx.xfel
#

if __name__=='__main__':

  from spotfinder.applications.xfel import cxi_phil
  horizons_phil = cxi_phil.cxi_versioned_extract()
  horizons_phil.persist.show()
