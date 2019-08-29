from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME xes.histograms
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import sys

from xfel.cxi.cspad_ana import xes_histograms

if __name__ == '__main__':
  xes_histograms.run(sys.argv[1:])
