from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.plot_spotfinder_stats_from_stats_pickle
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from xfel.ui.components.spotfinder_plotter import plot_spotfinder_stats

def easy_run_plot_spotfinder_stats(pickle):
  from libtbx import easy_pickle
  contents = easy_pickle.load(pickle)
  stats_tuple, spot_length_stats_tuple, run_tags, run_statuses, minimalist, interactive, xsize, ysize, high_vis, title = contents
  print(plot_spotfinder_stats(stats_tuple, spot_length_stats_tuple, run_tags=run_tags, run_statuses=run_statuses, minimalist=minimalist,
    interactive=interactive, xsize=xsize, ysize=ysize, high_vis=high_vis, title=title))

if __name__ == "__main__":
  import sys
  easy_run_plot_spotfinder_stats(sys.argv[1])
