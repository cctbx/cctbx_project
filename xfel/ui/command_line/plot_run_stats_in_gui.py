from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.plot_run_stats_from_stats_pickle
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from xfel.ui.components.run_stats_plotter import plot_run_stats

def easy_run_plot_multirun_stats(pickle):
  from libtbx import easy_pickle
  contents = easy_pickle.load(pickle)
  stats_tuple, d_min, n_multiples, run_tags, run_statuses, minimalist, interactive, xsize, ysize, high_vis, title = contents
  print(plot_run_stats(stats_tuple, d_min, n_multiples, run_tags=run_tags, run_statuses=run_statuses, minimalist=minimalist,
    interactive=interactive, xsize=xsize, ysize=ysize, high_vis=high_vis, title=title))
  #plot_run_stats(stats_tuple, d_min, run_tags=run_tags, run_statuses=run_statuses, minimalist=minimalist,
  #  interactive=True, xsize=xsize, ysize=ysize, high_vis=high_vis, title=title)

if __name__ == "__main__":
  import sys
  easy_run_plot_multirun_stats(sys.argv[1])

