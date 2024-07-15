# LIBTBX_SET_DISPATCHER_NAME libtbx.resource_monitor_plot
from __future__ import division

import argparse
import inspect

from libtbx.resource_monitor import plot_logs


if __name__ == '__main__':  # make the plot in case the original monitor failed
  parser = argparse.ArgumentParser(description=str(inspect.getdoc(plot_logs)))
  parser.add_argument('prefix', nargs='?', type=str, default='monitor*.log',
                      help='Glob matching all log files to be plotted')
  parser.add_argument('-o', '--output', type=str, default='monitor.png',
                      help='Filepath to save the final plot under')
  args = parser.parse_args()
  plot_logs(log_glob=args.prefix, save_path=args.output)
