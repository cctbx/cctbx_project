from __future__ import absolute_import, print_function, division
import matplotlib.pyplot as plt
import sys
import os
from iotbx.detectors.cspad_detector_formats import reverse_timestamp
from libtbx.phil import parse
from libtbx.utils import Sorry
from scitbx.array_family import flex
from scitbx.math import five_number_summary

message = """
A script to get a sense of the computational performance of every rank while
processing data. End product is a plot of wall time vs MPI rank number with
every data point being that of a frame processed by dials.stills_process.
The information is read in from the debug files created by dials.stills_process.
Example usage on cxic0415 processed demo data:
    cctbx.xfel.weather input_path=cxic0415/output/debug
""".strip()
phil_scope = parse('''
  input_path = .
    .type = str
    .help = path to where the processing results are. For example path to XXX_rgYYYY
  num_nodes = 1
    .type = int
    .help = Number of nodes used to do data processing. Used in timing information
  num_cores_per_node = 72
    .type = int
    .help = Number of cores per node in the machine (default is for Cori KNL)
  wall_time = None
    .type = int
    .help = total wall time (seconds) taken for job to finish. Used for plotting node-partitioning
  plot_title = Computational weather plot
    .type = str
    .help = title of the computational weather plot
  show_plot = True
    .type = bool
    .help = flag to indicate if plot should be displayed on screen
  pickle_plot = False
    .type = bool
    .help = If True, will pickle matplotlib session so that it can be opened later for analysis/viewing \
            https://stackoverflow.com/questions/29160177/matplotlib-save-file-to-be-reedited-later-in-ipython
  pickle_filename = fig_object.pickle
    .type = str
    .help = Default name of pickled matplotlib plot saved to disk
''')


def params_from_phil(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      user_phil.append(parse(file_name=arg))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s" % arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  return params


def timestamp_to_seconds(ts):
  sec, ms = reverse_timestamp(ts)
  return sec + ms * 0.001


def run(params):
  root = params.input_path
  fig_object = plt.figure()
  good_total = fail_total = 0
  fail_deltas = []
  good_deltas = []
  rank_walltimes = []
  for filename in os.listdir(root):
    if os.path.splitext(filename)[1] != '.txt': continue
    if 'debug' not in filename: continue
    reference = None
    fail_timepoints = []
    good_timepoints = []
    run_timepoints = []
    rank = int(filename.split('_')[1].split('.')[0])
    for line in open(os.path.join(root, filename)):
      try:
        hostname, psanats, ts, status, result = line.strip().split(',')
      except ValueError:
        continue
      if reference is None:
        reference = timestamp_to_seconds(ts)
        run_timepoints.append(0)
        assert status not in ['stop', 'done', 'fail']

      if status in ['stop', 'done', 'fail']:
        timepoint = timestamp_to_seconds(ts) - reference
        run_timepoints.append(timepoint)
        if status == 'done':
          good_timepoints.append(timepoint)
          good_deltas.append(good_timepoints[-1] - run_timepoints[-2])
        else:
          fail_timepoints.append(timepoint)
          fail_deltas.append(fail_timepoints[-1] - run_timepoints[-2])
        processing_of_most_recent_still_terminated = True
      else:
        processing_of_most_recent_still_terminated = False
    rank_walltimes.append(timestamp_to_seconds(ts) - reference)
    plt.plot(fail_timepoints, [rank]*len(fail_timepoints), 'b.')
    plt.plot(good_timepoints, [rank]*len(good_timepoints), 'g.')
    fail_total += len(fail_timepoints)
    good_total += len(good_timepoints)
    if not processing_of_most_recent_still_terminated:
      plt.plot([rank_walltimes[-1]], [rank], 'rx')

  msg = "Five number summary of {} (s): {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}, {:7.2f}"
  if fail_deltas:
    process = '{:5d} fail image processing times'.format(fail_total)
    print(msg.format(process, *five_number_summary(flex.double(fail_deltas))))
  if good_deltas:
    process = '{:5d} good image processing times'.format(good_total)
    print(msg.format(process, *five_number_summary(flex.double(good_deltas))))
  if rank_walltimes:
    process = "{:5d} individual ranks' walltimes".format(len(rank_walltimes))
    print(msg.format(process, *five_number_summary(flex.double(rank_walltimes))))

  if params.wall_time and params.num_nodes and params.num_cores_per_node:
    for i in range(params.num_nodes):
      plt.plot([0, params.wall_time], [i*params.num_cores_per_node-0.5, i*params.num_cores_per_node-0.5], 'r-')
  plt.xlabel('Wall time (sec)')
  plt.ylabel('MPI Rank Number')
  plt.title(params.plot_title)
  if params.pickle_plot:
    from libtbx.easy_pickle import dump
    dump('%s' % params.pickle_filename, fig_object)
  if params.show_plot:
    plt.show()


if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
  run(params)
