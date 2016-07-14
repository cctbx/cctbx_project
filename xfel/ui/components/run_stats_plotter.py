from __future__ import division

from dials.array_family import flex
from matplotlib import pyplot as plt
import random

# get_hitrate_stats takes a tuple (run, trial, rungroup)
# and returns a tuple of flex arrays as follows:
# time (s) -- flex.double, timestamp of the shot,
# n_strong -- flex.int, number of strong spots identified by hitfinder,
# intensities -- flex.double, the average intensity to 10 Angstroms of each shot, if it indexed,
# sigmas -- flex.double, the average sigma to 10 Angstroms of each shot, if it indexed,
# I_sig_I -- flex.double, the average I/sig(I) to 10 Angstroms of each shot, if it indexed

# # test data
# time_zero = 19483.3
# test_runs = []
# for r in xrange(5): # test on some mock runs
#   # timestamps in seconds
#   t = flex.double()
#   for i in xrange(1000):
#     t.append(time_zero)
#     time_zero += i
#   # number strong spots identified
#   n = flex.int()
#   for i in xrange(1000):
#     n.append(max(random.randint(-200, 100), 0))
#   # intensities, sigmas and avg I/sig(I) of spots, if indexed, cut to 10 Angstroms
#   I = flex.double() # not used
#   S = flex.double() # not used
#   I_sig_I = flex.double()
#   for i in xrange(1000):
#     I.append(random.random()*10000 - 3000)
#     S.append(random.random()*1000 - 300)
#     if n[i] > 0:
#       I_sig_I.append(random.random()*6)
#     else:
#       I_sig_I.append(0)
#   test_runs.append((t, n, I, S, I_sig_I))


def plot_run_stats(timestamps, n_shots, i_sig_i, tuple_of_timestamp_boundaries):
  iterator = xrange(len(i_sig_i))
  # indexing rate in a sliding window
  half_idx_rate_window = int(len(i_sig_i)//20)
  idx_bool = flex.bool()
  for i in iterator:
    idx_bool.append(i_sig_i[i] > 0)
  idx_rate = flex.double()
  for i in iterator:
    idx_min = max(0, i - half_idx_rate_window)
    idx_max = min(i + half_idx_rate_window, len(i_sig_i))
    idx_span = idx_max - idx_min
    idx_sele = idx_bool[idx_min:idx_max]
    idx_local_rate = idx_sele.count(True)/idx_span
    idx_rate.append(idx_local_rate)
  # <I/sig(I)> per shot
  # I_over_sig_I = flex.double()
  # for i in iterator:
  #   I_over_sig_I.append(I[i]/S[i])
  # avg_I_over_sig_I = flex.double()
  # for i in iterator:
  #   avg = I_over_sig_I[i].min_max_mean().mean
  #   if avg is not None:
  #     avg_I_over_sig_I.append(avg)
  #   else:
  #     avg_I_over_sig_I.append(0)
  # plot hitrate (as n_strong), indexing rate and <I/sig(I)> on same axes against timestamp, per run
  t, hitrate = timestamps, n_shots
  f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
  ax1.scatter(t, hitrate, edgecolors="none")
  ax1.axis('tight')
  ax1.set_ylabel("strong spots")
  ax2.plot(t, idx_rate)
  ax2.axis('tight')
  ax2.set_ylabel("indexing rate per %d" % (2*half_idx_rate_window))
  ax3.scatter(t, i_sig_i, edgecolors="none")
  ax3.axis('tight')
  ax3.set_ylabel("signal-to-noise")
  f.subplots_adjust(hspace=0)
  ax3.set_xlabel("timestamp (s)")
  plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
  # add lines at the timestamp boundaries
  for boundary in tuple_of_timestamp_boundaries:
    for a in (ax1, ax2, ax3):
      a.axvline(x=boundary, ymin=0, ymax=3, linewidth=1, color='k')
  plt.show()

def plot_multirun_stats(runs):
  tset = flex.double()
  nset = flex.int()
  Iset = flex.double()
  Sset = flex.double()
  I_sig_I_set = flex.double()
  boundaries = []
  for r in runs:
    tset.extend(r[0])
    nset.extend(r[1])
    Iset.extend(r[2])
    Sset.extend(r[3])
    I_sig_I_set.extend(r[4])
    boundaries.append(r[0][0])
    boundaries.append(r[0][-1])
  plot_run_stats(tset, nset, I_sig_I_set, tuple(boundaries))

if __name__ == "__main__":
  # plot_multirun_stats(test_runs)
  import sys
  plot_multirun_stats(sys.argv[1])
