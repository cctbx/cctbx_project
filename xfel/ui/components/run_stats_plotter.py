from __future__ import division

from dials.array_family import flex
from matplotlib import pyplot as plt

# get_hitrate_stats takes a tuple (run, trial, rungroup)
# and returns a tuple of flex arrays as follows:
# time (s) -- flex.double, timestamp of the shot,
# n_strong -- flex.int, number of strong spots identified by hitfinder,
# I_sig_I_low -- flex.double, the average I/sig(I) in the low res bin of each shot, if it indexed
# I_sig_I_high -- flex.double, the average I/sig(I) in the high res bin of each shot, if it indexed

def plot_run_stats(timestamps,
                   n_strong,
                   i_sig_i_low,
                   i_sig_i_high,
                   tuple_of_timestamp_boundaries,
                   lengths,
                   run_numbers):
  iterator = xrange(len(i_sig_i_low))
  # indexing rate in a sliding window
  half_idx_rate_window = min(50, int(len(i_sig_i_low)//20))
  idx_bool = i_sig_i_low > 0
  idx_rate = flex.double()
  for i in iterator:
    idx_min = max(0, i - half_idx_rate_window)
    idx_max = min(i + half_idx_rate_window, len(i_sig_i_low))
    idx_span = idx_max - idx_min
    idx_sele = idx_bool[idx_min:idx_max]
    idx_local_rate = idx_sele.count(True)/idx_span
    idx_rate.append(idx_local_rate)
  idxrate_sel = i_sig_i_low > 0
  t, hits = timestamps, n_strong >= 40
  f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=False)
  ax1.scatter(t, n_strong, edgecolors="none", color ='orange')
  ax1.scatter(t.select(idxrate_sel), n_strong.select(idxrate_sel), edgecolors="none", color='blue')
  ax1.axis('tight')
  ax1.set_ylabel("strong spots")
  ax2.plot(t, idx_rate*100)
  ax2.axis('tight')
  ax2.set_ylabel("indexing rate per %d frames (%%)" % (2*half_idx_rate_window))
  ax3.scatter(t, i_sig_i_low, edgecolors="none", color='red')
  ax3.scatter(t, i_sig_i_high, edgecolors="none", color='orange')
  ax3.axis('tight')
  ax3.set_ylabel("signal-to-noise")
  f.subplots_adjust(hspace=0)
  ax3.set_xlabel("timestamp (s)")
  # add lines and text summaries at the timestamp boundaries
  for boundary in tuple_of_timestamp_boundaries:
    for a in (ax1, ax2, ax3):
      a.axvline(x=boundary, ymin=0, ymax=3, linewidth=1, color='k')
  run_starts = tuple_of_timestamp_boundaries[0::2]
  run_ends = tuple_of_timestamp_boundaries[1::2]
  start = 0
  end = -1
  for idx in xrange(len(run_numbers)):
    start_t = run_starts[idx]
    end_t = run_ends[idx]
    end += lengths[idx]
    slice_t = t[start:end+1]
    slice_hits = hits[start:end+1]
    n_hits = slice_hits.count(True)
    slice_idx_bool = idx_bool[start:end+1]
    n_idx = slice_idx_bool.count(True)
    ax4.text(start_t, .9, "run %d" % run_numbers[idx])
    ax4.text(start_t, .7, "%d hits" % lengths[idx])
    ax4.text(start_t, .5, "%d idx" % n_idx)
    ax4.text(start_t, .3, "%-5.1f%% hit" % (100*n_hits/lengths[idx]))
    ax4.text(start_t, .1, "%-5.1f%% idx" % (100*n_idx/lengths[idx]))
    #ax4.gca().yaxis.set_major_locator(plt.NullLocator())
    start += lengths[idx]
  plt.show()

def plot_multirun_stats(runs, run_numbers):
  tset = flex.double()
  nset = flex.int()
  I_sig_I_low_set = flex.double()
  I_sig_I_high_set = flex.double()
  boundaries = []
  lengths = []
  for r in runs:
    if len(r[0]) > 0:
      tset.extend(r[0])
      nset.extend(r[1])
      I_sig_I_low_set.extend(r[2])
      I_sig_I_high_set.extend(r[3])
      boundaries.append(r[0][0])
      boundaries.append(r[0][-1])
      lengths.append(len(r[0]))
  plot_run_stats(tset,
                 nset,
                 I_sig_I_low_set,
                 I_sig_I_high_set,
                 tuple(boundaries),
                 tuple(lengths),
                 run_numbers)

if __name__ == "__main__":
  # plot_multirun_stats(test_runs)
  import sys
  plot_multirun_stats(sys.argv[1])
