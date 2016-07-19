from __future__ import division

from dials.array_family import flex
from matplotlib import pyplot as plt

# get_hitrate_stats takes a tuple (run, trial, rungroup, d_min)
# and returns a tuple of flex arrays as follows:
# time (s) -- flex.double, timestamp of the shot,
# n_strong -- flex.int, number of strong spots identified by hitfinder,
# I_sig_I_low -- flex.double, the average I/sig(I) in the low res bin of each shot, if it indexed
# I_sig_I_high -- flex.double, the average I/sig(I) in the high res bin of each shot, if it indexed

def get_run_stats(timestamps,
                   n_strong,
                   isigi_low,
                   isigi_high,
                   tuple_of_timestamp_boundaries,
                   lengths,
                   run_numbers,
                   n_strong_cutoff):
  iterator = xrange(len(isigi_low))
  # indexing rate in a sliding window
  half_idx_rate_window = min(50, int(len(isigi_low)//20))
  idx_low_sel = (isigi_low > 0) & (n_strong >= n_strong_cutoff)
  idx_high_sel = (isigi_high > 0) & (n_strong >= n_strong_cutoff)
  idx_rate = flex.double()
  for i in iterator:
    idx_min = max(0, i - half_idx_rate_window)
    idx_max = min(i + half_idx_rate_window, len(isigi_low))
    idx_span = idx_max - idx_min
    idx_sel = idx_low_sel[idx_min:idx_max]
    idx_local_rate = idx_sel.count(True)/idx_span
    idx_rate.append(idx_local_rate)
  # hit rate as dependent on n_strong_cutoff
  hits = n_strong >= n_strong_cutoff
  return (timestamps,
          n_strong,
          hits,
          idx_rate,
          idx_low_sel,
          idx_high_sel,
          isigi_low,
          isigi_high,
          half_idx_rate_window*2,
          lengths,
          tuple_of_timestamp_boundaries,
          run_numbers)

def plot_run_stats(stats, d_min, interactive=True, xsize=100, ysize=100):
  t, n_strong, hits, idx_rate, idx_low_sel, idx_high_sel, isigi_low, isigi_high, \
  window, lengths, boundaries, run_numbers = stats
  if len(t) == 0:
    return None
  f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=False)
  ax1.scatter(t.select(~idx_low_sel), n_strong.select(~idx_low_sel), edgecolors="none", color ='grey')
  ax1.scatter(t.select(idx_low_sel), n_strong.select(idx_low_sel), edgecolors="none", color='blue')
  ax1.axis('tight')
  ax1.set_ylabel("strong spots\nblue: indexed\ngray: did not index").set_fontsize(4)
  ax2.plot(t, idx_rate*100)
  ax2.axis('tight')
  ax2.set_ylabel("indexing rate (%%)\nper %d frames" % window).set_fontsize(4)
  ax3.scatter(t, isigi_low, edgecolors="none", color='red')
  ax3.scatter(t, isigi_high, edgecolors="none", color='orange')
  ax3.axis('tight')
  ax3.set_ylabel("signal-to-noise\nred: low res\nyellow: %3.1f Angstroms" % d_min).set_fontsize(4)
  for a in [ax1, ax2, ax3, ax4]:
    xlab = a.get_xticklabels()
    ylab = a.get_yticklabels()
    for l in xlab + ylab:
      l.set_fontsize(4)
  f.subplots_adjust(hspace=0)
  # add lines and text summaries at the timestamp boundaries
  for boundary in boundaries:
    for a in (ax1, ax2, ax3):
      a.axvline(x=boundary, ymin=0, ymax=3, linewidth=1, color='k')
  run_starts = boundaries[0::2]
  run_ends = boundaries[1::2]
  start = 0
  end = -1
  for idx in xrange(len(run_numbers)):
    start_t = run_starts[idx]
    end_t = run_ends[idx]
    end += lengths[idx]
    slice_t = t[start:end+1]
    slice_hits = hits[start:end+1]
    n_hits = slice_hits.count(True)
    slice_idx_low_sel = idx_low_sel[start:end+1]
    slice_idx_high_sel = idx_high_sel[start:end+1]
    n_idx_low = slice_idx_low_sel.count(True)
    n_idx_high = slice_idx_high_sel.count(True)
    ax4.text(start_t, .9, "run %d" % run_numbers[idx]).set_fontsize(4)
    ax4.text(start_t, .7, "%d f/%d h" % (lengths[idx], n_hits)).set_fontsize(4)
    ax4.text(start_t, .5, "%d (%d) idx" % (n_idx_low, n_idx_high)).set_fontsize(4)
    ax4.text(start_t, .3, "%-3.1f%% hit" % (100*n_hits/lengths[idx],)).set_fontsize(4)
    ax4.text(start_t, .1, "%-3.1f (%-3.1f)%% idx" % \
      (100*n_idx_low/lengths[idx], 100*n_idx_high/lengths[idx])).set_fontsize(4)
    ax4.set_xlabel("timestamp (s)\n# images shown as all (%3.1f Angstroms)" % d_min).set_fontsize(4)
    ax4.set_yticks([])
    for item in [ax1, ax2, ax3, ax4]:
      item.tick_params(labelsize=4)
    start += lengths[idx]
  if interactive:
    plt.show()
  else:
    f.savefig("runstats_tmp.png", bbox_inches='tight', figsize=(xsize, ysize), dpi=300)
    plt.close(f)
    return "runstats_tmp.png"

def plot_multirun_stats(runs,
                        run_numbers,
                        d_min,
                        n_strong_cutoff=40,
                        interactive=False,
                        compress_runs=True,
                        xsize=100,
                        ysize=100):
  tset = flex.double()
  nset = flex.int()
  I_sig_I_low_set = flex.double()
  I_sig_I_high_set = flex.double()
  boundaries = []
  lengths = []
  runs_with_data = []
  offset = 0
  for idx in xrange(len(runs)):
    r = runs[idx]
    if len(r[0]) > 0:
      if compress_runs:
        tslice = r[0] - r[0][0] + offset
        offset += (r[0][-1] - r[0][0])
      else:
        tslice = r[0]
      last_end = r[0][-1]
      tset.extend(tslice)
      nset.extend(r[1])
      I_sig_I_low_set.extend(r[2])
      I_sig_I_high_set.extend(r[3])
      boundaries.append(tslice[0])
      boundaries.append(tslice[-1])
      lengths.append(len(tslice))
      runs_with_data.append(run_numbers[idx])
  stats_tuple = get_run_stats(tset,
                              nset,
                              I_sig_I_low_set,
                              I_sig_I_high_set,
                              tuple(boundaries),
                              tuple(lengths),
                              runs_with_data,
                              n_strong_cutoff=n_strong_cutoff)
  png = plot_run_stats(stats_tuple, d_min, interactive=interactive)
  return png

if __name__ == "__main__":
  import sys
  plot_multirun_stats(sys.argv[1])
