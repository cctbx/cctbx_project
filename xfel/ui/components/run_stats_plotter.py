from __future__ import division

from dials.array_family import flex
from matplotlib import pyplot as plt
from xfel.ui.components.timeit import duration
import time

# get_hitrate_stats takes a tuple (run, trial, rungroup, d_min)
# and returns a tuple of flex arrays as follows:
# time (s) -- flex.double, timestamp of the shot,
# ratio -- flex.double, ratio of intensities at two angles in the radial average
# n_strong -- flex.int, number of strong spots identified by hitfinder,
# I_sig_I_low -- flex.double, the average I/sig(I) in the low res bin of each shot, if it indexed
# I_sig_I_high -- flex.double, the average I/sig(I) in the high res bin of each shot, if it indexed

def get_should_have_indexed_timestamps(timestamps,
                                       n_strong,
                                       isigi_low,
                                       n_strong_cutoff,
                                       indexed=False):
  if indexed:
    should_have_indexed_sel = (isigi_low > 0) & (n_strong >= n_strong_cutoff)
  else:
    should_have_indexed_sel = (isigi_low == 0) & (n_strong >= n_strong_cutoff) # isigi = 0 if not indexed
  should_have_indexed = timestamps.select(should_have_indexed_sel)
  return should_have_indexed

def get_multirun_should_have_indexed_timestamps(stats_by_run,
                                                run_numbers,
                                                d_min,
                                                n_strong_cutoff=40,
                                                indexed=False):
  timestamps = []
  for idx in xrange(len(stats_by_run)):
    r = stats_by_run[idx]
    if len(r[0]) > 0:
      timestamps.append(
        get_should_have_indexed_timestamps(r[0], r[3], r[4], n_strong_cutoff, indexed=indexed))
    else:
      timestamps.append(flex.double())
  return (run_numbers, timestamps)

def get_string_from_timestamp(ts, long_form=False):
  import time, math
  time_seconds = int(math.floor(ts))
  time_milliseconds = int(round((ts - time_seconds)*1000))
  time_obj = time.gmtime(time_seconds)
  if long_form:
    string = "%04d-%02d-%02dT%02d:%02dZ%02d.%03d" % (
      time_obj.tm_year,
      time_obj.tm_mon,
      time_obj.tm_mday,
      time_obj.tm_hour,
      time_obj.tm_min,
      time_obj.tm_sec,
      time_milliseconds)
  else:
    string = "%04d%02d%02d%02d%02d%02d%03d" % (
      time_obj.tm_year,
      time_obj.tm_mon,
      time_obj.tm_mday,
      time_obj.tm_hour,
      time_obj.tm_min,
      time_obj.tm_sec,
      time_milliseconds)
  return string

def get_strings_from_timestamps(timestamps, long_form=False):
  import os
  get_strings = lambda ts: get_string_from_timestamp(ts, long_form=long_form)
  names = map(get_strings, timestamps)
  return names

def get_paths_from_timestamps(timestamps,
                              prepend="",
                              tag="idx",
                              ext="cbf",
                              long_form=False):
  import os
  def convert(s):
    timestamp_string = get_string_from_timestamp(s, long_form=long_form)
    name = "%s-%s.%s" % (
      tag,
      timestamp_string,
      ext)
    return name
  names = map(convert, timestamps)
  paths = [os.path.join(prepend, name) for name in names]
  return paths

def get_run_stats(timestamps,
                   two_theta_low,
                   two_theta_high,
                   n_strong,
                   isigi_low,
                   isigi_high,
                   tuple_of_timestamp_boundaries,
                   lengths,
                   run_numbers,
                   ratio_cutoff=1,
                   n_strong_cutoff=40,
                   i_sigi_cutoff=1,
                   ):
  iterator = xrange(len(isigi_low))
  # hit rate of drops (observe solvent) or crystals (observe strong spots)
  # since -1 is used as a flag for "did not store this value", and we want a quotient,
  # set the numerator value to 0 whenever either the numerator or denominator is -1
  invalid = (two_theta_low <= 0) or (two_theta_high < 0) # <= to prevent /0
  numerator = two_theta_high.set_selected(invalid, 0)
  denominator = two_theta_low.set_selected(two_theta_low == 0, 1) # prevent /0
  drop_ratios = numerator/denominator
  drop_hits = drop_ratios >= ratio_cutoff
  xtal_hits = n_strong >= n_strong_cutoff
  # indexing and droplet hit rate in a sliding window
  half_idx_rate_window = min(50, int(len(isigi_low)//20))
  half_hq_rate_window = 500
  low_sel = isigi_low > 0
  high_sel = isigi_high > i_sigi_cutoff
  idx_rate = flex.double()
  hq_rate = flex.double()
  drop_hit_rate = flex.double()
  for i in iterator:
    idx_min = max(0, i - half_idx_rate_window)
    idx_max = min(i + half_idx_rate_window, len(isigi_low))
    idx_span = idx_max - idx_min
    idx_sel = low_sel[idx_min:idx_max]
    idx_local_rate = idx_sel.count(True)/idx_span
    idx_rate.append(idx_local_rate)
    hq_min = max(0, i - half_hq_rate_window)
    hq_max = min(i + half_hq_rate_window, len(isigi_low))
    hq_span = hq_max - hq_min
    hq_sel = low_sel[hq_min:hq_max]
    hq_local_sel = high_sel[hq_min:hq_max]
    if hq_sel.count(True) > 0:
      hq_rate.append(hq_local_sel.count(True)/hq_sel.count(True))
    else:
      hq_rate.append(0)
    drop_sel = drop_hits[idx_min:idx_max]
    drop_local_rate = drop_sel.count(True)/idx_span
    drop_hit_rate.append(drop_local_rate)
  return (timestamps,
          drop_ratios,
          drop_hits,
          drop_hit_rate,
          n_strong,
          xtal_hits,
          idx_rate,
          hq_rate,
          low_sel,
          high_sel,
          isigi_low,
          isigi_high,
          half_idx_rate_window*2,
          lengths,
          tuple_of_timestamp_boundaries,
          run_numbers)

def plot_run_stats(stats,
                   d_min,
                   run_tags=[],
                   run_statuses=[],
                   minimalist=False,
                   interactive=True,
                   xsize=30,
                   ysize=10,
                   high_vis=False,
                   title=None,
                   ext='cbf',
                   ):
  t1 = time.time()
  plot_ratio = max(min(xsize, ysize)/2.5, 3)
  if high_vis:
    spot_ratio = plot_ratio*4
    text_ratio = plot_ratio*4.5
  else:
    spot_ratio = plot_ratio*2
    text_ratio = plot_ratio*3
  t, drop_ratios, drop_hits, drop_hit_rate, n_strong, xtal_hits, \
  idx_rate, hq_rate, low_sel, high_sel, isigi_low, isigi_high, \
  window, lengths, boundaries, run_numbers = stats
  if len(t) == 0:
    return None
  n_runs = len(boundaries)//2
  if len(run_tags) != n_runs:
    run_tags = [[] for i in xrange(n_runs)]
  if len(run_statuses) != n_runs:
    run_statuses = [None for i in xrange(n_runs)]
  if minimalist:
    print "Minimalist mode activated."
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
    axset = (ax1, ax2, ax3)
  else:
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=False)
    axset = (ax1, ax2, ax3, ax4)
  for a in axset:
    a.tick_params(axis='x', which='both', bottom='off', top='off')
  ax1.scatter(t.select(~low_sel), n_strong.select(~low_sel), edgecolors="none", color ='#d9d9d9', s=spot_ratio)
  ax1.scatter(t.select(low_sel), n_strong.select(low_sel), edgecolors="none", color='blue', s=spot_ratio)
  ax1.set_ylim(ymin=0)
  ax1.axis('tight')
  ax1.set_ylabel("strong spots\nblue: indexed\ngray: did not index", fontsize=text_ratio)
  ax2.plot(t, idx_rate*100)
  ax2_twin = ax2.twinx()
  ax2_twin.plot(t, drop_hit_rate*100, color='green')
  ax2_twin.set_ylim(ymin=0)
  ax2.axis('tight')
  ax2.set_ylabel("blue:% indexed", fontsize=text_ratio)
  ax2_twin.set_ylabel("green: % solvent", fontsize=text_ratio)
  #ax3.scatter(t, isigi_low, edgecolors="none", color='red', s=spot_ratio)
  gtz = isigi_high > 0
  ax3.scatter(t.select(gtz), isigi_high.select(gtz), edgecolors="none", color='orange', s=spot_ratio)
  ax3.set_ylim(ymin=0)
  ax3_twin = ax3.twinx()
  ax3_twin.plot(t, hq_rate*100, color='orange')
  ax3_twin.set_ylim(ymin=0)
  ax3.axis('tight')
  ax3.set_ylabel("I/sig(I)\nred: low\nyellow: %3.1f Ang" % d_min, fontsize=text_ratio)
  ax3_twin.set_ylabel("line:% HQ", fontsize=text_ratio)
  axset_with_twins = list(axset) + [ax2_twin, ax3_twin]
  for a in axset_with_twins:
    xlab = a.get_xticklabels()
    ylab = a.get_yticklabels()
    for l in xlab + ylab:
      l.set_fontsize(text_ratio)
  f.subplots_adjust(hspace=0)
  # add lines and text summaries at the timestamp boundaries
  if not minimalist:
    for boundary in boundaries:
      if boundary is not None:
        for a in (ax1, ax2, ax3):
          a.axvline(x=boundary, ymin=0, ymax=3, linewidth=1, color='k')
  run_starts = boundaries[0::2]
  run_ends = boundaries[1::2]
  start = 0
  end = -1
  for idx in xrange(len(run_numbers)):
    start_t = run_starts[idx]
    end_t = run_ends[idx]
    if start_t is None or end_t is None: continue
    end += lengths[idx]
    slice_t = t[start:end+1]
    slice_hits = xtal_hits[start:end+1]
    n_hits = slice_hits.count(True)
    slice_drops = drop_hits[start:end+1]
    n_drops = slice_drops.count(True)
    slice_low_sel = low_sel[start:end+1]
    slice_high_sel = high_sel[start:end+1]
    n_idx_low = slice_low_sel.count(True)
    n_idx_high = slice_high_sel.count(True)
    tags = run_tags[idx]
    status = run_statuses[idx]
    if status == "DONE":
      status_color = 'blue'
    elif status in ["RUN", "PEND", "SUBMITTED"]:
      status_color = 'green'
    elif status is None:
      status_color = 'black'
    else:
      status_color = 'red'
    if minimalist:
      ax3.set_xlabel("timestamp (s)\n# images shown as all (%3.1f Angstroms)" % d_min, fontsize=text_ratio)
      ax3.set_yticks([])
    else:
      ax4.text(start_t, 3.85, " " + ", ".join(tags) + " [%s]" % status, fontsize=text_ratio, color=status_color)
      ax4.text(start_t, .85, "run %d" % run_numbers[idx], fontsize=text_ratio)
      ax4.text(start_t, .65, "%d img/%d hit" % (lengths[idx], n_hits), fontsize=text_ratio)
      ax4.text(start_t, .45, "%d (%d) idx" % (n_idx_low, n_idx_high), fontsize=text_ratio)
      ax4.text(start_t, .25, "%-3.1f%% solv/%-3.1f%% xtal" % ((100*n_drops/lengths[idx]),(100*n_hits/lengths[idx])), fontsize=text_ratio)
      ax4.text(start_t, .05, "%-3.1f (%-3.1f)%% idx" % \
        (100*n_idx_low/lengths[idx], 100*n_idx_high/lengths[idx]), fontsize=text_ratio)
      ax4.set_xlabel("timestamp (s)\n# images shown as all (%3.1f Angstroms)" % d_min, fontsize=text_ratio)
      ax4.set_yticks([])
    for item in axset:
      item.tick_params(labelsize=text_ratio)
    start += lengths[idx]
  if title is not None:
    plt.title(title)
  if interactive:
    def onclick(event):
      import math
      ts = event.xdata
      diffs = flex.abs(t - ts)
      ts = t[flex.first_index(diffs, flex.min(diffs))]
      print get_paths_from_timestamps([ts], tag="shot", ext=ext)[0]

    f.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
  else:
    f.set_size_inches(xsize, ysize)
    f.savefig("runstats_tmp.png", bbox_inches='tight', dpi=100)
    plt.close(f)
    t2 = time.time()
    # print "plot_run_stats took %s" % duration(t1, t2)
    return "runstats_tmp.png"

def plot_multirun_stats(runs,
                        run_numbers,
                        d_min,
                        ratio_cutoff=1,
                        n_strong_cutoff=40,
                        i_sigi_cutoff=1,
                        run_tags=[],
                        run_statuses=[],
                        minimalist=False,
                        interactive=False,
                        easy_run=False,
                        compress_runs=True,
                        xsize=30,
                        ysize=10,
                        high_vis=False,
                        title=None):
  tset = flex.double()
  two_theta_low_set = flex.double()
  two_theta_high_set = flex.double()
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
      two_theta_low_set.extend(r[1])
      two_theta_high_set.extend(r[2])
      nset.extend(r[3])
      I_sig_I_low_set.extend(r[4])
      I_sig_I_high_set.extend(r[5])
      boundaries.append(tslice[0])
      boundaries.append(tslice[-1])
      lengths.append(len(tslice))
      runs_with_data.append(run_numbers[idx])
    else:
      boundaries.extend([None]*2)
  stats_tuple = get_run_stats(tset,
                              two_theta_low_set,
                              two_theta_high_set,
                              nset,
                              I_sig_I_low_set,
                              I_sig_I_high_set,
                              tuple(boundaries),
                              tuple(lengths),
                              runs_with_data,
                              ratio_cutoff=ratio_cutoff,
                              n_strong_cutoff=n_strong_cutoff,
                              i_sigi_cutoff=i_sigi_cutoff)
  if easy_run:
    from libtbx import easy_run, easy_pickle
    easy_pickle.dump("plot_run_stats_tmp.pickle", (stats_tuple, d_min, run_tags, run_statuses, minimalist, interactive, xsize, ysize, high_vis, title))
    result = easy_run.fully_buffered(command="cctbx.xfel.plot_run_stats_from_stats_pickle plot_run_stats_tmp.pickle")
    #try:
    #  png = result.stdout_lines[0]
    #  if png == "None":
    #    return None
    #except Exception:
    #  return None
    return None
  else:
    png = plot_run_stats(stats_tuple, d_min, run_tags=run_tags, run_statuses=run_statuses, minimalist=minimalist,
      interactive=interactive, xsize=xsize, ysize=ysize, high_vis=high_vis, title=title)
  return png

if __name__ == "__main__":
  import sys
  easy_run_plot_multirun_stats(sys.argv[1])

