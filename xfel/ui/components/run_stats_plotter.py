from __future__ import absolute_import, division, print_function
from six.moves import range

from dials.array_family import flex
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
from xfel.ui.components.timeit import duration
import time, math
from six.moves import map

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
                                                n_strong_cutoff=16,
                                                indexed=False):
  timestamps = []
  for idx in range(len(stats_by_run)):
    r = stats_by_run[idx]
    if len(r[0]) > 0:
      timestamps.append(
        get_should_have_indexed_timestamps(r[0], r[3], r[4], n_strong_cutoff, indexed=indexed))
    else:
      timestamps.append(flex.double())
  return (run_numbers, timestamps)

def get_string_from_timestamp(ts, long_form=False):
  import time
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
  names = [get_strings(ts) for ts in  timestamps]
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

def drop_hit_selection(two_theta_low, two_theta_high, ratio_cutoff=1):
  """ Images whose radial average at the high two theta, relative to the low
      one, reaches ratio_cutoff: solvent is present. Returns (ratios, hits).
      -1 flags a value that was not stored; such images are not hits. """
  invalid = (two_theta_low <= 0) | (two_theta_high < 0) # <= to prevent /0
  numerator = two_theta_high.set_selected(invalid, 0)
  denominator = two_theta_low.set_selected(two_theta_low == 0, 1) # prevent /0
  drop_ratios = numerator/denominator
  return drop_ratios, drop_ratios >= ratio_cutoff

def high_quality_selection(resolutions, d_min):
  """ Indexed images whose best resolution bin passing the I/sig(I) cutoff
      reaches d_min. Unindexed images have resolution 0. """
  return (resolutions > 0) & (resolutions <= d_min)

def get_run_stats(timestamps,
                   two_theta_low,
                   two_theta_high,
                   n_strong,
                   resolutions,
                   n_lattices,
                   wavelengths,
                   tuple_of_timestamp_boundaries,
                   lengths,
                   run_numbers,
                   n_multiples=2,
                   ratio_cutoff=1,
                   n_strong_cutoff=16,
                   i_sigi_cutoff=1,
                   d_min=2,
                   ):
  print("")
  print("%d shots" % len(timestamps))
  print("%d first lattices" % (n_lattices >= 1).count(True))
  print("%d multiple lattices" % (n_lattices >= 2).count(True))
  print("%d total lattices" % (flex.sum(n_lattices)))
  iterator = range(len(resolutions))
  # hit rate of drops (observe solvent) or crystals (observe strong spots)
  drop_ratios, drop_hits = drop_hit_selection(two_theta_low, two_theta_high, ratio_cutoff)
  xtal_hits = n_strong >= n_strong_cutoff
  half_idx_rate_window = min(50, max(int(len(timestamps)//20), 1))
  half_hq_rate_window = 500
  indexed_sel = n_lattices > 0
  hq_sel = high_quality_selection(resolutions, d_min)
  # indexing and droplet hit rate in a sliding window
  idx_rate = flex.double()
  multiples_rate = flex.double()
  hq_rate = flex.double()
  drop_hit_rate = flex.double()
  for i in iterator:
    idx_min = max(0, i - half_idx_rate_window)
    idx_max = min(i + half_idx_rate_window, len(resolutions))
    n_lattices_local = n_lattices[idx_min:idx_max]
    shots_this_span = len(n_lattices_local)
    first_lattices_local = n_lattices_local >= 1
    idx_local_rate = first_lattices_local.count(True)/shots_this_span
    idx_rate.append(idx_local_rate)
    multiples_sel = n_lattices_local >= n_multiples
    multiples_local_rate = multiples_sel.count(True)/shots_this_span
    multiples_rate.append(multiples_local_rate)
    drop_sel = drop_hits[idx_min:idx_max]
    drop_local_rate = drop_sel.count(True)/shots_this_span
    drop_hit_rate.append(drop_local_rate)
    # different sliding window for "high quality" xtals
    hq_min = max(0, i - half_hq_rate_window)
    hq_max = min(i + half_hq_rate_window, len(resolutions))
    n_lattices_local_hq = n_lattices[hq_min:hq_max]
    first_lattices_local_hq = n_lattices_local_hq >= 1
    hq_high_sel = hq_sel[hq_min:hq_max].select(first_lattices_local_hq)
    n_first_lattices_local_hq = first_lattices_local_hq.count(True)
    if n_first_lattices_local_hq > 0:
      hq_rate.append(hq_high_sel.count(True)/n_first_lattices_local_hq)
    else:
      hq_rate.append(0)
  return (timestamps,
          drop_ratios,
          drop_hits,
          drop_hit_rate,
          n_strong,
          xtal_hits,
          idx_rate,
          multiples_rate,
          hq_rate,
          indexed_sel,
          hq_sel,
          resolutions,
          half_idx_rate_window*2,
          wavelengths,
          lengths,
          tuple_of_timestamp_boundaries,
          run_numbers)

def plot_run_stats(stats,
                   d_min,
                   n_multiples=2,
                   run_tags=[],
                   run_statuses=[],
                   minimalist=False,
                   interactive=True,
                   xsize=30,
                   ysize=10,
                   high_vis=False,
                   title=None,
                   ext='cbf',
                   figure=None,
                   ):
  t1 = time.time()
  if figure:
    f = figure
  else:
    f = plt.figure()
  plot_ratio = max(min(xsize, ysize)/2.5, 3)
  if high_vis:
    spot_ratio = plot_ratio*4
    text_ratio = plot_ratio*4.5
  else:
    spot_ratio = plot_ratio*2
    text_ratio = plot_ratio*3
  t, drop_ratios, drop_hits, drop_hit_rate, n_strong, xtal_hits, \
  idx_rate, multiples_rate, hq_rate, indexed_sel, hq_sel, resolutions, \
  window, wavelengths, lengths, boundaries, run_numbers = stats
  if len(t) == 0:
    return None
  n_runs = len(boundaries)//2
  if len(run_tags) != n_runs:
    run_tags = [[] for i in range(n_runs)]
  if len(run_statuses) != n_runs:
    run_statuses = [None for i in range(n_runs)]
  if minimalist:
    print("Minimalist mode activated.")
    ax1, ax2, ax3 = f.subplots(3, sharex=True, sharey=False)
    axset = (ax1, ax2, ax3)
  else:
    ax1, ax2, ax3, ax4 = f.subplots(4, sharex=True, sharey=False)
    axset = (ax1, ax2, ax3, ax4)
  ax1.scatter(t.select(~indexed_sel), n_strong.select(~indexed_sel), edgecolors="none", color ='#d9d9d9', s=spot_ratio)
  ax1.scatter(t.select(indexed_sel), n_strong.select(indexed_sel), edgecolors="none", color='blue', s=spot_ratio)
  ax1.set_ylim(ymin=0)
  ax1.axis('tight')
  ax1.set_ylabel("strong spots\nblue: idx\ngray: did not idx", fontsize=text_ratio)
  ax1.set_ylim(0, min(flex.max(n_strong), 5000))
  ax2.plot(t, idx_rate*100)
  ax2.plot(t, multiples_rate*100, color='magenta')
  ax2_twin = ax2.twinx()
  ax2_twin.plot(t, drop_hit_rate*100, color='green')
  ax2_twin.set_ylim(ymin=0)
  ax2.axis('tight')
  ax2.set_ylabel("blue:%% idx\npink:%% %d lattices" % n_multiples, fontsize=text_ratio)
  ax2_twin.set_ylabel("green:\n% solvent", fontsize=text_ratio)
  gtz = resolutions > 0 # no indexing solution stored as resolution of 0
  ax3.scatter(t.select(gtz), 1/(resolutions.select(gtz)**2), edgecolors="none", color='orange', s=spot_ratio)
  def resolution(x,pos):
    if x <= 0:
      return '-'
    return "%.1f"%(1/math.sqrt(x))
  formatter = FuncFormatter(resolution)
  ax3.yaxis.set_major_formatter(formatter)
  ax3_twin = ax3.twinx()
  ax3.set_zorder(0.1)
  ax3_twin.plot(t, hq_rate*100, color='orange')
  ax3_twin.set_ylim(ymin=0)
  ax3.axis('tight')
  ax3.set_ylabel("high res\nbin (Ang)", fontsize=text_ratio)
  ax3_twin.set_ylabel("line:%% <=\n%.2f Ang" % d_min, fontsize=text_ratio)
  axset_with_twins = list(axset) + [ax2_twin, ax3_twin]
  for a in axset_with_twins:
    xlab = a.get_xticklabels()
    ylab = a.get_yticklabels()
    for l in xlab + ylab:
      l.set_fontsize(text_ratio)
  f.subplots_adjust(hspace=0)
  if minimalist:
    ax3.set_xlabel("timestamp (s)\n# images shown as all (%3.1f Angstroms)" % d_min, fontsize=text_ratio)
    ax3.set_yticks([])
  else:
    sel = wavelengths > 0
    energies = 12398.4/wavelengths.select(sel)
    ax4.scatter(t.select(sel), energies, c='lightgray', s=spot_ratio, alpha=0.5)
    ax4.set_xlabel("timestamp (s)\n# images shown as all (%3.1f Angstroms)" % d_min, fontsize=text_ratio)
    ax4.set_ylabel("energy (eV)", fontsize=text_ratio)
  for item in axset:
    item.tick_params(labelsize=text_ratio)
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
  for idx in range(len(run_numbers)):
    start_t = run_starts[idx]
    end_t = run_ends[idx]
    if start_t is None or end_t is None: continue
    end += lengths[idx]
    slice_t = t[start:end+1]
    slice_hits = xtal_hits[start:end+1]
    n_hits = slice_hits.count(True)
    slice_drops = drop_hits[start:end+1]
    n_drops = slice_drops.count(True)
    slice_indexed_sel = indexed_sel[start:end+1]
    slice_hq_sel = hq_sel[start:end+1]
    n_idx_low = slice_indexed_sel.count(True)
    n_idx_high = slice_hq_sel.count(True)
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
    if not minimalist:
      def clipit(ax, x, y, text, max_width, fontsize):
        # Draw text at (x, y), shortened with "..." if it is wider than
        # max_width data units.
        t = ax.text(x, y, text, fontsize=fontsize, transform=ax.get_xaxis_transform())
        width = t.get_window_extent().width
        ax_width = ax.get_window_extent().width
        xmin, xmax = ax.get_xlim()
        max_width_in_pixels = (max_width / (xmax-xmin)) * (ax_width)
        if width <= max_width_in_pixels:
          return
        # Start from the proportional estimate and shorten until it fits.
        n = min(len(text) - 1, int(len(text) * max_width_in_pixels / width))
        for i in range(max(n, 1), 0, -1):
          t.set_text(text[:i] + "...")
          if t.get_window_extent().width <= max_width_in_pixels:
            break

      ax1.text(start_t, 1.05, " " + ", ".join(tags) + " [%s]" % status, fontsize=text_ratio, color=status_color, rotation=10, transform=ax1.get_xaxis_transform(), clip_on=False, zorder=10)

      if end_t:
        width = end_t - start_t
        clipit(ax4, start_t, .85, "run %s" % str(run_numbers[idx]), width, fontsize=text_ratio)
        clipit(ax4, start_t, .65, "%d img/%d hit" % (lengths[idx], n_hits), width, fontsize=text_ratio)
        clipit(ax4, start_t, .45, "%d (%d) idx" % (n_idx_low, n_idx_high), width, fontsize=text_ratio)
        clipit(ax4, start_t, .25, "%-3.1f%% solv/%-3.1f%% xtal" % ((100*n_drops/lengths[idx]),(100*n_hits/lengths[idx])), width, fontsize=text_ratio)
        clipit(ax4, start_t, .05, "%-3.1f (%-3.1f)%% idx" % \
          (100*n_idx_low/lengths[idx], 100*n_idx_high/lengths[idx]), width, fontsize=text_ratio)
      else:
        ax4.text(start_t, .85, "run %s" % str(run_numbers[idx]), fontsize=text_ratio, transform=ax4.get_xaxis_transform())
        ax4.text(start_t, .65, "%d img/%d hit" % (lengths[idx], n_hits), fontsize=text_ratio, transform=ax4.get_xaxis_transform())
        ax4.text(start_t, .45, "%d (%d) idx" % (n_idx_low, n_idx_high), fontsize=text_ratio, transform=ax4.get_xaxis_transform())
        ax4.text(start_t, .25, "%-3.1f%% solv/%-3.1f%% xtal" % ((100*n_drops/lengths[idx]),(100*n_hits/lengths[idx])), fontsize=text_ratio, transform=ax4.get_xaxis_transform())
        ax4.text(start_t, .05, "%-3.1f (%-3.1f)%% idx" % \
          (100*n_idx_low/lengths[idx], 100*n_idx_high/lengths[idx]), fontsize=text_ratio, transform=ax4.get_xaxis_transform())
    start += lengths[idx]
  if title is not None:
    plt.title(title)
  if interactive:
    if not figure:
      plt.show()
  else:
    f.set_size_inches(xsize, ysize)
    f.savefig("runstats_tmp.png", bbox_inches='tight', dpi=100)
    plt.close(f)
    t2 = time.time()
    # print "plot_run_stats took %s" % duration(t1, t2)
    return "runstats_tmp.png"

# Per-run view colours: the series keep the identities of the per-shot plot
# (indexed blue, solvent green, high quality orange); job statuses use
# reserved status colours.
SUMMARY_SERIES = (('hits', 'n_hits', '#eda100'),
                  ('indexed', 'n_indexed', '#2a78d6'),
                  ('high quality', 'n_high_quality', '#eb6834'),
                  ('solvent', 'n_drop_hits', '#1baf7a'))
SUMMARY_STATUS_COLORS = {'DONE': '#898781', 'RUN': '#fab219', 'PEND': '#fab219', 'SUBMITTED': '#fab219'}

def plot_run_summaries(summaries, statuses, figure, high_vis=False):
  """ The per-run view: rates as a percentage of frames over runs, and the
      frames per run as bars coloured by job status. summaries is a list of
      RunSummary in run order and statuses the job status of each. """
  fontsize = 12 if high_vis else 9
  ax_rates, ax_frames = figure.subplots(2, sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.08))
  n = len(summaries)
  x = list(range(n))
  def pct(a, b):
    return 100. * a / b if b else 0.
  for name, key, color in SUMMARY_SERIES:
    ax_rates.plot(x, [pct(getattr(s, key), s.n_frames) for s in summaries], '-',
                  color=color, linewidth=2, marker='o', markersize=3.5, label=name)
  ax_rates.set_ylabel('% of frames', fontsize=fontsize)
  ax_rates.set_ylim(0, 100)
  ax_rates.grid(True, color='#e1e0d9', linewidth=0.8)
  ax_rates.set_axisbelow(True)
  # Legend above the axes so it never covers a line
  ax_rates.legend(loc='lower left', bbox_to_anchor=(0, 1), ncol=4, frameon=False, fontsize=fontsize)
  for ax in (ax_rates, ax_frames):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=fontsize)
  colors = [SUMMARY_STATUS_COLORS.get(status, '#d03b3b') for status in statuses]
  frames = [s.n_frames / 1000. for s in summaries]
  ax_frames.bar(x, frames, width=0.8, color=colors)
  ax_frames.set_ylabel('frames (k)', fontsize=fontsize)
  ax_frames.grid(True, axis='y', color='#e1e0d9', linewidth=0.8)
  ax_frames.set_axisbelow(True)
  # At least ten slots so a few runs do not become one wide bar
  ax_frames.set_xlim(-0.5, max(n, 10) - 0.5)
  step = max(1, n // 20)
  ax_frames.set_xticks(x[::step])
  ax_frames.set_xticklabels([summaries[i].run for i in x[::step]], fontsize=fontsize,
                            rotation=45 if n > 20 else 0, ha='right' if n > 20 else 'center')
  ax_frames.set_xlabel('run', fontsize=fontsize)
  # Name the status of runs not done, once per stretch of the same status
  for i, status in enumerate(statuses):
    if status != 'DONE' and (i == 0 or statuses[i - 1] != status):
      ax_frames.annotate(status, (i, frames[i]), textcoords='offset points', xytext=(0, 3),
                         ha='center', fontsize=fontsize - 2, color='#52514e')
  figure.subplots_adjust(left=0.07, right=0.99, top=0.92, bottom=0.22)

def get_multirun_stats(runs,
                       run_numbers,
                       d_min,
                       n_multiples=2,
                       ratio_cutoff=1,
                       n_strong_cutoff=40,
                       i_sigi_cutoff=1,
                       run_tags=[],
                       run_statuses=[],
                       compress_runs=True):
  """ Concatenate the per-run stats and compute the sliding-window rates.
      Returns (stats_tuple, run_tags, run_statuses) for the runs that have
      data, ready for plot_run_stats. Does no plotting, so it can run off
      the GUI thread. """
  tset = flex.double()
  two_theta_low_set = flex.double()
  two_theta_high_set = flex.double()
  nset = flex.int()
  resolutions_set = flex.double()
  n_lattices = flex.int()
  wavelengths = flex.double()
  boundaries = []
  lengths = []
  runs_with_data = []
  run_tags_with_data = []
  run_statuses_with_data =[]
  offset = 0
  for idx in range(len(runs)):
    r = runs[idx]
    if len(r[0]) > 0:
      if compress_runs:
        tslice = r[0] - r[0][0] + offset
        offset += (r[0][-1] - r[0][0] + 1)
      else:
        tslice = r[0]
      tset.extend(tslice)
      two_theta_low_set.extend(r[1])
      two_theta_high_set.extend(r[2])
      nset.extend(r[3])
      resolutions_set.extend(r[4])
      n_lattices.extend(r[5])
      wavelengths.extend(r[6])
      boundaries.append(tslice[0])
      boundaries.append(tslice[-1])
      lengths.append(len(tslice))
      runs_with_data.append(run_numbers[idx])
      if run_tags:
        run_tags_with_data.append(run_tags[idx])
      if run_statuses:
        run_statuses_with_data.append(run_statuses[idx])
  stats_tuple = get_run_stats(tset,
                              two_theta_low_set,
                              two_theta_high_set,
                              nset,
                              resolutions_set,
                              n_lattices,
                              wavelengths,
                              tuple(boundaries),
                              tuple(lengths),
                              runs_with_data,
                              n_multiples=n_multiples,
                              ratio_cutoff=ratio_cutoff,
                              n_strong_cutoff=n_strong_cutoff,
                              i_sigi_cutoff=i_sigi_cutoff,
                              d_min=d_min)
  return stats_tuple, run_tags_with_data, run_statuses_with_data

def plot_multirun_stats(runs,
                        run_numbers,
                        d_min,
                        n_multiples=2,
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
                        title=None,
                        figure=None):
  stats_tuple, run_tags_with_data, run_statuses_with_data = get_multirun_stats(
    runs, run_numbers, d_min, n_multiples=n_multiples, ratio_cutoff=ratio_cutoff,
    n_strong_cutoff=n_strong_cutoff, i_sigi_cutoff=i_sigi_cutoff,
    run_tags=run_tags, run_statuses=run_statuses, compress_runs=compress_runs)
  if easy_run:
    from libtbx import easy_run, easy_pickle
    easy_pickle.dump("plot_run_stats_tmp.pickle", (stats_tuple, d_min, n_multiples, run_tags_with_data, run_statuses_with_data, minimalist, interactive, xsize, ysize, high_vis, title))
    result = easy_run.fully_buffered(command="cctbx.xfel.plot_run_stats_from_stats_pickle plot_run_stats_tmp.pickle")
    try:
      png = result.stdout_lines[-1]
      if png == "None":
        return None
    except Exception:
      return None
  else:
    png = plot_run_stats(stats_tuple, d_min, n_multiples=n_multiples, run_tags=run_tags_with_data, run_statuses=run_statuses_with_data, minimalist=minimalist,
      interactive=interactive, xsize=xsize, ysize=ysize, high_vis=high_vis, title=title, figure=figure)
  return png
