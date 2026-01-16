from __future__ import absolute_import, division, print_function
import os, glob
from matplotlib import pyplot as plt

"""
Searches the cctbx.xfel.merge log files for statistics tables.

Makes lots of assumptions regarding formatting.
"""

# Name of stat: (Line signal, offset to first line of table, index of value for a bin line, index of value for an 'All' line)
types = {
  "% accepted": ("Lattices resolution", 6, 4, 1),
  "Multiplicity": ("Intensity Statistics (all accepted experiments)", 14, 6, 3),
  "Completeness": ("Intensity Statistics (all accepted experiments)", 14, 5, 2),
  "CC1/2": ("Table of Scaling Results", 6, 5, 2),
  "Merged I/sigI": ("Intensity Statistics (all accepted experiments)", 14, 11, 7),
}

class Scraper(object):
  def __init__(self, output_path, accepted):
    self.output_path = output_path
    assert accepted in ['%','#']
    self.accepted = accepted

  def scrape(self):
    results = {}
    try:
      path = glob.glob(os.path.join(self.output_path, "*main*"))[0]
    except IndexError: return results
    lines = open(path).readlines()
    for i, line in enumerate(lines):
      for name in types:
        signal, offset, _, _ = types[name]
        if signal in line:
          table = lines[i+offset:]
          for l, dontjudgeme in enumerate(table):
            if dontjudgeme.startswith("All"):
              results[name] = table[:l+1]
              break
    for t in results:
      parsed = []
      for line in results[t]:
        s = line.replace('[', ' ').replace(']', ' ').split()
        if not s: continue
        if s[0] == 'All':
          value = float(s[types[t][3]].rstrip('%'))
          parsed.append((s[0], value))
        else:
          bin_id = int(s[0])
          d_max = float(s[1])
          d_min = float(s[3])
          value = float(s[types[t][2]].rstrip('%'))
          if d_max < 0: d_max = 100
          parsed.append((bin_id, d_max, d_min, value))
      results[t] = parsed
    if '% accepted' in results:
      if self.accepted == '%':
        self._num_to_percent(results, '% accepted')
      else:
        results['# accepted'] = results['% accepted']
        del results['% accepted']
      return results

  def _num_to_percent(self, results, key):
    data = results[key]
    for line in data:
      if line[0] == 'All':
        denom = line[1]
    new_data = []
    for line in data:
      if line[0] == 'All':
        line = line[0], 100*line[1]/denom
      else:
        line = line[0], line[1], line[2], 100*line[3]/denom
      new_data.append(line)
    results[key] = new_data

  def plot_single_results(self, results, title, xsize=30, ysize=10, interactive = True):
    from matplotlib.ticker import FuncFormatter
    import numpy as np
    import math
    fig = plt.figure()
    ax = ax1 = fig.gca()
    ax2 = ax1.twinx()

    colors = {
      "% accepted": 'orange',
      "Multiplicity": 'red',
      "Completeness": 'green',
      "CC1/2": 'blue'
    }

    for name in results:
      if name == 'Multiplicity':
        ax = ax2
      else:
        ax = ax1

      x = []; y = []
      for data in results[name]:
        if data[0] == 'All':
          continue
        bin_num, d_max, d_min, value = data
        x.append((d_max+d_min)/2)
        y.append(value)
      ax.plot(1/(np.array(x)**2), y, '-', label = name, color = colors[name])

    def resolution(x, pos):
      if x <= 0:
        return '-'
      return "%.1f"%(1/math.sqrt(x))
    formatter = FuncFormatter(resolution)
    ax1.xaxis.set_major_formatter(formatter)
    ax1.set_xlabel(r'Resolution ${\AA}$')
    ax1.set_ylabel('%')
    ax2.set_ylabel('Multiplicity')
    handles, labels = ax1.get_legend_handles_labels()
    handles.extend(ax2.get_legend_handles_labels()[0])
    labels.extend(ax2.get_legend_handles_labels()[1])
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
    plt.title(title)

    if interactive:
      plt.show()
    else:
      fig.set_size_inches(xsize, ysize)
      fig.savefig("datasets_tmp.png", bbox_inches='tight', dpi=100)
      plt.close(fig)
      return "datasets_tmp.png"

  def plot_many_results(self, all_results, title, xsize=30, ysize=10, interactive = True):
    from matplotlib.ticker import FuncFormatter
    import numpy as np
    import math
    fig, (ax1a, ax2) = plt.subplots(2,1, sharex=True)
    ax1b = ax1a.twinx()

    xvals = []
    overall_cc = []
    overall_mult = []
    overall_comp = []
    cc_cutoff = []
    mult_cutoff = []
    comp_cutoff = []
    for r in all_results:
      if r is None or '# accepted' not in r or r['# accepted'] is None: continue
      name, value = r['# accepted'][-1]
      assert name == 'All'
      #if xvals: assert value >= xvals[-1]
      xvals.append(value)
      for key, array in zip(['CC1/2', 'Multiplicity', 'Completeness'], [overall_cc, overall_mult, overall_comp]):
        if key in r:
          name, value = r[key][-1]
          assert name == 'All'
          array.append(value)
        else:
          array.append(0)

      if 'CC1/2' in r:
        last = last_res = 0
        for row in r['CC1/2']:
          if 'All' in row:
            break
          row_n, d_max, d_min, cc = row
          if cc > 0 and not last:
            last = cc; last_res = (d_max + d_min) / 2
          elif cc > 0 and cc <= last:
            last = cc; last_res = (d_max + d_min) / 2
          else:
            break
        cc_cutoff.append(last_res)
      else:
        cc_cutoff.append(0)

      if 'Multiplicity' in r:
        last = last_res = 0
        for row in r['Multiplicity']:
          if row[0]=='All':
            break
          row_n, d_max, d_min, mult = row
          if mult >= 10:
            last = mult; last_res = (d_max + d_min) / 2
          else:
            break
        mult_cutoff.append(last_res)
      else:
        mult_cutoff.append(0)

      if 'Completeness' in r:
        last = last_res = 0
        for row in r['Completeness']:
          if row[0]=='All':
            break
          row_n, d_max, d_min, comp = row
          if comp >= 90:
            last = comp; last_res = (d_max + d_min) / 2
          else:
            break
        comp_cutoff.append(last_res)
      else:
        comp_cutoff.append(0)


    ax1a.plot(xvals, overall_cc, 'o-', color='blue')
    ax1b.plot(xvals, overall_mult, 'o-', color='red')
    ax2.plot(xvals, 1/(np.array(cc_cutoff)**2), 'o-', color='blue')
    ax2.plot(xvals, 1/(np.array(mult_cutoff)**2), 'o-', color='red')
    ax2.plot(xvals, 1/(np.array(comp_cutoff)**2), 'o-', color='green')

    ax2.legend(["CC1/2 (monotonic)", "Multiplicity (10x)", "Completeness (90%)"])

    ax2.set_xlabel("N images")
    ax1a.set_ylabel("Overall CC1/2 (%)")
    ax1b.set_ylabel("Overall multiplicity")
    ax2.set_ylabel(r"Resolution ($\AA$")
    ax1a.set_title(title)

    def resolution(y, pos):
      if y <= 0:
        return '-'
      return "%.2f"%(1/math.sqrt(y))
    formatter = FuncFormatter(resolution)
    ax2.yaxis.set_major_formatter(formatter)

    #ax1a.set_xlabel

    if interactive:
      plt.show()
    else:
      fig.set_size_inches(xsize, ysize)
      fig.savefig("datasets_tmp.png", bbox_inches='tight', dpi=100)
      plt.close(fig)
      return "datasets_tmp.png"

if __name__ == "__main__":
  import sys
  args = sys.argv[1:]
  if len(args) > 1:
    all_results = []
    for folder in args:
      scraper = Scraper(folder, '#')
      all_results.append(scraper.scrape())
    scraper.plot_many_results(all_results, sys.argv[1])

    xvals = []
    all_y = {}
    for r in all_results:
      if r is None or '# accepted' not in r or r['# accepted'] is None: continue
      name, value = r['# accepted'][-1]
      assert name == 'All'
      xvals.append(value)
      key = 'Merged I/sigI'
      if key in r:
        for row in r[key]:
          if 'All' not in row:
            row_n, d_max, d_min, isigi = row
            d = (d_max + d_min) / 2
            if d not in all_y:
              all_y[d] = []
            all_y[d].append(isigi)
    for d in all_y:
      plt.plot(xvals, all_y[d])
      if d <= 5:
        plt.text(xvals[-1], all_y[d][-1], "%.2f"%d)
      else:
        plt.text(xvals[-1], all_y[d][-1], "-")
    plt.xlabel('N images')
    plt.ylabel('Merged I/sigI')
    plt.plot([xvals[0], xvals[-1]], [2.0, 2.0], 'k:')
    plt.show()
  else:
    scraper = Scraper(args[0], '%')
    results = scraper.scrape()
    for t in results:
      for j in results[t]:
        print (t, j)

    scraper.plot_single_results(results, sys.argv[1])
