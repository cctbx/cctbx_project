from __future__ import absolute_import, division, print_function
import os, glob

"""
Searches the cctbx.xfel.merge log files for statistics tables.

Makes lots of assumptions regarding formatting so this is a work in progres.
"""

# Name of stat: (Line signal, offset to first line of table, index of value for a bin line, index of value for an 'All' line)
types = {
  "% accepted": ("Lattices resolution", 6, 4, 1),
  "Multiplicity": ("Intensity Statistics (all accepted experiments)", 14, 6, 3),
  "CC1/2": ("Table of Scaling Results", 6, 5, 2)
}

class Scraper(object):
  def __init__(self, output_path):
    self.output_path = output_path

  def scrape(self):
    path = glob.glob(os.path.join(self.output_path, "*main*"))[0]
    results = {}
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
    self._num_to_percent(results, '% accepted')
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

  def plot_single_results(self, fig, results):
    from matplotlib.ticker import FuncFormatter
    import numpy as np
    import math
    ax1 = fig.gca()
    ax2 = ax1.twinx()

    for name, c in zip(results, ['b', 'r', 'g']):
      if name == 'Multiplicity':
        ax = ax2
      else:
        ax = ax1

      x = []; y = []
      #x = np.array([]); y = np.array([])
      for data in results[name]:
        if data[0] == 'All':
          continue
        bin_num, d_max, d_min, value = data
        x.append((d_max+d_min)/2)
        y.append(value)
      ax.plot(1/(np.array(x)**2), y, '-', label = name, color = c)

    def resolution(x, pos):
      if x <= 0:
        return '-'
      return "%.1f"%(1/math.sqrt(x))
    formatter = FuncFormatter(resolution)
    ax1.xaxis.set_major_formatter(formatter)
    ax1.set_xlabel(u'Resolution ${\AA}$')
    ax1.set_ylabel('%')
    ax2.set_ylabel('Multiplicity')
    fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)

  def plot_many_results(self, all_results):
    pass

if __name__ == "__main__":
  from matplotlib import pyplot as plt
  import sys
  fig = plt.figure()

  args = sys.argv[1:]
  if len(args) > 1:
    all_results = []
    for folder in args:
      scraper = Scraper(arg])
      all_results.append(scraper.scrape)
    scraper.plot_many_results(all_results)
  else:
    results = scraper.scrape()
    for t in results:
      for j in results[t]:
        print (t, j)

    scraper.plot_single_results(fig, results)
  plt.show()
