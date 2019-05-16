"""
Rolling rotamericity metric for EM-Ringer.

Reference:
  Barad BA, Echols N, Wang RYR, Cheng YC, DiMaio F, Adams PD, Fraser JS. (2015)
  Side-chain-directed model and map validation for 3D Electron Cryomicroscopy.
  Nature Methods, in press.

"""

########################################################################
# Package imports
from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
from collections import defaultdict
import argparse
import os
import sys

import six
from six.moves import range

# from matplotlib import rcParams
# rcParams['figure.autolayout'] = True
# rcParams['xtick.labelsize'] = 16
# rcParams['ytick.labelsize'] = 16
# rcParams['axes.labelsize'] = 24
# rcParams['axes.titlesize'] = 24

Residue_codes = ["ARG","ASN","ASP","CYS","GLU","GLN","HIS",
"LEU","LYS","MET","PHE","SER","TRP","TYR","SEC","PYL"]
Branched_residues = ["THR","VAL","ILE"]
No_c_gamma = ["ALA", "GLY"]
Weird = ["PRO"]


class RingerDict(object):
  '''Ringerdict: A dictionary accessible form of the output of ringer'''
  def __init__(self, resultlist, offset):
    self.dict = {}
    for residue in resultlist:
      if residue.resname in Residue_codes:
        residue.resid = int(residue.resid)+offset
        self.add_residue(residue)

  def add_residue(self, residue):
    if residue.chain_id not in self.dict:
      self.dict[residue.chain_id] = {}
    if 1 in residue._angles:  # TODO: verify this is a dict
      self.dict[residue.chain_id][residue.resid] = residue._angles[1]

  def get_peak(self, chain_id, residue_id):
    if chain_id in self.dict and residue_id in self.dict[chain_id]:
      return self.dict[chain_id][residue_id]
    else:
      return None

  def get_chains(self):
    return list(self.dict.keys())

  def get_residues(self, chain_id):
    return sorted(self.dict[chain_id].keys())

def ranges(p):
    q = sorted(p)
    i = 0
    for j in range(1,len(q)):
        if q[j] > 1+q[j-1]:
            yield (q[i],q[j-1])
            i = j
    yield (q[i], q[-1])

def identify_regions(results,
      thresholded_cutoff=0.8,
      rotamer_cutoff=0.5,
      extension=10,
      out=sys.stdout):
  import numpy as np
  for chain, chain_out in six.iteritems(results):
    outliers = []
    print("For Chain %s:" % chain, file=out)
    for k in chain_out:
      if ((np.divide(k[2],k[1]) > thresholded_cutoff) and
          (np.divide(k[3],k[2]) < rotamer_cutoff)):
        for i in range(k[0]-extension, k[0]+extension):
          outliers.append(i)
    if len(outliers) > 0:
      print(list(ranges(outliers)), file=out)
      print("", file=out)
    else:
      print("No outliers at this threshold \n", file=out)

def make_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)

class main(object):
  def __init__(self,
      ringer_results,
      dir_name=None,
      threshold=0,
      extension=10,
      thresholded_cutoff=0.8,
      rotamer_cutoff=0.5,
      graph=False,
      save=True,
      rel=False,
      out=sys.stdout):
    self.threshold = threshold
    self.extension = extension
    self.threshold_cutoff = thresholded_cutoff
    self.rotamer_cutoff = rotamer_cutoff
    if (dir_name is None) and (save):
      dir_name = os.getcwd()
    hierarchy = RingerDict(ringer_results, 0)
    self.results_a = defaultdict(list)
    for chain in hierarchy.get_chains():
      # Results will be a list of tuples of the form residue number,
      # number checked in window, number passing threshold in window,
      # number deviating in window.
      for i in hierarchy.get_residues(chain):
        total_n = 0.0
        threshold_n = 0.0
        # threshold_deviation = 0
        n_deviate = 0.0
        for j in range(-extension, extension+1):
          chi = hierarchy.get_peak(chain, int(i)+j)
          if chi:
            total_n += 1
            if rel:
              if chi.relrho > threshold:
                threshold_n += 1
                if chi.deviation <= 30:
                  n_deviate += 1
            else:
              if chi.peak_rho > threshold:
                threshold_n += 1
                if chi.deviation <= 30:
                  n_deviate += 1
        self.results_a[chain].append((i, total_n, threshold_n, n_deviate))
    print("====Low-scoring, high-signal regions====", file=out)
    identify_regions(self.results_a, out=out)
    if graph or save:
      plot_results(self.results_a,
        dir_name=dir_name,
        threshold=threshold,
        graph=graph,
        save=save)

  def draw_wx_plot(self, plot, chain_id, threshold=0):
    plot.figure.clear()
    _plot_results_for_chain(
      figure=plot.figure,
      results_a=self.results_a,
      chain_id=chain_id,
      threshold=self.threshold,
      extension=self.extension)

  @property
  def chain_ids(self):
    return sorted(self.results_a.keys())

def _plot_results_for_chain(figure, results_a, chain_id, extension=10,
    threshold=0):
  import numpy as np
  ax = figure.add_subplot(111)
  ax.set_title("Rolling window - Chain %s, Threshold %f" %
    (chain_id, threshold), y=1.05)
  x_a = [k[0] for k in results_a[chain_id]]
  y_a = [np.divide(k[3],k[2]) for k in results_a[chain_id]]
  y_b = [np.divide(k[2],k[1]) for k in results_a[chain_id]]
  ax.plot(x_a, y_a, linewidth=3.0, alpha=0.9, label="Fraction Rotameric (Passing Threshold)")
  ax.plot(x_a, y_b, linewidth=3.0, alpha=0.9, label="Fraction above Threshold")
  ax.set_xlabel("Center Residue of %d-Residue Window" % (2*extension+1), labelpad=10)
  ax.set_ylabel("Fraction of Residues", labelpad=10)
  ax.set_ylim(0,1)
  ax.legend(loc=4)
  # ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
  #     ncol=2, mode="expand", borderaxespad=0., fontsize="small")
  ax.yaxis.set_ticks_position('left') # this one is optional but I still recommend it...
  ax.xaxis.set_ticks_position('bottom')

def plot_results(results_a, dir_name, threshold, graph=False, save=True):
  import matplotlib.pyplot as plt
  for chain in results_a.keys():
    fig = plt.figure(1)
    _plot_results_for_chain(
      figure=fig,
      results_a=results_a,
      chain_id=chain,
      threshold=threshold)
    if graph:
      fig.show()
    if save:
      output = os.path.join(dir_name, chain + "_rolling.png")
      fig.savefig(output)
    plt.close()

def run(args, out=sys.stdout):
  parser = argparse.ArgumentParser()
  parser.add_argument("file",nargs="?")
  parser.add_argument("-o", dest="offset", type=int, default=0)
  parser.add_argument("-t", "--threshold", dest="threshold",
    help='Threshold cutoff for rho density',
    nargs='?', type = float, default=0)
  parser.add_argument("-w", "--extension_around_center", dest = "extension",
    help='Number of amino acids to extend around the center in both directions. \
    The total window will therefore be twice this number plus one for the center.'
    , nargs="?", type=int, default=10)
  parser.add_argument("--percent_passing_cutoff", dest = "thresholded_cutoff",
    help='Minimum %% passing threshold to flag as a bad region...'
    , nargs="?", type=float, default=0.8)
  parser.add_argument("--rotamericity_cutoff", dest = "rotamer_cutoff",
    help='Maximum rotamericity to be flagged.'
    , nargs="?", type=float, default=0.5)
  parser.add_argument("--graph", dest = "graph", action='store_true')
  parser.add_argument("--save", dest = "save", action='store_true')
  parser.add_argument("-r", "--rel", dest = "rel", action='store_true')
  parser.set_defaults(rel=False, graph=False, save=True)
  options = parser.parse_args(args)
  import matplotlib
  matplotlib.use("Agg")
  dir_name = None
  if (options.save):
    dir_name = os.path.splitext(options.file)[0] + ".output"
    make_dir(dir_name)
  result = easy_pickle.load(options.file)
  return main(
    ringer_results=result,
    dir_name=dir_name,
    threshold=options.threshold,
    extension=options.extension,
    thresholded_cutoff=options.thresholded_cutoff,
    rotamer_cutoff=options.rotamer_cutoff,
    graph=options.graph,
    save=options.save,
    rel=options.rel,
    out=out)

if __name__ == "__main__":
  run(sys.argv[1:])
