# LIBTBX_SET_DISPATCHER_NAME cluster.intensity_statistics
from __future__ import division
__author__ = 'zeldin'

import logging
from xfel.clustering.cluster import Cluster
from xfel.clustering.cluster_groups import unit_cell_info
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)

def run(_args):
  if _args < 2:
    raise IOError("Must give at least one path to folder of pickles")

  ucs = Cluster.from_directories(_args.folders, "cluster_42")
  logging.info("Data imported.")

  fig = plt.figure(figsize=(17, 7))

  errors = [i.wilson_err['Standard Error'] for i in ucs.members]
  plt.subplot(2,3,1)
  plt.hist(errors, 50, range=[0,200])
  plt.title("Distribution of Standard Errors on the Wilson fit")

  plt.subplot(2,3,2)
  rs = [-1 * i.minus_2B / 2 for i in ucs.members]
  plt.hist(rs, 50, range=[-50,200])
  plt.title("Distribution of B values for the Wilson plot")

  plt.subplot(2,3,3)
  plt.plot([i.G for i in ucs.members],
           [-1 * i.minus_2B / 2 for i in ucs.members], 'x')
  plt.xlabel("G")
  plt.ylabel("B")
  plt.title("G and B for all members")

  ax_allb = plt.subplot(2,1,2)
  for frame in ucs.members:
    plt.plot([0, -1 * frame.G / frame.minus_2B], [frame.G, 0], 'y--', lw=1)
  plt.xlim([0, max([-1 * frame.G / frame.minus_2B for frame in ucs.members])])
  plt.xlabel("1/d**2")
  plt.ylabel("ln(I)")
  plt.title("All (positive B) fitted lines of G and B on one plot")
  plt.tight_layout()
  plt.show()


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description=('Shows overall unmerged I vs. 1'
                                                '/d**2 statistics'))
  parser.add_argument('folders', type=str, nargs='+',
                      help='One or more folers containing integration pickles.')
  args = parser.parse_args()
  run(args)
