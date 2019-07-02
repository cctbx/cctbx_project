# LIBTBX_SET_DISPATCHER_NAME cluster.intensity_statistics
from __future__ import absolute_import, division, print_function
__author__ = 'zeldin'

import logging
from xfel.clustering.cluster import Cluster
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)


def run(_args):
  if _args < 2:
    raise IOError("Must give at least one path to folder of pickles")

  ucs = Cluster.from_directories(_args.folders, "cluster_intensity_stats")
  logging.info("Data imported.")
  plt.figure(figsize=(20,10))
  gs = gridspec.GridSpec(3, 2, width_ratios=[1, 3])
  inten_axes = [plt.subplot(gs[0,0]),
                plt.subplot(gs[1,0]),
                plt.subplot(gs[2,0])]
  big_axes = plt.subplot(gs[:,1])

  ucs.intensity_statistics(ax=inten_axes)
  ucs.all_frames_intensity_stats(ax=big_axes)
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
