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
  ucs.intensity_statistics()


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description=('Shows overall unmerged I vs. 1'
                                                '/d**2 statistics'))
  parser.add_argument('folders', type=str, nargs='+',
                      help='One or more folers containing integration pickles.')
  args = parser.parse_args()
  run(args)
