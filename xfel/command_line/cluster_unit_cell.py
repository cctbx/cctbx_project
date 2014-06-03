# LIBTBX_SET_DISPATCHER_NAME cluster.unit_cell
from __future__ import division
import logging
from xfel.clustering.cluster import Cluster
from xfel.clustering.cluster_groups import unit_cell_info
import matplotlib.pyplot as plt

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)

def run(_args):
  if _args < 2:
    raise IOError("Must give at least one path to folder of pickles")

  ucs = Cluster.from_directories(_args.folders, "cxi_targt_uc")

  fig = plt.figure("Andrews-Bernstein distance dendogram", figsize=(12,8))
  ax = plt.gca()
  clusters, cluster_axes = ucs.ab_cluster(_args.t, log=_args.log, ax=ax)

  print unit_cell_info(clusters)
  plt.tight_layout()
  plt.show()

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description=('Find the best target cell from'
                                                'a set of indexing pickles.'))
  parser.add_argument('folders', type=str, nargs='+',
                      help='One or more folers containing integration pickles.')
  parser.add_argument('-t', type=float, default=5000,
                      help='threshold value for the clustering. Default = 5000')
  parser.add_argument('-o', type=str, default='clustering',
                      help='output file name for unit cells.')
  parser.add_argument('--noplot', action='store_true',
                      help="Do not display plots")
  parser.add_argument('--log', action='store_true',
                    help="Display the dendrogram with a log scale")
  args = parser.parse_args()
  run(args)
#  parser.add_argument('-m', type=str, default='distance',
#                      help='Clustering method for numpy clustering.')
#  parser.add_argument('-l', type=str, default='single',
#                      help='Linkage method for clustering. Default single.')
