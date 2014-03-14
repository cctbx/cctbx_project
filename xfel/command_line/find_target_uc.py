# LIBTBX_SET_DISPATCHER_NAME cxi.find_target_uc
from __future__ import division
#from libtbx.utils import multi_out

def run(args):
  if args < 2: raise IOError("Must give at least one path to folder of pickles")
  from cctbx.uctbx.determine_unit_cell.target_uc import target
  ucs = target(args.folders)
  ucs.cluster(args.t, args.m, args.l)
  if not args.noplot:
    ucs.plot_clusters(ucs.clusters, args.log)


if (__name__ == "__main__"):
  import argparse
  parser = argparse.ArgumentParser(description=
      'Find the best target cell from a collection of indexing pickles.')
  parser.add_argument('folders',  type=str, nargs='+',
                       help='One or more folers containing integration pickles.')
  parser.add_argument('-t', type=float, default=5000,
                       help='threshold value for the clustering. Default = 5000')
  parser.add_argument('-m', type=str, default='distance',
      help='Clustering method for numpy clustering. Options are: inconsistent' +
           ', distance, maxclust, monocri, maxclust_monocrit')
  parser.add_argument('-l', type=str, default='single',
      help='Linkage method for clustering. Default is single. Other options are' +
            'complete, average, weighted, centroid, median, ward.')
  parser.add_argument('--log',  action='store_true',
                      help="Display the dendrogram with a log scale")
  parser.add_argument('--noplot',  action='store_true',
                      help="Do not display plots")
  args = parser.parse_args()
  result = run(args)
