# LIBTBX_SET_DISPATCHER_NAME cxi.find_target_uc
from __future__ import division
from libtbx.utils import multi_out
import sys

def run(args):
  if args < 2: raise IOError("Must give at least one path to folder of pickles") 
  from cctbx.uctbx.determine_unit_cell.target_uc import target
  ucs = target(args.folders) 
  ucs.cluster(args.t, args.m, args.l) 
  ucs.plot_clusters(ucs.clusters) 

if (__name__ == "__main__"):
  import sys
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
  args = parser.parse_args()
  result = run(args)
 
