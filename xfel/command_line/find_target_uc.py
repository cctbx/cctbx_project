# LIBTBX_SET_DISPATCHER_NAME cxi.find_target_uc
from __future__ import division
from libtbx.utils import multi_out
import sys

def run(args):
  if args < 2: raise IOError("Must give at least one path to folder of pickles") 
  from cctbx.uctbx.determine_unit_cell.target_uc import target
  ucs = target(args.folders) 
  ucs.cluster(args.t)
  ucs.plot_clusters(ucs.clusters)
  

if (__name__ == "__main__"):
  import sys
  import argparse
  parser = argparse.ArgumentParser(description='Process some integers.')
  parser.add_argument('folders',  type=str, nargs='+',
                       help='One or more folers containing integration pickles.')
  parser.add_argument('-t', type=float, default=50,
                       help='threshold value for the clustering. Note that this is given in A, but will be squared for clustering, since the clustering occurs in G6 space. Default = 50')

  args = parser.parse_args()
  result = run(args)

