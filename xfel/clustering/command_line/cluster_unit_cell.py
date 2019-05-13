# LIBTBX_SET_DISPATCHER_NAME cluster.unit_cell
from __future__ import division
from __future__ import print_function
import logging
from xfel.clustering.cluster import Cluster
from xfel.clustering.cluster_groups import unit_cell_info
import matplotlib.pyplot as plt

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)

def run(_args):
  if _args < 2:
    raise IOError("Must provide location(s) of pickles")
  if _args.paths:
    ucs = Cluster.from_files(raw_input=_args.dirs, n_images=_args.n, dials=_args.dials, json=_args.json)
  elif _args.text:
    assert len(_args.dirs)==1 # one file
    one_file = _args.dirs[0]
    ucs = Cluster.from_list(one_file)
  else:
    ucs = Cluster.from_directories(_args.dirs, n_images=_args.n, dials=_args.dials)

  if not _args.noplot:
    clusters, _ = ucs.ab_cluster(_args.t, log=_args.log,
                               write_file_lists=_args.nofiles,
                               schnell=_args.schnell,
                               doplot=_args.noplot)
    print(unit_cell_info(clusters))
  else:
    plt.figure("Andrews-Bernstein distance dendogram", figsize=(12, 8))
    ax = plt.gca()
    clusters, cluster_axes = ucs.ab_cluster(_args.t, log=_args.log, ax=ax,
                                            write_file_lists=_args.nofiles,
                                            schnell=_args.schnell,
                                            doplot=_args.noplot)
    print(unit_cell_info(clusters))
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description=('Find the best target cell from'
                                                'a set of indexing pickles.'))
  parser.add_argument('dirs', type=str, nargs='+',
                      help='One or more paths to directories, integration pickles or DIALS experiment lists.')
  parser.add_argument('--text', action='store_true',
                      help='Interpret the argument as a text file containing unit cells and space group type.')
  parser.add_argument('--paths', action='store_true',
                      help='Interpret the arguments as complete paths to pickles or tarred pickles, not directories.')
  parser.add_argument('--dials', action='store_true',
                      help='Interpret the arguments as DIALS-format pickles and jsons.')
  parser.add_argument('--json', action='store_true',
                      help='Interpret the arguments as DIALS-format jsons only, no reflection tables.'
                           "specifically implemented for unit cell clustering without reflections.")
  parser.add_argument('-t', type=float, default=5000,
                      help='threshold value for the clustering. Default = 5000')
  parser.add_argument('--noplot', action='store_false',
                      help="Do not display plots")
  parser.add_argument('--log', action='store_true',
                    help="Display the dendrogram with a log scale")
  parser.add_argument('--schnell', action='store_true',
                    help="Use euclidian distance only for increased speed."\
                    "Risky!")
  parser.add_argument('--nofiles', action='store_false',
                      help="Write files with lists of the images making up "
                           "each cluster")
  parser.add_argument('-n', type=int, default=None, help="Maximum number of files to use")
  args = parser.parse_args()
  run(args)
  #import cProfile
  #cProfile.run('run(args)')
#  parser.add_argument('-m', type=str, default='distance',
#                      help='Clustering method for numpy clustering.')
#  parser.add_argument('-l', type=str, default='single',
#                      help='Linkage method for clustering. Default single.')
