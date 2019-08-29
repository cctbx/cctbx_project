# LIBTBX_SET_DISPATCHER_NAME cluster.visualize_orientations
from __future__ import absolute_import, division, print_function
__author__ = 'zeldin'
#from libtbx.utils import multi_out


def run(_args):
  if _args < 2:
    raise IOError("Must give at least one path to folder of pickles")
  import logging
  from xfel.clustering.cluster import Cluster
  FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
  logging.basicConfig(level=logging.WARNING, format=FORMAT)

  cluster = Cluster.from_directories(_args.folders,
                                          'Command line visualisation')
  logging.info("data imported")
  cluster.visualise_orientational_distribution()

if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser(description=('''Visualise the orientational
  distribution of a set of integration pickles'''))
  parser.add_argument('folders', type=str, nargs='+',
                      help='One or more folers containing integration pickles.')
  args = parser.parse_args()
  result = run(args)

