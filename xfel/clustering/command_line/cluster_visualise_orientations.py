# LIBTBX_SET_DISPATCHER_NAME cluster.visualize_orientations
from __future__ import absolute_import, division, print_function
__author__ = 'zeldin'
#from libtbx.utils import multi_out


def run(_args):
  import logging
  from xfel.clustering.cluster import Cluster
  FORMAT = '%(message)s'
  logging.basicConfig(level=logging.WARNING, format=FORMAT)

  if _args.paths:
    cluster = Cluster.from_files(raw_input=_args.folders, dials=_args.dials)
  else:
    cluster = Cluster.from_directories(_args.folders, 'Command line visualisation', dials=_args.dials)

  logging.info("data imported")
  cluster.visualise_orientational_distribution()

if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser(description=('''Visualise the orientational
  distribution of a set of integration pickles'''))
  parser.add_argument('folders', type=str, nargs='+',
                      help='One or more folers containing integration pickles.')
  parser.add_argument('--dials', action='store_true',
                      help='Interpret the arguments as DIALS-format pickles and jsons.')
  parser.add_argument('--paths', action='store_true',
                      help='Interpret the arguments as complete paths to pickles or tarred pickles, not directories.')
  args = parser.parse_args()
  result = run(args)

