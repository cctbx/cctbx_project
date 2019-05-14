# LIBTBX_SET_DISPATCHER_NAME iota.double_gauss_visualizer
from __future__ import absolute_import, division, print_function
import logging
import prime.iota.double_gauss as dg

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)

def run(args):
  if args < 2:
    raise IOError("Must specify an integration pickle image")
  img_data = dg.get_img_data(args.image)
  spot_man  = dg.SpotManager(img_data, args.t, args.b, args.s,
                            args.threshold_range, args.background_range,
                            args.smooth_range)
  spot_man.plot_spots()

if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser(description=('Visually optimize double-Gauss spotfinder'))
  parser.add_argument('image', type=str,
                      help='Path to an image pickle.')
  parser.add_argument('-t', type=float, default=1.5,
                      help='Initial threshold value for spotfinding.')
  parser.add_argument('-b', type=float, default=10,
                      help='Initial threshold value for background Gaussian kernel.')
  parser.add_argument('-s', type=float, default=1.5,
                      help='Initial threshold value for spot Gaussian kernel.')
  parser.add_argument('--threshold-range', type=float, nargs=2, default=(1, 4),
                      help="Range for threshold")
  parser.add_argument('--background-range', type=float, nargs=2, default=(0,100),
                      help="Range for background slider")
  parser.add_argument('--smooth-range', type=float, nargs=2, default=(0, 15),
                      help="Range for smoothing slider")
  args = parser.parse_args()
  run(args)
