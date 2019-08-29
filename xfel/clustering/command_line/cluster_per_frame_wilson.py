# LIBTBX_SET_DISPATCHER_NAME cluster.individual_frame_intensity
from __future__ import absolute_import, division, print_function
__author__ = 'zeldin'

import logging
from xfel.clustering.cluster import Cluster
import matplotlib.pyplot as plt

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)

class Key_event:
  def __init__(self, ax, members, fig):
    self.image_index = 0
    self.ax = ax
    self.members = members
    self.fig = fig

  def key_event(self, e):
    if e.key == "right":
      self.image_index += 1
    elif e.key == "left":
      self.image_index -= 1
    else:
      return
    self.ax.cla()
    self.image_index %= len(self.members)
    self.ax = self.members[self.image_index].plot_wilson(ax=self.ax)
    self.fig.canvas.draw()

def run(_args):
  if _args < 2:
    raise IOError("Must give at least one path to folder of pickles")

  ucs = Cluster.from_directories(_args.folders, "Per-frame-Wilson")
  logging.info("Data imported.")


  fig = plt.figure(figsize=(10,10))
  ax = plt.gca()
  ucs.members[0].plot_wilson(ax=ax)

  browser = Key_event(ax, ucs.members, fig)

  fig.canvas.mpl_connect('key_press_event', browser.key_event)
  plt.show()

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description="Create image-by-image Wilson "
               "plots that can be clicked through using left and right arrows.")
  parser.add_argument('folders', type=str, nargs='+',
                      help='One or more folers containing integration pickles.')
  args = parser.parse_args()
  run(args)
