from __future__ import absolute_import, division, print_function
from matplotlib.pyplot import *

def plot_pairs(xy, *a, **k):
  x, y = zip(*xy)
  plot(x, y, *a, **k)
