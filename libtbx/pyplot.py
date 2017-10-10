from __future__ import division
from builtins import zip
from matplotlib.pyplot import *

def plot_pairs(xy, *a, **k):
  x, y = list(zip(*xy))
  plot(x, y, *a, **k)
