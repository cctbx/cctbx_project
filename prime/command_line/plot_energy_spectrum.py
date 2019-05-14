# LIBTBX_SET_DISPATCHER_NAME prime.plot_energy_spectrum
"""
Author      : Uervirojnangkoorn, M.
Created     : 5/16/2016
Description : read a single energy spectrum file (csv) and plot the spectrum
"""
from __future__ import absolute_import, division, print_function
import sys
from cctbx.array_family import flex
import matplotlib.pyplot as plt

if (__name__ == "__main__"):
  if len(sys.argv)==1:
    print('Usage: prime.plot_energy_spectrum energy_as_comma_separated.file')
    exit()
  energy_file = sys.argv[1]

  pf = open(energy_file,'r')
  data = pf.read().split('\n')
  x = flex.double()
  y = flex.double()
  for row in data:
    cols = row.split(',')
    if len(cols) == 2:
      x.append(float(cols[1]))
      y.append(float(cols[0]))
  plt.plot(x,y)
  plt.show()
