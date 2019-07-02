#
# spot_xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Authors: James Parkhurst, Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip

class reader(object):
  '''Class to read the SPOT.XDS file.'''

  def __init__(self):
    ''' Init the reader.'''
    pass

  def read_file(self, filename):
    '''Read the spot file.'''

    # Setup the list
    self.centroid = []
    self.intensity = []
    self.miller_index = []

    # Open the file
    with open(filename, 'r') as handle:
      for line in handle:
        tokens = line.split()
        if len(tokens) < 4:
          raise IndexError('Not enough tokens')
        if len(tokens) >= 4:
          self.centroid.append(tuple([float(t) for t in tokens[0:3]]))
          self.intensity.append(float(tokens[3]))
        if len(tokens) >= 7:
          self.miller_index.append(tuple([int(t) for t in tokens[4:7]]))

class writer(object):
  '''Class to write the SPOT.XDS file.'''

  def __init__(self, centroids, intensities, miller_indices=None):
    self.centroids = centroids
    self.intensities = intensities
    self.miller_indices = miller_indices

    assert len(self.centroids) == len(self.intensities)
    if self.miller_indices is not None:
      assert len(self.centroids) == len(self.miller_indices)

  def write_file(self, filename=None):
    '''Write the spot file.'''

    with open(filename, 'w') as f:
      for i in range(len(self.centroids)):
        print(" %.2f"*3 %self.centroids[i], end=' ', file=f)
        print("%.2f" %self.intensities[i], end=' ', file=f)
        if self.miller_indices is not None:
          print(" %i"*3 %self.miller_indices[i], end=' ', file=f)
        print("\n", end='', file=f)

if __name__ == '__main__':
  import sys

  handle = reader()
  handle.read_file(sys.argv[1])

  for (xc, yc, zc), i, (h, k, l) in zip(handle.centroid,
                                        handle.intensity,
                                        handle.miller_index):
    print((xc, yc, zc), i, (h, k, l))
