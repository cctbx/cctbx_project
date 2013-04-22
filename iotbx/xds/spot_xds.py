#
# spot_xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

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
                self.centroid.append(map(float, tokens[0:3]))
                self.intensity.append(float(tokens[3]))
            if len(tokens) >= 7:
                self.miller_index.append(map(int, tokens[4:7]))


if __name__ == '__main__':
  import sys

  handle = reader()
  handle.read_file(sys.argv[1])

  for (xc, yc, zc), i, (h, k, l) in zip(handle.centroid,
                                        handle.intensity,
                                        handle.miller_index):
    print (xc, yc, zc), i, (h, k, l)
