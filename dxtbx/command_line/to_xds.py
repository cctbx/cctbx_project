from __future__ import division
#!/usr/bin/env python
# to_xds.py
#
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Print out the contents of the dxtbx understanding of a bunch of images to
# an example XDS.INP file. This should illustrate the usage of the dxtbx
# classes.

import sys
from dxtbx.serialize import xds

def run(file_names):
  if len(file_names) == 1 and file_names[0].endswith('json'):
    from dxtbx.serialize import load
    sweep = load.imageset(file_names[0])
  else:
    from dxtbx.imageset import ImageSetFactory
    sweep = ImageSetFactory.new(file_names)[0]
  xsx = xds.to_xds(sweep)
  xsx.XDS_INP()


if __name__ == '__main__':
  run(sys.argv[1:])
