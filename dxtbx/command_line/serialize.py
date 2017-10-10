from __future__ import absolute_import, division
from __future__ import print_function
#!/usr/bin/env python
#
# dxtbx.serialize.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

if __name__ == '__main__':

  from optparse import OptionParser
  from dxtbx.serialize import dump
  from dxtbx.imageset import ImageSetFactory

  # Specify the command line options
  usage  = "usage: %prog [options] /path/to/image/files.ext"
  parser = OptionParser(usage)

  # Add a verbose option (False by default)
  parser.add_option('-o', '--output-file',
                    dest='output_file', type="string",
                    default="imageset.json",
                    help='Enter a destination filename for serialization')

  # Parse the arguments
  (options, args) = parser.parse_args()

  # Print help if no arguments specified, otherwise call spot prediction
  if len(args) == 0:
    print(parser.print_help())

  else:
    imagesets = ImageSetFactory.new(args)
    if len(imagesets) == 0:
      print("Error: no imagesets to serialize.")
    elif len(imagesets) > 1:
      print("Error: more than 1 imageset has been specified")
    else:
      dump.imageset(imagesets[0], options.output_file)
      print("Serialized imageset to {0}".format(options.output_file))
