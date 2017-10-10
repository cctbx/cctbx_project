#!/usr/bin/env python
#
# import.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division
from __future__ import print_function
from builtins import str
from dxtbx.datablock import DataBlockFactory, DataBlockDumper

if __name__ == '__main__':

  from optparse import OptionParser
  usage = "usage: %prog [options] /path/to/image/files"
  parser = OptionParser(usage)

  # Print verbose output
  parser.add_option(
    "-v", "--verbose",
    dest = "verbose",
    action = "count", default = 0,
    help = "Set the verbosity level (-vv gives a verbosity level of 2)")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-o", "--output",
    dest = "output",
    type = "string", default = None,
    help = "The output JSON or pickle file (filename.json | filename.pickle)")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-c", "--compact",
    dest = "compact",
    action = "store_true", default = False,
    help = "For JSON output, use compact representation")

  # Don't sort input filenames
  parser.add_option(
    "-n", "--no-sort",
    dest = "sort",
    action = "store_false", default = True,
    help = "Don't sort input files (default is True)")

  # Parse the command line arguments
  (options, args) = parser.parse_args()
  if len(args) == 0:
    parser.print_help()

  # Sort arguments
  if options.sort:
    args = sorted(args)

  # Get the data blocks from the input files
  # We've set verbose to print out files as they're tested.
  unhandled = []
  datablocks = DataBlockFactory.from_args(args,
    verbose=options.verbose, unhandled=unhandled)

  # Print out any unhandled files
  if len(unhandled) > 0:
    print('-' * 80)
    print('The following command line arguments were not handled:')
    for filename in unhandled:
      print('  %s' % filename)

  # Loop through the data blocks
  for i, datablock in enumerate(datablocks):

    # Extract any sweeps
    sweeps = datablock.extract_sweeps()

    # Extract any stills
    stills = datablock.extract_stills()
    if not stills:
      num_stills = 0
    else:
      num_stills = len(stills)

    # Print some data block info
    print("-" * 80)
    print("DataBlock %d" % i)
    print("  format: %s" % str(datablock.format_class()))
    print("  num images: %d" % datablock.num_images())
    print("  num sweeps: %d" % len(sweeps))
    print("  num stills: %d" % num_stills)

    # Loop through all the sweeps
    if options.verbose > 1:
      for j, sweep in enumerate(sweeps):
        print("")
        print("Sweep %d" % j)
        print("  length %d" % len(sweep))
        print(sweep.get_beam())
        print(sweep.get_goniometer())
        print(sweep.get_detector())
        print(sweep.get_scan())

  # Write the datablock to a JSON or pickle file
  if options.output:
    print("-" * 80)
    print('Writing datablocks to %s' % options.output)
    dump = DataBlockDumper(datablocks)
    dump.as_file(options.output, compact=options.compact)
