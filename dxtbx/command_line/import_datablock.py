from dxtbx.datablock import DataBlockFactory

if __name__ == '__main__':

  from optparse import OptionParser
  usage = "usage: %prog [options] /path/to/image/files"
  parser = OptionParser(usage)

  # Print verbose output
  parser.add_option(
    "-v", "--verbose",
    dest = "verbose",
    action = "store_true", default = False,
    help = "Write extra output")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-o", "--output",
    dest = "output",
    type = "string", default = None,
    help = "The output JSON or pickle file (filename.json | filename.pickle)")

  # Parse the command line arguments
  (options, args) = parser.parse_args()
  if len(args) == 0:
    parser.print_help()

  # Get the data blocks from the input files
  # We've set verbose to print out files as they're tested.
  datablocks = DataBlockFactory.from_filenames(args, verbose=options.verbose)

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
    print "-" * 80
    print "DataBlock %d" % i
    print "  format: %s" % str(datablock.format_class())
    print "  num images: %d" % len(datablock)
    print "  num sweeps: %d" % len(sweeps)
    print "  num stills: %d" % num_stills

    # Loop through all the sweeps
    if options.verbose:
      for j, sweep in enumerate(sweeps):
        print ""
        print "Sweep %d" % j
        print "  length %d" % len(sweep)
        print sweep.get_beam()
        print sweep.get_goniometer()
        print sweep.get_detector()
        print sweep.get_scan()

  # Write the datablock to a JSON or pickle file
  if options.output:
    print "-" * 80
    print 'Writing datablocks to %s' % options.output
    import os
    import json
    import cPickle as pickle
    ext = os.path.splitext(options.output)[1]
    if ext == '.json':
      dictionary = [db.to_dict() for db in datablocks]
      json.dump(dictionary, open(options.output, "w"),
        indent=2, ensure_ascii=True)
    elif ext == '.pickle':
      pickle.dump(datablocks, open(options.output, "w"),
        protocol=pickle.HIGHEST_PROTOCOL)
    else:
      raise RuntimeError('expected extension .json or .pickle, got %s' % ext)
