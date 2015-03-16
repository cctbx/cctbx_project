#!/usr/bin/env python
#
# install_format.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from libtbx.phil import parse

phil_scope = parse(
'''
  local = True
    .type = bool
    .help = "Install the format locally. If True, the format class will be"
            "installed in the users home directory under .dxtbx/filename.py."
            "If False, the format class will be installed in the system"
            "path under ./cctbx/sources/dxtbx/format/filename.py."

''')

if __name__ == '__main__':
  import sys
  from os.path import expanduser, join, exists
  from os import makedirs
  import libtbx.load_env
  from shutil import copyfile

  # Create the command line interpretro
  interpretor = phil_scope.command_line_argument_interpreter()

  # Process the phil arguments
  good = []
  bad = []
  for arg in sys.argv[1:]:
    try:
      interpretor.process_arg(arg)
      good.append(arg)
    except Exception:
      bad.append(arg)
  processed = interpretor.process(args=good)

  # Get the phil parameters
  phil = phil_scope.fetch(sources=processed)

  # Get the parameters
  params = phil.extract()

  # Make sure we just have the input filenames
  assert(len(bad) == 1)
  filename = bad[0]

  # Get the target directory
  if params.local == True:
    target = join(expanduser('~'), '.dxtbx')
    if not exists(target):
      makedirs(target)
  else:
    dxtbx = libtbx.env.dist_path('dxtbx')
    target = join(dxtbx, "format")
  target = join(target, filename)

  # Copy the file
  print "Copying %s -> %s" % (filename, target)
  copyfile(filename, target)
