from __future__ import division
#!/usr/bin/env python
# TestFormat.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Code to give the Format static methods a workout.

import sys

from Toolkit.ImageFormat.Format import Format
from Toolkit.ImageFormat.FormatCBF import FormatCBF

def workout(args):
    for arg in args:
        if Format.is_bz2(arg):
            print 'bzip2 %s' % arg
        if Format.is_gzip(arg):
            print 'gzip %s' % arg

    if FormatCBF.understand(arg):
        print FormatCBF.get_cbf_header(arg)

if __name__ == '__main__':
    workout(sys.argv[1:])
