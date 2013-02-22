#!/usr/bin/env python
# test_beam.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Tests for the beam class.

import os
import sys

from dxtbx.model.beam import beam
from dxtbx.model.beam import beam_factory

def test_beam():
    '''A test class for the beam class.'''

    cbf = beam_factory.imgCIF('phi_scan_001.cbf')

    print cbf

if __name__ == '__main__':

    test_beam()
