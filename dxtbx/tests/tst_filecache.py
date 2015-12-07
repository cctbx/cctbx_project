#!/usr/bin/env python
#
#   Copyright (C) 2015 Diamond Light Source, Markus Gerstel
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Tests for the filecache classes. This compares behaviour against StringIO.

from __future__ import division

def test_filecache():
  import dxtbx.filecache
  import libtbx.load_env
  import os
  from StringIO import StringIO

  dxtbx_dir = libtbx.env.dist_path('dxtbx')
  image = os.path.join(dxtbx_dir, 'tests', 'phi_scan_001.cbf')

  with open(image, 'rb') as fh:
    correct_data = fh.read()

  # Create a caching object
  cache = dxtbx.filecache.lazy_file_cache(open(image, 'rb'))

  # read 100 bytes
  sh = StringIO(correct_data)
  with cache.open() as fh:
    actual   = fh.read(100)
    expected = sh.read(100)
    assert(actual == expected)
    actual   = fh.read(0)
    expected = sh.read(0)
    assert(actual == expected)
    actual   = fh.read(5000)
    expected = sh.read(5000)
    assert(actual == expected)

  # readlines
  sh = StringIO(correct_data)
  with cache.open() as fh:
    actual   = fh.readlines()
    expected = sh.readlines()
    assert(actual == expected)

  # 5x readline
  sh = StringIO(correct_data)
  with cache.open() as fh:
    actual   = [ fh.readline() for n in range(5) ]
    expected = [ sh.readline() for n in range(5) ]
    assert(actual == expected)

  # Get a new cache object
  cache.close()
  cache = dxtbx.filecache.lazy_file_cache(open(image, 'rb'))

  sh = StringIO(correct_data)
  fh = dxtbx.filecache.pseudo_file(cache)
  import random

  random_a, random_b = random.randint(0, 10000), random.randint(0,150000)
  print "Running test for parameters %d %d" % (random_a, random_b)

  actual   = (fh.readline(), \
              fh.read(random_a), \
              fh.readline(), \
              fh.read(random_b))
  expected = (sh.readline(), \
              sh.read(random_a), \
              sh.readline(), \
              sh.read(random_b))
  assert(actual == expected)

if __name__ == '__main__':
  test_filecache()
