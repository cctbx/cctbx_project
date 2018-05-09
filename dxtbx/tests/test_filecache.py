from __future__ import absolute_import, division, print_function

import os

def test_filecache():
  import dxtbx.filecache
  import libtbx.load_env
  from StringIO import StringIO # this is not cStringIO on purpose!

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
  cache._page_size = 5
  fh = dxtbx.filecache.pseudo_file(cache)

  # readline stress test
  sh = StringIO(correct_data)
  with cache.open() as fh:
    actual   = fh.readline()
    expected = sh.readline()
    assert(actual == expected)

    actual   = fh.read(68)
    expected = sh.read(68)
    assert(actual == expected)

    actual   = fh.readline()
    expected = sh.readline()
    assert(actual == expected)

    actual   = fh.read(1)
    expected = sh.read(1)
    assert(actual == expected)

  # Get a new cache object
  cache.close()
  cache = dxtbx.filecache.lazy_file_cache(open(image, 'rb'))

  sh = StringIO(correct_data)
  fh = dxtbx.filecache.pseudo_file(cache)
  import random

  random_a, random_b = random.randint(0, 10000), random.randint(0,150000)
  print("Running test for parameters %d %d" % (random_a, random_b))

  actual   = (fh.readline(), \
              fh.read(random_a), \
              fh.readline(), \
              fh.read(random_b))
  expected = (sh.readline(), \
              sh.read(random_a), \
              sh.readline(), \
              sh.read(random_b))
  assert(actual == expected)

def test_filecache_more(dials_regression):
  import dxtbx.filecache

  filename = os.path.join(dials_regression, 'image_examples', 'MacScience', 'reallysurprise_001.ipf')
  cache = dxtbx.filecache.lazy_file_cache(open(filename, 'rb'))
  with dxtbx.filecache.pseudo_file(cache) as fh:
    fh.seek(3000 * 3000 * 2)
    fh.read(1024)
    assert fh.read(1) == ''

  cache = dxtbx.filecache.lazy_file_cache(open(filename, 'rb'))
  with dxtbx.filecache.pseudo_file(cache) as fh:
    data = fh.read(3000 * 3000 * 2)
    assert len(data) == 3000 * 3000 * 2
    assert len(fh.read(1024)) == 1024
    assert fh.read(1) == ''

  cache = dxtbx.filecache.lazy_file_cache(open(filename, 'rb'))
  with dxtbx.filecache.pseudo_file(cache) as fh:
    fh.read(1024)
    data = fh.read()
    assert len(data) == 3000 * 3000 * 2
