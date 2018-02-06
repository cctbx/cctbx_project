# coding: utf-8
from __future__ import absolute_import, division, print_function

import dxtbx.filecache
import dxtbx.filecache_controller as fcc
from mock import Mock, create_autospec


def test_invalid_cache(monkeypatch):
  """Tests a condition where the filecache controller can give out an invalid cache"""

  mocklazy = create_autospec(dxtbx.filecache.lazy_file_cache)
  monkeypatch.setattr(fcc.dxtbx.filecache, "lazy_file_cache", mocklazy)

  # Create the cache
  cache = fcc.simple_controller()

  # Set up an intial, working cache
  good_file_opener = Mock()
  cache.check("working", lambda: good_file_opener)
  # This should have used the lazy_file_cache to set up
  mocklazy.assert_called_with(good_file_opener)
  mocklazy.return_value.open.assert_called()

  # Now, pass the cache an opener that fails
  try:
    badfile = Mock(side_effect=IOError("Testing bad file"))
    cache.check("not_working", badfile)
    assert False, "Failed to raise IOError in file cache"
  except IOError:
    pass

  # The bug: Calling check with the same tag shouldn't use the invalid cache
  # To Test: Do the working test, but with the failed tag. If the invalid
  #          cache is used, then we won't have asked for reconstruction
  mocklazy.reset_mock()
  cache.check("not_working", lambda: good_file_opener)
  mocklazy.assert_called_with(good_file_opener)
  mocklazy.return_value.open.assert_called()
