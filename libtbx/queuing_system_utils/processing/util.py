from __future__ import absolute_import, division, print_function
from functools import reduce

def get_regex(pattern, flags):

  import re
  import operator
  return re.compile(
    pattern,
    reduce( operator.or_, [ getattr( re, f ) for f in flags ], 0 ),
    )


def get_lazy_initialized_regex(pattern, flags = []):

  from libtbx.object_oriented_patterns import lazy_initialization

  return lazy_initialization( func = get_regex, pattern = pattern, flags = flags )
