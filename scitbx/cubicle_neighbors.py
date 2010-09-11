import scitbx.array_family.flex
import scitbx.stl.map

import boost.python
ext = boost.python.import_ext("scitbx_cubicle_neighbors_ext")
from scitbx_cubicle_neighbors_ext import *

def __cubicles_max_memory_allocation_set():
  from libtbx.introspection import machine_memory_info
  m = machine_memory_info().memory_total()
  if (m is None): m = 1000000000
  l = 2**(8*boost.python.sizeof_void_ptr-1)-1
  if (m > l): m = l
  cubicles_max_memory_allocation_set(number_of_bytes=m//2)
__cubicles_max_memory_allocation_set()
