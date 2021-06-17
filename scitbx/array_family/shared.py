from __future__ import absolute_import, division, print_function
import scitbx.stl.vector # import dependency

import boost_adaptbx.boost.python as bp
bp.import_ext("scitbx_array_family_shared_ext")
from scitbx_array_family_shared_ext import *
import scitbx_array_family_shared_ext as ext

class pickle_import_trigger(object): pass

@bp.inject_into(ext.stl_set_unsigned)
class _():

  def __getstate__(self): # XXX slow, move to C++
    version = 2
    list_of_lists = []
    for elem in self:
      list_of_lists.append(list(elem))
    return (version, pickle_import_trigger(), list_of_lists)

  def __setstate__(self, state):
    assert len(state) >= 2
    version = state[0]
    if   (version == 1): assert len(state) == 2
    elif (version == 2): assert len(state) == 3
    else: raise RuntimeError("Unknown version of pickled state.")
    list_of_lists = state[-1]
    for elem in list_of_lists:
      self.append(scitbx.stl.set.unsigned(elem))
