from scitbx.python_utils import misc
ext = misc.import_ext("scitbx_boost.array_family.flex_scitbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

def to_list(array):
  """Workaround for C++ exception handling bugs
     (list(array) involves C++ exceptions)"""
  result = []
  for i in xrange(array.size()):
    result.append(array[i])
  return result
