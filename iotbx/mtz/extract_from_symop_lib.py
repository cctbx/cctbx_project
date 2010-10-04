from cctbx import sgtbx
import libtbx.load_env
import os.path as op

if (libtbx.env.has_module("ccp4io")):
  ccp4io_dist = libtbx.env.dist_path("ccp4io")
  ccp4io_symop_lib_path = op.normpath(op.join(
    ccp4io_dist, "lib/data/symop.lib"))
else:
  ccp4io_dist = None

_ccp4_symbol_cache = {}

def ccp4_symbol(space_group_info, require_at_least_one_symop_lib=True):
  lookup_symbol = space_group_info.type().lookup_symbol()
  result = _ccp4_symbol_cache.get(lookup_symbol, "..unknown..")
  if (result != "..unknown.."):
    return result
  result = None
  symop_lib_paths = []
  if (ccp4io_dist is not None):
    symop_lib_paths.append(ccp4io_symop_lib_path)
  symop_lib_paths.append(op.expandvars("$CCP4_LIB/data/symop.lib"))
  found_at_least_one_symop_lib = False
  for symop_lib_path in symop_lib_paths:
    if (op.isfile(symop_lib_path)):
      found_at_least_one_symop_lib = True
      file_iter = open(symop_lib_path)
      result = search_for_ccp4_symbol(space_group_info, file_iter)
      if (result is not None):
        break
  else:
    if (require_at_least_one_symop_lib):
      assert found_at_least_one_symop_lib
  _ccp4_symbol_cache[lookup_symbol] = result
  return result

def search_for_ccp4_symbol(space_group_info, file_iter):
  given_space_group_number = space_group_info.type().number()
  for line in file_iter:
    flds = line.split(None, 4)
    space_group_number = int(flds[0][-3:])
    order_z = int(flds[1])
    if (space_group_number != given_space_group_number):
      for i in xrange(order_z):
        file_iter.next()
    else:
      result = flds[3]
      group = collect_symops(file_iter=file_iter, order_z=order_z)
      if (space_group_info.group() == group):
        return result
  return None

def collect_symops(file_iter, order_z):
  result = sgtbx.space_group()
  for i in xrange(order_z):
    line = file_iter.next().strip()
    result.expand_smx(sgtbx.rt_mx(line))
  return result
