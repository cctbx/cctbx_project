from cctbx import sgtbx
import libtbx.load_env
import os.path as op

if (libtbx.env.has_module("ccp4io")):
  ccp4io_lib_data = libtbx.env.under_dist(
    module_name="ccp4io", path="lib/data")
else:
  ccp4io_lib_data = None

_ccp4_symbol_cache = {"symop.lib": {}, "syminfo.lib": {}}

def ccp4_symbol(space_group_info, lib_name, require_at_least_one_lib=True):
  assert lib_name in _ccp4_symbol_cache.keys()
  lookup_symbol = space_group_info.type().lookup_symbol()
  cache = _ccp4_symbol_cache[lib_name]
  result = cache.get(lookup_symbol, "..unknown..")
  if (result != "..unknown.."):
    return result
  result = None
  lib_paths = []
  if (ccp4io_lib_data is not None):
    lib_paths.append(op.join(ccp4io_lib_data, lib_name))
  lib_paths.append(op.expandvars("$CCP4_LIB/data/"+lib_name))
  found_at_least_one_lib = False
  for lib_path in lib_paths:
    if (op.isfile(lib_path)):
      found_at_least_one_lib = True
      file_iter = open(lib_path)
      if (lib_name == "symop.lib"):
        result = search_symop_lib_for_ccp4_symbol(
          space_group_info=space_group_info, file_iter=file_iter)
      else:
        result = search_syminfo_lib_for_ccp4_symbol(
          space_group_info=space_group_info, file_iter=file_iter)
      if (result is not None):
        break
  else:
    if (require_at_least_one_lib):
      assert found_at_least_one_lib
  cache[lookup_symbol] = result
  return result

def search_symop_lib_for_ccp4_symbol(space_group_info, file_iter):
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

def search_syminfo_lib_for_ccp4_symbol(space_group_info, file_iter):
  given_space_group_number = space_group_info.type().number()
  for line in file_iter:
    l = line.strip()
    if (l == "begin_spacegroup"):
      symbols = {}
      for line in file_iter:
        l = line.strip()
        if (l == "end_spacegroup"):
          assert len(symbols) == 3
          group = sgtbx.space_group(symbols["hall"])
          if (group == space_group_info.group()):
            def get_shortest(s_list):
              result = None
              for s in s_list:
                if (len(s) == 0): continue
                if (result is None or len(result) > len(s)):
                  result = s
              return result
            result = get_shortest(symbols["old"])
            if (result is None):
              if (len(symbols["xhm"]) != 0):
                result = symbols["xhm"]
              else:
                raise RuntimeError("Missing both xHM and old symbols")
            return result
          break
        if (l.startswith("number ")):
          flds = l.split()
          assert len(flds) == 2
          number = int(flds[1])
          if (number != given_space_group_number):
            break
        elif (l.startswith("symbol ")):
          flds = l.split(None, 2)
          assert len(flds) == 3
          stype = flds[1].lower()
          if (stype in ["hall", "xhm", "old"]):
            assert stype not in symbols
            symbol = flds[2].strip()
            assert len(symbol) >= 2
            assert symbol.startswith("'")
            assert symbol.endswith("'")
            if (stype == "old"):
              symbols[stype] = " ".join(symbol[1:-1].split()).split("' '")
            else:
              symbols[stype] = symbol[1:-1]
      else:
        raise RuntimeError("Missing end_spacegroup")
  return None
