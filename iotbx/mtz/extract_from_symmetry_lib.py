from cctbx import sgtbx
import libtbx.load_env
import os.path as op

if (libtbx.env.has_module("ccp4io")):
  ccp4io_lib_data = libtbx.env.under_dist(
    module_name="ccp4io", path="lib/data")
else:
  ccp4io_lib_data = None

_ccp4_symbol_cache = {"symop.lib": {}, "syminfo.lib": {}}
_syminfo_lib_cache = []

syminfo_lib_bad_old = set("""
P 21/m 21/m 2/n a
""".splitlines())

def ccp4_symbol(space_group_info, lib_name, require_at_least_one_lib=True):
  assert lib_name in _ccp4_symbol_cache.keys()
  sg_type = space_group_info.type()
  lookup_symbol = sg_type.lookup_symbol()
  cache = _ccp4_symbol_cache[lib_name]
  result = cache.get(lookup_symbol, "..unknown..")
  if (result != "..unknown.."):
    return result
  if (lib_name != "syminfo.lib" or len(_syminfo_lib_cache) == 0):
    lib_paths = []
    if (ccp4io_lib_data is not None):
      lib_paths.append(op.join(ccp4io_lib_data, lib_name))
    lib_paths.append(op.expandvars("$CCP4_LIB/data/"+lib_name))
    found_at_least_one_lib = False
    for lib_path in lib_paths:
      if (op.isfile(lib_path)):
        found_at_least_one_lib = True
        if (lib_name == "symop.lib"):
          ccp4_symbol = search_symop_lib_for_ccp4_symbol(
            space_group_info=space_group_info,
            file_iter=open(lib_path))
          if (ccp4_symbol is not None):
            cache[lookup_symbol] = ccp4_symbol
            return ccp4_symbol
        else:
          build_syminfo_lib_cache(lib_path)
          break
    else:
      if (require_at_least_one_lib):
        assert found_at_least_one_lib
  if (lib_name == "syminfo.lib"):
    for hall,ccp4_symbol in _syminfo_lib_cache[sg_type.number()]:
      sgi = sgtbx.space_group_info(symbol="Hall: "+hall)
      lus = sgi.type().lookup_symbol()
      cache[lus] = ccp4_symbol
      if (lus == lookup_symbol):
        return ccp4_symbol
  return None

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

def build_syminfo_lib_cache(lib_path):
  _syminfo_lib_cache.append(None)
  for number in xrange(230):
    _syminfo_lib_cache.append([])
  file_iter = open(lib_path)
  for line in file_iter:
    l = line.strip()
    if (l == "begin_spacegroup"):
      number = None
      symbols = {}
      for line in file_iter:
        l = line.strip()
        if (l == "end_spacegroup"):
          assert number is not None
          assert len(symbols) == 3
          def get_shortest(s_list):
            result = None
            for s in s_list:
              if (len(s) == 0): continue
              if (result is None or len(result) > len(s)):
                result = s
            return result
          ccp4_symbol = get_shortest(symbols["old"])
          if (   ccp4_symbol is None
              or ccp4_symbol in syminfo_lib_bad_old):
            if (len(symbols["xhm"]) != 0):
              ccp4_symbol = symbols["xhm"]
            else:
              raise RuntimeError("Missing both xHM and old symbols")
          _syminfo_lib_cache[number].append((symbols["hall"], ccp4_symbol))
          break
        if (l.startswith("number ")):
          flds = l.split()
          assert len(flds) == 2
          number = int(flds[1])
          assert number >= 1
          assert number <= 230
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
