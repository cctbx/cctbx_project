from cctbx import sgtbx
import libtbx.env
import os

try:
  ccp4io_dist = libtbx.env.dist_path("ccp4io")
except KeyError:
  ccp4io_dist = None

def ccp4_symbol(space_group_info):
  symop_lib_paths = [os.path.expandvars("$CCP4_LIB/data/symop.lib")]
  if (ccp4io_dist is not None):
    symop_lib_paths.append(os.path.normpath(os.path.join(
      ccp4io_dist, "lib/data/symop.lib")))
  found_at_least_one_symop_lib = False
  for symop_lib_path in symop_lib_paths:
    if (os.path.isfile(symop_lib_path)):
      found_at_least_one_symop_lib = True
      file_iter = iter(open(symop_lib_path, "r"))
      symbol = search_for_ccp4_symbol(space_group_info, file_iter)
      if (symbol is not None):
        return symbol
  assert found_at_least_one_symop_lib
  return None

def search_for_ccp4_symbol(space_group_info, file_iter):
  given_space_group_number = space_group_info.type().number()
  while 1:
    try: line = file_iter.next()
    except StopIteration: break
    flds = line.split(None, 4)
    space_group_number = int(flds[0])
    order_z = int(flds[1])
    if (space_group_number != given_space_group_number):
      for i in xrange(order_z):
        file_iter.next()
    else:
      space_group_symbol = flds[3]
      group = sgtbx.space_group()
      for i in xrange(order_z):
        line = file_iter.next().strip()
        group.expand_smx(sgtbx.rt_mx(line))
      if (space_group_info.group() == group):
        return space_group_symbol
  return None

def exercise():
  for space_group_number in xrange(1,231):
    space_group_info = sgtbx.space_group_info(
      number=space_group_number,
      table_id="A1983")
    symbol = ccp4_symbol(space_group_info=space_group_info)
    if (symbol[0] == "H"):
      symbol = "R" + symbol[1:] + ":H"
    assert sgtbx.space_group_info(
      symbol=symbol,
      table_id="A1983").group() == space_group_info.group()

if (__name__ == "__main__"):
  exercise()
  print "OK"
