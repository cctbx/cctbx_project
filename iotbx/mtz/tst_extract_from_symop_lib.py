from iotbx.mtz import extract_from_symop_lib
from cctbx import sgtbx
from libtbx.utils import format_cpu_times
import sys

def exercise_230():
  for space_group_number in xrange(1,231):
    space_group_info = sgtbx.space_group_info(
      number=space_group_number,
      table_id="A1983")
    symbol = extract_from_symop_lib.ccp4_symbol(
      space_group_info=space_group_info)
    if (symbol[0] == "H"):
      symbol = "R" + symbol[1:] + ":H"
    assert sgtbx.space_group_info(
      symbol=symbol,
      table_id="A1983").group() == space_group_info.group()

def exercise_symop_lib_recycling():
  file_iter = open(extract_from_symop_lib.ccp4io_symop_lib_path)
  for line in file_iter:
    flds = line.split(None, 4)
    ccp4_id = flds[0]
    space_group_number = int(ccp4_id[-3:])
    order_z = int(flds[1])
    given_ccp4_symbol = flds[3]
    group = extract_from_symop_lib.collect_symops(
      file_iter=file_iter, order_z=order_z)
    assert group.order_z() == order_z
    space_group_info = sgtbx.space_group_info(group=group)
    retrieved_ccp4_symbol = extract_from_symop_lib.ccp4_symbol(
      space_group_info=space_group_info)
    assert retrieved_ccp4_symbol == given_ccp4_symbol
    assert space_group_info.type().number() == space_group_number

def exercise(args):
  assert len(args) == 0
  if (extract_from_symop_lib.ccp4io_dist is None):
    print "Skipping iotbx/mtz/tst_extract_from_symop_lib.py:" \
      " ccp4io not available"
    return
  exercise_230()
  exercise_symop_lib_recycling()
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(args=sys.argv[1:])
