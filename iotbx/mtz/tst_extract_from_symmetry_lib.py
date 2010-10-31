from iotbx.mtz import extract_from_symmetry_lib
from cctbx import sgtbx
from libtbx.utils import format_cpu_times
import libtbx
import sys, os
op = os.path

try:
  import ccp4io_adaptbx
except ImportError:
  ccp4io_adaptbx = None

def exercise_230():
  for space_group_number in xrange(1,231):
    space_group_info = sgtbx.space_group_info(
      number=space_group_number,
      table_id="A1983")
    symbol = extract_from_symmetry_lib.ccp4_symbol(
      space_group_info=space_group_info,
      lib_name="symop.lib")
    if (symbol[0] == "H"):
      symbol = "R" + symbol[1:] + ":H"
    assert sgtbx.space_group_info(
      symbol=symbol,
      table_id="A1983").group() == space_group_info.group()
    #
    symbol = extract_from_symmetry_lib.ccp4_symbol(
      space_group_info=space_group_info,
      lib_name="syminfo.lib")
    if (symbol[0] == "H"):
      symbol = "R" + symbol[1:] + ":H"
    assert sgtbx.space_group_info(
      symbol=symbol,
      table_id="A1983").group() == space_group_info.group()

def exercise_symop_lib_recycling():
  file_iter = open(op.join(
    extract_from_symmetry_lib.ccp4io_lib_data, "symop.lib"))
  ccp4_id_counts = libtbx.dict_with_default_0()
  ccp4_symbol_counts = libtbx.dict_with_default_0()
  for line in file_iter:
    flds = line.split(None, 4)
    ccp4_id = flds[0]
    ccp4_id_counts[ccp4_id] += 1
    space_group_number = int(ccp4_id[-3:])
    order_z = int(flds[1])
    given_ccp4_symbol = flds[3]
    ccp4_symbol_counts[given_ccp4_symbol] += 1
    group = extract_from_symmetry_lib.collect_symops(
      file_iter=file_iter, order_z=order_z)
    assert group.order_z() == order_z
    space_group_info = sgtbx.space_group_info(group=group)
    retrieved_ccp4_symbol = extract_from_symmetry_lib.ccp4_symbol(
      space_group_info=space_group_info,
      lib_name="symop.lib")
    assert retrieved_ccp4_symbol == given_ccp4_symbol
    assert space_group_info.type().number() == space_group_number
    if (1):
      from iotbx.pdb import format_cryst1_sgroup
      sgroup = format_cryst1_sgroup(space_group_info=space_group_info)
      if (len(sgroup) > 11):
        print "ccp4 symop.lib setting leads to pdb CRYST1 overflow:",\
          ccp4_id, given_ccp4_symbol, sgroup
  for ccp4_id,count in ccp4_id_counts.items():
    if (count != 1):
      raise RuntimeError(
        'ccp4 id "%s" appears %d times (should be unique).'
          % (ccp4_id, count))
  ccp4_symbol_counts = libtbx.dict_with_default_0()
  for ccp4_symbol,count in ccp4_symbol_counts.items():
    if (count != 1):
      raise RuntimeError(
        'ccp4 symbol "%s" appears %d times (should be unique).'
          % (ccp4_symbol, count))

def exercise_mmdb_cryst1_interpretation(sgi_hall, pdb_str):
  print >> open("tmp.pdb", "w"), pdb_str
  mgr = ccp4io_adaptbx.mmdb.Manager()
  mgr.ReadPDBASCII(fileName="tmp.pdb", gzipMode=0)
  if (mgr.isSpaceGroup() == 0):
    print "MMDB does not recognize space group symbol:"
    print pdb_str
    print
  else:
    sg = sgtbx.space_group()
    for i in xrange(mgr.GetNumberOfSymOps()):
      s = mgr.GetSymOp(Nop=i)
      sg.expand_smx(s)
    sgi = sgtbx.space_group_info(group=sg)
    if (sgi.group() != sgi_hall.group()):
      print "MMDB symmetry mismatch:"
      print pdb_str
      print sgi_hall
      print sgi
      print

def exercise_syminfo_lib_pdb_cryst1_recycling():
  # this call is to build _syminfo_lib_cache
  assert extract_from_symmetry_lib.ccp4_symbol(
    space_group_info=sgtbx.space_group_info("P 1"),
    lib_name="syminfo.lib") == "P 1"
  #
  import iotbx.pdb.cryst1_interpretation
  n_need_more_special = 0
  for number in xrange(1,230+1):
    for hall,ccp4_symbol in \
          extract_from_symmetry_lib._syminfo_lib_cache[number]:
      sgi_hall = sgtbx.space_group_info(symbol="Hall: "+hall)
      if (sgi_hall.group().is_centric()):
        continue
      sgroup = iotbx.pdb.format_cryst1_sgroup(space_group_info=sgi_hall)
      if (len(sgroup) > 11):
        print "ccp4 syminfo.lib setting leads to pdb CRYST1 overflow:",\
          ccp4_symbol, sgroup
      cs = sgi_hall.any_compatible_crystal_symmetry(volume=1000)
      pdb_str = iotbx.pdb.format_cryst1_record(crystal_symmetry=cs)
      cs2 = iotbx.pdb.cryst1_interpretation.crystal_symmetry(
        cryst1_record=pdb_str)
      if (cs2.space_group_info() is None):
        if (n_need_more_special == 0): print
        print '"%s": "Hall: %s",' % (
          ccp4_symbol.replace(" ", "").upper(), hall)
        n_need_more_special += 1
      else:
        assert cs2.is_similar_symmetry(other=cs)
      if (ccp4io_adaptbx is not None):
        exercise_mmdb_cryst1_interpretation(
          sgi_hall=sgi_hall, pdb_str=pdb_str)
  if (n_need_more_special != 0):
    print
    from libtbx.utils import plural_s
    raise RuntimeError("""\
Please edit iotbx/pdb/cryst1_interpretation.py:
  Add the %d line%s above with "Hall:" to the "special" dictionary.
""" % plural_s(n_need_more_special))

def exercise(args):
  assert len(args) == 0
  if (extract_from_symmetry_lib.ccp4io_lib_data is None):
    print "Skipping iotbx/mtz/tst_extract_from_symmetry_lib.py:" \
      " ccp4io not available"
    return
  exercise_230()
  exercise_symop_lib_recycling()
  exercise_syminfo_lib_pdb_cryst1_recycling()
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(args=sys.argv[1:])
