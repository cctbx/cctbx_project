"""Check if sgtbx Hermann-Mauguin vs. Hall symbol table is fully
   consistent with ITVB 2001.
"""
from cctbx import sgtbx
import urllib
import os

html_file = "itvb_2001_table_a1427_hall_symbols.html"

def get_test_files():
  try:
    urllib.urlretrieve("http://cci.lbl.gov/sginfo/%s" % html_file, html_file)
  except IOError:
    return False
  return True

def get_and_check_file():
  if (   not get_test_files()
      or not os.path.isfile(html_file)
      or open(html_file).read().lower().find("not found") >= 0):
    print "Skipping exercise(): input file not available"
    return
  table_lines = open(html_file).read().splitlines()[11:-4]
  assert len(table_lines) == 530
  space_group_symbol_iterator = sgtbx.space_group_symbol_iterator()
  for line in table_lines:
    flds = line.split()
    assert len(flds) == 3
    nc, hm, hall = flds
    assert hall.lower() == hall
    symbols = sgtbx.space_group_symbols(symbol=hm)
    hm_sgtbx = symbols.universal_hermann_mauguin().replace(" ", "_") \
      .replace("_:",":") \
      .replace(":H",":h") \
      .replace(":R",":r")
    hall_sgtbx = symbols.hall().lower().replace(" ", "_")
    if (hall_sgtbx[0] == "_"): hall_sgtbx = hall_sgtbx[1:]
    assert hm_sgtbx == hm
    assert hall_sgtbx == hall
    symbols_i = space_group_symbol_iterator.next()
    assert symbols_i.universal_hermann_mauguin() \
        == symbols.universal_hermann_mauguin()

def exercise():
  get_and_check_file()
  print "OK"

if (__name__ == "__main__"):
  exercise()
