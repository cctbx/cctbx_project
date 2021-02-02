"""Check if sgtbx Hermann-Mauguin vs. Hall symbol table is fully
   consistent with ITVB 2001.
"""
from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
import os

html_file = "itvb_2001_table_a1427_hall_symbols.html"
html_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), html_file)

def get_and_check_file():
  not_found = False
  if os.path.isfile(html_file):
    with open(html_file) as f:
      not_found = f.read().lower().find("not found") >= 0
  if not os.path.isfile(html_file) or not_found:
    print("Skipping exercise(): input file not available")
    return
  with open(html_file) as f:
    table_lines = f.read().splitlines()[11:-4]
  assert len(table_lines) == 530, "%d != 530" % len(table_lines)
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
    symbols_i = next(space_group_symbol_iterator)
    assert symbols_i.universal_hermann_mauguin() \
        == symbols.universal_hermann_mauguin()

def exercise():
  get_and_check_file()
  print("OK")

if (__name__ == "__main__"):
  exercise()
