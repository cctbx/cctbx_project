# Generate SHELX LATT and SYMM cards for a given space group.

from cctbx import sgtbx
from cctbx.web import cgi_utils
from iotbx.shelx.write_ins import LATT_SYMM
import sys

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("sgsymbol", "P1"),
     ("convention", "")))
  return inp

def run(server_info, inp, status):
  print "<pre>"
  space_group_info = sgtbx.space_group_info(
    symbol=inp.sgsymbol,
    table_id=inp.convention)
  space_group_info.show_summary()
  print
  LATT_SYMM(sys.stdout, space_group_info.group())
  print "</pre>"
