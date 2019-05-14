from __future__ import absolute_import, division, print_function
from cctbx.examples import all_axes
from cctbx import sgtbx
from cctbx.web import cgi_utils

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("sgsymbol", "P1"),
     ("convention", "")))
  return inp

def run(server_info, inp, status):
  print("<pre>")
  space_group_info = sgtbx.space_group_info(
    symbol=inp.sgsymbol,
    table_id=inp.convention)
  all_axes.list_all_axes(space_group_info=space_group_info)
  print("</pre>")
