from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import sgtbx
import iotbx.command_line.lattice_symmetry
from cctbx.web import cgi_utils
import sys

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("ucparams", "1 1 1 90 90 90"),
     ("sgsymbol", "P1"),
     ("convention", ""),
     ("max_delta", "5")))
  return inp

def run(server_info, inp, status):
  sys.stdout.write("<pre>")
  z = inp.sgsymbol.strip().upper()
  if (z in ("P","A","B","C","I","R","F")):
    inp.sgsymbol = "Hall: %s 1" % z
    inp.convention = ""
  input_symmetry = crystal.symmetry(
    unit_cell=inp.ucparams,
    space_group_info=sgtbx.space_group_info(
      symbol=inp.sgsymbol,
      table_id=inp.convention))
  max_delta = float(inp.max_delta)
  Groups = iotbx.command_line.lattice_symmetry.metric_subgroups(
    input_symmetry=input_symmetry,
    max_delta=max_delta)
  Groups.show()
  sys.stdout.write("</pre>")
