from cctbx.examples import all_axes
from cctbx import sgtbx

class empty: pass

def interpret_form_data(form):
  inp = empty()
  for key in (("sgsymbol", "P1"),
              ("convention", "")):
    if (form.has_key(key[0])):
      inp.__dict__[key[0]] = form[key[0]].value.strip()
    else:
      inp.__dict__[key[0]] = key[1]
  return inp

def run(server_info, inp, status):
  print "<pre>"
  space_group_info = sgtbx.space_group_info(
    symbol=inp.sgsymbol,
    table_id=inp.convention)
  all_axes.list_all_axes(space_group_info=space_group_info)
  print "</pre>"
