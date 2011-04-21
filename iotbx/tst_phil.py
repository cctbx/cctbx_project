import iotbx.phil
from libtbx.test_utils import show_diff
from libtbx import Auto

def exercise():
  master = iotbx.phil.parse(input_string="""\
u=10,12 13 80,90 100
  .type=unit_cell
s=19
  .type=space_group
U=None
  .type=unit_cell
S=None
  .type=space_group
ua=Auto
  .type=unit_cell
sa=Auto
  .type=space_group
""")
  assert not show_diff(master.format(master.extract()).as_str(), """\
u = 10 12 13 80 90 100
s = "P 21 21 21"
U = None
S = None
ua = Auto
sa = Auto
""")
  custom = iotbx.phil.parse(input_string="""\
ua = None
s = Auto
U = 1,2,3
u = Auto
S = 33
sa = None
""")
  assert not show_diff(
    master.fetch(source=custom).extract_format().as_str(), """\
u = Auto
s = Auto
U = 1 2 3 90 90 90
S = "P n a 21"
ua = None
sa = None
""")
  #
  master = iotbx.phil.parse(input_string="""\
s=None
  .type=space_group
  .multiple=True
""")
  custom = iotbx.phil.parse(input_string="""\
s=p212121
s=19
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
s = 19
""")
  custom = iotbx.phil.parse(input_string="""\
s=19
s=p212121
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
s = p212121
""")
  custom = iotbx.phil.parse(input_string="""\
s=18
s=p212121
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
s = 18
s = p212121
""")
  #
  master = iotbx.phil.parse(input_string="""\
sel = None
  .type = atom_selection
  .multiple = True
""")
  clai = master.command_line_argument_interpreter()
  assert not show_diff(clai.process(
    "sel=altloc  ' ' or  name C\\* ").as_str(), """\
sel = altloc ' ' or name C\*
""")
  user = iotbx.phil.parse(input_string="""\
sel = "altloc ' ' or name C\*"
sel = altloc ' ' or name D\*
""")
  work = master.fetch(user)
  assert not show_diff(work.as_str(), """\
sel = "altloc ' ' or name C\\\\*"
sel = altloc ' ' or name D\\*
""")
  ex = work.extract()
  assert len(ex.sel) == 2
  assert ex.sel[0] == "altloc ' ' or name C\*"
  assert ex.sel[1] == "altloc ' ' or name D\*"
  #
  print "OK"

if (__name__ == "__main__"):
  exercise()
