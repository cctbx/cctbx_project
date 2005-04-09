import iotbx.phil
from libtbx.test_utils import show_diff
from cStringIO import StringIO

def exercise():
  parameters = iotbx.phil.parse(input_string="""\
u=10,12 13 80,90 100
  .type=unit_cell
s=19
  .type=space_group
U=none
  .type=unit_cell
S=none
  .type=space_group
""")
  out = StringIO()
  extracted = parameters.fetch(source=parameters).extract()
  parameters.format(extracted).show(out=out)
  assert out.getvalue() == """\
u = 10 12 13 80 90 100
s = "P 21 21 21"
U = None
S = None
"""
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
  print "OK"

if (__name__ == "__main__"):
  exercise()
