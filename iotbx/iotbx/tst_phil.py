import iotbx.phil
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
  print "OK"

if (__name__ == "__main__"):
  exercise()
