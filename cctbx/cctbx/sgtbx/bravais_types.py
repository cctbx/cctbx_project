acentric = (
"P 1",
"P 1 2 1",
"C 1 2 1",
"P 2 2 2",
"C 2 2 2",
"F 2 2 2",
"I 2 2 2",
"P 4 2 2",
"I 4 2 2",
"P 6 2 2",
"R 3 2 :H",
"P 4 3 2",
"I 4 3 2",
"F 4 3 2",
)

centric = (
"P -1",
"P 1 2/m 1",
"C 1 2/m 1",
"P m m m",
"C m m m",
"F m m m",
"I m m m",
"P 4/m m m",
"I 4/m m m",
"P 6/m m m",
"R -3 m :H",
"P m -3 m",
"I m -3 m",
"F m -3 m",
)

def exercise():
  from cctbx import sgtbx
  for symbol in acentric + centric:
    space_group_info = sgtbx.space_group_info(symbol=symbol)
    assert str(space_group_info) == symbol
    assert space_group_info.is_reference_setting()
  print "OK"

if (__name__ == "__main__"):
  exercise()
