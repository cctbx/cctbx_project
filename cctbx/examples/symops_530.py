"""\
Loop over 530 conventional settings of the 230 space groups,
show symmetry operations in various formats.

See also:
  List of 530 settings:
    Shmueli U, Hall SR, Grosse-Kunstleve RW:
    Space-Group Symbols for numeric and symbolic computations.
    In International Tables for Crystallography, Volume B:
    Reciprocal space, U. Shmueli, Ed.,
    Kluwer Academic Publishers (Dordrecht), 2001, 107-119.

  Universal Hermann-Mauguin symbols: section 2.1 of:
    Zwart P, Grosse-Kunstleve RW, Lebedev AA, Murshudov GN, Adams PD:
    Surprises and pitfalls arising from(pseudo)symmetry
    Acta Cryst. 2008, D64, 99-107.
    http://scripts.iucr.org/cgi-bin/paper?ba5111
"""
from __future__ import absolute_import, division, print_function

from cctbx import sgtbx

def run():
  for symbols in sgtbx.space_group_symbol_iterator():
    symbol = symbols.universal_hermann_mauguin()
    print(symbol)
    space_group_info = sgtbx.space_group_info(symbol=symbol)
    for s in space_group_info.group():
      print(s.as_xyz())
    for s in space_group_info.group():
      sr = s.as_rational()
      print(sr.r.elems, sr.t.elems)
    for s in space_group_info.group():
      print(s.r().num(), s.r().den(), s.t().num(), s.t().den())
    print()

if (__name__ == "__main__"):
  run()
