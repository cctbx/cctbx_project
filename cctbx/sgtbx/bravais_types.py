from __future__ import absolute_import, division, print_function
from six.moves import range
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

class bravais_lattice(object):
  """Reference:  Table 1.2 (page 4) and Table 2.1.1 (page 13) of
  International Tables for Crystallography, Volume A Space Group Symmetry,
  ed. Theo Hahn, Fourth, revised edition, Kluwer Academic Publishers, 1996.
  """

  def __init__(self,
               symbol=None,
               group=None,
               number=None):
    from cctbx import sgtbx
    self.space_group_info = sgtbx.space_group_info(symbol=symbol,
                                                   group=group,
                                                   number=number)
    self.space_group = self.space_group_info.group()
    self.crystal_system = self.space_group.crystal_system()
    self.centring_symbol = self.space_group.conventional_centring_type_symbol()

  crystal_family_symbols = {"Triclinic":'a',
                            "Monoclinic":'m',
                            "Orthorhombic":'o',
                            "Tetragonal":'t',
                            "Trigonal":'h',
                            "Hexagonal":'h',
                            "Cubic":'c'
                            }

  def __str__(self):
    return "%s%s"%(self.crystal_family_symbols[self.crystal_system],
                   self.centring_symbol)

  def __eq__(self, other):
    return str(self) == str(other)

  def __ne__(self, other):
    return str(self) != str(other)

def tst_bravais_types(verbose):
  from cctbx import sgtbx
  from cctbx.sgtbx.subgroups import subgroups

  bravais_types_reference = {"aP":2,"mP":8,"mC":5,"mA":0,"mI":0,"mB":0,"oP":30,"oC":11,
                             "oA":4,"oB":0,"oI":9,"oF":5,"tP":49,"tI":19,"hP":45,"hR":7,
                             "cP":15,"cI":10,"cF":11}

  bravais_types_tally = {}
  for key in bravais_types_reference: bravais_types_tally[key]=0

  crystal_systems_reference = {"Triclinic":2,"Monoclinic":13,"Orthorhombic":59,"Tetragonal":68,
                               "Trigonal":25,"Hexagonal":27,"Cubic":36}
  crystal_systems_tally = {}
  for key in crystal_systems_reference: crystal_systems_tally[key]=0

  for space_group_number in range(1,231):
    space_group_info = sgtbx.space_group_info(number=space_group_number)
    GC = bravais_lattice(number=space_group_number)
    GC_1 = bravais_lattice(group=space_group_info.group())
    GC_2 = bravais_lattice(symbol=str(space_group_info))
    assert GC == GC_1
    assert GC == GC_2
    bravais_types_tally[str(GC)]+=1
    crystal_systems_tally[GC.crystal_system]+=1
    if verbose:
      GC.space_group_info.show_summary()
      print(str(GC), GC.crystal_system)

    # idea, not fully implemented--more extensive testing by generating all subgroups
    if False:
      subgrs = subgroups(parent_group_info).groups_parent_setting()
      for subgroup in subgrs:
        subgroup_info = sgtbx.space_group_info(group=subgroup)
        subgroup_info.show_summary()
        bravais_lattice(subgroup_info.type().lookup_symbol())

  return (bravais_types_tally==bravais_types_reference and
          crystal_systems_tally==crystal_systems_reference)

def exercise():
  from cctbx import sgtbx
  import sys
  for symbol in acentric + centric:
    space_group_info = sgtbx.space_group_info(symbol=symbol)
    assert str(space_group_info) == symbol
    assert space_group_info.is_reference_setting()
  if ("--Verbose" in sys.argv[1:]):
    for symbol in centric:
      print("/* %s */ %d," % (
        symbol, sgtbx.space_group_info(symbol=symbol).type().number()))
  assert tst_bravais_types(verbose=("--Verbose" in sys.argv[1:]))
  print("OK")

if (__name__ == "__main__"):
  exercise()
