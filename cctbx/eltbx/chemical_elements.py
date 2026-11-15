"""
This module can be used to retrieve a list or set of all chemical symbols on the
periodic table, from hydrogen up to roentgenium, and the set of the metals among
them.
"""
from __future__ import absolute_import, division, print_function
import scitbx.stl.set # import dependency

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_eltbx_chemical_elements_ext")
from cctbx_eltbx_chemical_elements_ext import *

proper_caps_set.__doc__ = """
Retrieves element symbols with only their first letter capitalized.

Returns
-------
scitbx.stl.set.stl_string
"""

proper_caps_list.__doc__ = """
Retrieves element symbols with only their first letter capitalized, sorted by
atomic number.

Returns
-------
list of str

Examples
--------
>>> from cctbx.eltbx.chemical_elements import proper_caps_list
>>> proper_caps_list()[:5]
['H', 'He', 'Li', 'Be', 'B']
"""

proper_upper_set.__doc__ = """
Retrieves element symbols with all letters capitalized.

Returns
-------
scitbx.stl.set.stl_string

Examples
--------
>>> from cctbx.eltbx.chemical_elements import proper_upper_set
>>> list(proper_upper_set())[:5]
['AC', 'AG', 'AL', 'AM', 'AR']
"""

proper_upper_list.__doc__ = """
Retrieves element symbols with all letters capitalized, sorted by atomic
number.

Returns
-------
list of str

Examples
--------
>>> from cctbx.eltbx.chemical_elements import proper_upper_list
>>> proper_upper_list()[:5]
['H', 'HE', 'LI', 'BE', 'B']
"""

proper_and_isotopes_upper_set.__doc__ = """
Like proper_upper_set, but also includes symbols for deuterium and tritium.

Returns
-------
scitbx.stl.set.stl_string
"""

_metals_upper_set = frozenset((
  'LI NA K RB CS FR BE MG CA SR BA RA '
  'SC TI V CR MN FE CO NI CU ZN Y ZR NB MO TC RU RH PD AG CD '
  'HF TA W RE OS IR PT AU HG RF DB SG BH HS MT DS RG '
  'LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU '
  'AC TH PA U NP PU AM CM BK CF ES FM MD NO LR '
  'AL GA IN SN TL PB BI PO').split())

def metals_upper_set():
  """
  Retrieves the symbols of the metals, from lithium up to roentgenium.

  All letters are capitalized. The set holds the alkali and alkaline-earth
  metals, the elements of groups 3 to 12, the lanthanides, the actinides, and
  Al, Ga, In, Sn, Tl, Pb, Bi and Po. The metalloids B, Si, Ge, As, Sb, Te and
  At are not included; periodic tables differ on Po and At.

  Returns
  -------
  frozenset of str

  Examples
  --------
  >>> from cctbx.eltbx.chemical_elements import metals_upper_set
  >>> "ZN" in metals_upper_set(), "AS" in metals_upper_set()
  (True, False)
  """
  return _metals_upper_set
