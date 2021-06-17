"""
This module can be used to retrieve a list or set of all chemical symbols on the
periodic table, from hydrogen up to roentgenium.
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
