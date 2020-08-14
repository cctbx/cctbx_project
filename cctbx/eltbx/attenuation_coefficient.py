from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex # import dependency
import boost_adaptbx.boost.python as bp
from six.moves import range
ext = bp.import_ext("cctbx_eltbx_attenuation_coefficient_ext")
from cctbx_eltbx_attenuation_coefficient_ext import *


class nist_elements(object):
  ''' A table of nist elements and composite materials.
      Note that elements are only defined up to uranium. '''

  def __init__(self):
    ''' Initialise the table. '''
    self._elements = [
        '',  'H', 'He', 'Li', 'Be',  'B',  'C',  'N',  'O',  'F', 'Ne',
      'Na', 'Mg', 'Al', 'Si',  'P',  'S', 'Cl', 'Ar',  'K', 'Ca', 'Sc',
      'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
      'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',  'Y', 'Zr', 'Nb', 'Mo', 'Tc',
      'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',  'I', 'Xe',
      'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
      'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',  'W', 'Re', 'Os',
      'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr',
      'Ra', 'Ac', 'Th', 'Pa', 'U',
      'CdTe', 'GaAs']

  def __len__(self):
    ''' Return the number of elemets '''
    return len(self._elements) - 1

  def symbol_list(self):
    ''' Return the list of elements. '''
    return list(self._elements[1:])

  def atomic_number_list(self):
    ''' Return a list of atomic numbers. '''
    return list(range(1, len(self._elements)))

  def atomic_number(self, symbol):
    ''' Get the atomic number from the symbol.

    Params:
        symbol The element symbol

    Return:
        The atomic number

    '''
    return self._elements.index(symbol)

  def symbol(self, number):
    ''' Get the symbol from the atomic number

    Params:
        number The atomic number

    Return:
        The symbol

    '''
    assert(number > 0 and number < len(self._elements))
    return self._elements[number]


def chemlex(formula):
  ''' Function taking a chemical formula and getting the component
  elements and number.

  Parsing is simply done by regex and no checking for whether the
  symbol corresponds to anything real is done. Proper capitalisation
  must be used.

  Params:
      formula A chemical formula string

  Returns:
      A list of (symbol, number) tuples

  '''
  import re

  # Check that the formula is valid
  if re.match('^([A-Z][a-z]{0,2}[1-9]?)*$', formula):

    # Split the components by capital letter
    comp = re.findall('[A-Z][^A-Z]*', formula)

    # For each component split into element and number
    el_num = []
    for c in comp:

      r = re.compile('([A-Z][a-z]{0,2})([1-9]+)')
      m = r.match(c)
      if m:
        el = m.group(1)
        num = int(m.group(2))
      else:
        el = c
        num = 1

      # Append the element with number
      el_num.append((el, num))

    # Return the list of components
    return el_num

  else:
    raise ValueError('Unidentified formula')


def get_table(index):
    ''' Get the table for a given element or composite

    Params:
        index Either an atomic number of symbol

    Returns:
        The table for the requested element

    '''
    # Get the index as atomic number
    if isinstance(index, str):
      atomic_number = nist_elements().atomic_number(index)
    else:
      atomic_number = index

    # return the nist element table
    return table(atomic_number)
