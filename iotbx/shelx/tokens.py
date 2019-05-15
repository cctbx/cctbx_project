from __future__ import absolute_import, division, print_function
from six.moves import range

(atom_tok,
 element_tok,
 forward_range_tok,
 backward_range_tok,
 residue_class_tok,
 residue_number_tok,
 all_residues_tok,
 plus_tok,
 minus_tok) = range(9)


class atomname_token(object):
  """
  An atomname token that can be passed to a SHELXL instruction
  """
  __slots__ = ("name", "symmetry", "residue_number", "plus_minus")

  def __init__(self, name, symmetry=None, residue_number=None, plus_minus=None):
    self.name = name
    self.symmetry = symmetry
    self.residue_number = residue_number
    self.plus_minus = plus_minus

  def __eq__(self, other):
    return (self.name == other.name
            and self.symmetry == other.symmetry
            and self.residue_number == other.residue_number
            and self.plus_minus == other.plus_minus)


class element_token(object):
  """
  An element name may also be passed to a SHELXL instruction in place of
  an atom name
  """
  __slots__ = ("element", "residue_number", "plus_minus")

  def __init__(self, element, residue_number=None, plus_minus=None):
    self.element = element
    self.residue_number = residue_number
    self.plus_minus = plus_minus

  def __eq__(self, other):
    return (self.element == other.element
            and self.residue_number == other.residue_number
            and self.plus_minus == other.plus_minus)


class residue_token(object):
  __slots__ = ("class_", "number")

  def __init__(self, class_=None, number=None):
    self.class_ = class_
    self.number = number

  def __eq__(self, other):
    return (self.class_ == other.class_ and self.number == other.number)
