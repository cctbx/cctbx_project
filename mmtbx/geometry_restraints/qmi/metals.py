"""Canonical set of metal element symbols + lightweight helpers for
finding metal atoms in a model or atoms array.

Single source of truth for "what counts as a metal" across qmi and
endo_exo. Both modules import ``METALS`` from here instead of carrying
private copies that can drift apart.

Element symbols are in Title case (``'Mg'``, ``'Fe'``, ``'ZN'`` ->
``'Zn'``) to match the ``atom.element.strip().capitalize()`` lookup
convention used downstream.
"""

from __future__ import absolute_import, division, print_function


METALS = {
  'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr',
  'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
  'Sc', 'Y',  # lanthanides and actinides to be added later
  'Ti', 'Zr', 'Hf', 'Rf',
  'V', 'Nb', 'Ta', 'Db',
  'Cr', 'Mo', 'W', 'Sg',
  'Mn', 'Tc', 'Re', 'Bh',
  'Fe', 'Ru', 'Os', 'Hs',
  'Co', 'Rh', 'Ir', 'Mt',
  'Ni', 'Pd', 'Pt', 'Ds',
  'Cu', 'Ag', 'Au', 'Rg',
  'Zn', 'Cd', 'Hg', 'Cn',
  'Al', 'Ga', 'In', 'Tl', 'Nh',
  'Sn', 'Pb', 'Fl',
}


def _atoms_from(atoms_or_model):
  """Accept either an atoms array or an mmtbx.model.manager; return atoms."""
  if hasattr(atoms_or_model, 'get_hierarchy'):
    return atoms_or_model.get_hierarchy().atoms()
  return atoms_or_model


def normalize_element(element_symbol):
  """Element symbol in the Title-case convention qmi uses for lookups
  (``' ZN '`` / ``'zn'`` -> ``'Zn'``). Single source of the strip/capitalize
  idiom used across qmi."""
  return element_symbol.strip().capitalize()


def is_metal(element_symbol, metals=None):
  """True iff ``element_symbol`` (any case, with/without whitespace) is
  a metal per ``metals`` (default ``METALS``)."""
  if metals is None:
    metals = METALS
  return normalize_element(element_symbol) in metals


def metal_atoms(atoms_or_model, metals=None):
  """Return the metal atom objects in ``atoms_or_model``.

  Accepts an iotbx.pdb.hierarchy atoms array or an mmtbx.model.manager.
  """
  return [a for a in _atoms_from(atoms_or_model) if is_metal(a.element, metals)]


def metal_iseqs(atoms_or_model, metals=None):
  """Return iseqs (positional indices) of metal atoms in ``atoms_or_model``."""
  return [i for i, a in enumerate(_atoms_from(atoms_or_model))
          if is_metal(a.element, metals)]


def count_metals(atoms_or_model, metals=None):
  """Return the number of metal atoms in ``atoms_or_model``."""
  return sum(1 for a in _atoms_from(atoms_or_model) if is_metal(a.element, metals))
