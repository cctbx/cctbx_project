"""Seed-atom discovery for the endo_exo QM-region builder."""

from __future__ import absolute_import, division, print_function

from mmtbx.geometry_restraints.qmi import metals as qmi_metals


class SeedFinder:
  """Locate seed atoms within a CCTBX model hierarchy.

  Seeds are either all metal atoms in the structure (default) or the atoms
  matched by one or more user-supplied CCTBX selection strings.
  """

  def find_by_element(self, model, element_filter=None):
    """Return atom objects whose element matches *element_filter*.

    Delegates to ``mmtbx.geometry_restraints.qmi.metals.metal_atoms``. When
    *element_filter* is given, its element symbols are passed straight through
    as the ``metals=`` argument, restricting the scan to those element(s),
    which need not be metals. When ``None``, the canonical ``METALS``
    recognition list is used as the default filter.

    Parameters
    ----------
    model : mmtbx.model.manager
    element_filter : iterable of str or None, optional
        Element symbols (case-/whitespace-tolerant, matched against
        ``element.strip().capitalize()``).

    Returns
    -------
    list of iotbx.pdb.hierarchy.atom
    """
    if element_filter:
      wanted = {
        str(e).strip().capitalize()
        for e in element_filter if str(e).strip()
      }
      if wanted:
        return qmi_metals.metal_atoms(model, metals=wanted)
    return qmi_metals.metal_atoms(model)

  def find_by_selection(self, model, selection_str):
    """Return atom objects matched by *selection_str*.

    Parameters
    ----------
    model : mmtbx.model.manager
    selection_str : str
        CCTBX atom selection string.

    Returns
    -------
    list of iotbx.pdb.hierarchy.atom
    """
    mask = model.selection(selection_str)
    atoms = model.get_hierarchy().atoms()
    return [atoms[i] for i in mask.iselection()]

  def find(self, model, selection_strings=None, element_filter=None):
    """Return seed groups for the given model.

    Each group is a ``(label, atoms)`` tuple where *label* is the selection
    string that produced the group (or ``None`` for metal-scan groups) and
    *atoms* is the list of seed atoms for that group.

    When *selection_strings* is a non-empty list each entry produces one
    group, and *element_filter* is ignored.  When it is empty or ``None``
    every metal atom in the model becomes its own group, optionally
    restricted to the element(s) listed in *element_filter*.

    Parameters
    ----------
    model : mmtbx.model.manager
    selection_strings : list of str or None, optional
    element_filter : iterable of str or None, optional
        Element symbols (case-insensitive) restricting the metal-scan.
        Only consulted when *selection_strings* is empty.

    Returns
    -------
    list of tuple
        Each element is ``(str or None, list of iotbx.pdb.hierarchy.atom)``.
    """
    if selection_strings:
      return [
        (sel_str, self.find_by_selection(model, sel_str))
        for sel_str in selection_strings
      ]
    return [
      (None, [m])
      for m in self.find_by_element(model, element_filter=element_filter)
    ]

  def count_metals(self, model):
    """Return the number of metal atoms in *model*.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    int
    """
    return qmi_metals.count_metals(model)
