# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
Utility functions used within this module.
"""
from __future__ import division

def iterate_sites(pdb_hierarchy, split_sites=False, res_filter=None):
  """
  Returns a generator iterating over all atoms in pdb_hierarchy. Optionally
  skips sites with alternate conformations and can filter by residue name.

  Parameters
  ----------
  pdb_hierarchy: iotbx.pdb_hierarchy
  split_sites: bool, optional
      Indicates whether to iterate over sites with alternate conformations, by
      default they are not included.
  res_filter: list of str, optional
      List of residue names to include, by default, all residues are examined.
  """

  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for atom_group in residue_group.atom_groups():
          resname = atom_group.resname.strip().upper()
          if res_filter is None or resname in res_filter:
            atoms = atom_group.atoms()
            if len(atoms) == 1 or split_sites:
              for atom in atoms:
                element = atom.element.strip().upper()
                if element in ["H", "D"]:
                  continue
                yield atom
