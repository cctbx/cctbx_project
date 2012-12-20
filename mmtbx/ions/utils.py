
from __future__ import division
import sys

def anonymize_ions (hierarchy, log=sys.stdout) :
  """
  Convert any elemental ions in the PDB hierarchy to water, resetting the
  occupancy and scaling the B-factor.  The atom segids will be set to the old
  resname.  NOTE: this does not change the corresponding scatterer in the xray
  structure, but a new xray structure can be obtained by calling
  hierarchy.extract_xray_structure(crystal_symmetry).
  """
  from cctbx.eltbx import sasaki
  from cctbx.eltbx import chemical_elements
  elements = chemical_elements.proper_upper_list()
  n_converted = 0
  hierarchy = hierarchy.deep_copy()
  for model in hierarchy.models() :
    for chain in model.chains() :
      for residue_group in chain.residue_groups() :
        for atom_group in residue_group.atom_groups() :
          if (atom_group.resname.strip() in elements) :
            atoms = atom_group.atoms()
            id_strs = []
            for atom in atoms :
              elem = atom.element.strip()
              atomic_number = sasaki.table(elem).atomic_number()
              id_strs.append(atom.id_str())
              atom.segid = atom_group.resname
              atom.name = " O  "
              atom.element = "O"
              atom.charge = ""
              atom.occupancy = 1.0
              atom.b = atom.b * (10 / atomic_number)
            atom_group.resname = "HOH"
            for atom, id_str in zip(atoms, id_strs) :
              print >> log, "%s --> %s, B-iso = %.2f" % (id_str, atom.id_str(),
                atom.b)
              n_converted += 1
  return hierarchy, n_converted
