
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

def sort_atoms_permutation (pdb_atoms, xray_structure) :
  assert (pdb_atoms.size() == xray_structure.scatterers().size())
  from scitbx.array_family import flex
  pdb_atoms.reset_i_seq()
  atoms_sorted = sorted(pdb_atoms, cmp_atom)
  sele = flex.size_t([ atom.i_seq for atom in atoms_sorted ])
  return sele

# compare mass, then occupancy, then B_iso
def cmp_atom (a, b) :
  from cctbx.eltbx import sasaki
  mass_a = sasaki.table(a.element.strip().upper()).atomic_number()
  mass_b = sasaki.table(b.element.strip().upper()).atomic_number()
  if (mass_a == mass_b) :
    if (a.occ == b.occ) :
      return cmp(b.b, a.b)
    else :
      return cmp(b.occ, a.occ)
  else :
    return cmp(mass_b, mass_a)

def collect_ions (pdb_hierarchy) :
  from cctbx.eltbx import chemical_elements
  elements = chemical_elements.proper_upper_list()
  ions = []
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for residue_group in chain.residue_groups() :
        for atom_group in residue_group.atom_groups() :
          if (atom_group.resname.strip() in elements) :
            atoms = atom_group.atoms()
            assert (len(atoms) == 1)
            for atom in atoms:
              ions.append(atom)
  return ions
