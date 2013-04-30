
from __future__ import division
from libtbx.utils import Sorry
import sys

def compare_ligands (ligand_code,
    hierarchy_1=None,
    hierarchy_2=None,
    pdb_file_1=None,
    pdb_file_2=None,
    max_distance_between_centers_of_mass=8.0,
    exclude_hydrogens=True,
    out=sys.stdout) :
  from scitbx.array_family import flex
  from scitbx.matrix import col
  assert (ligand_code.isalnum())
  assert (([hierarchy_1, hierarchy_2] == [None, None]) or
          ([pdb_file_1, pdb_file_2] == [None, None]))
  if (hierarchy_1 is None) :
    assert (hierarchy_2 is None)
    from iotbx import file_reader
    pdb_1 = file_reader.any_file(pdb_file_1, force_type="pdb")
    pdb_1.check_file_type("pdb")
    hierarchy_1 = pdb_1.file_object.construct_hierarchy()
    pdb_2 = file_reader.any_file(pdb_file_2, force_type="pdb")
    pdb_2.check_file_type("pdb")
    hierarchy_2 = pdb_2.file_object.construct_hierarchy()
  def extract_ligand (hierarchy) :
    copies = []
    for chain in hierarchy.only_model().chains() :
      for residue_group in chain.residue_groups() :
        for atom_group in residue_group.atom_groups() :
          if (atom_group.resname == ligand_code) :
            copies.append(atom_group)
    return copies
  rmsds = []
  ligands_1 = extract_ligand(hierarchy_1)
  ligands_2 = extract_ligand(hierarchy_2)
  if (len(ligands_1) == 0) or (len(ligands_2) == 0) :
    raise Sorry("One or both models missing residue '%s'!" % ligand_code)
  print >> out, "%d copies in 1st model, %d copies in 2nd model" % \
    (len(ligands_1), len(ligands_2))
  for ligand_1 in ligands_1 :
    matching = []
    sites_1 = ligand_1.atoms().extract_xyz()
    xyz_mean_1 = sites_1.mean()
    for ligand_2 in ligands_2 :
      sites_2 = ligand_2.atoms().extract_xyz()
      xyz_mean_2 = sites_2.mean()
      dxyz = abs(col(xyz_mean_1) - col(xyz_mean_2))
      if (dxyz < max_distance_between_centers_of_mass) :
        matching.append(ligand_2)
    if (len(matching) == 1) :
      ligand_2 = matching[0]
      isel_1 = flex.size_t()
      isel_2 = flex.size_t()
      for i_seq, atom_1 in enumerate(ligand_1.atoms()) :
        if (atom_1.element.strip() in ["H","D"]) and (exclude_hydrogens) :
          continue
        for j_seq, atom_2 in enumerate(ligand_2.atoms()) :
          if (atom_1.name == atom_2.name) :
            isel_1.append(i_seq)
            isel_2.append(j_seq)
            break
      if (len(isel_1) == 0) :
        raise Sorry("No matching atoms found!")
      sites_1 = sites_1.select(isel_1)
      sites_2 = ligand_2.atoms().extract_xyz().select(isel_2)
      rmsd = sites_1.rms_difference(sites_2)
      print >> out, "  '%s' matches '%s': rmsd=%.3f" % (ligand_1.id_str(),
        ligand_2.id_str(), rmsd)
      rmsds.append(rmsd)
  return rmsds
