from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import server
import sys

def run(args):
  assert len(args) == 0
  srv = server.server()
  standard_amino_acids = [
    "GLY", "VAL", "ALA", "LEU", "ILE", "PRO", "MET", "PHE", "TRP", "SER",
    "THR", "TYR", "CYS", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"]
  for comp_id in standard_amino_acids:
    comp_comp_id = srv.get_comp_comp_id_direct(comp_id)
    print(comp_comp_id.chem_comp.id.strip(), end=' ')
    print(comp_comp_id.chem_comp.name.strip(), end=' ')
    print(comp_comp_id.chem_comp.group.strip())
    hydrogens = {}
    for atom in comp_comp_id.atom_list:
      if (atom.type_symbol == "H"):
        hydrogens[atom.atom_id] = 0
    bond_counts = dict(hydrogens)
    for bond in comp_comp_id.bond_list:
      for atom_id in [bond.atom_id_1, bond.atom_id_2]:
        if (atom_id in hydrogens):
          bond_counts[atom_id] += 1
    if (list(bond_counts.values()) != [1] * len(bond_counts)):
      print("bad bond counts:", bond_counts)
      raise AssertionError
    angle_counts = dict(hydrogens)
    for angle in comp_comp_id.angle_list:
      for atom_id in [angle.atom_id_1, angle.atom_id_2, angle.atom_id_3]:
        if (atom_id in hydrogens):
          angle_counts[atom_id] += 1
    #print angle_counts
    assert min(angle_counts.values()) > 0
    assert max(angle_counts.values()) <= 3
    for atom_id,count in angle_counts.items():
      if (count == 3):
        print("three angles:", atom_id)
        for bond in comp_comp_id.bond_list:
          atom_ids = [bond.atom_id_1, bond.atom_id_2]
          if (atom_id in atom_ids):
            print("  bond:", atom_ids)
        for angle in comp_comp_id.angle_list:
          atom_ids = [angle.atom_id_1, angle.atom_id_2, angle.atom_id_3]
          if (atom_id in atom_ids):
            print("  angle:", atom_ids, angle.value_angle)
    tor_counts = dict(hydrogens)
    for tor in comp_comp_id.tor_list:
      for atom_id in [tor.atom_id_1,
                      tor.atom_id_2,
                      tor.atom_id_3,
                      tor.atom_id_4]:
        if (atom_id in hydrogens):
          tor_counts[atom_id] += 1
    #print tor_counts
    assert max(tor_counts.values()) <= 1
    print("no tor:", list(tor_counts.values()).count(0))
    chir_counts = dict(hydrogens)
    for chir in comp_comp_id.chir_list:
      for atom_id in [chir.atom_id_centre,
                      chir.atom_id_1,
                      chir.atom_id_2,
                      chir.atom_id_3]:
        if (atom_id in hydrogens):
          chir_counts[atom_id] += 1
    #print chir_counts
    assert max(chir_counts.values()) == 0
    plane_counts = dict(hydrogens)
    for plane in comp_comp_id.get_planes():
      for plane_atom in plane.plane_atoms:
        if (plane_atom.atom_id in hydrogens):
          plane_counts[plane_atom.atom_id] += 1
    #print plane_counts
    assert max(plane_counts.values()) <= 1
    for atom_id,count in plane_counts.items():
      if (count == 0): continue
      assert angle_counts[atom_id] == 2

if (__name__ == "__main__"):
  run(sys.argv[1:])
