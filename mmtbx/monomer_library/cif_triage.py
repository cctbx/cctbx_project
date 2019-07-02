from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library.server
from libtbx.str_utils import show_string, show_sorted_by_counts
from libtbx import dict_with_default_0

def check_comp(file_name):
  result = 0
  print("file name:", file_name)
  cif_object = mmtbx.monomer_library.server.read_cif(file_name)
  for comp_comp_id in mmtbx.monomer_library.server.convert_comp_list(
                        source_info=file_name, cif_object=cif_object):
    result += 1
    atom_names = set()
    for atom in comp_comp_id.atom_list:
      atom_name = atom.atom_id
      assert atom_name.count(" ") == 0
      if (atom_name in atom_names):
        raise RuntimeError(
          "Duplicate atom name: %s" % show_string(atom_name))
      atom_names.add(atom_name)
    print("  number of atoms:", len(atom_names))
    #
    bond_atom_ids = set()
    for bond in comp_comp_id.bond_list:
      atom_ids = [bond.atom_id_1, bond.atom_id_2]
      for atom_name in atom_ids:
        if (atom_name not in atom_names):
          raise RuntimeError(
            "Unknown bond atom name: %s" % show_string(atom_name))
      atom_ids = tuple(sorted(atom_ids))
      if (atom_ids in bond_atom_ids):
        raise RuntimeError(
          "Duplicate bond: %s - %s" % tuple([show_string(s)
            for s in atom_ids]))
      bond_atom_ids.add(atom_ids)
    print("  number of bonds:", len(bond_atom_ids))
    #
    angle_atom_ids = set()
    for angle in comp_comp_id.angle_list:
      atom_ids = [angle.atom_id_1, angle.atom_id_2, angle.atom_id_3]
      for atom_name in atom_ids:
        if (atom_name not in atom_names):
          raise RuntimeError(
            "Unknown angle atom name: %s" % show_string(atom_name))
      atom_ids = tuple(sorted(atom_ids))
      if (atom_ids in angle_atom_ids):
        raise RuntimeError(
          "Duplicate angle: %s - %s - %s" % tuple([show_string(s)
            for s in atom_ids]))
      angle_atom_ids.add(atom_ids)
    print("  number of angles:", len(angle_atom_ids))
    #
    tor_atom_ids = set()
    for tor in comp_comp_id.tor_list:
      atom_ids = [tor.atom_id_1, tor.atom_id_2, tor.atom_id_3, tor.atom_id_4]
      for atom_name in atom_ids:
        if (atom_name not in atom_names):
          raise RuntimeError(
            "Unknown tor atom name: %s" % show_string(atom_name))
      atom_ids = tuple(sorted(atom_ids))
      if (atom_ids in tor_atom_ids):
        raise RuntimeError(
          "Duplicate tor: %s - %s - %s - %s" % tuple([show_string(s)
            for s in atom_ids]))
      tor_atom_ids.add(atom_ids)
    print("  number of tors:", len(tor_atom_ids))
    tor_atom_ids = {}
    for tor in comp_comp_id.tor_list:
      atom_ids = tuple(sorted([tor.atom_id_2, tor.atom_id_3]))
      tor_atom_ids.setdefault(atom_ids, []).append(tor)
    for atom_ids,tors in tor_atom_ids.items():
      if (len(tors) != 1):
        print("    redundant tors:", ", ".join([tor.id for tor in tors]))
    #
    chir_atom_ids = set()
    for chir in comp_comp_id.chir_list:
      atom_ids = [
        chir.atom_id_1, chir.atom_id_2, chir.atom_id_3, chir.atom_id_centre]
      for atom_name in atom_ids:
        if (atom_name not in atom_names):
          raise RuntimeError(
            "Unknown chir atom name: %s" % show_string(atom_name))
      atom_ids = tuple(sorted(atom_ids))
      if (atom_ids in chir_atom_ids):
        raise RuntimeError(
          "Duplicate chir: %s - %s - %s - %s" % tuple([show_string(s)
            for s in atom_ids]))
      chir_atom_ids.add(atom_ids)
    print("  number of chirs:", len(chir_atom_ids))
    #
    plane_atom_counts = dict_with_default_0()
    for plane_atom in comp_comp_id.plane_atom_list:
      if (plane_atom.atom_id not in atom_names):
        raise RuntimeError(
          "Unknown plane atom name: %s" % show_string(plane_atom.atom_id))
      plane_atom_counts[plane_atom.plane_id] += 1
    print("  number of planes:", len(plane_atom_counts))
    if (len(plane_atom_counts) != 0):
      show_sorted_by_counts(
        label_count_pairs=list(plane_atom_counts.items()),
        prefix="    ")
      assert min(plane_atom_counts.values()) >= 3
    #
    rotamer_info = comp_comp_id.rotamer_info()
    if (rotamer_info is not None):
      print("  rotamer_info.tor_ids:", rotamer_info.tor_ids)
      for tor_id in rotamer_info.tor_ids:
        assert tor_id.strip() == tor_id
        assert tor_id.split() == [tor_id]
      for tor_atom_ids in rotamer_info.tor_atom_ids:
        assert len(tor_atom_ids) == 5
        assert tor_atom_ids[0] in rotamer_info.tor_ids
        for atom_id in tor_atom_ids[1:]:
          assert atom_id.strip() == atom_id
          assert atom_id.split() == [atom_id]
      atom_ids = rotamer_info.atom_ids_not_handled
      if (atom_ids is not None):
        for atom_id in atom_ids:
          assert atom_id.strip() == atom_id
          assert atom_id.split() == [atom_id]
      assert (
           rotamer_info.constrain_dihedrals_with_sigma_less_than_or_equal_to
             is None
        or rotamer_info.constrain_dihedrals_with_sigma_less_than_or_equal_to
             > 0)
      print("  number of rotamers:", len(rotamer_info.rotamer))
      n_missing_frequencies = 0
      for rotamer in rotamer_info.rotamer:
        assert rotamer.id is not None
        assert len(rotamer.id.strip()) == len(rotamer.id)
        assert len(rotamer.id.split()) == 1
        if (rotamer.frequency is None):
          if (rotamer.frequency_annotation != "for more uniform sampling"):
            n_missing_frequencies += 1
        else:
          assert rotamer.frequency > 0
          assert rotamer.frequency < 1
        assert rotamer.angles is not None
        assert len(rotamer.angles) == len(rotamer_info.tor_ids)
        for angle in rotamer.angles:
          assert angle is None or -180 < angle <= 180
      if (n_missing_frequencies != 0):
        print("  WARNING: number of missing frequencies:", \
          n_missing_frequencies)
  return result
