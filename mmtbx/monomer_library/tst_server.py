from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import server
from mmtbx.monomer_library import cif_types
from cctbx.array_family import flex
from scitbx.python_utils import dicts
from libtbx import easy_pickle
import random
import sys, os
from six.moves import zip

def check_chem_comp(chem_comp_from_list_file, comp_comp_id):
  if (comp_comp_id is None):
    print("Missing file:", chem_comp_from_list_file.id)
    return
  for key in chem_comp_from_list_file.cif_keywords():
    a = str(getattr(chem_comp_from_list_file, key)).strip()
    b = str(getattr(comp_comp_id.chem_comp, key)).strip()
    if (a != "" and b != "" and a != b):
      print("%s: %s: %s, %s" % (
        chem_comp_from_list_file.id, key, str(a), str(b)))
  a = chem_comp_from_list_file.number_atoms_all
  b = len(comp_comp_id.atom_list)
  if (a != b):
    print("%s: number_atoms_all: %s, %d counted" % (
      chem_comp_from_list_file.id, str(a), b))
  a = chem_comp_from_list_file.number_atoms_nh
  b = 0
  for atom in comp_comp_id.atom_list:
    if (atom.type_symbol != "H"): b += 1
  if (a != b):
    print("%s: number_atoms_nh: %s, %d counted" % (
      chem_comp_from_list_file.id, str(a), b))

def exercise():
  verbose = "--verbose" in sys.argv[1:]
  quick = "--quick" in sys.argv[1:]
  list_cif = server.mon_lib_list_cif()
  srv = server.server(list_cif=list_cif)
  print("srv.root_path:", srv.root_path)
  default_switch = "--default_off" not in sys.argv[1:]
  if (False or default_switch):
    monomers_with_commas = {}
    atom_id_counts = dicts.with_default_value(0)
    for row_id in list_cif.cif["comp_list"]["_chem_comp.id"]:
      if (quick and random.random() < 0.95): continue
      if (verbose): print("id:", row_id)
      comp_comp_id = srv.get_comp_comp_id_direct(comp_id=row_id)
      if (comp_comp_id is None):
        print("Could not instantiating comp_comp_id(%s)" % row_id)
      else:
        has_primes = False
        has_commas = False
        for atom in comp_comp_id.atom_list:
          atom_id_counts[atom.atom_id] += 1
          if (atom.atom_id.find("'") >= 0):
            has_primes = True
          if (atom.atom_id.find(",") >= 0):
            has_commas = True
        if (has_commas):
          monomers_with_commas[comp_comp_id.chem_comp.id] = has_primes
    print(monomers_with_commas)
    atom_ids = flex.std_string(list(atom_id_counts.keys()))
    counts = flex.size_t(list(atom_id_counts.values()))
    perm = flex.sort_permutation(data=counts, reverse=True)
    atom_ids = atom_ids.select(perm)
    counts = counts.select(perm)
    for atom_id,count in zip(atom_ids, counts):
      print(atom_id, count)
  if (False or default_switch):
    for row in list_cif.cif["comp_list"]["_chem_comp"].iterrows():
      if (quick and random.random() < 0.95): continue
      if (verbose): print("id:", row["_chem_comp.id"])
      comp_comp_id = srv.get_comp_comp_id_direct(comp_id=row["_chem_comp.id"])
      check_chem_comp(cif_types.chem_comp(**row), comp_comp_id)
    if ("--pickle" in sys.argv[1:]):
      easy_pickle.dump("mon_lib.pickle", srv)
  if (False or default_switch):
    comp = srv.get_comp_comp_id_direct("GLY")
    comp.show()
    mod = srv.mod_mod_id_dict["COO"]
    comp.apply_mod(mod).show()
  if (False or default_switch):
    comp = srv.get_comp_comp_id_direct("LYS")
    comp.show()
    mod = srv.mod_mod_id_dict["B2C"]
    comp.apply_mod(mod).show()
  if (False or default_switch):
    for row in list_cif.cif["comp_list"]["_chem_comp"].iterrows():
      if (quick and random.random() < 0.95): continue
      comp_comp_id = srv.get_comp_comp_id_direct(row["_chem_comp.id"])
      if (comp_comp_id is not None):
        if (comp_comp_id.classification == "peptide"):
          print(comp_comp_id.chem_comp.id, comp_comp_id.chem_comp.name, end=' ')
          print(row["_chem_comp.group"], end=' ')
          grp = row["_chem_comp.group"].lower().strip()
          if (grp not in ("l-peptide", "d-peptide", "polymer")):
            print("LOOK", end=' ')
            #if (not os.path.isdir("look")): os.makedirs("look")
            #open("look/%s.cif" % row["_chem_comp.id"], "w").write(
              #open(comp_comp_id.file_name).read())
          print()
        elif (row["_chem_comp.group"].lower().find("peptide") >= 0
              or comp_comp_id.chem_comp.group.lower().find("peptide") >= 0):
          print(comp_comp_id.chem_comp.id, comp_comp_id.chem_comp.name, end=' ')
          print(row["_chem_comp.group"], "MISMATCH")
        if (comp_comp_id.classification in ("RNA", "DNA")):
          print(comp_comp_id.chem_comp.id, comp_comp_id.chem_comp.name, end=' ')
          print(row["_chem_comp.group"], end=' ')
          if (comp_comp_id.classification != row["_chem_comp.group"].strip()):
            print(comp_comp_id.classification, "MISMATCH", end=' ')
          print()
        elif (row["_chem_comp.group"].lower().find("NA") >= 0
              or comp_comp_id.chem_comp.group.lower().find("NA") >= 0):
          print(comp_comp_id.chem_comp.id, comp_comp_id.chem_comp.name, end=' ')
          print(row["_chem_comp.group"], "MISMATCH")
  if (False or default_switch):
    for row in list_cif.cif["comp_list"]["_chem_comp"].iterrows():
      if (quick and random.random() < 0.95): continue
      comp_comp_id = srv.get_comp_comp_id_direct(row["_chem_comp.id"])
      if (comp_comp_id is not None):
        planes = comp_comp_id.get_planes()
        for plane in planes:
          dist_esd_dict = {}
          for plane_atom in plane.plane_atoms:
            dist_esd_dict[str(plane_atom.dist_esd)] = 0
          # FIXME: might break compat for py2/3 because indexing a values call
          if (len(dist_esd_dict) != 1 or list(dist_esd_dict.keys())[0] != "0.02"):
            print(comp_comp_id.chem_comp.id, plane.plane_id, end=' ')
            print(list(dist_esd_dict.keys()))
  if (False or default_switch):
    standard_amino_acids = [
      "GLY", "VAL", "ALA", "LEU", "ILE", "PRO", "MET", "PHE", "TRP", "SER",
      "THR", "TYR", "CYS", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"]
    for row in list_cif.cif["comp_list"]["_chem_comp"].iterrows():
      is_standard_aa = row["_chem_comp.id"] in standard_amino_acids
      if (1 and not is_standard_aa):
        continue
      comp_comp_id = srv.get_comp_comp_id_direct(row["_chem_comp.id"])
      if (is_standard_aa):
        assert comp_comp_id is not None
        assert comp_comp_id.chem_comp.group.strip() == "L-peptide", '%s %s' % (row["_chem_comp.id"], comp_comp_id.chem_comp.group)
      if (comp_comp_id is not None):
        print(comp_comp_id.chem_comp.id.strip(), end=' ')
        print(comp_comp_id.chem_comp.name.strip(), end=' ')
        print(comp_comp_id.chem_comp.group.strip())
        for tor in comp_comp_id.tor_list:
          print("  tor:", tor.atom_id_1, tor.atom_id_2, end=' ')
          print(tor.atom_id_3, tor.atom_id_4, tor.value_angle, end=' ')
          print(tor.value_angle_esd, tor.period)
        for chir in comp_comp_id.chir_list:
          print("  chir:", chir.atom_id_centre, chir.atom_id_1, end=' ')
          print(chir.atom_id_2, chir.atom_id_3, chir.volume_sign)
  if (False or default_switch):
    elib = server.ener_lib()
    if (False or default_switch):
      for syn in elib.lib_synonym.items():
        print(syn)
    if (False or default_switch):
      for vdw in elib.lib_vdw:
        vdw.show()
  print("OK")

if (__name__ == "__main__"):
  exercise()
