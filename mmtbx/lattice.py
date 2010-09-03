
from libtbx import group_args
import sys, os

def find_lattice_contacts (xray_structure,
                           pdb_atoms,
                           selected_atoms=None,
                           distance_cutoff=3.5,
                           ignore_same_asu=True,
                           ignore_waters=True) :
  from scitbx.array_family import flex
  sites_frac = xray_structure.sites_frac()
  unit_cell = xray_structure.unit_cell()
  pair_asu_table = xray_structure.pair_asu_table(
    distance_cutoff=distance_cutoff)
  pair_sym_table = pair_asu_table.extract_pair_sym_table()
  contacts = []
  if (selected_atoms is None) :
    selected_atoms = flex.bool(pdb_atoms.size(), True)
  for i_seq,pair_sym_dict in enumerate(pair_sym_table):
    if (not selected_atoms[i]) :
      continue
    site_i = sites_frac[i_seq]
    atom_i = pdb_atoms[i_seq]
    resname_i = atom_i.resname
    atmname_i = atom_i.name
    chainid_i = atom_i.chain_id
    for j_seq,sym_ops in pair_sym_dict.items():
      site_j = sites_frac[j_seq]
      atom_j = pdb_atoms[j_seq]
      resname_j = atom_j.resname
      atmname_j = atom_j.name
      chainid_j = atom_j.chain_id
      for sym_op in sym_ops:
        if sym_op.is_unit_mx() :
          if ignore_same_asu :
            continue
          elif (chainid_i == chainid_j) :
            continue
        if (resname_j in ["HOH","WAT"] and ignore_waters) :
          continue
        site_ji = sym_op * site_j
        distance = unit_cell.distance(site_i, site_ji)
        contacts.append((i_seq, j_seq, sym_op, distance))
        #print resname_i, atmname_i, resname_j, atmname_j, str(sym_op), distance
  return contacts

def show_contacts (contacts, pdb_atoms) :
  for contact in contacts :
    (i_seq, j_seq, sym_op, distance) = contact
    atom_i = pdb_atoms[i_seq]
    atom_j = pdb_atoms[j_seq]
    fmt_i = atom_i.id_str()[5:20]
    fmt_j = atom_j.id_str()[5:20]
    #fmt_i = "%-2s %4s %3s %4s" % (atom_i.chain_id, atom_i.resid(),
    #  atom_i.resname, atom_i.name)
    #fmt_j = "%-2s %4s %3s %4s" % (atom_j.chain_id, atom_j.resid(),
    #  atom_j.resname, atom_j.name)
    print "%s %s %5.2f %s" % (fmt_i,fmt_j,distance,str(sym_op))

def show_contacts_for_pymol (contacts, pdb_atoms, object_name,
    distance_cutoff=3.5) :
  for contact in contacts :
    (i_seq, j_seq, sym_op, distance) = contact
    atom_i = pdb_atoms[i_seq]
    atom_j = pdb_atoms[j_seq]
    s1 = "(%s and chain '%s' and resi %d and name %s)" % (object_name,
      atom_i.chain_id, atom_i.resseq_as_int(), atom_i.name)
    s2 = "((not %s) and (chain '%s' and resi %d and name %s))" % (
      object_name, atom_j.chain_id, atom_j.resseq_as_int(), atom_j.name)
    print "dist %s, %s within %.1f of %s" % (s1, s2, distance_cutoff+0.1, s1)

def apply_sym_op_to_pdb (pdb_hierarchy, sym_op, unit_cell) :
  #import scitbx.matrix
  r = sym_op.r()
  t = sym_op.t()
  #rt = scitbx.matrix.rt((r.as_double(), t.as_double()))
  new_hierarchy = pdb_hierarchy.deep_copy()
  atoms = pdb_hierarchy.atoms()
  sites_frac = unit_cell.fractionalize(sites_cart=atoms.extract_xyz())
  new_sites = sites_frac * r.as_double() + t.as_double()
  atoms.set_xyz(unit_cell.orthogonalize(sites_frac=new_sites))
  return new_hierarchy

def apply_biological_unit (pdb_in) :
  atoms = pdb_in.atoms()
  remark = pdb_in.remark_section()
  if (remark.size() == 0) :
    raise Sorry("No REMARK records in this PDB file.")
  return pdb_out
