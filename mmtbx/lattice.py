
from libtbx import group_args
import sys, os

def find_crystal_contacts (xray_structure,
                           pdb_atoms, # atom_with_labels, not atom!
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
    selected_atoms = flex.bool(len(pdb_atoms), True)
  for i_seq,pair_sym_dict in enumerate(pair_sym_table):
    if (not selected_atoms[i_seq]) :
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

def find_crystal_contacts_by_residue (xray_structure,
                                      pdb_hierarchy,
                                      **kwds) :
  contacts_by_residue = {}
  atoms = list(pdb_hierarchy.atoms_with_labels())
  contacts = find_crystal_contacts(xray_structure, atoms, **kwds)
  for (i_seq, j_seq, sym_op, distance) in contacts :
    atom_rec = atoms[i_seq].fetch_labels()
    residue_key = (atom_rec.chain_id, atom_rec.resname, atom_rec.resid(),
      atom_rec.altloc)
    if (not residue_key in contacts_by_residue) :
      contacts_by_residue[residue_key] = []
    contacts_by_residue[residue_key].append((j_seq, sym_op, distance))
  all_residues = []
  for chain in pdb_hierarchy.models()[0].chains() :
    chain_id = chain.id
    for residue_group in chain.residue_groups() :
      resid = residue_group.resid()
      for atom_group in residue_group.atom_groups() :
        resname = atom_group.resname
        altloc = atom_group.altloc
        residue_key = (chain_id, resname, resid, altloc)
        residue_contacts = contacts_by_residue.get(residue_key, [])
        all_residues.append((residue_key, residue_contacts))
  return all_residues

def extract_closest_contacting_residues (residue_contacts,
                                         pdb_atoms) :
  reduced_contacts = []
  for (residue_key, contacts) in residue_contacts :
    if (len(contacts) == 0) :
      reduced_contacts.append((residue_key, None, None, None))
    else :
      contacts.sort(lambda x,y: cmp(x[2], y[2]))
      (j_seq, sym_op, distance) = contacts[0]
      atom_rec = pdb_atoms[j_seq].fetch_labels()
      contact_key = (atom_rec.chain_id, atom_rec.resname, atom_rec.resid(),
        atom_rec.altloc)
      reduced_contacts.append((residue_key, contact_key, sym_op, distance))
  return reduced_contacts

def summarize_contacts_by_residue (residue_contacts,
                                   pdb_hierarchy,
                                   out=sys.stdout) :
  from mmtbx.refinement.print_statistics import make_header
  summary = extract_closest_contacting_residues(residue_contacts,
    pdb_hierarchy.atoms())
  make_header("Crystal contacts by residue", out=out)
  print >> out, "  %-16s %-16s %-16s %-16s" % ("residue", "closest contact",
    "symop", "distance (A)")
  print >> out, "-"*72
  for (residue_key, contact_key, sym_op, distance) in summary :
    (chain_id, resname, resid, altloc) = residue_key
    id_str = "%s%5s %3s %s" % (chain_id, resid, resname, altloc)
    if (contact_key is None) :
      print >> out, "  %-16s %-16s %-16s %-4s" % (id_str, "*","*","*")
    else :
      (chain_id, resname, resid, altloc) = contact_key
      id_str_2 = "%s%5s %3s %s" % (chain_id, resid, resname, altloc)
      print >> out, "  %-16s %-16s %-16s %-4.2f" % (id_str, id_str_2, sym_op,
        distance)

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

if __name__ == "__main__" :
  pdb_file = sys.argv[1]
  from iotbx import file_reader
  pdb_in = file_reader.any_file(pdb_file).file_object
  pdb_hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structure_simple()
  residue_contacts = find_crystal_contacts_by_residue(xrs, pdb_hierarchy)
  summarize_contacts_by_residue(residue_contacts, pdb_hierarchy)
