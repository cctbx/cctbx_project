from __future__ import absolute_import, division, print_function

def hydrogen_bond_from_selection_pair(donor_sele, acceptor_sele,
    selection_cache):
  isel = selection_cache.iselection
  donor_i_seqs = isel(donor_sele)
  acceptor_i_seqs = isel(acceptor_sele)
  n_donor_sel = donor_i_seqs.size()
  n_acceptor_sel = acceptor_i_seqs.size()
  if n_donor_sel == 0 or n_acceptor_sel == 0 :
    raise RuntimeError("""\
analyze_h_bonds: one or more atoms missing
  %s (%d atoms)
  %s (%d atoms)""" % (donor_sele, n_donor_sel, acceptor_sele, n_acceptor_sel))
  elif n_donor_sel > 1 or n_acceptor_sel > 1 :
    raise RuntimeError("""\
analyze_h_bonds: multiple atoms matching a selection
  %s (%d atoms)
  %s (%d atoms)""" % (donor_sele, n_donor_sel, acceptor_sele, n_acceptor_sel))
  return (donor_i_seqs[0], acceptor_i_seqs[0])

def get_residue_name_from_selection(resi_sele, selection_cache, atoms):
  i_seqs = selection_cache.iselection(resi_sele)
  if len(i_seqs) == 0 :
    raise RuntimeError("Empty selection '%s'" % resi_sele)
  resnames = []
  for i_seq in i_seqs :
    resnames.append(atoms[i_seq].fetch_labels().resname)
  unique_resnames = set(resnames)
  assert len(unique_resnames) != 0
  if len(unique_resnames) > 1 :
    raise RuntimeError("More than one residue name in selection '%s' (%s)" %
      (resi_sele, ", ".join(unique_resnames)))
  else :
    return resnames[0]
