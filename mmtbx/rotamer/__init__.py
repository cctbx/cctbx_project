from __future__ import absolute_import, division, print_function

def iterator(mon_lib_srv, residue, atom_selection_bool):
  atoms = residue.atoms()
  if (atom_selection_bool is not None):
    if (atom_selection_bool.select(
          indices=residue.atoms().extract_i_seq()).all_eq(False)):
      return None
  rotamer_iterator = mon_lib_srv.rotamer_iterator(
    comp_id=residue.resname,
    atom_names=residue.atoms().extract_name(),
    sites_cart=residue.atoms().extract_xyz())
  if (rotamer_iterator is None):
    return None
  if (rotamer_iterator.problem_message is not None):
    return None
  if (rotamer_iterator.rotamer_info is None):
    return None
  return rotamer_iterator

def improper_ncab_from_atoms(thisN, thisC, thisCA, thisCB):
  assert (not None in [thisN, thisC, thisCA, thisCB])
  return phi_from_sites(thisN.xyz, thisC.xyz, thisCA.xyz, thisCB.xyz)

def improper_cnab_from_atoms(thisC, thisN, thisCA, thisCB):
  assert (not None in [thisC, thisN, thisCA, thisCB])
  return phi_from_sites(thisC.xyz, thisN.xyz, thisCA.xyz, thisCB.xyz)

def omega_from_atoms(prevCA, prevC, thisN, thisCA):
  assert (not None in [prevCA, prevC, thisN, thisCA])
  return phi_from_sites(prevCA.xyz, prevC.xyz, thisN.xyz, thisCA.xyz)

def get_phi_psi_atoms(prev_res, residue, next_res):
  c1, n2, ca2, c2, n3 = None, None, None, None, None
  for atom in prev_res.atoms():
    if (atom.name == " C  "):
      c1 = atom
      break
  for atom in residue.atoms():
    if (atom.name == " N  "):
      n2 = atom
    elif (atom.name == " CA "):
      ca2 = atom
    elif (atom.name == " C  "):
      c2 = atom
  for atom in next_res.atoms():
    if (atom.name == " N  "):
      n3 = atom
  return (c1, n2, ca2, c2, n3)

def get_omega_atoms(prev_atoms, atoms):
  prevCA, prevC, thisN, thisCA = None, None, None, None
  if (prev_atoms is not None):
    for atom in prev_atoms:
      if (atom.name == " CA "): prevCA = atom
      if (atom.name == " C  "): prevC = atom
  if (atoms is not None):
    for atom in atoms:
      if (atom.name == " N  "): thisN = atom
      if (atom.name == " CA "): thisCA = atom
  return prevCA, prevC, thisN, thisCA

def is_cis_peptide(prev_atoms, atoms):
  prevCA, prevC, thisN, thisCA = get_omega_atoms(prev_atoms, atoms)
  if (not None in [ prevCA, prevC, thisN, thisCA]):
    omega = omega_from_atoms(prevCA, prevC, thisN, thisCA)
    if(omega > -30 and omega < 30):
      return True
  return False

def get_phi_psi_indices(prev_res, residue, next_res):
  (c1, n2, ca2, c2, n3) = get_phi_psi_atoms(
    prev_res=prev_res,
    residue=residue,
    next_res=next_res)
  if (not None in [c1, n2, ca2, c2, n3]):
    return [c1.i_seq, n2.i_seq, ca2.i_seq, c2.i_seq, n3.i_seq]
  return [None] * 5

def phi_from_atoms(prevC, resN, resCA, resC):
  assert (not None in [prevC, resN, resCA, resC])
  return phi_from_sites(prevC.xyz, resN.xyz, resCA.xyz, resC.xyz)

def phi_from_sites(prevC, resN, resCA, resC):
  from cctbx import geometry_restraints
  b = geometry_restraints.bond(
    sites=[prevC,resN],
    distance_ideal=1,
    weight=1)
  # check to see if residues are actually bonded.
  if (b.distance_model > 4): return None
  d = geometry_restraints.dihedral(
    sites=[prevC,resN,resCA,resC],
    angle_ideal=-40,
    weight=1)
  return d.angle_model

def psi_from_atoms(resN, resCA, resC, nextN):
  assert (not None in [resN, resCA, resC, nextN])
  return psi_from_sites(resN.xyz, resCA.xyz, resC.xyz, nextN.xyz)

def psi_from_sites(resN, resCA, resC, nextN):
  from cctbx import geometry_restraints
  b = geometry_restraints.bond(
    sites=[resC,nextN],
    distance_ideal=1,
    weight=1)
  if (b.distance_model > 4): return None
  d = geometry_restraints.dihedral(
    sites=[resN,resCA,resC,nextN],
    angle_ideal=-40,
    weight=1)
  return d.angle_model

def phi_psi_from_sites(i_seqs, sites_cart):
  assert (not None in i_seqs) and (len(i_seqs) == 5)
  phi = phi_from_sites(
    prevC=sites_cart[i_seqs[0]],
    resN=sites_cart[i_seqs[1]],
    resCA=sites_cart[i_seqs[2]],
    resC=sites_cart[i_seqs[3]])
  psi = psi_from_sites(
    resN=sites_cart[i_seqs[1]],
    resCA=sites_cart[i_seqs[2]],
    resC=sites_cart[i_seqs[3]],
    nextN=sites_cart[i_seqs[4]])
  return (phi, psi)
