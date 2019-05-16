from __future__ import absolute_import, division, print_function

from scitbx.math import dihedral_angle
from mmtbx.conformation_dependent_library.cdl_setup import \
  before_pro_groups, not_before_pro_groups
from six.moves import range

def distance2(a,b):
  d2 = 0
  for i in range(3):
    d2 += (a.xyz[i]-b.xyz[i])**2
  return d2

def get_c_ca_n(atom_group,
               atom_name_list=[" C  ", " CA ", " N  "],
               return_subset=False):
  assert atom_group
  tmp = []
  outl = []
  for name in atom_name_list:
    atom = atom_group.find_atom_by(name=name)
    if atom:
      tmp.append(atom)
    else:
      outl.append('    missing atom "%s %s %s"' % (
        name,
        atom_group.resname,
        atom_group.resseq,
      ))
      if return_subset:
        tmp.append(None)
      else:
        tmp = None
        break
  return tmp, outl

def round_to_int(d, n=10, wrap=True):
  t = int(round((float(d))/int(n)))*int(n)
  if wrap:
    if t==180: return -180
  return t

def round_to_ten(d):
  return round_to_int(d, 10)

def get_res_type_group(resname1, resname2):
  resname1=resname1.strip()
  resname2=resname2.strip()
  if resname2=="PRO":
    lookup = before_pro_groups
  else:
    lookup = not_before_pro_groups
  for key in lookup:
    if resname1 in lookup[key]:
      return key
  return None

def get_omega_value(residue1, residue2, verbose=False):
  ccn1, outl1 = get_c_ca_n(residue1, return_subset=True)
  ccn2, outl2 = get_c_ca_n(residue2, return_subset=True)
  ca1 = ccn1[1]
  n = ccn1[2]
  c = ccn2[0]
  ca2 = ccn2[1]
  omega_atoms = [ca1, n, c, ca2]
  if verbose:
    for atom in omega_atoms: print(atom.quote())
  if None in omega_atoms: return None
  omega = dihedral_angle(sites=[atom.xyz for atom in omega_atoms], deg=True)
  return omega

def get_phi_psi_atoms(residue1, residue2, residue3, verbose=False):
  backbone_i_minus_1, outl1 = get_c_ca_n(residue1, return_subset=True)
  backbone_i, outl2         = get_c_ca_n(residue2, return_subset=True)
  backbone_i_plus_1, outl3  = get_c_ca_n(residue3, return_subset=True)
  phi_atoms = [
    backbone_i_minus_1[0],
    backbone_i[2],
    backbone_i[1],
    backbone_i[0],
    ]
  psi_atoms = [
    backbone_i[2],
    backbone_i[1],
    backbone_i[0],
    backbone_i_plus_1[2],
    ]
  atoms = [phi_atoms, psi_atoms]
  if len(filter(None, atoms[0]))!=4: return None
  if len(filter(None, atoms[1]))!=4: return None
  if verbose:
    print(atoms)
    for group in atoms:
      for atom in group: print(atom.quote())
  return atoms

def get_phi_psi_angles(residues, verbose=False):
  assert len(residues)>=3
  dihedrals = []
  for i in range(len(residues)):
    if i<2: continue
    atoms = get_phi_psi_atoms(*tuple(residues[i-2:i+1]), verbose=verbose)
    if atoms is None: return None
    for dihedral in atoms:
      phi_or_psi=dihedral_angle(sites=[atom.xyz for atom in dihedral], deg=True)
      dihedrals.append(phi_or_psi)
  if verbose:
    print('dihedrals')
    for phi_or_psi in dihedrals:
      print('phi_or_psi',phi_or_psi)
  return dihedrals

def get_ca_dihedrals(residues, verbose=False):
  assert len(residues)>=4
  dihedrals = []
  atoms = []
  for residue in residues:
    atoms.append(residue.find_atom_by(name=' CA '))
    if len(atoms)==4:
      if verbose:
        print('CAs')
        for atom in atoms: print(atom.quote())
      dihedrals.append(dihedral_angle(sites=[atom.xyz for atom in atoms], deg=True))
      del atoms[0]
  return dihedrals

