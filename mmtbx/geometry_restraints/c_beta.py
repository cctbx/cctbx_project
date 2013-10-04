from __future__ import division
import cctbx.geometry_restraints
from cctbx.array_family import flex
from iotbx.pdb.amino_acid_codes import three_letter_l_given_three_letter_d

def get_c_beta_torsion_proxies(pdb_hierarchy,
                               selection=None,
                               sigma=2.5):
  if (selection is not None):
    if (isinstance(selection, flex.bool)):
      selection = selection.iselection()
  if selection is None:
    selection = flex.bool(
      len(pdb_hierarchy.atoms()),
      True).iselection()
  c_beta_dihedral_proxies = \
      cctbx.geometry_restraints.shared_dihedral_proxy()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for conformer in rg.conformers():
          if conformer.is_protein():
            for residue in conformer.residues():
              if residue.resname in three_letter_l_given_three_letter_d:
                continue
              N_atom = None
              CA_atom = None
              C_atom = None
              CB_atom = None
              for atom in residue.atoms():
                if atom.name.strip() == "N":
                  N_atom = atom
                elif atom.name.strip() == "CA":
                  CA_atom = atom
                elif atom.name.strip() == "C":
                  C_atom = atom
                elif atom.name.strip() == "CB":
                  CB_atom = atom
              if ( (N_atom is not None) and
                   (CA_atom is not None) and
                   (C_atom is not None) and
                   (CB_atom is not None) ):
                if ( (N_atom.i_seq not in selection) or
                     (CA_atom.i_seq not in selection) or
                     (C_atom.i_seq not in selection) or
                     (CB_atom.i_seq not in selection) ):
                  continue
                dihedralNCAB, dihedralCNAB = get_cb_target_angle_pair(
                                               resname=residue.resname)
                #NCAB
                i_seqs = [N_atom.i_seq,
                          C_atom.i_seq,
                          CA_atom.i_seq,
                          CB_atom.i_seq]
                dp_add = cctbx.geometry_restraints.dihedral_proxy(
                  i_seqs=i_seqs,
                  angle_ideal=dihedralNCAB,
                  weight=1/sigma**2)
                c_beta_dihedral_proxies.append(dp_add)
                #CNAB
                i_seqs = [C_atom.i_seq,
                          N_atom.i_seq,
                          CA_atom.i_seq,
                          CB_atom.i_seq]
                dp_add = cctbx.geometry_restraints.dihedral_proxy(
                  i_seqs=i_seqs,
                  angle_ideal=dihedralCNAB,
                  weight=1/sigma**2)
                c_beta_dihedral_proxies.append(dp_add)
  return c_beta_dihedral_proxies

def target_and_gradients(
      sites_cart,
      c_beta_dihedral_proxies,
      gradient_array,
      unit_cell=None):
  target = 0.0
  if unit_cell is None:
    target += cctbx.geometry_restraints.dihedral_residual_sum(
                sites_cart=sites_cart,
                proxies=c_beta_dihedral_proxies,
                gradient_array=gradient_array)
  else:
    target += cctbx.geometry_restraints.dihedral_residual_sum(
                unit_cell=unit_cell,
                sites_cart=sites_cart,
                proxies=c_beta_dihedral_proxies,
                gradient_array=gradient_array)
  return target

def get_cb_target_angle_pair(resname):
  if(resname == "ALA"):
    dihedralNCAB = 122.9
    dihedralCNAB = -122.6
  elif(resname == "PRO"):
    dihedralNCAB = 115.1
    dihedralCNAB = -120.7
  elif( (resname == "VAL") or
        (resname == "THR") or
        (resname == "ILE") ):
    dihedralNCAB = 123.4
    dihedralCNAB = -122.0
  elif(resname == "GLY"):
    dihedralNCAB = 121.6
    dihedralCNAB = -121.6
  else:
    dihedralNCAB = 122.8
    dihedralCNAB = -122.6
  return dihedralNCAB, dihedralCNAB
