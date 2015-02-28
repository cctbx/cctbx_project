from __future__ import division
import cctbx.geometry_restraints
from cctbx.array_family import flex
from iotbx.pdb.amino_acid_codes import three_letter_l_given_three_letter_d
from libtbx.utils import Sorry

def get_c_beta_torsion_proxies(pdb_hierarchy,
                               selection=None,
                               sigma=2.5):
  if (selection is not None):
    if (isinstance(selection, flex.bool)):
      actual_iselection = selection.iselection()
    elif (isinstance(selection, flex.size_t)):
      actual_iselection = selection
    else:
      raise Sorry("Bad selection supplied for c_beta restraints")
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
                if selection is not None:
                  if ( (N_atom.i_seq not in actual_iselection) or
                       (CA_atom.i_seq not in actual_iselection) or
                       (C_atom.i_seq not in actual_iselection) or
                       (CB_atom.i_seq not in actual_iselection) ):
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
  target_angle_dict = {
    "ALA" : (122.9, -122.6),
    "PRO" : (115.1, -120.7),
    "VAL" : (123.4, -122.0),
    "THR" : (123.4, -122.0),
    "ILE" : (123.4, -122.0),
    "GLY" : (121.6, -121.6)
  }
  dihedralNCAB, dihedralCNAB = target_angle_dict.get(resname, (122.8, -122.6))
  return dihedralNCAB, dihedralCNAB
