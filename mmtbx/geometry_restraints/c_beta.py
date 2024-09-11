from __future__ import absolute_import, division, print_function
import cctbx.geometry_restraints
from cctbx.array_family import flex
from iotbx.pdb.amino_acid_codes import three_letter_l_given_three_letter_d
from libtbx.utils import Sorry


def get_c_beta_torsion_proxies(pdb_hierarchy,
                               selection=None,
                               sigma=2.5):
  origin_ids = cctbx.geometry_restraints.linking_class.linking_class()
  c_beta_origin_id = origin_ids.get_origin_id('C-beta')
  if (selection is not None):
    if (isinstance(selection, flex.bool)):
      actual_bselection = selection
    elif (isinstance(selection, flex.size_t)):
      actual_bselection = flex.bool(pdb_hierarchy.atoms_size(), False)
      actual_bselection.set_selected(selection, True)
    else:
      raise Sorry("Bad selection supplied for c_beta restraints")
  if selection is None:
    actual_bselection = flex.bool(pdb_hierarchy.atoms_size(), True)
  cache = pdb_hierarchy.atom_selection_cache()
  sel = cache.selection("name N or name CA or name C or name CB")
  c_beta_dihedral_proxies = \
      cctbx.geometry_restraints.shared_dihedral_proxy()
  c_beta_residues_skipped = {}
  s0 = set([' N  ', ' CA ', ' C  ', ' CB '])
  for model in pdb_hierarchy.select(sel).models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for conformer in rg.conformers():
          for residue in conformer.residues():
            if(not s0.issubset(set(residue.atoms().extract_name()))): continue
            CB_atom = residue.find_atom_by(name=" CB ")
            if residue.resname in three_letter_l_given_three_letter_d:
              c_beta_residues_skipped.setdefault("d-peptide", [])
              c_beta_residues_skipped['d-peptide'].append(CB_atom)
              continue
            CA_atom = residue.find_atom_by(name=" CA ")
            N_atom  = residue.find_atom_by(name=" N  ")
            C_atom  = residue.find_atom_by(name=" C  ")
            if(N_atom is not None and CA_atom is not None
               and C_atom is not None and CB_atom is not None):
              if not (actual_bselection[N_atom.i_seq] and
                  actual_bselection[CA_atom.i_seq] and
                  actual_bselection[C_atom.i_seq] and
                  actual_bselection[CB_atom.i_seq] ):
                continue
              # check for correct chiral volume
              sites_cart = [N_atom.xyz, C_atom.xyz, CA_atom.xyz, CB_atom.xyz]
              CAC_dist = C_atom.distance(CA_atom)
              if CAC_dist>2.:
                c_beta_residues_skipped.setdefault('CA---C', [])
                c_beta_residues_skipped['CA---C'].append(CB_atom)
                continue
              chiral = cctbx.geometry_restraints.chirality(
                sites_cart,
                volume_ideal=0.,
                both_signs=True,
                weight=1.,
                )
              if(chiral.volume_model<0):
                c_beta_residues_skipped.setdefault("-ve", [])
                c_beta_residues_skipped["-ve"].append(CB_atom)
                continue
              dihedralNCAB, dihedralCNAB = get_cb_target_angle_pair(
                                             resname=residue.resname)
              #NCAB
              i_seqs = [N_atom.i_seq,C_atom.i_seq,CA_atom.i_seq,CB_atom.i_seq]
              dp_add = cctbx.geometry_restraints.dihedral_proxy(
                i_seqs=i_seqs,
                angle_ideal=dihedralNCAB,
                weight=1/sigma**2,
                origin_id=c_beta_origin_id)
              c_beta_dihedral_proxies.append(dp_add)
              #CNAB
              i_seqs = [C_atom.i_seq,N_atom.i_seq,CA_atom.i_seq,CB_atom.i_seq]
              dp_add = cctbx.geometry_restraints.dihedral_proxy(
                i_seqs=i_seqs,
                angle_ideal=dihedralCNAB,
                weight=1/sigma**2,
                origin_id=c_beta_origin_id)
              c_beta_dihedral_proxies.append(dp_add)
  return c_beta_dihedral_proxies, c_beta_residues_skipped # BAD

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
