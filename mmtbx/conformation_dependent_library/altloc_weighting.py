from __future__ import absolute_import, division, print_function
import sys
import math

import iotbx.pdb
from libtbx.utils import Sorry

def update_restraints(hierarchy,
                      restraints_manager,
                      altloc_weights_object={},
                      min_occ = 0.5,
                      sqrt_func=False,
                      func_factor=1.,
                      bond_weighting=True,
                      angle_weighting=True,
                      log=None,
                      verbose=False,
                      ):
  #raise Sorry("altloc weighting still being tested")
  if verbose:
    print('min_occ',min_occ)
    print('sqrt_func',sqrt_func)
    print('func_factor',func_factor)
    print('bond_weighting',bond_weighting)
    print('angle_weighting',angle_weighting)
  bond_params_table = restraints_manager.geometry.bond_params_table
  altloc_bonds = {}
  previous_residue_group = None
  for model in hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if verbose: print('chain: "%s"' % chain.id)
      for residue_group in chain.residue_groups():
        if verbose: print('  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode))
        # have residue_group and previous_residue_group
        for atom_group_i, atom_group1 in enumerate(residue_group.atom_groups()):
          #
          if atom_group1.resname in ["HOH"]: continue
          #
          if verbose: print('    %03d atom_group: altloc="%s" resname="%s"' % (
            atom_group_i, atom_group1.altloc, atom_group1.resname))
          for atom_group_j, atom_group2 in enumerate(residue_group.atom_groups()):
            #
            if atom_group_i<atom_group_j: continue
            #
            if atom_group2.resname in ["HOH"]: continue
            #
            if verbose: print('    %03d atom_group: altloc="%s" resname="%s"' % (
              atom_group_j, atom_group2.altloc, atom_group2.resname))
            #################################
            # both blank altloc
            #################################
            if not atom_group1.altloc.strip() and not atom_group2.altloc.strip():
              continue
            #################################
            # both altloc
            #################################
            if atom_group1.altloc.strip() and atom_group2.altloc.strip():
              if atom_group1.altloc.strip()!=atom_group2.altloc.strip():
                continue
            #################################
            # one blank altloc, one not
            #################################
            else:
              pass
            # intra-atom_group bonds
            for i, atom1 in enumerate(atom_group1.atoms()):
              for j, atom2 in enumerate(atom_group1.atoms()):
                if i<=j: continue
                bond = bond_params_table.lookup(atom1.i_seq, atom2.i_seq)
                if bond is None: continue
                altloc_bonds[tuple(sorted([atom1.i_seq, atom2.i_seq]))] = [bond, atom1.occ]
                #assert atom1.is_in_same_conformer_as(atom2)
                if verbose:
                  print('intra-atom_group')
                  print(atom1.format_atom_record())
                  print(atom2.format_atom_record())
                  print(bond.distance_ideal, bond.weight)

            # inter-atom_group bonds
            if atom_group_i==atom_group_j: continue
            for i, atom1 in enumerate(atom_group1.atoms()):
              for j, atom2 in enumerate(atom_group2.atoms()):
                bond = bond_params_table.lookup(atom1.i_seq, atom2.i_seq)
                if bond is None: continue
                altloc_bonds[tuple(sorted([atom1.i_seq, atom2.i_seq]))] = [bond, atom1.occ]
                #altloc_bonds.append(bond)
                #assert atom1.is_in_same_conformer_as(atom2)
                if verbose:
                  print('inter-atom_group')
                  print(atom1.format_atom_record())
                  print(atom2.format_atom_record())
                  print(bond.distance_ideal, bond.weight)

        # inter-residue_groups
        if previous_residue_group is None:
          previous_residue_group = residue_group
          continue
        for atom_group in residue_group.atom_groups():
          for atom_group_k, previous_atom_group in enumerate(previous_residue_group.atom_groups()):
            #################################
            # both blank altloc
            #################################
            if (not atom_group.altloc.strip() and
                not previous_atom_group.altloc.strip()
                ):
              continue
            #################################
            # both altloc
            #################################
            if atom_group.altloc.strip() and previous_atom_group.altloc.strip():
              if atom_group.altloc.strip()!=previous_atom_group.altloc.strip():
                continue
            #################################
            # one blank altloc, one not
            #################################
            else:
              pass
            if verbose: print('    %03d previous_atom_group: altloc="%s" resname="%s"' % (
                atom_group_k, previous_atom_group.altloc, previous_atom_group.resname))
            for i, atom1 in enumerate(atom_group.atoms()):
              for j, atom2 in enumerate(previous_atom_group.atoms()):
                bond = bond_params_table.lookup(atom1.i_seq, atom2.i_seq)
                if bond is None: continue
                altloc_bonds[tuple(sorted([atom1.i_seq, atom2.i_seq]))] = [bond, atom1.occ]
                if verbose:
                  print('inter-residue_group')
                  print(atom1.format_atom_record())
                  print(atom2.format_atom_record())
                  print(bond.distance_ideal, bond.weight)

        previous_residue_group = residue_group

  i_seqs = {}
  for key in altloc_bonds:
    bond, occ = altloc_bonds[key]
    #print 'altloc_bonds',key,bond,occ
    if key in altloc_weights_object:
      bond.weight = altloc_weights_object[key]
    else:
      altloc_weights_object[key] = bond.weight
    occ = max(min_occ, occ)
    if bond_weighting:
      #print 'weight bond',key,bond.weight,occ,
      if sqrt_func:
        bond.weight = bond.weight/math.sqrt(occ)
      else:
        bond.weight = bond.weight/occ*func_factor
      print(bond.weight,occ*func_factor,sqrt_func)
    i_seqs[key[0]]=occ
    i_seqs[key[1]]=occ
    if verbose:
      print('-'*80)
      for atom in hierarchy.atoms():
        if atom.i_seq==key[0]: print(atom.quote())
        if atom.i_seq==key[1]: print(atom.quote())


  if angle_weighting:
    for angle_proxy in restraints_manager.geometry.angle_proxies:
      if angle_proxy.i_seqs[0] not in i_seqs: continue
      if angle_proxy.i_seqs[1] not in i_seqs: continue
      if angle_proxy.i_seqs[2] not in i_seqs: continue
      occ = i_seqs[angle_proxy.i_seqs[1]]
      occ = max(min_occ, occ)
      if sqrt_func:
        angle_proxy.weight = angle_proxy.weight/math.sqrt(occ)
      else:
        angle_proxy.weight = angle_proxy.weight/occ*func_factor

  restraints_manager.geometry.reset_internals()
  #assert 0
  return altloc_weights_object

def run(filename):
  pdb_inp = iotbx.pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_serial()
  update_restraints(hierarchy,
                    restraints_manager,
                    verbose=True,
    )

if __name__=="__main__":
  run(sys.argv[1])
