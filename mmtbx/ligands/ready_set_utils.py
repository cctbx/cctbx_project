from __future__ import absolute_import, division, print_function

import iotbx
from iotbx.pdb import amino_acid_codes as aac
from scitbx.math import dihedral_angle
from mmtbx.ligands.ready_set_basics import construct_xyz
from mmtbx.ligands.ready_set_basics import generate_atom_group_atom_names
from mmtbx.ligands.ready_set_basics import get_hierarchy_h_atom
from mmtbx.ligands.ready_set_basics import get_proton_info
from mmtbx.ligands.hierarchy_utils import _add_atom_to_chain
from mmtbx.hydrogens.specialised_hydrogen_atoms import process_disulphide_hydrogen_atoms
from six.moves import range

get_class = iotbx.pdb.common_residue_names_get_class

def distance2(xyz1, xyz2):
  d2=0
  for k in range(3):
    d2+=(xyz2[k]-xyz1[k])**2
  return d2

def is_n_terminal_residue(residue_group):
  residues = []
  for atom_group in residue_group.atom_groups():
    if atom_group.resname not in residues: residues.append(atom_group.resname)
  assert len(residues)==1
  #if residues[0] in n_terminal_amino_acid_codes: return True
  return False

def _add_atoms_from_chains_to_end_of_hierarchy(hierarchy, chains):
  lookup = {}
  for chain in chains:
    lookup.setdefault(chain.id, [])
    lookup[chain.id].append(chain)
  model = hierarchy.models()[0]
  for i, chain_group in sorted(lookup.items()):
    tc = iotbx.pdb.hierarchy.chain()
    tc.id = i
    for chain in chain_group:
      for rg in chain.residue_groups():
        tc.append_residue_group(rg.detached_copy())
    model.append_chain(tc)

def validate_c_ca_n_for_n_terminal(c,ca,n,bonds):
  if bonds is None: return True
  nb = bonds.get(n.i_seq, None)
  if not nb: return False
  if len(nb) not in [1, 2]: return False
  return True

def validate_c_ca_n_for_c_terminal(c,ca,n,bonds):
  if bonds is None: return True
  cb = bonds.get(c.i_seq, None)
  if not cb: return False
  if len(cb) not in [2]: return False
  return True

def add_n_terminal_hydrogens_to_atom_group(ag,
                                           bonds=None,
                                           use_capping_hydrogens=False,
                                           append_to_end_of_model=False,
                                           retain_original_hydrogens=True,
                                           n_ca_c=None,
                                           verbose=False,
                                          ):
  if verbose and bonds is None:
    print('\n\t%s\n' % 'add_n_terminal_hydrogens_to_atom_group has not been given bonds')
  rc=[]
  if n_ca_c is not None:
    n, ca, c = n_ca_c
  else:
    n = ag.get_atom("N")
    if n is None: return 'no N'
    ca = ag.get_atom("CA")
    if ca is None: return 'no CA'
    c = ag.get_atom("C")
    if c is None: return 'no C'

  if not validate_c_ca_n_for_n_terminal(c, ca, n, bonds): return []

  proton_element, proton_name = get_proton_info(ag)
  atom = ag.get_atom(proton_element) # just so happens that the atom is named H/D
  dihedral=120.
  if atom:
    dihedral = dihedral_angle(sites=[atom.xyz,
                                     n.xyz,
                                     ca.xyz,
                                     c.xyz,
                                   ],
                              deg=True)
  if retain_original_hydrogens: pass
  else:
    if ag.get_atom(proton_name): # maybe needs to be smarter or actually work
      ag.remove_atom(ag.get_atom(proton_name))

  rh3 = construct_xyz(n, 1.0,
                      ca, 109.5,
                      c, dihedral,
                     )
  # this could be smarter
  if proton_element=='H':
    possible = ['H', 'H1', 'H2', 'H3', 'HT1', 'HT2']
  elif proton_element=='D':
    possible = ['D', 'D1', 'D2', 'D3'] #, 'HT1', 'HT2']
  h_count = 0
  d_count = 0
  for h in possible:
    if ag.get_atom(h): h_count+=1
    if ag.get_atom(h.replace('H', 'D')): d_count+=1
  number_of_hydrogens=3
  if verbose: print('starting number_of_hydrogens',number_of_hydrogens)
  if use_capping_hydrogens:
    if verbose: print('use_capping_hydrogens...')
    number_of_hydrogens-=1
    #if ag.atoms()[0].parent().resname=='PRO':
    #  number_of_hydrogens=-1
    #  # should name the hydrogens correctly
  if bonds:
    def _get_atom_with_i_seq(ag, i_seq):
      for atom in ag.atoms():
        if atom.i_seq==i_seq: return atom
      return None
    atoms=ag.atoms()
    i_seqs = atoms.extract_i_seq()
    number_of_heavy=0
    for i_seq in bonds.get(n.i_seq, []):
      if i_seq in i_seqs:
        ba = _get_atom_with_i_seq(ag, i_seq)
        if not ba.element_is_hydrogen():
          number_of_heavy+=1
      # else:
      #   number_of_heavy+=1
    if number_of_heavy==2 and ag.resname!='PRO': # need central PRO-like...
      number_of_hydrogens-=2
    if verbose: print('number_of_heavy', number_of_heavy)
  if verbose: print('number_of_hydrogens',number_of_hydrogens)
  if h_count+d_count>=number_of_hydrogens: return []
  j=0
  for i in range(0, number_of_hydrogens):
    name = " %s%d " % (proton_element, i+1)
    if number_of_hydrogens==1: name = ' %s  ' % proton_element
    if retain_original_hydrogens:
      if i==0 and ag.get_atom(proton_name): continue
      if i==1 and ag.get_atom(proton_name):
        retained = ag.get_atom(proton_name)
        d2 = distance2(retained.xyz, rh3[j])
        if d2<0.5: j+=1
    if ag.get_atom(name.strip()): continue
    if ag.resname=='PRO':
      if i==0:
        continue
    # _new_atom in heirarchy_utils.py
    atom = iotbx.pdb.hierarchy.atom()
    atom.name = name
    atom.element = proton_element
    atom.xyz = rh3[j]
    atom.occ = n.occ
    atom.b = n.b
    atom.segid = ' '*4
    if verbose: print('adding', atom.quote())
    if append_to_end_of_model and i+1==number_of_hydrogens:
      rg = _add_atom_to_chain(atom,
                              ag,
                              icode=n.parent().parent().icode)
      rc.append(rg)
    else:
      ag.append_atom(atom)
    j+=1
    if j==number_of_hydrogens: j=0
  return rc

def add_n_terminal_hydrogens_to_residue_group(residue_group,
                                              bonds=None,
                                              use_capping_hydrogens=False,
                                              append_to_end_of_model=False,
                                              retain_original_hydrogens=True,
                                             ):
  rc=[]
  for atom_group, atoms in generate_atom_group_atom_names(residue_group,
                                                       ['N', 'CA', 'C'],
                                                       verbose=False
                                                       ):
    if None not in [atom_group, atoms]:
      tmp = add_n_terminal_hydrogens_to_atom_group(
        atom_group,
        bonds=bonds,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
        retain_original_hydrogens=retain_original_hydrogens,
        n_ca_c=atoms,
      )
      assert type(tmp)!=type(''), 'not string "%s" %s' % (tmp, type(tmp))
      rc += tmp
  return rc

def add_c_terminal_oxygens_to_atom_group(ag,
                                         bonds=None,
                                         use_capping_hydrogens=False,
                                         append_to_end_of_model=False,
                                         c_ca_n=None,
                                        ):
  #
  # do we need ANISOU
  #
  proton_element, proton_name = get_proton_info(ag)
  rc = []
  atom_name=' OXT'
  atom_element = 'O'
  bond_length=1.231
  if use_capping_hydrogens:
    if ag.get_atom(atom_name.strip()): return []
    atom_name=" %sC " % proton_element
    atom_element=proton_element
    bond_length=1.
  if ag.get_atom(atom_name.strip()): return []
  if c_ca_n is not None:
    c, ca, n = c_ca_n
  else:
    c = ag.get_atom("C")
    if c is None: return
    ca = ag.get_atom("CA")
    if ca is None: return
    n = ag.get_atom("N")
    if n is None: return
  if not validate_c_ca_n_for_c_terminal(c, ca, n, bonds): return []
  atom = ag.get_atom('O')
  dihedral = dihedral_angle(sites=[atom.xyz,
                                   c.xyz,
                                   ca.xyz,
                                   n.xyz,
                                 ],
                            deg=True)
  ro2 = construct_xyz(c, bond_length,
                      ca, 120.,
                      n, dihedral,
                      period=2,
                     )
  oxys = [' O  ', atom_name]
  for i in range(0,2):
    name = oxys[i]
    atom = ag.get_atom(name.strip())
    if atom:
      pass #atom.xyz = ro2[i]
    else:
      atom = iotbx.pdb.hierarchy.atom()
      atom.name = name
      atom.element = atom_element
      atom.occ = c.occ
      atom.b = c.b
      atom.segid = ' '*4
      atom.xyz = ro2[i]
      if append_to_end_of_model:
        chain = _add_atom_to_chain(atom,
                                   ag,
                                   icode=c.parent().parent().icode)
        rc.append(chain)
      else:
        # add the atom to the hierarchy
        ag.append_atom(atom)
  return rc

def add_c_terminal_oxygens_to_residue_group(residue_group,
                                            bonds=None,
                                            use_capping_hydrogens=False,
                                            append_to_end_of_model=False,
                                           ):
  result=[]
  for ag, rc in generate_atom_group_atom_names(residue_group,
                                                       ['C', 'CA', 'N'],
                                                       ):
    if rc is None: continue
    c, ca, n = rc
    tmp = add_c_terminal_oxygens_to_atom_group(
      ag,
      bonds=bonds,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
      c_ca_n = [c, ca, n],
    )
    result += tmp
  return result

def add_main_chain_o_to_atom_group(ag, c_ca_n=None):
  # cetral functuon
  if c_ca_n is not None:
    c, ca, n = c_ca_n
  else:
    c = ag.get_atom("C")
    if c is None: return
    ca = ag.get_atom("CA")
    if ca is None: return
    n = ag.get_atom("N")
    if n is None: return
  atom = ag.get_atom('O')
  dihedral = dihedral_angle(sites=[atom.xyz,
                                   c.xyz,
                                   ca.xyz,
                                   n.xyz,
                                 ],
                            deg=True)
  ro2 = construct_xyz(c, bond_length,
                      ca, 120.,
                      n, dihedral,
                      period=2,
                     )
  oxys = [' O  ', atom_name]
  for i in range(0,2):
    name = oxys[i]
    atom = ag.get_atom(name.strip())
    if atom:
      pass #atom.xyz = ro2[i]
    else:
      atom = iotbx.pdb.hierarchy.atom()
      atom.name = name
      atom.element = atom_element
      atom.occ = c.occ
      atom.b = c.b
      atom.segid = ' '*4
      atom.xyz = ro2[i]
      if append_to_end_of_model:
        chain = _add_atom_to_chain(atom, ag)
        rc.append(chain)
      else:
        # add the atom to the hierarchy
        ag.append_atom(atom)

  #for
  assert atom

def add_main_chain_atoms_to_residue_group(residue_group):
  # I think this needs a three
  assert 0
  for ag, (c, ca, n, o) in generate_atom_group_atom_names(residue_group,
                                                          ['C', 'CA', 'N', 'O'],
                                                          return_Nones=True,
                                                          ):
    if o is None:
      add_main_chain_o_to_atom_group(ag,
                                     c_ca_n = [c, ca, n],
                                    )

def add_main_chain_atoms_to_protein_three(three):
  for ag, (c, ca, n, o) in generate_atom_group_atom_names(three[1],
                                                          ['C', 'CA', 'N', 'O'],
                                                          return_Nones=True,
                                                          ):
    if o is None:
      add_main_chain_o_to_atom_group(ag,
                                     c_ca_n = [c, ca, n],
                                    )

  assert 0

def generate_residues_via_conformer(hierarchy,
                                    backbone_only=False,
                                    verbose=False,
                                    ):
  backbone_asc = hierarchy.atom_selection_cache()
  backbone_sel = backbone_asc.selection("name ca or name c or name n or name o or name cb")
  backbone_hierarchy = hierarchy.select(backbone_sel)
  get_class = iotbx.pdb.common_residue_names_get_class
  loop_hierarchy=hierarchy
  if backbone_only: loop_hierarchy=backbone_hierarchy
  for model in loop_hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if verbose: print('chain: "%s"' % chain.id)
      for conformer in chain.conformers():
        if verbose: print('  conformer: altloc="%s"' % (conformer.altloc))
#        while threes: del threes[0]
#        threes.start=None
#        threes.end=None
#        list_of_threes = []
        for residue in conformer.residues():
          if verbose:
            if residue.resname not in ["HOH"]:
              print('    residue: resname="%s" resid="%s"' % (
                residue.resname, residue.resid()))
          if verbose: print('      residue class : %s' % get_class(residue.resname))
          if get_class(residue.resname) not in ["common_amino_acid",
                                                'modified_amino_acid',
                                              ]:
            # this needs to be moved to cctbx get_class
            #'ETA', # COOH terminal - not in modified
            if residue.resname not in aac.three_letter_l_given_three_letter_d:
              continue
          yield residue

def generate_protein_fragments(hierarchy,
                               geometry,
                               backbone_only=False,
                               use_capping_hydrogens=False,
                               verbose=False,
                               ):
  '''
  Called by qrefine. This should be replaced by the function in hierarchy_utils
  after much testing.
  '''
  from mmtbx.conformation_dependent_library.multi_residue_class import \
    ThreeProteinResidues, RestraintsRegistry
  registry = RestraintsRegistry()
  threes = ThreeProteinResidues(geometry, registry=registry)
  for residue in generate_residues_via_conformer(hierarchy,
                                                 backbone_only=backbone_only,
                                                 verbose=verbose,
                                                 ):
    list.append(threes, residue)
    if verbose: print('THREE',threes)
    sub_unit = threes.provide_second_sub_unit_if_unlinked()
    if verbose: print('THREE, SUBUNIT',threes, sub_unit)
    if sub_unit:
      threes.start = True
      threes.end = True
      yield threes
      threes = sub_unit
  threes.start = True
  threes.end = True
  yield threes

def _hierarchy_into_slots(hierarchy,
                          geometry_restraints_manager,
                          verbose=False,
                          ):
  def _is_linked(residue_group1, residue_group2, bpt):
    if residue_group2 is None: return False
    if residue_group2 is False: return False
    atom_group1 = residue_group1.only_atom_group()
    atom_group2 = residue_group2.only_atom_group()
    c_atom = atom_group1.get_atom('N') # order important
    n_atom = atom_group2.get_atom('C')
    if c_atom and n_atom:
      key = [c_atom.i_seq, n_atom.i_seq]
      bond = bpt.lookup(*tuple(key))
      if bond: return bond
    for i, atom1 in enumerate(atom_group1.atoms()):
      for j, atom2 in enumerate(atom_group2.atoms()):
        key = [atom2.i_seq, atom1.i_seq]
        bond = bpt.lookup(*tuple(key))
        if bond: return bond
    return False

  if verbose: print('_hierarchy_into_slots')
  terminal_entities = ['ACE']
  slots=[]
  #start=18
  assert len(hierarchy.models())==1
  for chain in hierarchy.chains():
    for j, residue_group in enumerate(chain.residue_groups()):
      if verbose: print(j, chain.id, residue_group.id_str())
      protein = True
      for atom_group in residue_group.atom_groups():
        gc_aa = get_class(atom_group.resname) in ["common_amino_acid",
                                                  "modified_amino_acid",
                                                 ]
        # not completely sure about this
        d_aa = atom_group.resname in aac.three_letter_l_given_three_letter_d
        t_aa = atom_group.resname in terminal_entities
        if verbose: print('get_class', gc_aa, 'd-amino', d_aa, 'term. ent', t_aa)
        if(not gc_aa and not d_aa and not t_aa):
          protein=False
          break
      if not protein:
        slots.append(False)
        continue
      if slots:
        bonded = _is_linked(residue_group,
                            slots[-1],
                            geometry_restraints_manager.bond_params_table)
        if getattr(bonded, 'origin_id', None)!=0: # default for peptide links...
          bonded=False
        if bonded:
          pass
        else:
          slots.append(None)
      slots.append(residue_group)
    slots.append(None)
  return slots

def generate_residue_group_with_start_and_end(hierarchy,
                                              geometry_restraints_manager,
                                              # ideal_hierarchy=None,
                                              ):
  # assert not ideal_hierarchy
  slots = _hierarchy_into_slots(hierarchy, geometry_restraints_manager)
  for i in range(len(slots)):
    start=False
    end=False
    if slots[i]:
      if i==0: start=True
      elif not slots[i-1]: start=True
      if i==len(slots)-1: end=True
      elif not slots[i+1]: end=True
      # does not work for chain ends
    else: continue
    residue_group = slots[i]
    yield residue_group, start, end

def add_terminal_hydrogens_old(hierarchy,
                           geometry_restraints_manager,
                           terminate_all_N_terminals=False,
                           terminate_all_C_terminals=False,
                           use_capping_hydrogens=False,
                           append_to_end_of_model=False,
                           verbose=False,
                           ):
  assert 0
  ptr=0 # belt and braces
  additional_hydrogens = []
  for residue_group, start, end in generate_residue_group_with_start_and_end(
    hierarchy,
    geometry_restraints_manager,
    verbose=verbose,
    ):
    if use_capping_hydrogens:
      conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                           residue_group)
    if start:
      ptr+=1
      assert ptr==1
      if is_n_terminal_residue(residue_group):
        rc = None
      else:
        rc = add_n_terminal_hydrogens_to_residue_group(
          residue_group,
          use_capping_hydrogens=use_capping_hydrogens,
          append_to_end_of_model=append_to_end_of_model,
        )
      if rc: additional_hydrogens.append(rc)
    if end:
      ptr-=1
      assert ptr==0
      rc = add_c_terminal_oxygens_to_residue_group(
        residue_group,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
    else:
      pass

def add_terminal_hydrogens_threes(hierarchy,
                                  geometry_restraints_manager,
                                  terminate_all_N_terminals=False,
                                  terminate_all_C_terminals=False,
                                  use_capping_hydrogens=False,
                                  append_to_end_of_model=False,
                                  verbose=False,
                                  ):
  assert 0
  from mmtbx.conformation_dependent_library import generate_protein_threes
  additional_hydrogens= [] #hierarchy_utils.smart_add_atoms()
  hierarchy.show()
  for three in generate_protein_threes(hierarchy,
                                       geometry_restraints_manager,
                                       include_non_linked=True,
                                       backbone_only=False,
                                       include_linked_via_restraints_manager=True,
                                       verbose=verbose,
                                       ):
    bond_params_table = geometry_restraints_manager.bond_params_table

    def get_bonds():
      bonds = {}
      for i, a1 in enumerate(residue_group.atoms()):
        for j, a2 in enumerate(residue_group.atoms()):
          if i>=j: continue
          bond=three.bond_params_table.lookup(a1.i_seq, a2.i_seq)
          if bond:
            bonds[(a1.i_seq, a2.i_seq)] = True
            bonds[(a2.i_seq, a1.i_seq)] = True
      return bonds

    if use_capping_hydrogens:
      for i in range(len(three)):
        residue_group=three.get_residue_group_from_hierarchy(hierarchy, i)
        rc = conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                                  residue_group)
      #assert not rc, '%s' % rc
    if three.start:
      residue_group=three.get_residue_group_from_hierarchy(hierarchy, 0)
      rc = add_n_terminal_hydrogens_to_residue_group(
        residue_group,
        bonds=get_bonds(),
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
    if three.end:
      residue_group=three.get_residue_group_from_hierarchy(hierarchy, 2)
      rc = add_c_terminal_oxygens_to_residue_group(
        residue_group,
        bonds=get_bonds(),
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
  return additional_hydrogens

def _check_for(atom_holder, resname, name):
  for atom in atom_holder.atoms():
    if atom.name.strip()!=name.strip(): continue
    p1 = atom.parent()
    if p1.resname.strip()!=resname.strip(): continue
    print(atom.quote(), resname, name)
    return True
  return False

def add_terminal_hydrogens_via_residue_groups(hierarchy,
                                              geometry_restraints_manager,
                                              terminate_all_N_terminals=False,
                                              terminate_all_C_terminals=False,
                                              use_capping_hydrogens=False,
                                              append_to_end_of_model=False,
                                              retain_original_hydrogens=True,
                                              verbose=False,
                                              ):
  # assert not _check_for(hierarchy, '0QS', 'HC')
  from mmtbx.ligands.hierarchy_utils import get_bonds_as_dict
  bonds=get_bonds_as_dict(geometry_restraints_manager)
  additional_hydrogens = []
  for residue_group, start, end in generate_residue_group_with_start_and_end(
    hierarchy,
    geometry_restraints_manager,
    ):
    if use_capping_hydrogens:
      # conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
      #                                      residue_group)
      process_disulphide_hydrogen_atoms(geometry_restraints_manager,
                                        residue_group)
    if start:
      # ptr+=1
      # assert ptr==1
      if is_n_terminal_residue(residue_group):
        rc = None
      else:
        rc = add_n_terminal_hydrogens_to_residue_group(
          residue_group,
          bonds=bonds,
          use_capping_hydrogens=use_capping_hydrogens,
          append_to_end_of_model=append_to_end_of_model,
          retain_original_hydrogens=retain_original_hydrogens,
        )
      if rc: additional_hydrogens.append(rc)
    if end:
      # ptr-=1
      # assert ptr==0
      rc = add_c_terminal_oxygens_to_residue_group(
        residue_group,
        bonds=bonds,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
    else:
      pass
  return additional_hydrogens

def add_terminal_hydrogens(hierarchy,
                           geometry_restraints_manager,
                           terminate_all_N_terminals=False,
                           terminate_all_C_terminals=False,
                           use_capping_hydrogens=False,
                           append_to_end_of_model=False,
                           retain_original_hydrogens=True,
                           verbose=False,
                           ):
  """Adds hydrogen atoms to a macromolecular model

  Args:
      hierarchy (pdb_hierarch): The model hierarchy
      geometry_restraints_manager (geometry restraints manager): Matching
        geometry restaints manager for hierarch
      terminate_all_N_terminals (bool, optional): Description
      terminate_all_C_terminals (bool, optional): Description
      use_capping_hydrogens (bool, optional): Description
      append_to_end_of_model (bool, optional): Description
      verbose (bool, optional): Description
  """
  additional_hydrogens = add_terminal_hydrogens_via_residue_groups(
    hierarchy,
    geometry_restraints_manager,
    terminate_all_N_terminals=terminate_all_N_terminals,
    terminate_all_C_terminals=terminate_all_C_terminals,
    use_capping_hydrogens=use_capping_hydrogens,
    append_to_end_of_model=append_to_end_of_model,
    retain_original_hydrogens=retain_original_hydrogens,
    verbose=verbose,
    )

  if append_to_end_of_model and additional_hydrogens:
    tmp = []
    for group in additional_hydrogens:
      for chain in group:
        tmp.append(chain)
    _add_atoms_from_chains_to_end_of_hierarchy(hierarchy, tmp)

def delete_charged_n_terminal_hydrogens(hierarchy):
  for ag in hierarchy.atom_groups():
    if get_class(ag.resname) not in ['common_amino_acid',
                                     'modified_amino_acid']: continue
    h1 = ag.get_atom('H1')
    if not h1: continue
    h2 = ag.get_atom('H2')
    if not h2: continue
    h3 = ag.get_atom('H3')
    if not h3: continue
    ag.remove_atom(h3)

def add_water_hydrogen_atoms_simple(hierarchy, log=None):
  displacements = [
    [ 1.0,  0.0,  0.0],
    [-0.7, -0.7,  0.0],
    [ 0.0,  1.0,  0.0],
    [-0.7, -0.7,  0.0],
    [ 0.0,  0.0,  1.0],
    [ 0.0, -0.7, -0.7],
    ]
  for atom_group in hierarchy.atom_groups():
    if atom_group.resname in ['HOH', 'DOD']:
      proton_element='H'
      if atom_group.resname=='DOD': proton_element='D'
      if len(atom_group.atoms())==3: continue
      names = [atom.name for atom in atom_group.atoms()]
      o_atom = atom_group.atoms()[0]
      for name in [' %s1 ' % proton_element, ' %s2 ' % proton_element]:
        if name not in names:
          j=0
          min_d2=0
          while min_d2<.9:
            xyz = (o_atom.xyz[0]+displacements[j][0],
                   o_atom.xyz[1]+displacements[j][1],
                   o_atom.xyz[2]+displacements[j][2])
            min_d2 = 1e9
            for atom in atom_group.atoms():
              d2 = distance2(xyz, atom.xyz)
              min_d2=min(d2, min_d2)
            j+=1
          h = get_hierarchy_h_atom(name, xyz, o_atom, proton_element=proton_element)
          atom_group.append_atom(h)
      # for atom in atom_group.atoms(): print(atom.format_atom_record())

def add_main_chain_atoms(hierarchy,
                         geometry_restraints_manager,
                         verbose=False,
                         ):
  from mmtbx.conformation_dependent_library import generate_protein_threes
  for three in generate_protein_threes(hierarchy,
                                       geometry_restraints_manager,
                                       verbose=verbose,
                                       ):
    print(three)
    add_main_chain_atoms_to_protein_three(three)
  assert 0

def perdeuterate_model_ligands(hierarchy,
                               verbose=False,
                               ):
  if verbose: hierarchy.show()
  for model in hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if verbose: print('chain: "%s"' % chain.id)
      for residue_group in chain.residue_groups():
        if verbose: print('  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode))
        altlocs = []
        for atom_group in residue_group.atom_groups():
          altlocs.append(atom_group.altloc)
        for atom_group in residue_group.atom_groups():
          if verbose: print('  atom_group: resname="%s" altloc="%s"' % (
            atom_group.resname, atom_group.altloc))
          for atom in atom_group.atoms():
            if atom.element.strip() in ["H"]:
              atom.element = " D"
              atom.name = atom.name.replace("H","D", 1)

  hierarchy.atoms().reset_serial()
  if verbose: hierarchy.show()
  return hierarchy

def main_hydrogen(model,
                  terminate_all_N_terminals=False,
                  terminate_all_C_terminals=False,
                  cap_all_terminals=False,
                  append_to_end_of_model=False,
                  verbose=False):
  """Main loop for ReadySet!

    Add specialised hydrogens not done by Reduce
    Curate ligands

  Args:
      model (TYPE): Model class with all the model info
      terminate_all_N_terminals (bool, optional): Description
      terminate_all_C_terminals (bool, optional): Description
      cap_all_terminals (bool, optional): Description
      append_to_end_of_model (bool, optional): Description
      verbose (bool, optional): Description

  Returns:
      TYPE: Description

  """
  #
  # make sure all pdb_interpretation parameters have been done
  #
  hierarchy = model.get_hierarchy()
  # should be automatic
  model.process(make_restraints=True)
  geometry_restraints_manager = model.get_restraints_manager().geometry
  atoms = hierarchy.atoms()

  #ready_set_utils.add_main_chain_atoms(hierarchy, geometry_restraints_manager)
  if (terminate_all_N_terminals or
      terminate_all_C_terminals or
      cap_all_terminals
      ):
    add_terminal_hydrogens(
      hierarchy,
      geometry_restraints_manager,
      terminate_all_N_terminals=terminate_all_N_terminals,
      terminate_all_C_terminals=terminate_all_C_terminals,
      use_capping_hydrogens=cap_all_terminals,
      append_to_end_of_model=append_to_end_of_model,
      verbose=False,
      )
  return
