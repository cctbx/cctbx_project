from __future__ import absolute_import, division, print_function

import iotbx
from iotbx.pdb import amino_acid_codes as aac
from scitbx.math import dihedral_angle
from mmtbx.ligands.ready_set_basics import construct_xyz
from mmtbx.ligands.ready_set_basics import generate_atom_group_atom_names
from mmtbx.hydrogens.specialised_hydrogen_atoms import conditional_add_cys_hg_to_atom_group
from six.moves import range

get_class = iotbx.pdb.common_residue_names_get_class

def is_n_terminal_residue(residue_group):
  residues = []
  for atom_group in residue_group.atom_groups():
    if atom_group.resname not in residues: residues.append(atom_group.resname)
  assert len(residues)==1
  #if residues[0] in n_terminal_amino_acid_codes: return True
  return False

def _add_atom_to_chain(atom, ag):
  rg = _add_atom_to_residue_group(atom, ag)
  chain = ag.parent().parent()
  tc = iotbx.pdb.hierarchy.chain()
  tc.id = chain.id
  tc.append_residue_group(rg)
  return tc

def _add_atom_to_residue_group(atom, ag):
  tag = iotbx.pdb.hierarchy.atom_group()
  tag.resname = ag.resname
  tag.append_atom(atom)
  rg = iotbx.pdb.hierarchy.residue_group()
  rg.resseq = ag.parent().resseq
  rg.append_atom_group(tag)
  for i, c in enumerate(letters):
    if c==ag.parent().parent().id:
      break
  atom.tmp = i
  return rg

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

def add_n_terminal_hydrogens_to_atom_group(ag,
                                           bonds=None,
                                           use_capping_hydrogens=False,
                                           append_to_end_of_model=False,
                                           retain_original_hydrogens=True,
                                           n_ca_c=None,
                                          ):
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

  if bonds:
    ns = []
    for key in bonds:
      if n.i_seq == key[0]:
        ns.append(key)
    if len(ns)>=3: return rc

  atom = ag.get_atom('H')
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
    if ag.get_atom("H"): # maybe needs to be smarter or actually work
      ag.remove_atom(ag.get_atom('H'))
  #if use_capping_hydrogens and 0:
  #  for i, atom in enumerate(ag.atoms()):
  #    if atom.name == ' H3 ':
  #      ag.remove_atom(i)
  #      break
  # add H1
  rh3 = construct_xyz(n, 0.9,
                      ca, 109.5,
                      c, dihedral,
                     )
  # this could be smarter
  possible = ['H', 'H1', 'H2', 'H3', 'HT1', 'HT2']
  h_count = 0
  for h in possible:
    if ag.get_atom(h): h_count+=1
  number_of_hydrogens=3
  if use_capping_hydrogens:
    number_of_hydrogens-=1
    #if ag.atoms()[0].parent().resname=='PRO':
    #  number_of_hydrogens=-1
    #  # should name the hydrogens correctly
  if h_count>=number_of_hydrogens: return []
  for i in range(0, number_of_hydrogens):
    name = " H%d " % (i+1)
    if retain_original_hydrogens:
      if i==0 and ag.get_atom('H'): continue
    if ag.get_atom(name.strip()): continue
    if ag.resname=='PRO':
      if i==0:
        continue
    atom = iotbx.pdb.hierarchy.atom()
    atom.name = name
    atom.element = "H"
    atom.xyz = rh3[i]
    atom.occ = n.occ
    atom.b = n.b
    atom.segid = ' '*4
    if append_to_end_of_model and i+1==number_of_hydrogens:
      rg = _add_atom_to_chain(atom, ag)
      rc.append(rg)
    else:
      ag.append_atom(atom)
  return rc

def add_n_terminal_hydrogens_to_residue_group(residue_group,
                                              bonds=None,
                                              use_capping_hydrogens=False,
                                              append_to_end_of_model=False,
                                             ):
  rc=[]
  for atom_group, atoms in generate_atom_group_atom_names(residue_group,
                                                       ['N', 'CA', 'C'],
                                                       verbose=False
                                                       ):
    if None not in [atom_group, atoms]:
#  for ag, (n, ca, c) in generate_atom_group_atom_names(residue_group,
#                                                       ['N', 'CA', 'C'],
#                                                       ):
      tmp = add_n_terminal_hydrogens_to_atom_group(
        atom_group,
        bonds=bonds,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
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
  rc = []
  atom_name=' OXT'
  atom_element = 'O'
  bond_length=1.231
  if use_capping_hydrogens:
    if ag.get_atom(atom_name.strip()): return []
    atom_name=" HC "
    atom_element="H"
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
  return rc

def add_c_terminal_oxygens_to_residue_group(residue_group,
                                            bonds=None,
                                            use_capping_hydrogens=False,
                                            append_to_end_of_model=False,
                                           ):
  rc=[]
  for ag, (c, ca, n) in generate_atom_group_atom_names(residue_group,
                                                       ['C', 'CA', 'N'],
                                                       ):
    tmp = add_c_terminal_oxygens_to_atom_group(
      ag,
      bonds=bonds,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
      c_ca_n = [c, ca, n],
    )
    rc += tmp
  return rc

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
    # print ag.id_str()
    # try: print c.quote(),ca.quote(),n.quote(),o.quote()
    # except: print c,ca,n,o
    if o is None:
      add_main_chain_o_to_atom_group(ag,
                                     c_ca_n = [c, ca, n],
                                    )

def add_main_chain_atoms_to_protein_three(three):
  for ag, (c, ca, n, o) in generate_atom_group_atom_names(three[1],
                                                          ['C', 'CA', 'N', 'O'],
                                                          return_Nones=True,
                                                          ):
    print(ag.id_str())
    if o is None:
      add_main_chain_o_to_atom_group(ag,
                                     c_ca_n = [c, ca, n],
                                    )

  assert 0

# def generate_residues_via_conformer(hierarchy,
#                                     backbone_only=False,
#                                     verbose=False,
#                                     ):
#   assert 0
#   backbone_asc = hierarchy.atom_selection_cache()
#   backbone_sel = backbone_asc.selection("name ca or name c or name n or name o or name cb")
#   backbone_hierarchy = hierarchy.select(backbone_sel)
#   get_class = iotbx.pdb.common_residue_names_get_class
#   loop_hierarchy=hierarchy
#   if backbone_only: loop_hierarchy=backbone_hierarchy
#   for model in loop_hierarchy.models():
#     if verbose: print 'model: "%s"' % model.id
#     for chain in model.chains():
#       if verbose: print 'chain: "%s"' % chain.id
#       for conformer in chain.conformers():
#         if verbose: print '  conformer: altloc="%s"' % (
#           conformer.altloc)
# #        while threes: del threes[0]
# #        threes.start=None
# #        threes.end=None
# #        list_of_threes = []
#         for residue in conformer.residues():
#           if verbose:
#             if residue.resname not in ["HOH"]:
#               print '    residue: resname="%s" resid="%s"' % (
#                 residue.resname, residue.resid())
#           if verbose: print '      residue class : %s' % get_class(residue.resname)
#           if get_class(residue.resname) not in ["common_amino_acid",
#                                                 'modified_amino_acid',
#                                               ]:
#             # this needs to be moved to cctbx get_class
#             #'ETA', # COOH terminal - not in modified
#             if residue.resname not in aac.three_letter_l_given_three_letter_d:
#               continue
#           yield residue

# def generate_protein_fragments(hierarchy,
#                                geometry,
#                                backbone_only=False,
#                                use_capping_hydrogens=False,
#                                verbose=False,
#                                ):
#   assert 0
#   from mmtbx.conformation_dependent_library.multi_residue_class import \
#     ThreeProteinResidues, RestraintsRegistry
#   registry = RestraintsRegistry()
#   threes = ThreeProteinResidues(geometry, registry=registry)
#   for residue in generate_residues_via_conformer(hierarchy,
#                                                  backbone_only=backbone_only,
#                                                  verbose=verbose,
#                                                  ):
#     list.append(threes, residue)
#     if verbose: print 'THREE',threes
#     sub_unit = threes.provide_second_sub_unit_if_unlinked()
#     if verbose: print 'THREE, SUBUNIT',threes, sub_unit
#     if sub_unit:
#       threes.start = True
#       threes.end = True
#       yield threes
#       threes = sub_unit
#   threes.start = True
#   threes.end = True
#   yield threes

def _hierarchy_into_slots(hierarchy,
                          verbose=False,
                          ):
  if verbose: print('_hierarchy_into_slots')
  slots=[]
  #start=18
  assert len(hierarchy.models())==1
  for chain in hierarchy.chains():
    for residue_group in chain.residue_groups():
      protein = True
      for atom_group in residue_group.atom_groups():
        if(get_class(atom_group.resname) not in ["common_amino_acid",
                                                 "modified_amino_acid",
                                                ] and
            atom_group.resname not in aac.three_letter_l_given_three_letter_d):
                                      # not completely sure about this
          protein=False
          break
      if not protein:
        slots.append(False)
        continue
      #residue_group_atom1_quote = residue_group.atoms()[0].quote()[start:]
      #for residue_group in hierarchy.residue_groups():
      #  if residue_group.atoms()[0].quote()[start:]==oresidue_group_atom1_quote:
      slots.append(residue_group)
      #    break
      #else:
      #  slots.append(0)
    slots.append(None)
  return slots

def generate_residue_group_with_start_and_end(hierarchy,
                                              geometry_restraints_manager,
                                              ideal_hierarchy=None,
                                              verbose=False,
                                              ):
  assert not ideal_hierarchy
  slots = _hierarchy_into_slots(hierarchy)
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
  print(additional_hydrogens)

def add_terminal_hydrogens_threes(hierarchy,
                                  geometry_restraints_manager,
                                  terminate_all_N_terminals=False,
                                  terminate_all_C_terminals=False,
                                  use_capping_hydrogens=False,
                                  append_to_end_of_model=False,
                                  verbose=False,
                                  ):
  from mmtbx.conformation_dependent_library import generate_protein_threes
  additional_hydrogens= [] #hierarchy_utils.smart_add_atoms()
  for three in generate_protein_threes(hierarchy,
                                       geometry_restraints_manager,
                                       #include_non_linked=False,
                                       backbone_only=False,
                                       include_linked_via_restraints_manager=True,
                                       verbose=verbose,
                                       ):
    # print three
    #print dir(three)
    # print geometry_restraints_manager
    #print dir(geometry_restraints_manager)
    bond_params_table = geometry_restraints_manager.bond_params_table
    # print bond_params_table
    #print dir(bond_params_table)
    # print 'use_capping_hydrogens',use_capping_hydrogens

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

def add_terminal_hydrogens(hierarchy,
                           geometry_restraints_manager,
                           terminate_all_N_terminals=False,
                           terminate_all_C_terminals=False,
                           use_capping_hydrogens=False,
                           append_to_end_of_model=False,
                           verbose=False,
                           ):
  additional_hydrogens = add_terminal_hydrogens_threes(
    hierarchy,
    geometry_restraints_manager,
    terminate_all_N_terminals=terminate_all_N_terminals,
    terminate_all_C_terminals=terminate_all_C_terminals,
    use_capping_hydrogens=use_capping_hydrogens,
    append_to_end_of_model=append_to_end_of_model,
    verbose=verbose,
    )

  if append_to_end_of_model and additional_hydrogens:
    tmp = []
    for group in additional_hydrogens:
      for chain in group:
        tmp.append(chain)
    _add_atoms_from_chains_to_end_of_hierarchy(hierarchy, tmp)

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

# def junk():
#   from mmtbx.conformation_dependent_library import generate_protein_threes
#   threes = []
#   for three in generate_protein_threes(hierarchy,
#                                        geometry_restraints_manager,
#                                        backbone_only=False,
#                                        #use_capping_hydrogens=use_capping_hydrogens,
#                                        ):
#     if verbose: print '...',three
#     if not len(three): continue
#     assert three.are_linked()
#     print '---',three
#     if three.start: yield three[0], three.start, False
#     yield three[1], False, False
#     if three.end: yield three[2], False, three.end
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

  assert 0


  n_done = []
  for three in hierarchy_utils.generate_protein_fragments(
    hierarchy,
    geometry_restraints_manager,
    backbone_only=False,
    #use_capping_hydrogens=use_capping_hydrogens,
    ):
    if verbose: print(three)
    if len(three)==1: continue
    for i, residue in enumerate(three):
      if not i: continue
      # this may not be necessary with the new threes
      residue = hierarchy_utils.get_residue_group(residue, atoms)
      h = hierarchy_utils.get_atom_from_residue_group(residue, 'H')
      if h is None:
        assert 0
        for ag, (n, ca, c) in ready_set_basics.generate_atom_group_atom_names(
            residue,
            ['N', 'CA', 'C'],
        ):
          if ag.resname in ['PRO']: continue
          if n in n_done: continue
          n_done.append(n)
          dihedral = 0
          rh3 = general_utils.construct_xyz(n, 0.9,
                                            ca, 109.5,
                                            c, dihedral,
          )
          atom = create_atom(' H  ', 'H', rh3[0], n)
          # adding to atom_group
          # need to add to geometry_restraints_manager
          ag.append_atom(atom)
          if verbose: print(atom.quote())
          assert ag.resname!='PRO'

