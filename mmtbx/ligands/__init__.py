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
  geometry_restraints_manager = model.get_restraints_manager().geometry
  atoms = hierarchy.atoms()

  from mmtbx.ligands import ready_set_utils
  #ready_set_utils.add_main_chain_atoms(hierarchy, geometry_restraints_manager)
  if (terminate_all_N_terminals or
      terminate_all_C_terminals or
      cap_all_terminals
      ):
    ready_set_utils.add_terminal_hydrogens(
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
        for ag, (n, ca, c) in hierarchy_utils.generate_atom_group_atom_names(
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

