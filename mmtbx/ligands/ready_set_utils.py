from __future__ import absolute_import, division, print_function
import iotbx
from iotbx.pdb import amino_acid_codes as aac
from scitbx.math import dihedral_angle
import six
from mmtbx.ligands.ready_set_basics import construct_xyz

get_class = iotbx.pdb.common_residue_names_get_class

def is_n_terminal_residue(residue_group):
  residues = []
  for atom_group in residue_group.atom_groups():
    if atom_group.resname not in residues: residues.append(atom_group.resname)
  assert len(residues)==1
  #if residues[0] in n_terminal_amino_acid_codes: return True
  return False

def generate_atom_group_atom_names(rg, names):
  '''
  Generate all alt. loc. groups of names
  '''
  atom_groups = rg.atom_groups()
  atom_altlocs = {}
  for ag in atom_groups:
    for atom in ag.atoms():
      atom_altlocs.setdefault(atom.parent().altloc, [])
      atom_altlocs[atom.parent().altloc].append(atom)
  if len(atom_altlocs)>1 and '' in atom_altlocs:
    for key in atom_altlocs:
      if key=='': continue
      for atom in atom_altlocs['']:
        atom_altlocs[key].append(atom)
    del atom_altlocs['']
  for key, value in six.iteritems(atom_altlocs):
    atoms=[]
    for name in names:
      for atom in value:
        if atom.name.strip()==name.strip():
          atoms.append(atom)
          break
      else:
        assert 0
    yield atoms[0].parent(), atoms

def _add_hydrogens_to_atom_group_using_bad(ag,
                                           atom_name,
                                           atom_element,
                                           bond_atom,
                                           angle_atom,
                                           dihedral_atom,
                                           bond_length,
                                           angle,
                                           dihedral,
                                           append_to_end_of_model=False,
                                           ):
  rc = []
  if ag.get_atom(atom_name.strip()): return []
  if type(bond_atom)==type(''):
    ba = ag.get_atom(bond_atom.strip())
    #print bond_atom,ba.quote()
    if ba is None: return
  else: ba = bond_atom
  if type(angle_atom)==type(''):
    aa = ag.get_atom(angle_atom.strip())
    #print angle_atom,aa.quote()
    if aa is None: return
  else: aa = angle_atom
  if type(dihedral_atom)==type(''):
    da = ag.get_atom(dihedral_atom.strip())
    #print dihedral_atom, da.quote()
    if da is None: return
  else: da = dihedral_atom
  ro2 = construct_xyz(ba, bond_length,
                      aa, angle,
                      da, dihedral,
                      period=1,
                     )
  atom = iotbx.pdb.hierarchy.atom()
  atom.name = atom_name
  atom.element = atom_element
  atom.occ = ba.occ
  atom.b = ba.b
  # altloc???
  atom.hetero = ba.hetero
  atom.segid = ' '*4
  atom.xyz = ro2[0]
  if append_to_end_of_model:
    chain = _add_atom_to_chain(atom, ag)
    rc.append(chain)
  else:
    ag.append_atom(atom)
  return rc

def add_cys_hg_to_atom_group(ag,
                             append_to_end_of_model=False,
                            ):
  #
  # do we need ANISOU
  #
  rc = _add_hydrogens_to_atom_group_using_bad(
    ag,
    ' HG ',
    'H',
    'SG',
    'CB',
    'CA',
    1.2,
    120.,
    160.,
    append_to_end_of_model=append_to_end_of_model,
   )
  return rc

def add_cys_hg_to_residue_group(rg,
                                append_to_end_of_model=False,
                               ):
  rc=[]
  for ag in rg.atom_groups():
    if ag.resname not in ['CYS']: continue
    rc += add_cys_hg_to_atom_group(
      ag,
      append_to_end_of_model=append_to_end_of_model,
    )
  return rc

def conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                         rg,
                                         ):
  # could be more general to include other disulphide amino acids
  resnames = []
  for ag in rg.atom_groups():
    resnames.append(ag.resname)
  if 'CYS' not in resnames: return -1
  sgs = []
  for atom in rg.atoms():
    if atom.name.strip()=='SG' and atom.parent().resname=='CYS':
      sgs.append(atom.i_seq)
  assert len(sgs) in [0, 1]
  sg_bonds = []
  if sgs:
    for bond in geometry_restraints_manager.get_all_bond_proxies():
      if not hasattr(bond, 'get_proxies_with_origin_id'): continue
      for p in bond.get_proxies_with_origin_id():
        assert p.origin_id==0
        if sgs[0] in p.i_seqs:
          sg_bonds.append(p.i_seqs)
  if len(sg_bonds)==1:
    add_cys_hg_to_residue_group(rg)

def add_n_terminal_hydrogens_to_atom_group(ag,
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
                                              use_capping_hydrogens=False,
                                              append_to_end_of_model=False,
                                             ):
  rc=[]
  for ag, (n, ca, c) in generate_atom_group_atom_names(residue_group,
                                                       ['N', 'CA', 'C'],
                                                       ):
    tmp = add_n_terminal_hydrogens_to_atom_group(
      ag,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
      n_ca_c=[n,ca,c],
    )
    assert type(tmp)!=type(''), 'not string "%s" %s' % (tmp, type(tmp))
    rc += tmp
  return rc

def add_c_terminal_oxygens_to_atom_group(ag,
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
                                            use_capping_hydrogens=False,
                                            append_to_end_of_model=False,
                                          ):
  rc=[]
  for ag, (c, ca, n) in generate_atom_group_atom_names(residue_group,
                                                       ['C', 'CA', 'N'],
                                                       ):
    tmp = add_c_terminal_oxygens_to_atom_group(
      ag,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
      c_ca_n = [c, ca, n],
    )
    rc += tmp
  return rc

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

def add_terminal_hydrogens(hierarchy,
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
