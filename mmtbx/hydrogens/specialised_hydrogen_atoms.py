from __future__ import absolute_import, division, print_function
from scitbx.math import dihedral_angle
from mmtbx.ligands.ready_set_basics import construct_xyz
from mmtbx.ligands.ready_set_basics import generate_atom_group_atom_names
from mmtbx.ligands.ready_set_basics import get_proton_info
from mmtbx.ligands.hierarchy_utils import new_atom_with_inheritance
from mmtbx.ligands.hierarchy_utils import add_hydrogens_to_atom_group_using_bad

from cctbx.geometry_restraints.linking_class import linking_class
origin_ids = linking_class()

def _generate_bonds_with_origin_ids_in_list(bond_proxies, specific_origin_ids=None):
  assert specific_origin_ids
  for specific_origin_id in specific_origin_ids:
    for p in bond_proxies.get_proxies_with_origin_id(specific_origin_id):
      yield p

def add_side_chain_acid_hydrogens_to_atom_group(atom_group,
                                                anchors=None,
                                                configuration_index=0,
                                                bond_length=0.95,
                                                element='H',
                                                ):
  """Add hydrogen atoms to side-chain acid in place

  Args:
      atom_group (TYPE): Atom group
      anchors (None, optional): Atoms that specify the acids moeity
      configuration_index (int, optional): Configuration to return

  """
  assert element in ['H', 'D']
  c, o1, o2 = anchors
  if configuration_index>=2:
    tmp = o1.name
    o1.name = o2.name
    o2.name = tmp
    tmp = o1
    o1 = o2
    o2 = tmp
    configuration_index=configuration_index%2
  if o2.name==' OD2':
    name = ' HD2'
    atom = atom_group.get_atom('CB')
  elif o2.name==' OE2':
    name = ' HE2'
    atom = atom_group.get_atom('CG')
  else: assert 0
  if element=='D': name = name.replace('H', 'D')
  dihedral = dihedral_angle(sites=[atom.xyz,
                                   c.xyz,
                                   o1.xyz,
                                   o2.xyz,
                                 ],
                            deg=True)
  ro2 = construct_xyz(o2, bond_length,
                      c, 120.,
                      o1, dihedral,
                      period=2,
                     )
  i = configuration_index
  atom = atom_group.get_atom(name.strip())
  if atom:
    pass #atom.xyz = ro2[i]
  else:
    atom = new_atom_with_inheritance(name, element, ro2[i], o2)
    atom_group.append_atom(atom)

def add_side_chain_acid_hydrogens_to_residue_group(residue_group,
                                                   configuration_index=0,
                                                   element='H',
                                                   ):
  """Adds hydrogen atoms to side-chain acid.

  Args:
      residue_group (TYPE): Specific residue group
  """
  def _get_atom_names(residue_group):
    assert len(residue_group.atom_groups())==1
    atom_group = residue_group.atom_groups()[0]
    lookup = {'ASP' : ['CG', 'OD1', 'OD2'],
              'GLU' : ['CD', 'OE1', 'OE2'],
    }
    return lookup.get(atom_group.resname, [])
  #
  if element=='H': bond_length=0.95
  elif element=='D': bond_length=1.00
  else: assert 0
  atoms = _get_atom_names(residue_group)
  for atom_group, atoms in generate_atom_group_atom_names(residue_group,
                                                          atoms,
                                                          ):
    if atom_group is None: continue
    tmp = add_side_chain_acid_hydrogens_to_atom_group(
      atom_group,
      # append_to_end_of_model=append_to_end_of_model,
      anchors = atoms,
      configuration_index=configuration_index,
      bond_length=bond_length,
      element=element,
    )

def add_side_chain_acid_hydrogens(hierarchy,
                                  configuration_index=0,
                                  element='H',
                                  ):
  """Add hydrogen atoms to every side-chain acid (ASP and GLU). Not very
  useful as adding to a single residue group (below) would be more prectical.

  Args:
      hierarchy (TYPE): Model hierarchy
      configuration_index (int, optional): Defaults to zero. Determines which
        of the four configurations the added hydrogen will be:
          0 - Current Ox2 gets Hx2 (x=D,E) pointing out
          1 - Current Ox2 gets Hx2 (x=D,E) pointing in
          2 - Current Ox1 gets swapped with Ox2, gets Hx2 (x=D,E) pointing out
          3 - Current Ox1 gets swapped with Ox2, gets Hx2 (x=D,E) pointing in
  """
  for residue_group in hierarchy.residue_groups():
    for atom_group in residue_group.atom_groups():
      if atom_group.resname in ['ASP', 'GLU']:
        add_side_chain_acid_hydrogens_to_residue_group(
          residue_group,
          configuration_index=configuration_index,
          element=element,
          )

#
# isolated CYS need HG
#
def add_cys_hg_to_atom_group(atom_group,
                             append_to_end_of_model=False,
                             element='H',
                             ):
  """Adds hydrogen to CYS

  Args:
      atom_group (TYPE): atom_group in hirarchy
      append_to_end_of_model (bool, optional): Some programs like the additional
        atoms added at end of PDB

  Returns:
      TYPE: New chains, if any
  """
  assert element in ['H', 'D']
  rc = add_hydrogens_to_atom_group_using_bad(
    atom_group,
    ' HG '.replace('H', element),
    element,
    'SG',
    'CB',
    'CA',
    1.2,
    120.,
    160.,
    append_to_end_of_model=append_to_end_of_model,
   )
  return rc

def add_cys_hg_to_residue_group(residue_group,
                                append_to_end_of_model=False,
                                element='H',
                               ):
  rc=[]
  for atom_group in residue_group.atom_groups():
    if atom_group.resname not in ['CYS']: continue
    rc += add_cys_hg_to_atom_group(
      atom_group,
      append_to_end_of_model=append_to_end_of_model,
      element=element,
    )
  return rc

def conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                         residue_group,
                                         element='H',
                                         append_to_end_of_model=False,
                                         ):
  """Adds HG atom to CYS if no disulfur bridge

  Args:
      geometry_restraints_manager (TYPE): GRM
      residue_group (TYPE): CYS residue group
  """
  # could be more general to include other disulphide amino acids
  resnames = []
  for atom_group in residue_group.atom_groups():
    resnames.append(atom_group.resname)
  if 'CYS' not in resnames: return None
  sgs = []
  for atom in residue_group.atoms():
    if atom.name.strip()=='SG' and atom.parent().resname=='CYS':
      sgs.append(atom.i_seq)
  assert len(sgs) in [0, 1]
  sg_bonds = []
  if sgs:
    specific_origin_ids = [origin_ids.get_origin_id('SS BOND'),
                           origin_ids.get_origin_id('metal coordination'),
                           origin_ids.get_origin_id('Misc. bond'),
      ]
    for bond in geometry_restraints_manager.get_all_bond_proxies():
      if not hasattr(bond, 'get_proxies_with_origin_id'): continue
      for p in _generate_bonds_with_origin_ids_in_list(bond, specific_origin_ids):
        assert p.origin_id in specific_origin_ids
        if sgs[0] in p.i_seqs:
          sg_bonds.append(p.i_seqs)
  rc = []
  if len(sg_bonds)==0:
    rc += add_cys_hg_to_residue_group(residue_group,
                                      element=element,
                                      append_to_end_of_model=append_to_end_of_model,
                                      )
  return rc

def add_disulfur_hydrogen_atoms(geometry_restraints_manager,
                                hierarchy,
                                element='H',
                                ):
  """Example of usage

  """
  for residue_group in hierarchy.residue_groups():
    resnames=[]
    for atom_group in residue_group.atom_groups():
      resnames.append(atom_group.resname)
    if 'CYS' in resnames:
      rc = conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                                residue_group,
                                                element=element,
                                                )

def remove_cys_hg_from_residue_group(rg):
  proton_element, proton_name = get_proton_info(rg)
  for ag in rg.atom_groups():
    if ag.resname not in ['CYS']: continue
    for atom in ag.atoms():
      if atom.name==' %sG ' % proton_element:
        ag.remove_atom(atom)
        break

def generate_bonded_i_seqs(geometry_restraints_manager, rg, j_seq):
  def _not_j_seq(j_seq, i_seqs):
    i_seqs.remove(j_seq)
    return i_seqs[0]
  bonds = []
  for bond in geometry_restraints_manager.pair_proxies().bond_proxies.simple:
    if j_seq in bond.i_seqs:
      bonds.append(_not_j_seq(j_seq, list(bond.i_seqs)))
  for bond in geometry_restraints_manager.pair_proxies().bond_proxies.asu:
    if j_seq in [bond.i_seq, bond.j_seq]:
      bonds.append(_not_j_seq(j_seq, [bond.i_seq, bond.j_seq]))
  return bonds

def conditional_remove_cys_hg_to_atom_group(geometry_restraints_manager,
                                            rg,
                                            ):
  sgs = None
  for atom in rg.atoms():
    if atom.name.strip()=='SG' and atom.parent().resname=='CYS':
      sgs = atom.i_seq
      break
  if sgs:
    sg_bonds = generate_bonded_i_seqs(geometry_restraints_manager, rg, sgs)
    if len(sg_bonds)>2:
      remove_cys_hg_from_residue_group(rg)

def process_disulphide_hydrogen_atoms(geometry_restraints_manager,
                                      residue_group,
                                      element='H',
                                      append_to_end_of_model=False,
                                      ):
  rc = conditional_add_cys_hg_to_atom_group(geometry_restraints_manager,
                                            residue_group,
                                            element=element,
                                            append_to_end_of_model=append_to_end_of_model)
  assert not rc, '%s' % rc
  conditional_remove_cys_hg_to_atom_group(geometry_restraints_manager,
                                          residue_group)
  return rc
