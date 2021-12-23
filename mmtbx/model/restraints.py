from __future__ import absolute_import,division, print_function

from math import sqrt

from libtbx.utils import Sorry

from mmtbx import monomer_library
import mmtbx.monomer_library.server

def atom_as_restraint(atom):
  resname = atom.parent().resname
  name = atom.name
  element = atom.element
  atom_type = element
  charge = atom.charge
  # outl = ' %(resname)s %(name)s %(element)s %(atom_type)s %(charge)s ' % locals()
  outl = ' %(resname)s %(name)s %(element)s' % locals()
  return outl

def bond_as_restraint(bond, atom1, atom2):
  resname = atom1.parent().resname
  name1 = atom1.name
  name2 = atom2.name
  ideal = bond.distance_ideal
  esd = 1/sqrt(bond.weight)
  outl = ' %(resname)s %(name1)s %(name2)s coval %(ideal).3f %(esd)s ' % locals()
  return outl

def get_bond_dictionary(ligand_model, ideal=True):
  ligand_grm = ligand_model.get_restraints_manager()
  bond_proxies_simple, asu = ligand_grm.geometry.get_all_bond_proxies(
    sites_cart=ligand_model.get_sites_cart())
  sorted_table, n_not_shown = bond_proxies_simple.get_sorted(
    'delta',
    ligand_model.get_sites_cart())
  atoms = ligand_model.get_atoms()
  rc = {}
  for info in sorted_table:
    (i_seq, j_seq, i_seqs, distance_ideal, distance_model, slack, delta, sigma, weight, residual, sym_op_j, rt_mx) = info
    i_atom=atoms[i_seq]
    j_atom=atoms[j_seq]
    key = [i_atom.name.strip(), j_atom.name.strip()]
    key.sort()
    if ideal:
      rc[tuple(key)]=distance_ideal
    else:
      rc[tuple(key)]=distance_model
  return rc

def get_angle_dictionary(ligand_model, ideal=True):
  ligand_grm = ligand_model.get_restraints_manager()
  sorted_table, n_not_shown = ligand_grm.geometry.angle_proxies.get_sorted(
      'delta',
      ligand_model.get_sites_cart())
  atoms = ligand_model.get_atoms()
  rc={}
  for info in sorted_table:
    (i_seqs, angle_ideal, angle_model, delta, sigma, weight, residual) = info
    i_atom=atoms[int(i_seqs[0])]
    j_atom=atoms[int(i_seqs[1])]
    k_atom=atoms[int(i_seqs[2])]
    key = [i_atom.name.strip(), j_atom.name.strip(), k_atom.name.strip()]
    key.sort()
    if ideal:
      rc[tuple(key)]=angle_ideal
    else:
      rc[tuple(key)]=angle_model
  return rc

def get_torsion_dictionary(ligand_model, ideal=True):
  ligand_grm = ligand_model.get_restraints_manager()
  sorted_table, n_not_shown = ligand_grm.geometry.dihedral_proxies.get_sorted(
      'delta',
      ligand_model.get_sites_cart())
  atoms = ligand_model.get_atoms()
  rc = {}
  if sorted_table is None: sorted_table=[]
  i=0
  for info in sorted_table:
    (i_seqs, angle_ideal, angle_model, delta, period, sigma, weight, residual) = info
    i_atom=atoms[int(i_seqs[0])]
    j_atom=atoms[int(i_seqs[1])]
    k_atom=atoms[int(i_seqs[2])]
    l_atom=atoms[int(i_seqs[3])]
    key = [i_atom.name.strip(), j_atom.name.strip(), k_atom.name.strip(), l_atom.name.strip()]
    key.sort()
    if ideal:
      rc[tuple(key)]=angle_ideal
    else:
      rc[tuple(key)]=angle_model
  return rc

def get_restraints_from_model_via_grm(ligand_model, ideal=True, verbose=True):
  bond_lookup = get_bond_dictionary(ligand_model, ideal=ideal)
  angle_lookup = get_angle_dictionary(ligand_model, ideal=ideal)
  codes = []
  atoms = ligand_model.get_atoms()
  for atom in atoms:
    resname = atom.parent().resname
    if resname not in codes:
      codes.append(resname)
  if len(codes)>1:
    raise Sorry('more than one entity sent to restraints writer %s' % codes)
  resname = codes[0]
  mon_lib_srv = monomer_library.server.server()
  rc = mon_lib_srv.get_comp_comp_id_direct(resname)
  co = rc.cif_object
  for key, loop in co.loops.items():
    if key=='_chem_comp_bond':
      values = []
      for row in loop.iterrows():
        key = [row['_chem_comp_bond.atom_id_1'], row['_chem_comp_bond.atom_id_2']]
        key.sort()
        distance_ideal = bond_lookup.get(tuple(key), None)
        assert distance_ideal
        values.append(distance_ideal)
      loop.update_column('_chem_comp_bond.value_dist', values)
    elif key=='_chem_comp_angle':
      values = []
      for row in loop.iterrows():
        key = [row['_chem_comp_angle.atom_id_1'],
               row['_chem_comp_angle.atom_id_2'],
               row['_chem_comp_angle.atom_id_3'],
               ]
        key.sort()
        angle_ideal = angle_lookup.get(tuple(key), None)
        assert angle_ideal
        values.append(angle_ideal)
      loop.update_column('_chem_comp_angle.value_angle', values)
  assert 0
  return co
