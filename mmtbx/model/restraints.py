from __future__ import absolute_import,division, print_function

import copy

from math import sqrt

from libtbx.utils import Sorry

import iotbx
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

def get_bond_dictionary(ligand_model, ligand_grm=None, ideal=True):
  if ligand_grm is None:
    ligand_grm = ligand_model.get_restraints_manager()
  bond_proxies_simple, asu = ligand_grm.geometry.get_all_bond_proxies(
    sites_cart=ligand_model.get_sites_cart())
  sorted_table, n_not_shown = bond_proxies_simple.get_sorted(
    'delta',
    ligand_model.get_sites_cart())
  if sorted_table is None: sorted_table=[]
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

def get_angle_dictionary(ligand_model, ligand_grm=None, ideal=True):
  if ligand_grm is None:
    ligand_grm = ligand_model.get_restraints_manager()
  sorted_table, n_not_shown = ligand_grm.geometry.angle_proxies.get_sorted(
      'delta',
      ligand_model.get_sites_cart())
  if sorted_table is None: sorted_table=[]
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

def get_torsion_dictionary(ligand_model, ligand_grm=None, ideal=True):
  if ligand_grm is None:
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
    if ideal:
      rc[tuple(key)]=angle_ideal
    else:
      rc[tuple(key)]=angle_model
  return rc

def get_restraints_from_model_via_grm(ligand_model,
                                      ligand_grm=None,
                                      ideal=True,
                                      cartesian_coordinates=True,
                                      ):
  """Write the restraints from the geometry using the CIF object from the beginning

  Args:
      ligand_model (TYPE): Model object of one entity with no alt. loc.
      ideal (bool, optional): Use the ideal distance from the proxy rather than the acutal
      cartesian_coordinates (bool, optional): Update the atom loop with the hierachy coordinates

  Returns:
      TYPE: CIF object

  Raises:
      Sorry: Description
  """
  from libtbx.utils import null_out
  ligand_model = copy.deepcopy(ligand_model)
  ligand_model.unset_restraints_manager()
  ligand_model.log=null_out()
  ligand_model.process(make_restraints=True)

  bond_lookup = get_bond_dictionary(ligand_model, ligand_grm=ligand_grm, ideal=ideal)
  angle_lookup = get_angle_dictionary(ligand_model, ligand_grm=ligand_grm, ideal=ideal)
  torsion_lookup = get_torsion_dictionary(ligand_model, ligand_grm=ligand_grm, ideal=ideal)
  codes = []
  atoms = ligand_model.get_atoms()
  for atom in atoms:
    resname = atom.parent().resname
    if resname not in codes:
      codes.append(resname)
  if len(codes)>1:
    raise Sorry('more than one entity sent to restraints writer %s' % codes)
  #
  def get_ligand_code_from_cif_object(cif_object):
    code=None
    if 'comp_list' in cif_object.keys():
      for key, co in cif_object.items():
        if key=='comp_list':
          for field, loop in co.loops.items():
            for row in loop.iterrows():
              code = row['_chem_comp.id']
              break
    return code
  #
  # needs to be smarter
  #
  resname = codes[0]
  restraints = ligand_model.get_restraint_objects()
  cif_object=None
  for filename, cif_object in restraints:
    code = get_ligand_code_from_cif_object(cif_object)
    if code==resname:
      break
  else:
    mon_lib_srv = monomer_library.server.server()
    filename = mon_lib_srv.get_comp_comp_id_direct(resname, return_filename=True)
    cif_object = iotbx.cif.reader(filename).model()
  #
  # find the columns to replace
  #
  if 'comp_list' in cif_object.keys():
    for key, co in cif_object.items():
      if key=='comp_%s' % resname.strip():
        break
    else:
      outl = 'To write restraints, a restraints file needs to be provided to refinement.'
      raise Sorry(outl)
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
    elif key=='_chem_comp_tor':
      def _find_torsion_match(key):
        """Matches the torsion keys and flips the sign of the value if needed

        Args:
            key (List): List of atom names in torsion

        Returns:
            float: Value of torsion sign corrected
        """
        k1 = copy.copy(key)
        k1.sort()
        for lookup_key, ai in torsion_lookup.items():
          k2 = list(copy.copy(lookup_key))
          k2.sort()
          if k1==k2:
            if key[0]!=lookup_key[0]:
              ai=ai*-1
            if key[1]!=lookup_key[1]:
              ai=ai*-1
            return ai
        return None
      values = []
      for row in loop.iterrows():
        key = [row['_chem_comp_tor.atom_id_1'],
               row['_chem_comp_tor.atom_id_2'],
               row['_chem_comp_tor.atom_id_3'],
               row['_chem_comp_tor.atom_id_4'],
               ]
        angle_ideal = torsion_lookup.get(tuple(key), None)
        if angle_ideal is None: angle_ideal = _find_torsion_match(key)
        if angle_ideal:
          values.append(angle_ideal)
        else:
          tor_id = row.get('_chem_comp_tor.id')
          if tor_id.find('CONST')==-1:
            print('WARNING: torsion angle %s not transfered' % (row.get('_chem_comp_tor.id', None)))
          values.append(float(row.get('_chem_comp_tor.value_angle')))
      loop.update_column('_chem_comp_tor.value_angle', values)
    elif key=='_chem_comp_atom' and cartesian_coordinates:
      ags=[]
      for i, ag in enumerate(ligand_model.get_hierarchy().atom_groups()):
        ags.append(ag)
      assert len(ags)==1
      ag=ags[0]
      values=[[],[],[]]
      for row in loop.iterrows():
        name = row['_chem_comp_atom.atom_id']
        atom = ag.get_atom(name)
        if atom is None:
          print('  Failed to find atom %s' % name)
        else:
          for i in range(3):
            values[i].append('%.5f' % atom.xyz[i])
      loop.update_column('_chem_comp_atom.x', values[0])
      loop.update_column('_chem_comp_atom.y', values[1])
      loop.update_column('_chem_comp_atom.z', values[2])
  return cif_object
