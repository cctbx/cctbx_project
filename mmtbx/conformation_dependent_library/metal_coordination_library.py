from __future__ import division
import sys

from cctbx import geometry_restraints
from cctbx.geometry_restraints import linking_class
origin_ids = linking_class.linking_class()

headers = ['Zn-SG (A)',
           'SG-Zn-SG (degrees)',
           'Zn-N (A)',
           'N-Zn-N (degrees)',
           'ZnCysxHisy',
           ]
database = {'Zn2+ tetrahedral': {
    (4,0) : {
        ('ZN', 'SG') : (2.330, 0.029, 1033),
        ('SG', 'ZN', 'SG') : (109.45, 5.46, 1553),
        }, # n/a        n/a     Cys4
    (3,1) : {
        ('ZN', 'SG') : (2.318, 0.027, 912),
        ('SG', 'ZN', 'SG') : (112.15, 3.96, 912),
        ('ZN', 'ND1') : (2.074, 0.056, 303),
        }, #n/a Cys3His1
    (2,2) : {
        ('ZN', 'SG') : (2.306, 0.029, 76),
        ('SG', 'ZN', 'SG') : (116.23, 4.58, 38),
        ('ZN', 'ND1') : (2.040, 0.050, 65),
        ('ND1', 'ZN', 'ND1') : (102.38, 5.44, 38),
        }, #    Cys2His2
    (1,3) : {
        ('ZN', 'SG') : (2.298, 0.017, 12),
        ('ZN', 'ND1') : (2.002, 0.045, 36),
        ('ND1', 'ZN', 'ND1') : (107.23, 4.78, 36),
        }, #    Cys1His3
    (0,4) : {
        }, #n/a n/a     Insufficient data       Insufficient data       His4
    }
  }
database['ZN'] = database['Zn2+ tetrahedral']

for nums, restraints in database['Zn2+ tetrahedral'].items():
  for atoms, values in list(restraints.items()):
    if 'ND1' in atoms:
      key = list(atoms)
      for i, atom in enumerate(key):
        if atom=='ND1':
          key[i]='NE2'
          restraints[tuple(key)]=values

def print_restraints(db):
  print('-'*80)
  for coordination, library in sorted(db.items()):
    print(coordination)
    for nums, restraints in sorted(library.items()):
      print('  %s' % str(nums))
      for atoms, values in sorted(restraints.items()):
        print('    ',atoms, values)

def check_other_in_database(metal, other):
  sub = database.get(metal.name.strip(), None)
  if sub is None: return False
  for key, item in sub.items():
    for atoms in item:
      if other.name.strip() in atoms: return True
  return False

def get_metal_coordination_proxies(pdb_hierarchy,
                                   nonbonded_proxies,
                                   sorted_nb_proxies_res=None,
                                   prefix=None,
                                   params=None,
                                   log=sys.stdout,
                                   add_segid=None,
                                   verbose=False):
  def is_residue_already_linked_to_metal(linked, atom):
    for link in linked:
      if link.parent().id_str()==atom.parent().id_str():
        break
    else:
      return False
    return True
  hbond_distance_cutoff = 3
  if params is not None:
    hbond_distance_cutoff = params.hbond_distance_cutoff
  mbonds = {}
  atoms = pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  get_sorted_result = sorted_nb_proxies_res
  if get_sorted_result is None:
    get_sorted_result = nonbonded_proxies.get_sorted(
        by_value="delta",
        sites_cart=sites_cart)
    # its never None, it is always a tuple ([], int)
    if get_sorted_result is None: return mbonds
  sorted_nonb, n_not_shown = get_sorted_result
  n_nonb = len(sorted_nonb)
  i = 0
  while i < n_nonb and sorted_nonb[i][3] < hbond_distance_cutoff:
    (labels, i_seq, j_seq, dist, vdw_distance, sym_op_j, rt_mx) = sorted_nonb[i]
    a1 = atoms[i_seq]
    ag1 = a1.parent()
    a2 = atoms[j_seq]
    ag2 = a2.parent()
    metal = None
    for metal_element in database:
      if a1.name.strip()==metal_element:
        metal = a1
        other = a2
      elif a2.name.strip()==metal_element:
        metal = a2
        other = a1
      if metal:
        if not check_other_in_database(metal, other): continue
        mbonds.setdefault(metal.i_seq, {})
        mbonds[metal.i_seq]['metal'] = metal
        mbonds[metal.i_seq].setdefault('others', [])
        if not is_residue_already_linked_to_metal(mbonds[metal.i_seq]['others'],
                                                  other,
                                                  ):
          mbonds[metal.i_seq]['others'].append(other)
    i += 1
  pairs = []
  if verbose:
    for key, item in mbonds.items():
      for label, l in item.items():
        if type(l)==type([]):
          for atom in l:
            print('  ',atom.quote())
        else:
          print(l.quote())
  return mbonds

def get_proxies(coordination, verbose=False):
  #
  # TODO
  #   - check that only one link is made to each resiude
  #     e.g. 1a6y "2080 ZN    ZN B 451 .*." "1874  CB  CYS B 153 .*."
  #
  def _bond_generator(atoms):
    for atom in atoms['others']:
      yield atoms['metal'], atom
  def _angle_generator(atoms):
    for i, a1 in enumerate(atoms['others']):
     for j, a2 in enumerate(atoms['others']):
       if i==j: break
       yield a1, atoms['metal'], a2
  def _default_bonds(atoms):
    tmp = []
    for a1, a2 in _bond_generator(atoms):
      p = geometry_restraints.bond_simple_proxy(
        i_seqs=[a1.i_seq, a2.i_seq],
        distance_ideal=2.3,
        weight=1.0/0.03**2,
        slack=0,
        top_out=False,
        limit=1,
        origin_id=origin_ids.get_origin_id('metal coordination'))
      tmp.append(p)
    return tmp
  bonds = []
  angles = []
  if coordination is None: return bonds, angles
  atoms = None
  for metal_i_seq, atoms in coordination.items():
    if len(atoms['others'])<4:
      bonds += _default_bonds(atoms)
      continue
    cyss = []
    hiss = []
    for atom in atoms['others']:
      if atom.parent().resname=='CYS':
        cyss.append(atom)
      elif atom.parent().resname=='HIS':
        hiss.append(atom)
    key = (len(cyss), len(hiss))
    metal_name = atoms['metal'].name.strip()
    if key not in database[metal_name]:
      if verbose:
        outl = '\n'
        for atom in atoms['others']:
          outl += '    %s\n' % atom.quote()
        print('''  Metal %s has coordination not in MCL%s
              ''' % (atoms['metal'].quote(), outl)
              )
      continue
    ideals = database[metal_name][key]
    if not atoms: continue #return None, None
    for a1, a2 in _bond_generator(atoms):
      key = (a1.name.strip(), a2.name.strip())
      if key not in ideals: continue
      t = ideals[key]
      p = geometry_restraints.bond_simple_proxy(
        i_seqs=[a1.i_seq, a2.i_seq],
        distance_ideal=t[0],
        weight=1.0/t[1]**2,
        slack=0,
        top_out=False,
        limit=1,
        origin_id=origin_ids.get_origin_id('metal coordination'))
      bonds.append(p)
    for a1, a2, a3 in _angle_generator(atoms):
      key = (a1.name.strip(), a2.name.strip(), a3.name.strip())
      if key not in ideals: continue
      t = ideals[key]
      p = geometry_restraints.angle_proxy(
        i_seqs=[a1.i_seq, a2.i_seq, a3.i_seq],
        angle_ideal=t[0],
        weight=1.0/t[1]**2,
        origin_id=origin_ids.get_origin_id('metal coordination'))
      angles.append(p)
  return bonds, angles

def run(model_filename=None):
  import mmtbx.monomer_library.pdb_interpretation as pdb_inter
  #print_restraints(database)
  if model_filename is not None:
    from iotbx import pdb
    pdb_inp = pdb.input(model_filename)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    #pdb_hierarchy.show()
    pdb_processed_file = pdb_inter.run(
      args=[model_filename],
      assume_hydrogens_all_missing=False,
      hard_minimum_nonbonded_distance=0.0,
      nonbonded_distance_threshold=None,
      substitute_non_crystallographic_unit_cell_if_necessary=True,
      )
    grm = pdb_processed_file.geometry_restraints_manager()
    xrs = pdb_processed_file.xray_structure()
    sites_cart = xrs.sites_cart()
    site_labels = xrs.scatterers().extract_labels()
    pair_proxies=grm.pair_proxies(sites_cart=sites_cart,site_labels=site_labels)
    proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted(
        by_value="delta",
        sites_cart=sites_cart,
        site_labels=site_labels)
    rc = get_metal_coordination_proxies(pdb_hierarchy,
                                        pair_proxies.nonbonded_proxies,
    )
    bonds, angles = get_proxies(rc)
    print('\n\tbonds, angles : %d, %d\n\n' % (len(bonds), len(angles)))


if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
