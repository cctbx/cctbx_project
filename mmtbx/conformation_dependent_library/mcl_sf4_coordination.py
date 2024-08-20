from __future__ import absolute_import, division, print_function
import sys

from cctbx import geometry_restraints
from six.moves import range

origin_ids = geometry_restraints.linking_class.linking_class()

# defaults to CYS
sf4_coordination = {
  'CYS' : {
    ("FE", "S")      : [  2.268, 0.017*2],
    ("S", "FE", "S") : [114.24,  5.75*2],
  },
  'MET' : {
    ("FE", "S")      : [  2.311, 0.006*2],
    ("S", "FE", "S") : [113.97,  8.764*2],
  },
  'HIS' : {
    ('FE', 'N')        : [  2.04,  0.05],
  }
}
fes_coordination = {
  'CYS' : {
    ("FE", "S")        : [  2.305, 0.022*2],
    ("S", "FE", "S")   : [111.20,  4.05*2],
    ('SG', 'FE', 'SG') : [107.77,  4.08*2],
  },
  'HIS' : {
    ('FE', 'N')        : [  2.14,  0.05],
  },
}
# parent pair
fes_coordination['CYS'][('CYS', 'CYS')] = ('SG', 'FE', 'SG')
#
f3s_coordination = {
  'CYS' : {
    ('FE', 'S')      : [  2.318, 0.008*2],
    ('S', 'FE', 'S') : [112.23,  6.03*2],
  },
}
# not coodinated number FE !- S
f3s_naming = {
  1 : 4,
  3 : 2,
  4 : 1,
}
coordination_defaults = {
  'SF4' : sf4_coordination,
  'F3S' : f3s_coordination,
  'FES' : fes_coordination,
}

phil_str = '''
'''

sf_clusters = set(['SF4', 'F3S', 'FES'])

def get_cluster_name(a1, a2, a3=None, other=False):
  resname = [a1.parent().resname,a2.parent().resname]
  if a3: resname.append(a3.parent().resname)
  resname = set(resname)
  if other:
    resname = resname.difference(set(sf_clusters))
  else:
    resname = sf_clusters.intersection(set(resname))
  if len(resname)==1: resname = resname.pop()
  elif len(resname)>1: resname = None
  else: resname = None
  return resname

def get_lookup(a1, a2, a3=None):
  from mmtbx.monomer_library import bondlength_defaults
  resname = get_cluster_name(a1, a2, a3)
  ligand = get_cluster_name(a1, a2, a3, other=True)
  cluster_lookup = coordination_defaults.get(resname, None)
  assert cluster_lookup, 'library for %s not found' % (resname)
  ligand_lookup = cluster_lookup.get(ligand, None)
  if ligand_lookup is None:
    ligand_lookup = cluster_lookup.get('CYS', None)
  if ligand_lookup is None:
    ans = bondlength_defaults.run(a1, a2)
    ligand_lookup = {}
    ligand_lookup[(a1.element.strip().upper(), a2.element.strip().upper())]=[ans, 0.1]
  return ligand_lookup

def get_distance_ideal_and_weight(a1, a2):
  ligand_lookup = get_lookup(a1, a2)
  key = (a1.element.strip().upper(), a2.element.strip().upper())
  if key not in ligand_lookup:
    return None, ' Atom pair %s %s not found in MCL' % (a1.quote(), a2.quote())
  distance_ideal=ligand_lookup[key][0]
  weight=1.0/ligand_lookup[key][1]**2
  return distance_ideal, weight

def get_angle_ideal_and_weight(a1,a2,a3):
  ligand_lookup = get_lookup(a1, a2, a3)
  key = (a1.element.strip().upper(),
         a2.element.strip().upper(),
         a3.element.strip().upper(),
         )
  if key not in ligand_lookup: return None, None
  parent_pair = (a1.parent().resname, a3.parent().resname)
  if parent_pair in ligand_lookup:
    key = ligand_lookup[parent_pair]
    angle_ideal=ligand_lookup[key][0]
    weight=1.0/ligand_lookup[key][1]**2
  else:
    angle_ideal = ligand_lookup[key][0]
    weight = 1.0/ligand_lookup[key][1]**2
  return angle_ideal, weight

def get_sulfur_iron_cluster_coordination(pdb_hierarchy,
                                         nonbonded_proxies,
                                         sorted_nb_proxies_res=None,
                                         coordination_distance_cutoff=3.5,
                                         #params=None,
                                         log=sys.stdout,
                                         verbose=False,
                                       ):
  coordination = []
  done_aa = []
  atoms = pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  selection_string = " or ".join(['resname %s' %x for x in sf_clusters])
  s = pdb_hierarchy.atom_selection_cache().selection(selection_string)
  if s.all_eq(False):
    return coordination
  for item in nonbonded_proxies.sorted_value_proxies_generator(
      by_value="delta",
      sites_cart=sites_cart,
      cutoff=coordination_distance_cutoff):
    i_seq, j_seq, dist, sym_op_j, rt_mx, proxy = item
    a1 = atoms[i_seq]
    ag1 = a1.parent()
    a2 = atoms[j_seq]
    ag2 = a2.parent()
    current = set([ag1.resname, ag2.resname])
    intersection = sf_clusters.intersection(current)
    if len(intersection)==2:
      if ag1.id_str()!=ag2.id_str():
        print('Two residues (%s, %s) are close enough to coordinate! ODD!' % (
          ag1.id_str(),
          ag2.id_str()), file=log)
    elif len(intersection)==1:
      if rt_mx:
        coordination = []
        break
      resname = intersection.pop()
      sf4=a2
      sf4g=ag2
      aa=a1
      aag=ag1
      if ag1.resname==resname:
        sf4=a1
        sf4g=ag2
        aa=a2
        aag=ag2
      if aa.element.strip() not in ['S', 'N']: continue
      # if aa.element.strip() in ['H', 'D']: continue
      if verbose: print('%s-aa' % resname,sf4.quote(),aa.quote(),dist)
      if sf4.element.lower()=="fe":
        if aag.id_str() not in done_aa:
          #coordination.append((i_seq, j_seq))
          coordination.append((sf4, aa))
          done_aa.append(aag.id_str())
  return coordination

def get_bond_proxies(coordination, verbose=False):
  #
  bonds = []
  if coordination is None: return bonds
  for a1, a2 in coordination:
    distance_ideal, weight = get_distance_ideal_and_weight(a1, a2)
    if distance_ideal is None:
      if verbose: print('no distance_ideal %s %s' % (distance_ideal, weight))
      continue
    p = geometry_restraints.bond_simple_proxy(
      i_seqs=[a1.i_seq, a2.i_seq],
      distance_ideal=distance_ideal,
      weight=weight,
      slack=0,
      top_out=False,
      limit=1,
      origin_id=origin_ids.get_origin_id('metal coordination'))
    bonds.append(p)
  return bonds

def get_angle_proxies_for_bond(coordination):
  #
  def _get_angle_atoms(a1, a2, resname, second_residues):
    atoms = []
    ii=int(a1.name.strip()[-1])
    if resname=='F3S':
      for i in range(1,5):
        if i == f3s_naming.get(ii, -1): continue
        name = 'S%d' % i
        a3 = a1.parent().get_atom(name)
        if a3: atoms.append(a3)
    else:
      # SF4 has a special naming scheme
      for i in range(1,5):
        if i==ii: continue
        name = 'S%d' % i
        a3 = a1.parent().get_atom(name)
        if a3: atoms.append(a3)
    if resname in ['FES']:
      for ag in second_residues:
        if ag.id_str()==a2.parent().id_str(): continue
        for name in ['SG']:
          sg = ag.get_atom(name)
          if sg and sg.distance(a1)<3.5:
            atoms.append(sg)
    return atoms
  #
  angles = []
  if coordination is None: return angles
  second_residues = []
  for a1, a2 in coordination:
    second_residues.append(a2.parent())
  for a1, a2 in coordination:
    assert a1.name.find("FE")>-1
    resname = get_cluster_name(a1, a2)
    if resname in sf_clusters:
      atoms = _get_angle_atoms(a1, a2, resname, second_residues)
      for a3 in atoms:
        angle_ideal, weight = get_angle_ideal_and_weight(a3, a1, a2)
        if angle_ideal is None: continue
        p = geometry_restraints.angle_proxy(
          i_seqs=[a3.i_seq, a1.i_seq, a2.i_seq],
          angle_ideal=angle_ideal,
          weight=weight,
          origin_id=origin_ids.get_origin_id('metal coordination'))
        angles.append(p)
  return angles

def get_all_proxies(coordination, resname=None):
  return get_bond_proxies(coordination), \
      get_angle_proxies_for_bond(coordination)

def run(pdb_filename):
  print("run",pdb_filename)
  from mmtbx.command_line.geometry_minimization import \
    get_geometry_restraints_manager, master_params
  import mmtbx.monomer_library.pdb_interpretation
  from mmtbx import monomer_library

  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = pdb_filename,
    #force_symmetry = True,
  )
  xrs = processed_pdb_file.xray_structure()
  #work_params = master_params().extract()
  #work_params.reference_model.enabled=True
  #work_params.reference_model.use_starting_model_as_reference=True
  grm = get_geometry_restraints_manager(
    processed_pdb_file,
    xrs,
    #params=work_params,
    #log=null_out(),
  )
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  rc = get_sf4_coordination(
    pdb_hierarchy=pdb_hierarchy,
    nonbonded_proxies=grm.geometry.pair_proxies(
      sites_cart=pdb_hierarchy.atoms().extract_xyz()).nonbonded_proxies,
    #verbose=True,
  )
  bproxies, aproxies = get_all_proxies(rc)
  print(len(bproxies),len(aproxies))
  grm.geometry.add_new_bond_restraints_in_place(
    proxies=bproxies,
    sites_cart=pdb_hierarchy.atoms().extract_xyz(),
  )
  grm.geometry.add_angles_in_place(aproxies)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
