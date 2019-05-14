from __future__ import absolute_import, division, print_function
import sys

from cctbx import geometry_restraints

origin_ids = geometry_restraints.linking_class.linking_class()

sf4_coordination = {
  ("FE", "S")      : [  2.268, 0.017*2],
  ("S", "FE", "S") : [114.24,  5.75*2],
                    }

phil_str = '''
'''

def get_sulfur_iron_cluster_coordination(pdb_hierarchy,
                                         nonbonded_proxies,
                                         coordination_distance_cutoff=3.5,
                                         #params=None,
                                         log=sys.stdout,
                                         verbose=False,
                                       ):
  coordination = []
  done_aa = []
  atoms = pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  get_sorted_result = nonbonded_proxies.get_sorted(
      by_value="delta",
      sites_cart=sites_cart)
  if get_sorted_result is None:
    return None
  sorted_nonb, n_not_shown = get_sorted_result

  # Get potential hbonds
  n_nonb = len(sorted_nonb)
  i = 0
  while i < n_nonb and sorted_nonb[i][3] < coordination_distance_cutoff:
    (labels, i_seq, j_seq, dist, vdw_distance, sym_op_j, rt_mx) = sorted_nonb[i]
    a1 = atoms[i_seq]
    ag1 = a1.parent()
    a2 = atoms[j_seq]
    ag2 = a2.parent()
    for resname in ['SF4', 'F3S']:
      if (ag1.resname==resname and ag2.resname==resname):
        if ag1.id_str()!=ag2.id_str():
          print('Two %(resname)s residues are close enough to coordinate! ODD!' % locals(), file=log)
      elif (ag1.resname==resname or ag2.resname==resname):
        sf4=a2
        sf4g=ag2
        aa=a1
        aag=ag1
        if ag1.resname==resname:
          sf4=a1
          sf4g=ag2
          aa=a2
          aag=ag2
        if aa.element.strip() in ['H', 'D']: continue
        if aa.element.strip() in ['CA', 'C', 'O', 'N']: continue
        if verbose: print('%s-aa' % resname,sf4.quote(),aa.quote(),dist)
        if sf4.element.lower()=="fe":
          if aag.id_str() not in done_aa:
            #coordination.append((i_seq, j_seq))
            coordination.append((sf4, aa))
            done_aa.append(aag.id_str())
    i += 1
  return coordination

def get_bond_proxies(coordination):
  #
  # works for SF4, F3S
  #
  bonds = []
  if coordination is None: return bonds
  for a1, a2 in coordination:
    p = geometry_restraints.bond_simple_proxy(
      i_seqs=[a1.i_seq, a2.i_seq],
      distance_ideal=sf4_coordination[('FE', 'S')][0],
      weight=1.0/sf4_coordination[('FE', 'S')][1]**2,
      slack=0,
      top_out=False,
      limit=1,
      origin_id=origin_ids.get_origin_id('metal coordination'))
    bonds.append(p)
  return bonds

def get_angle_proxies_for_bond(coordination):
  #
  # works for SF4
  # does not work for F3S
  #
  angles = []
  if coordination is None: return angles
  for a1, a2 in coordination:
    assert a1.name.find("FE")>-1
    if a1.parent().resname in ['F3S']: break
    assert a1.parent().resname in ['SF4']
    ii=int(a1.name.strip()[-1])
    for i in range(1,5):
      if i==ii: continue
      name = 'S%d' % i
      a3 = a1.parent().get_atom(name)
      angle_ideal = sf4_coordination[('S', 'FE', 'S')][0]
      weight = sf4_coordination[('S', 'FE', 'S')][1]
      p = geometry_restraints.angle_proxy(
        i_seqs=[a3.i_seq, a1.i_seq, a2.i_seq],
        angle_ideal=angle_ideal,
        weight=1./weight**2,
        origin_id=origin_ids.get_origin_id('metal coordination'))
      angles.append(p)
  return angles

def get_all_proxies(coordination):
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
