import os, sys

from cctbx import geometry_restraints
origin_ids = geometry_restraints.linking_class.linking_class()

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
	}, # n/a	n/a	Cys4
    (3,1) : {
        ('ZN', 'SG') : (2.318, 0.027, 912),
	('SG', 'ZN', 'SG') : (112.15, 3.96, 912),
        ('ZN', 'ND1') : (2.074, 0.056, 303),
	}, #n/a	Cys3His1
    (2,2) : {
        ('ZN', 'SG') : (2.306, 0.029, 76),
	('SG', 'ZN', 'SG') : (116.23, 4.58, 38),
        ('ZN', 'ND1') : (2.040, 0.050, 65),
        ('ND1', 'ZN', 'ND1') : (102.38, 5.44, 38),
        }, #	Cys2His2
    (1,3) : {
        ('ZN', 'SG') : (2.298, 0.017, 12),
        ('ZN', 'ND1') : (2.002, 0.045, 36),
        ('ND1', 'ZN', 'ND1') : (107.23, 4.78, 36),
        }, #	Cys1His3
    (0,4) : {
        }, #n/a	n/a	Insufficient data	Insufficient data	His4
    }
  }
database['ZN'] = database['Zn2+ tetrahedral']

# in code already
#database["SF4"] = { 
#    (('FEn', 'SF4'), ('SG', 'CYS')) : (2.3, 0.03, -1),
#  }

for nums, restraints in database['Zn2+ tetrahedral'].items():
  for atoms, values in restraints.items():
    if 'ND1' in atoms:
      key = list(atoms)
      for i, atom in enumerate(key):
        if atom=='ND1':
          key[i]='NE2'
          restraints[tuple(key)]=values

def print_restraints(db):
  print '-'*80
  for coordination, library in sorted(db.items()):
    print coordination
    for nums, restraints in sorted(library.items()):
      print '  %s' % str(nums)
      for atoms, values in sorted(restraints.items()):
        print '    ',atoms, values

def get_metal_coordination_proxies(pdb_hierarchy,
                                   nonbonded_proxies,
                                   prefix=None,
                                   params=None,
                                   log=sys.stdout,
                                   add_segid=None,
                                   verbose=False):
  hbond_distance_cutoff = 3
  if params is not None:
    hbond_distance_cutoff = params.hbond_distance_cutoff
  mbonds = {}
  atoms = pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  get_sorted_result = nonbonded_proxies.get_sorted(
      by_value="delta",
      sites_cart=sites_cart)
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
        mbonds.setdefault(metal.i_seq, {})
        mbonds[metal.i_seq]['metal'] = metal
        mbonds[metal.i_seq].setdefault('others', [])
        mbonds[metal.i_seq]['others'].append(other)
    i += 1
    continue
    assert 0
    if (common_residue_names_get_class(ag1.resname, consider_ccp4_mon_lib_rna_dna=True) in \
          ["common_rna_dna", "ccp4_mon_lib_rna_dna"] and
        common_residue_names_get_class(ag2.resname, consider_ccp4_mon_lib_rna_dna=True) in \
          ["common_rna_dna", "ccp4_mon_lib_rna_dna"] and
        (a1.element in ["N", "O"] and a2.element in ["N", "O"]) and
        a1.name.find("P") < 0 and a2.name.find("P") < 0 and
        a1.name.find("'") < 0 and a2.name.find("'") < 0 and
        not consecutive_residues(a1, a2) and
        (ag1.altloc.strip() == ag2.altloc.strip()) and
         final_link_direction_check(a1, a2)):
      hbonds.append((i_seq, j_seq))
  # check and define basepairs
  pairs = []
  if 0: #for hb in hbonds:
    if verbose > 1:
      print >> log, "Making pair with", atoms[hb[0]].id_str(), atoms[hb[1]].id_str()
    new_hbonds, class_number = get_h_bonds_for_basepair(
        atoms[hb[0]],
        atoms[hb[1]],
        distance_cutoff=hbond_distance_cutoff,
        log=log,
        verbose=verbose)
    if verbose > 1:
      print >> log, "  Picked class: %d, number of h-bonds under cutoff:%d" % (class_number, len(new_hbonds)),
    if len(new_hbonds) > 1:
      p = make_phil_base_pair_record(atoms[hb[0]].parent(), atoms[hb[1]].parent(),
          params, saenger_class=class_number, add_segid=add_segid)
      if verbose > 1:
        print >> log, "  OK"
      pairs.append(p)
    else:
      if verbose > 0:
        s = " ".join(["Basepairing for residues '%s' and '%s'" % (
            atoms[hb[0]].id_str()[10:-1], atoms[hb[1]].id_str()[10:-1]),
          "was rejected because only 1 h-bond was found"])
        if verbose > 1:
          print >> log, "Rejected"
  return mbonds

def get_proxies(coordination):
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
  bonds = []
  angles = []
  if coordination is None: return bonds, angles
  atoms = None
  for metal_i_seq, atoms in coordination.items():
    if len(atoms['others'])!=4: continue
    cyss = []
    hiss = []
    for atom in atoms['others']:
      if atom.parent().resname=='CYS':
        cyss.append(atom)
      elif atom.parent().resname=='HIS':
        hiss.append(atom)
    key = (len(cyss), len(hiss))
    metal_name = atoms['metal'].name.strip()
    ideals = database[metal_name][key]
  if not atoms: return None, None
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
    rc = get_metal_coordination_bond_proxies(pdb_hierarchy,
                                             pair_proxies.nonbonded_proxies,
                                           )
    bonds, angles = get_proxies(rc)
    print bonds
    print angles
    assert 0

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
