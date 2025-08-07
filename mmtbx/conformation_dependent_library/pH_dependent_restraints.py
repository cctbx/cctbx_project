"""
Module to automatically attempt to load "low" pH range restraints from
restraints library - useful for various protonation states
"""
from __future__ import absolute_import, division, print_function
import time

from mmtbx.monomer_library import server
from cctbx import geometry_restraints

def process_bonds(gpr, bond, atom_dict, atom1, atom2, name1, name2, neutron1, neutron2, atoms=None, verbose=False):
  atoms_added = {}
  atoms_added_energy = {}
  bond_counters = [0,0]
  i_seqs = [atom1.i_seq, atom2.i_seq]
  i_seqs.sort()
  k=0
  l=1
  bond_table_entry = gpr.bond_simple.table[i_seqs[k]]
  if ( not bond_table_entry or
       i_seqs[l] not in gpr.bond_simple.table[i_seqs[k]]):
    k=1
    l=0
    bond_table_entry = gpr.bond_simple.table[i_seqs[k]]
  found_proxy=False
  for i, simple_proxy in enumerate(gpr.bond_simple.proxies):
    if simple_proxy.i_seqs==tuple(i_seqs):
      found_proxy=simple_proxy
      break
  if found_proxy:
    diff=abs(simple_proxy.distance_ideal-bond.value_dist)
    # if verbose: print('DIFF',simple_proxy.distance_ideal,bond.value_dist,diff)
    simple_proxy.distance_ideal = bond.value_dist
    simple_proxy.weight=1/bond.value_dist_esd**2
    if diff>bond.value_dist_esd: bond_counters[0]+=1 # changed
  else:
    if verbose: print('Adding bond for %s = %s' % (atom1.quote(), atom2.quote()))
    if neutron1 or neutron2:
      proxy = geometry_restraints.bond_simple_proxy(
        i_seqs=i_seqs,
        distance_ideal=bond.value_dist_neutron,
        weight=1/(bond.value_dist_esd**2),
      )
    else:
      proxy = geometry_restraints.bond_simple_proxy(
        i_seqs=i_seqs,
        distance_ideal=bond.value_dist,
        weight=1/(bond.value_dist_esd**2),
      )
    gpr.bond_simple.proxies.append(proxy)
    atoms_added_energy[atom1.i_seq] = atom_dict.get(name1.strip(), None)
    atoms_added_energy[atom2.i_seq] = atom_dict.get(name2.strip(), None)
    atoms_added[(atom1.i_seq, atom2.i_seq)] = (atom1, atom2)
    bond_counters[1]+=1
  return bond_counters, atoms_added_energy, atoms_added

def _get_atom_neutron(ag, name, bondlength=None):
  atom = ag.get_atom(name)
  if atom: return atom, atom.name, False
  if bondlength is not None and bondlength<1:
    atom = ag.get_atom(name.replace('H', 'D'))
    if atom: return atom, name, True
  return atom, None, False

def _generate_bond_atoms(rg, lookup_name1, lookup_name2, bondlength=None, verbose=False):
  for i, ag1 in enumerate(rg.atom_groups()):
    for j, ag2 in enumerate(rg.atom_groups()):
      if j>i: break
      atom1, name1, neutron1 = _get_atom_neutron(ag1,
                                                 lookup_name1,
                                                 bondlength=bondlength)
      if atom1 is None: continue
      if verbose: print('1',atom1.quote(), name1, neutron1)
      atom2, name2, neutron2 = _get_atom_neutron(ag2,
                                                 lookup_name2,
                                                 bondlength=bondlength)
      if atom2 is None: continue
      if verbose: print('2',atom2.quote(), name2, neutron2)
      yield atom1, atom2, name1, name2, neutron1, neutron2

def _generate_angle_atoms(rg, lookup_name1, lookup_name2, lookup_name3, verbose=False):
  for i, ag1 in enumerate(rg.atom_groups()):
    for j, ag2 in enumerate(rg.atom_groups()):
      if j>i: break
      atom1, name1, neutron1 = _get_atom_neutron(ag1,
                                                 lookup_name1,
                                                 bondlength=.5)
      if atom1 is None: continue
      if verbose: print('1',atom1.quote(), name1, neutron1)
      atom3, name3, neutron3 = _get_atom_neutron(ag2,
                                                 lookup_name3,
                                                 bondlength=.5)
      if atom3 is None: continue
      if verbose: print('3',atom3.quote(), name3, neutron3)
      for ag in [ag1, ag2]:
        atom2, name2, neutron2 = _get_atom_neutron(ag,
                                                   lookup_name2)
        if atom2: break
      if atom2 is None: continue
      if verbose: print('2',atom2.quote(), name2, neutron2)

      yield atom1, atom2, atom3, name1, name2, name3, neutron1, neutron2, neutron3

def _generate_dihedral_atoms(rg, lookup_name1, lookup_name2, lookup_name3, lookup_name4, verbose=False):
  for i, ag1 in enumerate(rg.atom_groups()):
    for j, ag2 in enumerate(rg.atom_groups()):
      # if j>i: break
      atom1, name1, neutron1 = _get_atom_neutron(ag1,
                                                 lookup_name1,
                                                 bondlength=.5)
      if atom1 is None: continue
      if verbose: print('1',atom1.quote(), name1, neutron1)
      atom4, name4, neutron4 = _get_atom_neutron(ag2,
                                                 lookup_name4,
                                                 bondlength=.5)
      if atom4 is None: continue
      if verbose: print('4',atom4.quote(), name4, neutron4)
      for ag in [ag1, ag2]:
        atom2, name2, neutron2 = _get_atom_neutron(ag,
                                                   lookup_name2)
        if atom2: break
      if atom2 is None: continue
      if verbose: print('2',atom2.quote(), name2, neutron2)

      for ag in [ag1, ag2]:
        atom3, name3, neutron3 = _get_atom_neutron(ag,
                                                   lookup_name3)
        if atom3: break
      if atom3 is None: continue
      if verbose: print('3',atom3.quote(), name3, neutron3)

      yield atom1, atom2, atom3, atom4, name1, name2, name3, name4, neutron1, neutron2, neutron3, neutron4

def adjust_geometry_proxies_registeries(hierarchy,
                                        gpr,
                                        error_i_seqs,
                                        log=None,
                                        verbose=False,
                                        ):
  t0=time.time()
  mon_lib_srv = server.server()
  pdb_atoms = hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  resnames=[]
  angle_counters = [0,0]
  checked=[]
  atoms_added={}
  atoms_added_energy={}
  for i_seq in error_i_seqs:
    atom = pdb_atoms[i_seq]
    ag = atom.parent()
    if ag.resname in checked: continue
    rg = ag.parent()
    # need to be able to check in user defined location
    for pH_range in ['neutron', 'low']:
      monomer_restraints = mon_lib_srv.get_comp_comp_id_direct(
        ag.resname,
        pH_range=pH_range,
        )
      if monomer_restraints: break
    checked.append(ag.resname)
    if monomer_restraints is None: continue
    atom_dict = monomer_restraints.atom_dict()
    resnames.append('"%s%s %s%5s"' % (' ',
                                      ag.resname,
                                      rg.parent().id,
                                      rg.resseq,
                                      ))
    bond_counters = [0,0]
    atoms_added = {}
    for bond in monomer_restraints.bond_list:
      for rc in _generate_bond_atoms(rg,
                                     bond.atom_id_1,
                                     bond.atom_id_2,
                                     bondlength=bond.value_dist,
                                     verbose=verbose):
        atom1, atom2, name1, name2, neutron1, neutron2 = rc
        bc, ade, ad = process_bonds(gpr,
                               bond,
                               atom_dict,
                               atom1,
                               atom2,
                               name1,
                               name2,
                               neutron1,
                               neutron2,
                               atoms=pdb_atoms,
                               verbose=verbose,
                               )
        for i in range(2):
          bond_counters[i]+=bc[i]
        if ade:
          atoms_added_energy.update(ade)
        if ad:
          atoms_added.update(ad)
    lookup={}
    for angle in monomer_restraints.angle_list:
      for rc in _generate_angle_atoms(rg,
                                      angle.atom_id_1,
                                      angle.atom_id_2,
                                      angle.atom_id_3,
                                      verbose=verbose):
        atom1, atom2, atom3, name1, name2, name3, neutron1, neutron2, neutron3 = rc
        i_seqs = (atom1.i_seq, atom2.i_seq, atom3.i_seq)
        lookup[i_seqs]=angle
        i_seqs = (atom3.i_seq, atom2.i_seq, atom1.i_seq)
        lookup[i_seqs]=angle
    for angle_proxy in gpr.angle.proxies:
      if angle_proxy.i_seqs in lookup:
        i_seqs = angle_proxy.i_seqs
        angle = lookup[i_seqs]
        diff=abs(angle.value_angle-angle_proxy.angle_ideal)
        if verbose: print('ANG',diff,angle_proxy.angle_ideal, angle.value_angle,angle.value_angle_esd)
        angle_proxy.angle_ideal = angle.value_angle
        angle_proxy.weight = 1/angle.value_angle_esd**2
        if diff>angle.value_angle_esd: angle_counters[0]+=1
        del lookup[i_seqs]
        i_seqs = list(i_seqs)
        i_seqs.reverse()
        del lookup[tuple(i_seqs)]
    if lookup:
      done = []
      for i_seqs in lookup:
        if i_seqs in done: continue
        angle = lookup[i_seqs]
        proxy =  geometry_restraints.angle_proxy(
          i_seqs=i_seqs,
          angle_ideal=angle.value_angle,
          weight=1/angle.value_angle_esd**2,
          )
        gpr.angle.add_if_not_duplicated(proxy)
        angle_counters[1]+=1
        i_seqs=list(i_seqs)
        i_seqs.reverse()
        done.append(tuple(i_seqs))
  if resnames:
    print("\n  Adjusted restraints in %d residue(s) for low pH or neutron in %0.1fs" % (
      len(resnames),
      time.time()-t0,
      ), file=log)
    print("    Residues changed", file=log)
    for resname in resnames:
      print("      %s" % resname, file=log)
    print("    Changed (significantly) %d bond restraint(s),  added %d bond restraint(s)" % (
      bond_counters[0],
      bond_counters[1],
      ), file=log)
    for i_seqs, atoms in atoms_added.items():
      print('      New atom/bond %s\n%s%s' % (atoms[0].quote(),' '*20, atoms[1].quote()), file=log)
    print("    Changed (significantly) %d angle restraint(s), added %d angle restraint(s)\n" % (
      angle_counters[0],
      angle_counters[1],
      ), file=log)
  #else:
  #  print >> log, "  Time to perform restraint checks: %0.1f" % (time.time()-t0)
  return atoms_added_energy

def adjust_geometry_restraints_manager(hierarchy,
                                       grm,
                                       error_i_seqs,
                                       log=None,
                                       ):
  # obsolete
  assert 0
  t0=time.time()
  mon_lib_srv = server.server()
  pdb_atoms = hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  resnames=[]
  bond_counters = [0,0]
  angle_counters = [0,0]
  checked=[]
  for i_seq in error_i_seqs:
    atom = pdb_atoms[i_seq]
    ag = atom.parent()
    if ag.resname in checked: continue
    rg = ag.parent()
    # need to be able to check in user defined location
    monomer_restraints = mon_lib_srv.get_comp_comp_id_direct(
      ag.resname,
      pH_range="low",
      )
    checked.append(ag.resname)
    if monomer_restraints is None: continue
    resnames.append('"%s%s %s%5s"' % (' ',
                                      ag.resname,
                                      rg.parent().id,
                                      rg.resseq,
                                      ))
    for bond in monomer_restraints.bond_list:
      bond.show()
      atom1 = ag.get_atom(bond.atom_id_1)
      atom2 = ag.get_atom(bond.atom_id_2)
      i_seqs = (atom1.i_seq, atom2.i_seq)
      bond_param = grm.bond_params_table.lookup(*list(i_seqs))
      if bond_param:
        bond_param.distance_ideal = bond.value_dist
        bond_counters[0]+=1
      else:
        proxy = geometry_restraints.bond_simple_proxy(
          i_seqs=i_seqs,
          distance_ideal=bond.value_dist,
          weight=1/(bond.value_dist_esd**2),
          )
        grm.add_new_bond_restraints_in_place([proxy], sites_cart)
        bond_counters[1]+=1
    lookup={}
    for angle in monomer_restraints.angle_list:
      atom1 = ag.get_atom(angle.atom_id_1)
      atom2 = ag.get_atom(angle.atom_id_2)
      atom3 = ag.get_atom(angle.atom_id_3)
      i_seqs = (atom1.i_seq, atom2.i_seq, atom3.i_seq)
      lookup[i_seqs]=angle
      i_seqs = (atom3.i_seq, atom2.i_seq, atom1.i_seq)
      lookup[i_seqs]=angle
    for angle_proxy in grm.angle_proxies:
      if angle_proxy.i_seqs in lookup:
        i_seqs = angle_proxy.i_seqs
        angle = lookup[i_seqs]
        angle_proxy.angle_ideal = angle.value_angle
        angle_proxy.weight = 1/angle.value_angle_esd**2
        angle_counters[0]+=1
        del lookup[i_seqs]
        i_seqs = list(i_seqs)
        i_seqs.reverse()
        del lookup[tuple(i_seqs)]
    if lookup:
      done = []
      for i_seqs in lookup:
        if i_seqs in done: continue
        proxy =  geometry_restraints.angle_proxy(
          i_seqs=i_seqs,
          angle_ideal=angle.value_angle,
          weight=1/angle.value_angle_esd**2,
          )
        grm.add_angles_in_place([proxy])
        angle_counters[1]+=1
        i_seqs=list(i_seqs)
        i_seqs.reverse()
        done.append(tuple(i_seqs))
