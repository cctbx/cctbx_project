"""
Module to automatically attempt to load "low" pH range restraints from
restraints library - useful for various protonation states
"""
from __future__ import division
from __future__ import print_function
import time

from mmtbx.monomer_library import server
from cctbx import geometry_restraints

def _get_atom_neutron(ag, name, bondlength=None):
  atom = ag.get_atom(name)
  if atom: return atom, atom.name, False
  if bondlength<1:
    atom = ag.get_atom(name.replace('H', 'D'))
    if atom: return atom, name, True
  return atom, None, False

def adjust_geometry_proxies_registeries(hierarchy,
                                        #bond_params_table,
                                        #bond_asu_table,
                                        gpr,
                                        error_i_seqs,
                                        log=None,
                                        ):
  t0=time.time()
  mon_lib_srv = server.server()
  pdb_atoms = hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  resnames=[]
  bond_counters = [0,0]
  angle_counters = [0,0]
  checked=[]
  atoms_added={}
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
    atom_dict = monomer_restraints.atom_dict()
    resnames.append('"%s%s %s%5s"' % (' ',
                                      ag.resname,
                                      rg.parent().id,
                                      rg.resseq,
                                      ))
    for bond in monomer_restraints.bond_list:
      atom1, name1, neutron1 = _get_atom_neutron(ag,
                                                 bond.atom_id_1,
                                                 bondlength=bond.value_dist)
      if atom1 is None: continue
      atom2, name2, neutron2 = _get_atom_neutron(ag,
                                                 bond.atom_id_2,
                                                 bondlength=bond.value_dist)
      if atom2 is None: continue
      i_seqs = [atom1.i_seq, atom2.i_seq]
      k=0
      l=1
      bond_table_entry = gpr.bond_simple.table[i_seqs[k]]
      if ( not bond_table_entry or
           i_seqs[l] not in gpr.bond_simple.table[i_seqs[k]]):
        k=1
        l=0
        bond_table_entry = gpr.bond_simple.table[i_seqs[k]]
      if i_seqs[l] in bond_table_entry:
        bond_simple = gpr.bond_simple.proxies[i_seqs[k]]
        bond_simple.distance_ideal = bond.value_dist
        bond_simple.weight=1/bond.value_dist_esd**2
        bond_counters[0]+=1 # changed
      else:
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
        atoms_added[atom1.i_seq] = atom_dict.get(name1.strip(), None)
        atoms_added[atom2.i_seq] = atom_dict.get(name2.strip(), None)
        bond_counters[1]+=1
    lookup={}
    for angle in monomer_restraints.angle_list:
      atom1, name1, neutron1 = _get_atom_neutron(ag,
                                                 angle.atom_id_1,
      )
      if atom1 is None: continue
      atom2, name2, neutron2 = _get_atom_neutron(ag,
                                                 angle.atom_id_2,
      )
      if atom2 is None: continue
      atom3, name3, neutron3 = _get_atom_neutron(ag,
                                                 angle.atom_id_3,
      )
      if atom3 is None: continue
      i_seqs = (atom1.i_seq, atom2.i_seq, atom3.i_seq)
      lookup[i_seqs]=angle
      i_seqs = (atom3.i_seq, atom2.i_seq, atom1.i_seq)
      lookup[i_seqs]=angle
    for angle_proxy in gpr.angle.proxies:
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
    print("\n  Adjusted restraints in %d residue(s) for low pH in %0.1fs" % (
      len(resnames),
      time.time()-t0,
      ), file=log)
    print("    Residues changed", file=log)
    for resname in resnames:
      print("      %s" % resname, file=log)
    print("    Changed %d bond restraint(s),  added %d bond restraint(s)" % (
      bond_counters[0],
      bond_counters[1],
      ), file=log)
    print("    Changed %d angle restraint(s), added %d angle restraint(s)\n" % (
      angle_counters[0],
      angle_counters[1],
      ), file=log)
  #else:
  #  print >> log, "  Time to perform restraint checks: %0.1f" % (time.time()-t0)
  return atoms_added

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
