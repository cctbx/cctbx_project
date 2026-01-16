"""
Module to automatically adjust for enol peptide bond
"""
from __future__ import absolute_import, division, print_function
import time

# from mmtbx.conformation_dependent_library import generate_protein_threes

from mmtbx.monomer_library import server
from cctbx import geometry_restraints

from mmtbx.conformation_dependent_library.pH_dependent_restraints import _generate_bond_atoms
from mmtbx.conformation_dependent_library.pH_dependent_restraints import _generate_angle_atoms
from mmtbx.conformation_dependent_library.pH_dependent_restraints import _generate_dihedral_atoms

class empty:
  def __repr__(self):
    outl='empty'
    for attr, item in self.__dict__.items():
      outl+='\n  %s : %s' % (attr, item)
    return outl

def get_bond(name1, name2):
  key=[name1,name2]
  key.sort()
  bond=empty()
  values=enol_peptide_values.get(tuple(key), None)
  if values:
    bond.value_dist=values[0]
    bond.value_dist_esd=values[1]
    bond.value_dist_neutron=values[2]
  return bond

def get_angle(name1, name2, name3):
  key=[name1,name2,name3]
  if name1<name3:
    key=[name3,name2,name1]
  values=enol_peptide_values.get(tuple(key), None)
  angle=None
  if values:
    angle=empty()
    angle.value_angle=values[0]
    angle.value_angle_esd=values[1]
  return angle

def get_dihedral(name1, name2, name3, name4):
  key=[name1,name2,name3,name4]
  if name1<name4:
    key=[name4,name3,name2,name1]
  values=enol_peptide_values.get(tuple(key), None)
  angle=None
  if values:
    angle=empty()
    angle.value_angle=values[0]
    angle.value_angle_esd=values[1]
    angle.period=values[2]
  return angle

class restraint_subsitute(dict):
  def __init__(self):
    # enol peptide
    self['HNO']='O' # missing atom is bonded to...
    self[('HNO', 'O')] = [0.97, 0.02, 1.05]
    self[('C', 'O')] =   [1.43, 0.02, 1.43]
    self[('HNO', 'O', 'C')] = [120., 3.]
    self[('HNO', 'O', 'C', 'CA')] = [180., 10., 18]

  def _internals(self, number, func):
    for key, item in self.items():
      if type(key)!=type(tuple([])): continue
      if len(key)!=number: continue
      bond=func(*key)
      yield key, bond

  def bonds(self):
    for key, bond in self._internals(2, get_bond):
      yield key, bond

  def angles(self):
    for key, angle in self._internals(3, get_angle):
      yield key, angle

  def dihedrals(self):
    for key, angle in self._internals(4, get_dihedral):
      yield key, angle

enol_peptide_values=restraint_subsitute()

def process_bonds(gpr, bond, atom1, atom2, neutron1, neutron2, verbose=False):
  name1=atom1.name.strip()
  name2=atom2.name.strip()
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
    if verbose: print('DIFF',simple_proxy.distance_ideal,bond.value_dist,diff)
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
    item=empty()
    item.type_energy='HOH1'
    atoms_added_energy[atom1.i_seq] = item
    atoms_added_energy[atom2.i_seq] = item
    atoms_added[(atom1.i_seq, atom2.i_seq)] = (atom1, atom2)
    bond_counters[1]+=1
  return bond_counters, atoms_added_energy, atoms_added

def adjust_geometry_proxies_registeries(hierarchy,
                                        gpr,
                                        error_i_seqs,
                                        log=None,
                                        verbose=False,
                                        ):
  t0=time.time()
  # mon_lib_srv = server.server()
  pdb_atoms = hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  resnames=[]
  bond_counters = [0,0]
  angle_counters = [0,0]
  # checked=[]
  atoms_added={}
  atoms_added_energy={}
  for i_seq in error_i_seqs:
    error_atom = pdb_atoms[i_seq]
    if error_atom.name.strip() not in enol_peptide_values: continue
    ag = error_atom.parent()
    # if ag.resname in checked: continue
    rg = ag.parent()
    # need to be able to check in user defined location
    # checked.append(ag.resname)
    # if monomer_restraints is None: continue
    # atom_dict = monomer_restraints.atom_dict()
    resnames.append('"%s%s %s%5s"' % (' ',
                                      ag.resname,
                                      rg.parent().id,
                                      rg.resseq,
                                      ))

    for names, bond in enol_peptide_values.bonds():
      atom_id_1, atom_id_2 = names
      for rc in _generate_bond_atoms(rg,
                                     atom_id_1,
                                     atom_id_2,
                                     bondlength=bond.value_dist,
                                     verbose=verbose):
        atom1, atom2, name1, name2, neutron1, neutron2 = rc
        bc, ade, ad = process_bonds(gpr,
                               bond,
                               # atom_dict,
                               atom1,
                               atom2,
                               # name1,
                               # name2,
                               neutron1,
                               neutron2,
                               # atoms=pdb_atoms,
                               verbose=verbose,
                               )
        for i in range(2):
          bond_counters[i]+=bc[i]
        if ade:
          atoms_added_energy.update(ade)
        if ad:
          atoms_added.update(ad)

    lookup={}
    for names, angle in enol_peptide_values.angles():
      atom_id_1, atom_id_2, atom_id_3 = names
      for rc in _generate_angle_atoms(rg,
                                      atom_id_1,
                                      atom_id_2,
                                      atom_id_3,
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

    lookup={}
    for names, angle in enol_peptide_values.dihedrals():
      atom_id_1, atom_id_2, atom_id_3, atom_id_4 = names
      for rc in _generate_dihedral_atoms(rg,
                                          atom_id_1,
                                          atom_id_2,
                                          atom_id_3,
                                          atom_id_4,
                                          verbose=False):
        atom1, atom2, atom3, atom4, name1, name2, name3, name4, neutron1, neutron2, neutron3, neutron4 = rc
        i_seqs = (atom1.i_seq, atom2.i_seq, atom3.i_seq, atom4.i_seq)
        lookup[i_seqs]=angle
        i_seqs = (atom4.i_seq, atom3.i_seq, atom2.i_seq, atom1.i_seq)
        lookup[i_seqs]=angle
    for angle_proxy in gpr.dihedral.proxies:
      if angle_proxy.i_seqs in lookup:
        i_seqs = angle_proxy.i_seqs
        angle = lookup[i_seqs]
        diff=abs(angle.value_angle-angle_proxy.angle_ideal)
        if 1: print('ANG',diff,angle_proxy.angle_ideal, angle.value_angle,angle.value_angle_esd)
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
        proxy =  geometry_restraints.dihedral_proxy(
          i_seqs=i_seqs,
          angle_ideal=angle.value_angle,
          weight=1/angle.value_angle_esd**2,
          periodicity=angle.period,
          )
        gpr.dihedral.add_if_not_duplicated(proxy)
        angle_counters[1]+=1
        i_seqs=list(i_seqs)
        i_seqs.reverse()
        done.append(tuple(i_seqs))

  if resnames :
    print("\n  Adjusted restraints in %d residue(s) for enol-peptides in %0.1fs" % (
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
    print("    Changed (significantly) %d angle/dihedral restraint(s), added %d angle/dihedral restraint(s)\n" % (
      angle_counters[0],
      angle_counters[1],
      ), file=log)
  #else:
  #  print >> log, "  Time to perform restraint checks: %0.1f" % (time.time()-t0)
  # assert 0
  return atoms_added_energy

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      current_geometry=None, # xray_structure!!
                      sites_cart=None,
                      # cdl_proxies=None,
                      # ideal=True,
                      # esd=True,
                      log=None,
                      verbose=False,
                      ):
  from mmtbx.conformation_dependent_library import generate_protein_tuples
  from mmtbx.conformation_dependent_library.cdl_utils import get_enol_atoms
  # global registry
  # registry = RestraintsRegistry()
  # if current_geometry:
  #   assert not sites_cart
  #   sites_cart = current_geometry.sites_cart()
  # if sites_cart:
  #   pdb_atoms = hierarchy.atoms()
  #   # XXX PDB_TRANSITION VERY SLOW
  #   for j_seq, atom in enumerate(pdb_atoms):
  #     atom.xyz = sites_cart[j_seq]

  twos = None
  average_updates = 0
  total_updates = 0
  enols=0
  for twos in generate_protein_tuples(hierarchy,
                                        geometry,
                                        2,
                                        # cdl_class=True,
                                        # omega_cdl=True,
                                        retain_selection="name ca or name c or name n or name o or name cb or name h or name HNO",
                                        #verbose=verbose,
                                        ):
    if twos.cis_group():
      continue

    if twos.enol_group():
      rc=get_enol_atoms(*twos)
      atoms1, atoms2 = rc
      c1,ca1,n1=atoms1
      c2,ca2,n2=atoms2
      i_seqs=[c1.i_seq, n2.i_seq]
      bond=twos.bond_params_table.lookup(*i_seqs)
      bond.distance_ideal=1.27
      enols+=1
  geometry.reset_internals()
  return geometry, enols

