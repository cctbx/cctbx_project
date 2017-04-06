import os, sys
import time

from cctbx.array_family import flex
from scitbx.math import superpose
from mmtbx.conformation_dependent_library import mcl_sf4_coordination

def update(grm,
           pdb_hierarchy,
           link_records=None,
           log=sys.stdout,
           verbose=False,
           ):
  # SF4
  rc = mcl_sf4_coordination.get_sulfur_iron_cluster_coordination(
    pdb_hierarchy=pdb_hierarchy,
    nonbonded_proxies=grm.pair_proxies(
      sites_cart=pdb_hierarchy.atoms().extract_xyz()).nonbonded_proxies,
    verbose=verbose,
  )
  if link_records is None: link_records={}
  link_records.setdefault('LINK', [])
  bproxies, aproxies = mcl_sf4_coordination.get_all_proxies(rc)
  if len(bproxies):
    print >> log, "  SF4 coordination"
    atoms = pdb_hierarchy.atoms()
    sf4_coordination = {}
    for bp in bproxies:
      sf4_ag = atoms[bp.i_seqs[0]].parent()
      sf4_coordination.setdefault(sf4_ag.id_str(), [])
      sf4_coordination[sf4_ag.id_str()].append((atoms[bp.i_seqs[0]], 
                                                atoms[bp.i_seqs[1]]))
      link = (atoms[bp.i_seqs[0]], atoms[bp.i_seqs[1]], 'x,y,z')
      if link not in link_records: link_records['LINK'].append(link)
    for sf4, aas in sorted(sf4_coordination.items()):
      print >> log, '    %s' % sf4
      for aa in sorted(aas):
        print >> log, '       %s - %s' % (aa[0].id_str(), aa[1].id_str())
    print >> log
  grm.add_new_bond_restraints_in_place(
    proxies=bproxies,
    sites_cart=pdb_hierarchy.atoms().extract_xyz(),
  )
  grm.add_angles_in_place(aproxies)

def _extract_sites_cart(ag, element=None):
  selection = []
  for atom in ag.atoms():
    if element and atom.element.upper()!=element.upper(): continue
    selection.append(atom.xyz)
  return flex.vec3_double(selection)

def generate_sites_fixed(pdb_hierarchy, resname, element=None):
  for ag in pdb_hierarchy.atom_groups():
    if ag.resname.strip().upper()==resname.upper():
      yield _extract_sites_cart(ag, element), ag

def superpose_ideal_sf4_coordinates(pdb_hierarchy, resname='SF4'):
  from iotbx import pdb
  t0=time.time()
  rmsd_list = {}
  if resname=='SF4':
    ideal = pdb.input(lines=mcl_sf4_coordination.ideal_sf4,
                      source_info='ideal',
                    )
  elif resname=='F3S':
    ideal = pdb.input(lines=mcl_sf4_coordination.ideal_f3s,
                      source_info='ideal',
                    )
  else:
    assert 0
  ideal_hierarchy = ideal.construct_hierarchy()
  sites_moving = _extract_sites_cart(ideal_hierarchy, 'Fe')
  for ideal_ag in ideal_hierarchy.atom_groups(): break
  for sites_fixed, ag in generate_sites_fixed(pdb_hierarchy, resname, 'Fe'):
    assert sites_fixed.size() == sites_moving.size(), '%(resname)s residue is missing atoms' % locals()
    lsq_fit = superpose.least_squares_fit(
      reference_sites = sites_fixed,
      other_sites     = sites_moving)
    new_atoms = ideal_ag.detached_copy().atoms()
    sites_new = new_atoms.extract_xyz()
    sites_new = lsq_fit.r.elems * sites_new + lsq_fit.t.elems
    rmsd = sites_fixed.rms_difference(lsq_fit.other_sites_best_fit())
    rmsd_list[ag.id_str()] = rmsd
    new_atoms.set_xyz(sites_new)
    for atom1 in ag.atoms():
      for atom2 in new_atoms:
        if atom1.name.strip()==atom2.name.strip():
          atom1.xyz=atom2.xyz
          break
      else:
        assert 0, 'not all atoms updated - missing %s' % atom1.quote()
  outl = ''
  if rmsd_list:
    outl = '  %(resname)s Regularisation' % locals()
    outl+= '\n    residue        rmsd'
    for id_str, rmsd in sorted(rmsd_list.items()):
      outl += '\n    "%s"   %0.1f' % (id_str, rmsd)
    outl += '\n  Time to superpose : %0.2fs\n' % (time.time()-t0)
  return outl

def superpose_ideal_residue_coordinates(pdb_hierarchy, resname='SF4'):
  if resname in ['SF4', 'F3S']:
    rc = superpose_ideal_sf4_coordinates(pdb_hierarchy, resname=resname)
  else:
    assert 0
  return rc

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
