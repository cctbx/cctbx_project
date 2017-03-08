import os, sys

from cctbx.array_family import flex
from scitbx.math import superpose
from mmtbx.conformation_dependent_library import mcl_sf4_coordination

def update(grm,
           pdb_hierarchy,
           log=sys.stdout,
           verbose=False,
           ):
  # SF4
  rc = mcl_sf4_coordination.get_sf4_coordination(
    pdb_hierarchy=pdb_hierarchy,
    nonbonded_proxies=grm.pair_proxies(
      sites_cart=pdb_hierarchy.atoms().extract_xyz()).nonbonded_proxies,
    verbose=verbose,
  )
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

def superpose_ideal_sf4_coordinates(pdb_hierarchy):
  from iotbx import pdb
  rmsd_list = {}
  ideal = pdb.input(lines=mcl_sf4_coordination.ideal_sf4,
                    source_info='ideal',
                  )
  ideal_hierarchy = ideal.construct_hierarchy()
  sites_moving = _extract_sites_cart(ideal_hierarchy, 'Fe')
  for ideal_ag in ideal_hierarchy.atom_groups(): break
  for sites_fixed, ag in generate_sites_fixed(pdb_hierarchy, 'SF4', 'Fe'):
    assert sites_fixed.size() == sites_moving.size(), 'SF4 residue is missing atoms'
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
        assert 0, 'not all atoms updated'
  outl = '\n  SF4 Regularisation'
  outl+= '\n    residue        rmsd'
  for id_str, rmsd in sorted(rmsd_list.items()):
    outl += '\n    "%s"   %0.1f' % (id_str, rmsd)
  outl += '\n'
  return outl

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
