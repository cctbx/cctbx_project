from __future__ import absolute_import, division, print_function
import sys

from cctbx.array_family import flex
from cctbx import geometry_restraints
from six.moves import range

def generate_atom_groups(pdb_hierarchy,
                         resnames=None,
                         ):
  assert resnames
  for atom_group in pdb_hierarchy.atom_groups():
    if atom_group.resname.strip() in resnames:
      print(atom_group)
      for atom in atom_group.atoms(): print(atom.quote())
      yield atom_group

def d2(atom1, atom2):
  from math import sqrt
  d2 = 0
  for i in range(3): d2 += (atom1.xyz[i]-atom2.xyz[i])**2
  return sqrt(d2)

def generate_proxies_from_geometry_restraints_manager_nonbonded_proxies(
        pdb_hierarchy,
        grm,
        log,
    ):
  #sites_cart = pdb_hierarchy.atoms().extract_xyz()
  nonbonded_proxies=grm.pair_proxies().nonbonded_proxies
  atoms = pdb_hierarchy.atoms()
  for atom in atoms: print(atom.format_atom_record())
  sites_cart = atoms.extract_xyz()
  get_sorted_result = nonbonded_proxies.get_sorted(
      by_value="delta",
      sites_cart=sites_cart)
  if get_sorted_result is None:
    return result
  sorted_nonb, n_not_shown = get_sorted_result
  print('1'*80)
  print(sorted_nonb)
  print(n_not_shown)
  # Get potential hbonds
  n_nonb = len(sorted_nonb)
  i = 0
  for tmp in sorted_nonb:
    print(tmp)
  while i < n_nonb: # and sorted_nonb[i][3] < 3.5:#hbond_distance_cutoff:
    (labels, i_seq, j_seq, dist, vdw_distance, sym_op_j, rt_mx) = sorted_nonb[i]
    a1 = atoms[i_seq]
    ag1 = a1.parent()
    a2 = atoms[j_seq]
    ag2 = a2.parent()
    if ag1.id_str()!=ag2.id_str():
      print(a1.quote(),a2.quote(),dist)
    i+=1
  assert 0

class minimum_distance_class(dict):
  def get_closest(self):
    assert 0

def generate_proxies_from_xray_structure(pdb_hierarchy,
                                         xray_structure,
                                         geometry_restraints_manager,
                                         log,
                                         ):
  assert 0
  for atom_group in generate_atom_groups(pdb_hierarchy, resnames=['CD']):
    selection = flex.bool(xray_structure.scatterers().size())
    for atom in atom_group.atoms(): selection[atom.i_seq]=True
    print(list(selection))
    selection_within = xray_structure.selection_within(3.5, selection)
    print(list(selection_within))
    #sel = selection.select(selection_within)
    sphere = pdb_hierarchy.select(selection_within)
    print(list(geometry_restraints_manager.nonbonded_charges))
    rc = geometry_restraints_manager.select(selection=selection_within)
    print(dir(rc))
    print(list(rc.nonbonded_charges))
    print(dir(rc))
    print(sphere)
    print(dir(sphere))
    sites_cart = sphere.atoms().extract_xyz()
    for atom in sphere.atoms(): print(atom.quote())
    #generate_proxies_from_geometry_restraints_manager_nonbonded_proxies(
    #    sphere, #pdb_hierarchy,
    #    rc,
    #    #sites_cart,
    #    log,
    #    )
    #assert 0
    #sphere.show()
    print('-'*80)
    for atom in atom_group.atoms():
      print(atom.format_atom_record())
      #if atom.name.strip().find("S"): continue
      min_dist = minimum_distance_class()
      for other in sphere.atoms():
        print(other.format_atom_record(), d2(atom, other))
        if atom.i_seq==other.i_seq: continue
        min_dist[d2(atom, other)] = other
        min_dist["%s %s %s" % (other.parent().resname,
                               other.parent().parent().resseq,
                               other.parent().altloc,
                               )] = d2(atom, other)
        #print min_dist
        #print min_dist.get_closest()
        #continue
        proxy = geometry_restraints.bond_simple_proxy(
          i_seqs=(atom.i_seq, other.i_seq),
          distance_ideal=1, #distance_ideal,
          weight=1, #weight/(sigma ** 2),
          #slack=slack,
          #top_out=top_out,
          #limit=limit,
        )
        #print dir(proxy)
        print(atom.i_seq, other.i_seq, proxy.rt_mx_ji)

def generate_proxies(pdb_hierarchy,
                     xray_structure,
                     geometry_restraints_manager,
                     log,
                     ):
  if 0:
    return generate_proxies_from_xray_structure(pdb_hierarchy,
                                                xray_structure,
                                                geometry_restraints_manager,
                                                log,
                                               )
  else:
    return generate_proxies_from_geometry_restraints_manager_nonbonded_proxies(
      pdb_hierarchy,
      geometry_restraints_manager,
      #pdb_hierarchy.atoms().extract_xyz(),
      log,
      )

def run(args):
  print("run",args)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
