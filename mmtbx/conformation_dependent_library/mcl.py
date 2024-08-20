from __future__ import absolute_import, division, print_function
import sys
import time

from cctbx.array_family import flex
from scitbx.math import superpose
from mmtbx.conformation_dependent_library import mcl_sf4_coordination
from six.moves import range
from mmtbx.conformation_dependent_library import metal_coordination_library

def get_pdb_hierarchy_from_restraints(code):
  from mmtbx.monomer_library import server
  from iotbx import pdb
  mon_lib_server = server.server()
  path = mon_lib_server.get_comp_comp_id_direct(code, return_filename=True)
  cif_obj = server.read_cif(path)
  ligand_inp=pdb.pdb_input(source_info="Model from %s" % path,
                          lines=flex.split_lines(""))
  ligand_hierarchy = ligand_inp.construct_hierarchy()
  model=pdb.hierarchy.model()
  chain=pdb.hierarchy.chain()
  chain.id='Z'
  rg=pdb.hierarchy.residue_group()
  ag=pdb.hierarchy.atom_group()
  for block, loops in cif_obj.blocks.items():
    if block=='comp_list': continue
    for loop in loops.iterloops():
      for row in loop.iterrows():
        if '_chem_comp_atom.comp_id' not in row: break
        ag.resname = row['_chem_comp_atom.comp_id']
        atom = pdb.hierarchy.atom()
        atom.name = row['_chem_comp_atom.atom_id']
        atom.element = '%2s' % row['_chem_comp_atom.type_symbol']
        atom.xyz = (
          float(row['_chem_comp_atom.x']),
          float(row['_chem_comp_atom.y']),
          float(row['_chem_comp_atom.z']),
                )
        ag.append_atom(atom)
  rg.append_atom_group(ag)
  chain.append_residue_group(rg)
  model.append_chain(chain)
  ligand_hierarchy.append_model(model)
  ligand_hierarchy.atoms().reset_i_seq()
  return ligand_hierarchy

def update(grm,
           pdb_hierarchy,
           link_records=None,
           log=sys.stdout,
           verbose=False,
           ):
  def _atom_id(a, show_i_seq=False):
    if show_i_seq:
      return '%s (%5d)' % (a.id_str(), a.i_seq)
    else:
      return '%s' % (a.id_str())
  if link_records is None: link_records={}
  link_records.setdefault('LINK', [])
  hooks = [
    ["Iron sulfur cluster coordination",
     mcl_sf4_coordination.get_sulfur_iron_cluster_coordination,
     mcl_sf4_coordination.get_all_proxies,
      ],
    ['Zn2+ tetrahedral coordination',
     metal_coordination_library.get_metal_coordination_proxies,
     metal_coordination_library.get_proxies_zn,
      ],
    # ['Mg2+ Nucleotide coordination',
    #  metal_coordination_library.get_metal_coordination_proxies,
    #  metal_coordination_library.get_proxies_mg_nuc,
    #   ],
    ]
  outl = ''
  outl_debug = ''

  sites_c = pdb_hierarchy.atoms().extract_xyz()
  nb_proxies = grm.pair_proxies(
        sites_cart=sites_c).nonbonded_proxies
  for label, get_coordination, get_all_proxies in hooks:
    rc = get_coordination(
      pdb_hierarchy=pdb_hierarchy,
      nonbonded_proxies=nb_proxies,
      verbose=verbose,
    )
    bproxies, aproxies = get_all_proxies(rc)
    if bproxies is None: continue
    if len(bproxies):
      outl += '    %s\n' % label
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
        outl += '%spdb="%s"\n' % (' '*6, sf4)
        outl_debug += '%spdb="%s"\n' % (' '*6, sf4)
        for aa in aas:
          outl += '%s%s - %s\n' % (' '*8, _atom_id(aa[0]), _atom_id(aa[1]))
          outl_debug += '%s%s - %s\n' % (' '*8,
                                         _atom_id(aa[0], True),
                                         _atom_id(aa[1], True))
    if bproxies:
      if verbose:
        atoms = pdb_hierarchy.atoms()
        for bp in bproxies:
          print(bp.i_seqs,
                atoms[bp.i_seqs[0]].quote(),
                atoms[bp.i_seqs[1]].quote(),
                bp.rt_mx_ji,
                )
      try:
        grm.add_new_bond_restraints_in_place(
          proxies=bproxies,
          sites_cart=pdb_hierarchy.atoms().extract_xyz(),
        )
      except RuntimeError as e:
        print('\n\n%s' % outl_debug)
        raise e
    #
    done = []
    remove = []
    for i, angle in enumerate(aproxies):
      i_seqs = list(angle.i_seqs)
      i_seqs.sort()
      if i_seqs in done:
        remove.append(i)
      else:
        done.append(i_seqs)
    if remove:
      remove.reverse()
      for r in remove:
        del aproxies[r]
    #
    if aproxies:
      outl += '%s%s' % (' '*6, 'Number of angles added : %d\n' % len(aproxies))
    grm.add_angles_in_place(aproxies)
  if outl:
    print('  Dynamic metal coordination', file=log)
    print(outl, file=log)

def _extract_sites_cart(ag, element=None):
  selection = []
  for atom in ag.atoms():
    if element and atom.element.upper().strip()!=element.upper().strip():
      continue
    selection.append(atom.xyz)
  return flex.vec3_double(selection)

def generate_sites_fixed(pdb_hierarchy, resname, element=None):
  for ag in pdb_hierarchy.atom_groups():
    if ag.resname.strip().upper()==resname.upper():
      yield _extract_sites_cart(ag, element), ag

def superpose_ideal_residue_coordinates(pdb_hierarchy,
                                        resname,
                                        superpose_element=None,
                                        ):
  element_lookup = {'SF4' : 'Fe',
                    'F3S' : 'S',
                    #'F4S' : 'S', # not done yet
                    #'CLF' : 'Fe', # too flexible
                    'DVT' : 'V',
                    }
  from mmtbx.monomer_library import pdb_interpretation
  t0=time.time()
  rmsd_list = {}
  if superpose_element is None:
    superpose_element = element_lookup.get(resname, None)
  if resname in pdb_interpretation.ideal_ligands:
    ideal_hierarchy = get_pdb_hierarchy_from_restraints(resname)
  else:
    assert 0
  sites_moving = _extract_sites_cart(ideal_hierarchy, superpose_element)
  assert len(sites_moving), 'No atoms %s found' % superpose_element
  for ideal_ag in ideal_hierarchy.atom_groups(): break
  for sites_fixed, ag in generate_sites_fixed(pdb_hierarchy,
                                              resname,
                                              superpose_element,
                                              ):
    assert sites_fixed.size() == sites_moving.size(), '%(resname)s residue is missing atoms' % locals()
    import random
    min_rmsd = 1e9
    min_sites_cart = None
    for i in range(100):
      random.shuffle(sites_moving)
      lsq_fit = superpose.least_squares_fit(
        reference_sites = sites_fixed,
        other_sites     = sites_moving)
      new_atoms = ideal_ag.detached_copy().atoms()
      sites_new = new_atoms.extract_xyz()
      sites_new = lsq_fit.r.elems * sites_new + lsq_fit.t.elems
      rmsd = sites_fixed.rms_difference(lsq_fit.other_sites_best_fit())
      if rmsd<min_rmsd:
        min_rmsd=rmsd
        min_sites_cart = sites_new
    rmsd_list[ag.id_str()] = min_rmsd
    sites_new = min_sites_cart
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
    outl = '\n  %(resname)s Regularisation' % locals()
    outl+= '\n    residue        rmsd'
    for id_str, rmsd in sorted(rmsd_list.items()):
      outl += '\n    "%s"   %0.1f' % (id_str, rmsd)
    outl += '\n  Time to superpose : %0.2fs\n' % (time.time()-t0)
  return outl

def superpose_ideal_ligand_on_poor_ligand(ideal_hierarchy,
                                          poor_hierarchy,
                                          ):
  """Function superpose an ideal ligand onto the mangled ligand from a
     ligand fitting procedure

  Args:
      ideal_hierarchy (pdb_hierarchy): Ideal ligand
      poor_hierarchy (pdb_hierarchy): Poor ligand with correct c.o.m. and same
        atom names in order. Could become more sophisticated.
  """
  sites_moving = flex.vec3_double()
  sites_fixed = flex.vec3_double()
  for atom1, atom2 in zip(ideal_hierarchy.atoms(), poor_hierarchy.atoms()):
    assert atom1.name==atom2.name, '%s!=%s' % (atom1.quote(),atom2.quote())
    sites_moving.append(atom1.xyz)
    sites_fixed.append(atom2.xyz)
  lsq_fit = superpose.least_squares_fit(
        reference_sites = sites_fixed,
        other_sites     = sites_moving)
  sites_new = ideal_hierarchy.atoms().extract_xyz()
  sites_new = lsq_fit.r.elems * sites_new + lsq_fit.t.elems
  # rmsd = sites_fixed.rms_difference(lsq_fit.other_sites_best_fit())
  ideal_hierarchy.atoms().set_xyz(sites_new)
  return ideal_hierarchy

if __name__=="__main__":
  from iotbx import pdb
  ideal_inp=pdb.pdb_input(sys.argv[1])
  ideal_hierarchy = ideal_inp.construct_hierarchy()
  poor_inp=pdb.pdb_input(sys.argv[2])
  poor_hierarchy = poor_inp.construct_hierarchy()
  ideal_hierarchy = superpose_ideal_ligand_on_poor_ligand(ideal_hierarchy, poor_hierarchy)
  ideal_hierarchy.write_pdb_file('new.pdb')
