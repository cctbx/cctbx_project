from __future__ import absolute_import, division, print_function
import sys
import time

from cctbx.array_family import flex
from scitbx.math import superpose
from mmtbx.conformation_dependent_library import mcl_sf4_coordination
from six.moves import range

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
    print("  SF4/F3S coordination", file=log)
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
      print('    %s' % sf4, file=log)
      for aa in sorted(aas):
        print('       %s - %s' % (aa[0].id_str(), aa[1].id_str()), file=log)
    print(file=log)
  grm.add_new_bond_restraints_in_place(
    proxies=bproxies,
    sites_cart=pdb_hierarchy.atoms().extract_xyz(),
  )
  grm.add_angles_in_place(aproxies)

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
  from iotbx import pdb
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

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
