from __future__ import absolute_import,division, print_function
from io import StringIO
import time
from math import sqrt

from libtbx import Auto
from libtbx.str_utils import make_header
from libtbx.utils import Sorry
from libtbx.utils import null_out

from cctbx.geometry_restraints.manager import manager as standard_manager
from mmtbx.geometry_restraints import quantum_interface
from mmtbx.geometry_restraints import qm_manager

class manager(standard_manager):
  def __init__(self,
               params,
               log=StringIO()
               ):
    pass

  def __repr__(self):
    return 'QM restraints manager'

  def cleanup(self, log=StringIO()):
    make_header('Cleaning up - QM restraints', out=log)
    assert 0

def digester(model,
             standard_geometry_restraints_manager,
             params,
             log=StringIO(),
             ):
  #
  # Digest
  #
  sgrm = standard_geometry_restraints_manager
  qm_grm = manager(params, log=log)
  for attr, value in list(vars(sgrm).items()):
    if attr.startswith('__'): continue
    setattr(qm_grm, attr, value)
  qm_grm.standard_geometry_restraints_manager = sgrm
  #
  make_header('QM Restraints Initialisation', out=log)
  qm_restraints_initialisation(params, log=log)
  run_program=True
  for i, qmr in enumerate(params.qi.qm_restraints):
    if qmr.run_in_macro_cycles=='test': break
  else: run_program=False
  update_restraints(model,
                    params,
                    run_program=run_program,
                    log=log,
                    )
  return qm_grm

def write_pdb_file(model, filename, log):
  outl = model.model_as_pdb()
  print('    Writing ligand : %s' % filename, file=log)
  f=open(filename, 'w')
  f.write(outl)
  del f

def add_hydrogen_atoms_to_model(model,
                                use_capping_hydrogens=False,
                               # use_neutron_distances=True,
                               # n_terminal_charge=True,
                               ):
  from mmtbx.ligands.ready_set_utils import add_terminal_hydrogens
  from mmtbx.ligands.ready_set_utils import add_water_hydrogen_atoms_simple
  rc = add_terminal_hydrogens( model.get_hierarchy(),
                               model.get_restraints_manager().geometry,
                               use_capping_hydrogens=use_capping_hydrogens,
                               )
  assert not rc
  rc = add_water_hydrogen_atoms_simple(model.get_hierarchy())

def select_and_reindex(model,
                       selection_str=None,
                       selection_array=None,
                       reindex=True,
                       verbose=False):
  from cctbx.array_family import flex
  from scitbx_array_family_flex_ext import reindexing_array
  def _reindexing(mod, sel, verbose=False):
    isel = sel.iselection()
    ra = reindexing_array(len(sel),
                          flex.size_t(isel).as_int())
    atoms = mod.get_atoms()
    for i_seq in isel:
      atoms[ra[i_seq]].tmp=i_seq
    if verbose:
      for atom in mod.get_atoms():
        print(atom.quote(), atom.i_seq, atom.tmp)
    return mod
  assert (selection_array, selection_str).count(None)==1
  if selection_str:
    selection_array = model.selection(selection_str)
  selected_model = model.select(selection_array)
  if reindex:
    rc = _reindexing(selected_model, selection_array, verbose=verbose)
  return selected_model

def super_cell_and_prune(buffer_model, ligand_model, buffer, prune_limit=5., write_steps=False):
  from cctbx.crystal import super_cell
  def dist2(r1,r2): return (r1[0]-r2[0])**2+(r1[1]-r2[1])**2+(r1[2]-r2[2])**2
  def min_dist2(residue_group1, residue_group2, limit=1e-2):
    min_d2=1e9
    min_ca_d2=1e9
    for atom1 in residue_group1.atoms():
      for atom2 in residue_group2.atoms(): # ligand but could be more than one rg!!!
        d2 = dist2(atom1.xyz, atom2.xyz)
        if d2<limit and min_ca_d2<1e9: return d2, min_ca_d2
        min_d2 = min(min_d2, d2)
        if atom1.name.strip()=='CA':
          min_ca_d2 = min(min_ca_d2, d2)
    return min_d2, min_ca_d2
  def find_movers(buffer_model, ligand_model, buffer, prune_limit=5.):
    buffer*=buffer
    prune_limit*=prune_limit
    movers = []
    close = []
    removers = []
    same = []
    prune_main = []
    for residue_group1 in buffer_model.get_hierarchy().residue_groups():
      for residue_group2 in ligand_model.get_hierarchy().residue_groups():
        if residue_group1.id_str()==residue_group2.id_str():
          same.append(residue_group1)
          continue
        if list(residue_group1.unique_resnames())==list(residue_group2.unique_resnames()):
          min_d2, min_ca_d2 = min_dist2(residue_group1, residue_group2, limit=.1)
          if min_d2<=.1:
            removers.append(residue_group1)
        else:
          min_d2, min_ca_d2 = min_dist2(residue_group1, residue_group2, limit=buffer)
        if min_d2<=buffer:
          close.append(residue_group1)
        else:
          movers.append(residue_group1)
        if min_ca_d2 is not None and min_ca_d2>prune_limit:
          prune_main.append(residue_group1)
      assert residue_group1 in movers or residue_group1 in close or residue_group1 in same
    for chain in buffer_model.get_hierarchy().chains():
      for rg in movers+removers:
        if rg in chain.residue_groups():
          chain.remove_residue_group(rg)
    # if prune_main: prune_main_chain(residue_group1, prune_main)
    return movers, removers
  def prune_main_chain(residue_group1, prune_main):
    mainchain = ['N', 'CA', 'C', 'O', 'OXT', 'H', 'HA', 'H1', 'H2', 'H3']
    for residue_group1 in prune_main:
      ur = residue_group1.unique_resnames()
      print(list(ur))
      if 'PRO' in ur or 'GLY' in ur:
        if list(ur)==['PRO']:
          assert 0
        continue
      cb_atom = None
      for atom in residue_group1.atoms():
        if atom.name.strip()=='CB':
          cb_atom=atom
          break
      if cb_atom is None: continue
      pruning=False
      for atom in residue_group1.atoms():
        if atom.name.strip()=='CA':
          atom.name=' HBC'
          atom.element='H'
          atom.xyz = ((atom.xyz[0]+cb_atom.xyz[0])/2,
                      (atom.xyz[1]+cb_atom.xyz[1])/2,
                      (atom.xyz[2]+cb_atom.xyz[2])/2,
                      )
          pruning=True
          break
      if pruning:
        for atom in residue_group1.atoms():
          if atom.name.strip() in mainchain:
              ag = atom.parent()
              ag.remove_atom(atom)

  complete_p1 = super_cell.run(
    pdb_hierarchy        = buffer_model.get_hierarchy(),
    crystal_symmetry     = buffer_model.crystal_symmetry(),
    select_within_radius = buffer,
    )
  if(write_steps):
    complete_p1.hierarchy.write_pdb_file(file_name="complete_p1.pdb",
      crystal_symmetry = complete_p1.crystal_symmetry)

  for atom1, atom2 in zip(buffer_model.get_atoms(), complete_p1.hierarchy.atoms()):
    atom2.tmp = atom1.tmp
  buffer_model._pdb_hierarchy = complete_p1.hierarchy
  movers, removers = find_movers(buffer_model, ligand_model, buffer)
  # write_pdb_file(buffer_model, 'test.pdb', None)

def get_ligand_buffer_models(model, qmr, verbose=False, write_steps=False):
  from cctbx.maptbx.box import shift_and_box_model
  ligand_model = select_and_reindex(model, qmr.selection)
  if len(ligand_model.get_atoms())==0:
    raise Sorry('selection "%s" results in empty model' % qmr.selection)
  buffer_selection_string = 'residues_within(%s, %s)' % (qmr.buffer,
                                                         qmr.selection)
  buffer_model = select_and_reindex(model, buffer_selection_string)
  buffer_model.remove_alternative_conformations(always_keep_one_conformer=True)
  for residue_group in buffer_model.get_hierarchy().residue_groups():
    if (residue_group.atom_groups_size() != 1):
      raise Sorry("Not implemented: cannot run QI on buffer "+
                  "molecules with alternate conformations")
  if write_steps: write_pdb_file(buffer_model, 'pre_super_cell.pdb', None)
  super_cell_and_prune(buffer_model, ligand_model, qmr.buffer, write_steps=write_steps)
  buffer_model_p1 = shift_and_box_model(model=buffer_model)
  for atom1, atom2 in zip(buffer_model_p1.get_atoms(), buffer_model.get_atoms()):
    atom1.tmp=atom2.tmp
  buffer_model = buffer_model_p1
  buffer_model.unset_restraints_manager()
  buffer_model.log=null_out()
  buffer_model.process(make_restraints=True)
  if write_steps: write_pdb_file(buffer_model, 'pre_add_terminii.pdb', None)
  add_hydrogen_atoms_to_model(buffer_model, use_capping_hydrogens=qmr.capping_groups)
  buffer_model.unset_restraints_manager()
  buffer_model.process(make_restraints=True)
  return ligand_model, buffer_model

def show_ligand_buffer_models(ligand_model, buffer_model):
  outl = '    Core atoms\n'
  ags = []
  for atom in ligand_model.get_atoms():
    outl += '      %s (%5d)\n' % (atom.id_str().replace('pdb=',''), atom.tmp)
  for atom in buffer_model.get_atoms():
    agi = atom.parent().id_str()
    if agi not in ags: ags.append(agi)
  outl += '    Buffer residues\n'
  for agi in ags:
    outl += '      %s\n' % agi
  return outl

def get_qm_manager(ligand_model, buffer_model, qmr, log=StringIO()):
  program = qmr.package.program
  if program=='test':
    qmm = qm_manager.base_manager
  elif program=='orca':
    qmm = qm_manager.orca_manager
    if qmr.package.method is Auto:
      qmr.package.method='AM1'
    if qmr.package.basis_set is Auto:
      qmr.package.basis_set=''
    if qmr.package.multiplicity is Auto:
      qmr.package.multiplicity=1
    if qmr.package.charge is Auto:
      qmr.package.charge=0
  else:
    assert 0

  total_charge = quantum_interface.electrons(buffer_model, log=log)
  if total_charge!=qmr.package.charge:
    print(u'update charge %s ~> %s' % (qmr.package.charge, total_charge),
          file=log)
    qmr.package.charge=total_charge
  qmm = qmm(buffer_model.get_atoms(),
            qmr.package.method,
            qmr.package.basis_set,
            '', #qmr.package.solvent_model,
            qmr.package.charge,
            qmr.package.multiplicity,
            # preamble='%02d' % (i+1),
            )
  qmm.program=program
  ligand_i_seqs = []
  for atom in ligand_model.get_atoms():
    ligand_i_seqs.append(atom.tmp)

  hd_selection = buffer_model.get_hd_selection()
  ligand_selection = []
  for i, (sel, atom) in enumerate(zip(hd_selection, buffer_model.get_atoms())):
    if atom.tmp in ligand_i_seqs:
      hd_selection[i]=True
      ligand_selection.append(True)
    else:
      ligand_selection.append(False)
  qmm.set_frozen_atoms(~hd_selection)
  qmm.set_interest_atoms(ligand_selection)
  return qmm

def get_all_xyz(inputs, nproc):
  from libtbx import easy_mp
  argss = []
  for i, (ligand_model, buffer_model, qmm, qmr) in enumerate(inputs):
    argss.append([qmm])
  print(argss)
  for args, res, err_str in easy_mp.multi_core_run(qm_manager.qm_runner,
                                                   tuple(argss),
                                                   nproc,
                                                   ):
    if res:
      print (args, res)
      # results.update(res)
      # results[args[0]]=res
    assert not err_str, 'Error %s' % err_str
  xyzs = []
  return xyzs

def qm_restraints_initialisation(params, log=StringIO()):
  for i, qmr in enumerate(params.qi.qm_restraints):
    if not i: print('  QM restraints selections', file=log)
    print('    %s - buffer %s (%s)' % (qmr.selection,
                                       qmr.buffer,
                                       qmr.capping_groups), file=log)
  print('',file=log)

def running_this_macro_cycle(qmr, macro_cycle):
  if qmr.run_in_macro_cycles=='all':
    return True
  elif qmr.run_in_macro_cycles=='first_only' and macro_cycle==1:
    return True
  elif qmr.run_in_macro_cycles=='test' and macro_cycle in [0, None]:
    return True
  else:
    return False

def write_restraints_from_model_via_grm(ligand_model, verbose=True):
  def atom_as_restraint(atom):
    resname = atom.parent().resname
    name = atom.name
    element = atom.element
    atom_type = element
    charge = atom.charge
    # outl = ' %(resname)s %(name)s %(element)s %(atom_type)s %(charge)s ' % locals()
    outl = ' %(resname)s %(name)s %(element)s' % locals()
    return outl
  def bond_as_restraint(bond, atom1, atom2):
    resname = atom.parent().resname
    name1 = atom1.name
    name2 = atom2.name
    ideal = bond.distance_ideal
    esd = 1/sqrt(bond.weight)
    outl = ' %(resname)s %(name1)s %(name2)s coval %(ideal).3f %(esd)s ' % locals()
    return outl
  ligand_grm = ligand_model.get_restraints_manager()
  atoms = ligand_model.get_atoms()
  outl = ''
  outl += '''
#
data_comp_%s
#''' % atoms[0].parent().resname
  outl += '''
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
'''
  for atom in atoms:
    outl += '%s\n' % atom_as_restraint(atom)
  bond_proxies_simple, asu = ligand_grm.geometry.get_all_bond_proxies(
    sites_cart=ligand_model.get_sites_cart())
  sorted_table, n_not_shown = bond_proxies_simple.get_sorted(
    'delta',
    ligand_model.get_sites_cart())
  outl += '''
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
'''
  for info in sorted_table:
    (i_seq, j_seq, i_seqs, distance_ideal, distance_model, slack, delta, sigma, weight, residual, sym_op_j, rt_mx) = info
    i_atom=atoms[i_seq]
    j_atom=atoms[j_seq]
    i_seqs=[i_seq, j_seq]
    bond=ligand_grm.geometry.bond_params_table.lookup(*list(i_seqs))
    outl += '%s\n' % bond_as_restraint(bond, i_atom, j_atom)
  # print(outl)

def update_restraints(model,
                      params,
                      macro_cycle=None,
                      run_program=True,
                      nproc=1,
                      log=StringIO(),
                      ):
  t0 = time.time()
  prefix = params.output.prefix
  if not model.restraints_manager_available():
    model.process(make_restraints=True)
  objects = []
  if quantum_interface.is_quantum_interface_active_this_macro_cycle(params,
                                                                    macro_cycle):
    print('  QM restraints calculations for macro cycle %s' % macro_cycle,
      file=log)
  for i, qmr in enumerate(params.qi.qm_restraints):
    if macro_cycle is not None and not running_this_macro_cycle(qmr, macro_cycle):
      print('    Skipping this selection in this macro_cycle : %s' % qmr.selection)
      continue
    preamble = quantum_interface.get_preamble(macro_cycle, i, qmr)
    #
    # get ligand and buffer region models
    #
    ligand_model, buffer_model = get_ligand_buffer_models(model, qmr)
    assert buffer_model.restraints_manager_available()
    if qmr.write_pdb_core:
      write_pdb_file(ligand_model, '%s_ligand_%s.pdb' % (prefix, preamble), log)
    if qmr.write_pdb_buffer:
      write_pdb_file(buffer_model, '%s_cluster_%s.pdb' % (prefix, preamble), log)
    #
    # get appropriate QM manager
    #
    qmm = get_qm_manager(ligand_model, buffer_model, qmr)
    qmm.preamble='%s_%s' % (prefix, preamble)
    objects.append([ligand_model, buffer_model, qmm, qmr])
  if not run_program: return
  print('',file=log)
  #
  # optimise
  #
  if nproc==1:
    xyzs = []
    xyzs_buffer = []
    for i, (ligand_model, buffer_model, qmm, qmr) in enumerate(objects):
      xyz, xyz_buffer = qm_manager.qm_runner(
        qmm,
        cleanup=qmr.cleanup,
        file_read=qmr.package.read_output_to_skip_opt_if_available,
        log=log,
        )
      print('  Time for calculation of "%s" using %s %s: %s' % (
        qmr.selection,
        qmr.package.method,
        qmr.package.basis_set,
        qmm.get_timings().split(':')[-1],
        ))
      xyzs.append(xyz)
      xyzs_buffer.append(xyz_buffer)
  else:
    print('nproc',nproc)
    assert 0
    xyzs, junk = get_all_xyz(objects, nproc)
  print('',file=log)
  #
  # update model restraints
  #
  for i, ((ligand_model, buffer_model, qmm, qmr), xyz, xyz_buffer) in enumerate(
                                                                    zip(objects,
                                                                        xyzs,
                                                                        xyzs_buffer,
                                                                        )):
    if qmr.package.view_output:
      qmm.view(qmr.package.view_output)
    if i: print(' ',file=log)
    print('  Updating QM restraints: "%s"' % qmr.selection, file=log)
    print(show_ligand_buffer_models(ligand_model, buffer_model), file=log)
    gs = ligand_model.geometry_statistics()
    print('  Starting stats: %s' % gs.show_bond_and_angle(), file=log)
    #
    # update coordinates of ligand
    #
    old = ligand_model.get_hierarchy().atoms().extract_xyz()
    rmsd = old.rms_difference(xyz)
    if rmsd>5:
      print('  QM minimisation has large rms difference in cartesian coordinates: %0.1f' % (rmsd),
            file=log)
      print('  Check the QM minimisation for errors or incorrect protonation.',
            file=log)
    ligand_model.get_hierarchy().atoms().set_xyz(xyz)
    old = buffer_model.get_hierarchy().atoms().extract_xyz()
    # rmsd = old.rms_difference(xyz_buffer)
    buffer_model.get_hierarchy().atoms().set_xyz(xyz_buffer)
    gs = ligand_model.geometry_statistics()
    print('  Interim stats : %s' % gs.show_bond_and_angle(), file=log)
    if qmr.write_final_pdb_core:
      write_pdb_file(ligand_model, '%s_ligand_final_%s.pdb' % (prefix, preamble), log)
    if qmr.write_final_pdb_buffer:
      write_pdb_file(buffer_model, '%s_cluster_final_%s.pdb' % (prefix, preamble), log)
    assert 0
    #
    # transfer geometry to proxies
    #  - bonds
    #
    model_grm = model.get_restraints_manager()
    ligand_grm = ligand_model.get_restraints_manager()
    buffer_grm = buffer_model.get_restraints_manager()
    atoms = ligand_model.get_atoms()
    #
    bond_proxies_simple, asu = ligand_grm.geometry.get_all_bond_proxies(
      sites_cart=ligand_model.get_sites_cart())
    sorted_table, n_not_shown = bond_proxies_simple.get_sorted(
      'delta',
      ligand_model.get_sites_cart())
    print('\n  Transfer', file=log)
    for info in sorted_table:
      (i_seq, j_seq, i_seqs, distance_ideal, distance_model, slack, delta, sigma, weight, residual, sym_op_j, rt_mx) = info
      i_atom=atoms[i_seq]
      j_atom=atoms[j_seq]
      i_seqs=[i_seq, j_seq]
      if i_atom.element_is_hydrogen(): continue
      if j_atom.element_is_hydrogen(): continue
      bond=ligand_grm.geometry.bond_params_table.lookup(*list(i_seqs))
      print('    %s - %s %5.2f ~> %5.2f' % (i_atom.id_str().replace('pdb=',''),
                                            j_atom.id_str().replace('pdb=',''),
                                            bond.distance_ideal,
                                            distance_model), file=log)
      bond.distance_ideal=distance_model
      i_seqs=[i_atom.tmp, j_atom.tmp]
      bond=model_grm.geometry.bond_params_table.lookup(*list(i_seqs))
      bond.distance_ideal=distance_model
    #
    #    - angles
    #
    sorted_table, n_not_shown = ligand_grm.geometry.angle_proxies.get_sorted(
      'delta',
      ligand_model.get_sites_cart())
    ligand_lookup = {}
    model_lookup = {}
    for info in sorted_table:
      (i_seqs, angle_ideal, angle_model, delta, sigma, weight, residual) = info
      i_atom=atoms[int(i_seqs[0])]
      j_atom=atoms[int(i_seqs[1])]
      k_atom=atoms[int(i_seqs[2])]
      if i_atom.element_is_hydrogen(): continue
      if j_atom.element_is_hydrogen(): continue
      if k_atom.element_is_hydrogen(): continue
      key = (int(i_seqs[0]), int(i_seqs[1]), int(i_seqs[2]))
      ligand_lookup[key]=angle_model
      print('    %s - %s - %s %5.2f ~> %5.2f' % (atoms[key[0]].id_str().replace('pdb=',''),
                                                 atoms[key[1]].id_str().replace('pdb=',''),
                                                 atoms[key[2]].id_str().replace('pdb=',''),
                                                 angle_ideal,
                                                 angle_model), file=log)
      key = (atoms[key[0]].tmp, atoms[key[1]].tmp, atoms[key[2]].tmp)
      model_lookup[key]=angle_model
      key = (int(i_seqs[1]), int(i_seqs[1]), int(i_seqs[0]))
      ligand_lookup[key]=angle_model
      key = (atoms[key[2]].tmp, atoms[key[1]].tmp, atoms[key[0]].tmp)
      model_lookup[key]=angle_model
    for angle_proxy in ligand_grm.geometry.angle_proxies:
      angle = ligand_lookup.get(angle_proxy.i_seqs, None)
      if angle is None: continue
      angle_proxy.angle_ideal=angle
    for angle_proxy in model_grm.geometry.angle_proxies:
      angle = model_lookup.get(angle_proxy.i_seqs, None)
      if angle is None: continue
      angle_proxy.angle_ideal=angle
    print('', file=log)
    # if qmr.write_restraints or 1:
    #   write_restraints_from_model_via_grm(ligand_model)
    #
    # final stats
    #
    gs = ligand_model.geometry_statistics()
    print('  Finished stats : %s' % gs.show_bond_and_angle(), file=log)
  print('  Total time for QM restaints: %0.1fs' % (time.time()-t0))
  print('='*80)

if __name__ == '__main__':
  print(quantum_chemistry_scope)
