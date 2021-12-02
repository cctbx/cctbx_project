from __future__ import absolute_import,division, print_function
from io import StringIO

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
  qm_restraints_initialisation(params.qi, log=log)
  update_restraints(model,
                    params.qi,
                    run_program=False,
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
  # model.log=StringIO()
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

def super_cell_and_prune(buffer_model, ligand_model, buffer, prune_limit=5.):
  from cctbx.crystal import super_cell
  def dist2(r1,r2):
    return (r1[0]-r2[0])**2+(r1[1]-r2[1])**2+(r1[2]-r2[2])**2
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
    if prune_main: prune_main_chain(residue_group1, prune_main)
    return movers, removers
  def prune_main_chain(residue_group1, prune_main):
    mainchain = ['N', 'CA', 'C', 'O', 'OXT', 'H', 'HA', 'H1', 'H2', 'H3']
    for residue_group1 in prune_main:
      ur = residue_group1.unique_resnames()
      if 'PRO' in ur or 'GLY' in ur:
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
  def move_movers(buffer_model, complete_p1_hierarchy, movers):
    assert 0
    def _get_resname_from_residue_group(residue_group): pass
    for residue_group1 in movers:
      for residue_group2 in complete_p1_hierarchy.residue_groups():
        # print(dir(residue_group1))
        # print(list(residue_group1.unique_resnames()),residue_group1.resseq,residue_group1.parent().id)
        # print(list(residue_group2.unique_resnames()),residue_group2.resseq,residue_group2.parent().id)
        if residue_group1.parent().id==residue_group2.parent().id: continue
        if residue_group1.resseq==residue_group2.resseq:
          if list(residue_group1.unique_resnames())==list(residue_group2.unique_resnames()):
            xyzs = residue_group2.atoms().extract_xyz()
            residue_group1.atoms().set_xyz(xyzs)
    return movers

  complete_p1 = super_cell.run(
    pdb_hierarchy        = buffer_model.get_hierarchy(),
    crystal_symmetry     = buffer_model.crystal_symmetry(),
    select_within_radius = buffer,
    )
  # if(0):
  #   complete_p1.hierarchy.write_pdb_file(file_name="complete_p1.pdb",
  #     crystal_symmetry = complete_p1.crystal_symmetry)

  for atom1, atom2 in zip(buffer_model.get_atoms(), complete_p1.hierarchy.atoms()):
    atom2.tmp = atom1.tmp
  buffer_model._pdb_hierarchy = complete_p1.hierarchy
  movers, removers = find_movers(buffer_model, ligand_model, buffer)
  # write_pdb_file(buffer_model, 'test.pdb', None)

def get_ligand_buffer_models(model, qmr, verbose=False):
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
  add_hydrogen_atoms_to_model(buffer_model)
  buffer_model.unset_restraints_manager()
  buffer_model.log=null_out()
  buffer_model.process(make_restraints=True)
  super_cell_and_prune(buffer_model, ligand_model, qmr.buffer)
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
  for i, qmr in enumerate(params.qm_restraints):
    if not i: print('  QM restraints calculations', file=log)
    print('    %s' % qmr.selection, file=log)

def update_restraints(model,
                      params, # just the qi scope
                      macro_cycle=None,
                      run_program=True,
                      nproc=1,
                      log=StringIO(),
                      ):
  if not model.restraints_manager_available():
    model.process(make_restraints=True)
  objects = []
  for i, qmr in enumerate(params.qm_restraints):
    if not i: print('  QM restraints calculations', file=log)
    preamble = quantum_interface.get_preamble(macro_cycle, i, qmr.selection)
    #
    # get ligand and buffer region models
    #
    ligand_model, buffer_model = get_ligand_buffer_models(model, qmr)
    assert buffer_model.restraints_manager_available()
    if qmr.write_pdb_core:
      write_pdb_file(ligand_model, 'ligand_%s.pdb' % preamble, log)
    if qmr.write_pdb_buffer:
      write_pdb_file(buffer_model, 'cluster_%s.pdb' % preamble, log)
    #
    # get appropriate QM manager
    #
    qmm = get_qm_manager(ligand_model, buffer_model, qmr)
    qmm.preamble=preamble
    objects.append([ligand_model, buffer_model, qmm, qmr])
  print('',file=log)
  if not run_program: return
  #
  # optimise
  #
  if nproc==1:
    xyzs = []
    for i, (ligand_model, buffer_model, qmm, qmr) in enumerate(objects):
      xyz = qm_manager.qm_runner(
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
  else:
    print('nproc',nproc)
    xyzs = get_all_xyz(objects, nproc)
  print('',file=log)
  #
  # update model restraints
  #
  for i, ((ligand_model, buffer_model, qmm, qmr), xyz) in enumerate(zip(objects, xyzs)):
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
    ligand_model.get_hierarchy().atoms().set_xyz(xyz)
    #need to update buffer_model also ???
    gs = ligand_model.geometry_statistics()
    print('  Interim stats : %s' % gs.show_bond_and_angle(), file=log)
    if qmr.write_final_pdb_core:
      write_pdb_file(ligand_model, 'ligand_final_%s.pdb' % preamble, log)
    if qmr.write_final_pdb_buffer:
      write_pdb_file(buffer_model, 'cluster_final_%s.pdb' % preamble, log)
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
    #
    # final stats
    #
    gs = ligand_model.geometry_statistics()
    print('  Finished stats : %s' % gs.show_bond_and_angle(), file=log)
  print('='*80)

if __name__ == '__main__':
  print(quantum_chemistry_scope)
