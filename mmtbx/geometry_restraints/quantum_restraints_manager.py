from __future__ import absolute_import,division, print_function
from io import StringIO
import copy

from libtbx import Auto
from libtbx.str_utils import make_header

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
  return qm_grm

def add_hydrogen_atoms_to_model(model,
                                use_capping_hydrogens=False,
                               # use_neutron_distances=True,
                               # n_terminal_charge=True,
                               ):
  from mmtbx.ligands.ready_set_utils import add_terminal_hydrogens
  from mmtbx.ligands.ready_set_utils import add_water_hydrogen_atoms_simple
  model.log=StringIO()
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

def get_ligand_buffer_models(model, qmr, verbose=False):
  ligand_model = select_and_reindex(model, qmr.selection)
  buffer_selection_string = 'residues_within(%s, %s)' % (qmr.buffer,
                                                         qmr.selection)
  buffer_model = select_and_reindex(model, buffer_selection_string)
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

  total_charge = quantum_interface.electrons(buffer_model)
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
  # for atom in buffer_model.get_atoms():
    # print(atom.quote(),atom.i_seq,atom.tmp)
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
  for i, (ligand_model, buffer_model, qmm) in enumerate(inputs):
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
    #
    # need to add H atoms including on water
    #
    # master_buffer_model = copy.deepcopy(buffer_model)
    # if not master_buffer_model.restraints_manager_available():
    #   master_buffer_model.process(make_restraints=True)
    # total_charge = quantum_interface.electrons(buffer_model)
    # print('total_charge',total_charge)
    add_hydrogen_atoms_to_model(buffer_model)
    buffer_model.unset_restraints_manager()
    buffer_model.process(make_restraints=True)
    # buffer_model.add_hydrogens(1., occupancy=1.)
    # buffer_model.unset_restraints_manager()
    # buffer_model.process(make_restraints=True)
    assert buffer_model.restraints_manager_available()
    # total_charge = quantum_interface.electrons(buffer_model)
    # print('total_charge',total_charge)
    if qmr.write_pdb_core:
      outl = ligand_model.model_as_pdb()
      lf = 'ligand_%s.pdb' % preamble
      print('    Writing ligand : %s' % lf, file=log)
      f=open(lf, 'w')
      f.write(outl)
      del f
    if qmr.write_pdb_buffer:
      outl = buffer_model.model_as_pdb()
      lf = 'cluster_%s.pdb' % preamble
      print('    Writing cluster: %s' % lf, file=log)
      f=open(lf, 'w')
      f.write(outl)
      del f
    #
    # get appropriate QM manager
    #
    qmm = get_qm_manager(ligand_model, buffer_model, qmr)
    qmm.preamble=preamble
    objects.append([ligand_model, buffer_model, qmm])
  print('',file=log)
  #
  # optimise
  #
  if nproc==1:
    xyzs = []
    for i, (ligand_model, buffer_model, qmm) in enumerate(objects):
      xyz = qm_manager.qm_runner(
        qmm,
        cleanup=qmr.cleanup,
        file_read=qmr.package.read_output_to_skip_opt_if_available,
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
  for i, ((ligand_model, buffer_model, qmm), xyz) in enumerate(zip(objects, xyzs)):
    if i: print(' ',file=log)
    print('  Updating QM restraints: "%s"' % qmr.selection, file=log)
    print(show_ligand_buffer_models(ligand_model, buffer_model), file=log)
    gs = ligand_model.geometry_statistics()
    print('  Starting stats: %s' % gs.show_bond_and_angle(), file=log)
    #
    # update coordinates of ligand
    #
    ligand_model.get_hierarchy().atoms().set_xyz(xyz)
    #need to update buffer_model also
    gs = ligand_model.geometry_statistics()
    print('  Interim stats : %s' % gs.show_bond_and_angle(), file=log)
    if qmr.write_final_pdb_core:
      outl = ligand_model.model_as_pdb()
      lf = 'ligand_%s.pdb' % preamble
      print('    Writing ligand : %s' % lf)
      f=open(lf, 'w')
      f.write(outl)
      del f
    if qmr.write_final_pdb_buffer:
      assert 0
      outl = buffer_model.model_as_pdb()
      lf = 'cluster_%s.pdb' % preamble
      print('    Writing cluster: %s' % lf)
      f=open(lf, 'w')
      f.write(outl)
      del f
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
