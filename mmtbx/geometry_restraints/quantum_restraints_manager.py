from __future__ import absolute_import,division, print_function
from io import StringIO
import time

from libtbx import Auto
from libtbx.str_utils import make_header
from libtbx.utils import Sorry
from libtbx.utils import null_out

from cctbx.geometry_restraints.manager import manager as standard_manager
from mmtbx.geometry_restraints import quantum_interface
from mmtbx.geometry_restraints import qm_manager
from mmtbx.geometry_restraints import mopac_manager
from mmtbx.geometry_restraints import orca_manager

from mmtbx.model.restraints import get_restraints_from_model_via_grm

WRITE_STEPS_GLOBAL=False

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
    if qmr.run_in_macro_cycles=='test': break # REMOVE AFTER TESTING!!!
  else: run_program=False
  run_energies(model,
               params,
               run_program=run_program,
               # transfer_internal_coordinates=transfer_internal_coordinates,
               log=log,
              )
  return qm_grm

def get_prefix(params):
  if hasattr(params, 'output') and hasattr(params.output, 'prefix'):
    prefix = params.output.prefix
  else:
    prefix = 'qmr'
    # print('changing prefix to %s' % prefix)
  return prefix

def write_pdb_file(model, filename, log=None):
  outl = model.model_as_pdb()
  print('    Writing PDB : %s' % filename, file=log)
  f=open(filename, 'w')
  f.write(outl)
  del f

def write_restraints(model, filename, header=None, log=None):
  """Write restraints from a model

  Args:
      model (model class): Model with one entity and no alt. loc.
      filename (TYPE): Description
      log (None, optional): Description

  """
  co = get_restraints_from_model_via_grm(model, ideal=False)
  print('\n  Writing restraints : %s' % filename, file=log)
  f=open(filename, 'w')
  if header:
    for line in header.splitlines():
      if not line.startswith('#'):
        line = '# %s' % line
      f.write('%s\n' % line)
  f.write(str(co))
  del f

def add_additional_hydrogen_atoms_to_model( model,
                                            use_capping_hydrogens=False,
                                            # use_neutron_distances=True,
                                            retain_original_hydrogens=True,
                                            append_to_end_of_model=False,
                                            n_terminal_charge=True,
                                            ):
  from mmtbx.ligands.ready_set_utils import add_terminal_hydrogens
  from mmtbx.ligands.ready_set_utils import add_water_hydrogen_atoms_simple
  from mmtbx.ligands.ready_set_utils import delete_charged_n_terminal_hydrogens
  rc = add_terminal_hydrogens( model.get_hierarchy(),
                               model.get_restraints_manager().geometry,
                               use_capping_hydrogens=use_capping_hydrogens,
                               retain_original_hydrogens=retain_original_hydrogens,
                               append_to_end_of_model=append_to_end_of_model,
                               )
  assert not rc
  rc = add_water_hydrogen_atoms_simple(model.get_hierarchy())
  if not n_terminal_charge:
    rc = delete_charged_n_terminal_hydrogens(model.get_hierarchy())

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
  assert len(buffer_model.get_atoms())>0

def use_neutron_distances_in_model_in_place(model):
  """Changes the X-H distances in place

  Args:
      model (model): model
  """
  params = model.get_current_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  model.log=null_out()
  model.process(make_restraints           = True,
                pdb_interpretation_params = params)
  model.set_hydrogen_bond_length(show=False)

def get_ligand_buffer_models(model, qmr, verbose=False, write_steps=False):
  from cctbx.maptbx.box import shift_and_box_model
  if WRITE_STEPS_GLOBAL: write_steps=True
  ligand_model = select_and_reindex(model, qmr.selection)
  #
  # check for to sparse selections like a ligand in two monomers
  #
  if len(ligand_model.get_atoms())==0:
    raise Sorry('selection "%s" results in empty model' % qmr.selection)
  if write_steps: write_pdb_file(ligand_model, 'model_selection.pdb', None)
  assert qmr.selection.find('within')==-1
  buffer_selection_string = 'residues_within(%s, %s)' % (qmr.buffer,
                                                         qmr.selection)
  buffer_model = select_and_reindex(model, buffer_selection_string)
  if write_steps: write_pdb_file(buffer_model, 'pre_remove_altloc.pdb', None)
  buffer_model.remove_alternative_conformations(always_keep_one_conformer=True)
  if write_steps: write_pdb_file(buffer_model, 'post_remove_altloc.pdb', None)
  for residue_group in buffer_model.get_hierarchy().residue_groups():
    if (residue_group.atom_groups_size() != 1):
      raise Sorry("Not implemented: cannot run QI on buffer "+
                  "molecules with alternate conformations")
  if write_steps: write_pdb_file(buffer_model, 'pre_super_cell.pdb', None)
  super_cell_and_prune(buffer_model, ligand_model, qmr.buffer, write_steps=write_steps)
  buffer_model_p1 = shift_and_box_model(model=buffer_model)
  for atom1, atom2 in zip(buffer_model_p1.get_atoms(), buffer_model.get_atoms()):
    atom1.tmp=atom2.tmp
  for atom1 in buffer_model.get_atoms():
    for atom2 in ligand_model.get_atoms():
      if atom1.id_str()==atom2.id_str():
        atom2.xyz=atom1.xyz
        break
  buffer_model = buffer_model_p1
  buffer_model.unset_restraints_manager()
  buffer_model.log=null_out()
  buffer_model.process(make_restraints=True)
  if write_steps: write_pdb_file(buffer_model, 'pre_add_terminii.pdb', None)
  add_additional_hydrogen_atoms_to_model(buffer_model,
                                         use_capping_hydrogens=qmr.capping_groups)
  buffer_model.unset_restraints_manager()
  buffer_model.log=null_out()
  buffer_model.process(make_restraints=True)
  if write_steps: write_pdb_file(buffer_model, 'post_add_terminii.pdb', None)
  ligand_atoms = ligand_model.get_atoms()
  buffer_atoms = buffer_model.get_atoms()
  for atom1 in ligand_atoms:
    for atom2 in buffer_atoms:
      if atom1.id_str()==atom2.id_str():
        break
    else:
      raise Sorry('''Bug alert
  Atom %s from ligand does not appear in buffer. Contact Phenix with input files.
  ''' % atom1.quote())
  use_neutron_distances_in_model_in_place(ligand_model)
  use_neutron_distances_in_model_in_place(buffer_model)
  if write_steps:
    write_pdb_file(buffer_model, 'post_neutron_cluster.pdb', None)
    write_pdb_file(ligand_model, 'post_neutron_ligand.pdb', None)
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

def get_qm_manager(ligand_model, buffer_model, qmr, program_goal, log=StringIO()):
  def default_defaults(qmr):
    if qmr.package.basis_set is Auto:
      qmr.package.basis_set=''
    if qmr.package.solvent_model is Auto:
      qmr.package.solvent_model=''
    if qmr.package.multiplicity is Auto:
      qmr.package.multiplicity=1
    if qmr.package.charge is Auto:
      qmr.package.charge=0
  program = qmr.package.program
  if program=='test':
    qmm = qm_manager.base_qm_manager.base_qm_manager
  elif program=='orca':
    qmm = orca_manager.orca_manager
    default_defaults(qmr)
    if qmr.package.method is Auto:
      qmr.package.method='AM1'
  elif program=='mopac':
    qmm = mopac_manager.mopac_manager
    default_defaults(qmr)
    if qmr.package.method is Auto:
      qmr.package.method='PM7'
  else:
    assert 0

  electron_model = None
  if program_goal in ['energy', 'strain']:
    electron_model = ligand_model
  elif program_goal in ['opt']:
    electron_model = buffer_model
  total_charge = quantum_interface.electrons(electron_model, log=log)
  if total_charge!=qmr.package.charge:
    print(u'  Update charge %s ~> %s' % (qmr.package.charge, total_charge),
          file=log)
    qmr.package.charge=total_charge
  qmm = qmm(electron_model.get_atoms(),
            qmr.package.method,
            qmr.package.basis_set,
            qmr.package.solvent_model,
            qmr.package.charge,
            qmr.package.multiplicity,
            # preamble='%02d' % (i+1),
            )
  qmm.program=program
  qmm.program_goal=program_goal
  ligand_i_seqs = []
  for atom in ligand_model.get_atoms():
    ligand_i_seqs.append(atom.tmp)
  hd_selection = electron_model.get_hd_selection()
  ligand_selection = []
  for i, (sel, atom) in enumerate(zip(hd_selection, electron_model.get_atoms())):
    if atom.tmp in ligand_i_seqs:
      ligand_selection.append(True)
    else:
      ligand_selection.append(False)
  if qmr.include_nearest_neighbours_in_optimisation:
    for i, (sel, atom) in enumerate(zip(ligand_selection, electron_model.get_atoms())):
      if atom.name.strip() in ['CA', 'C', 'N', 'O', 'OXT']: continue
      ligand_selection[i]=True
  qmm.set_ligand_atoms(ligand_selection)
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
  quantum_interface.get_qi_macro_cycle_array(params, verbose=True, log=log)

def running_this_macro_cycle(qmr,
                             macro_cycle,
                             number_of_macro_cycles=-1,
                             energy_only=False,
                             pre_refinement=True,
                             ):
  if qmr.run_in_macro_cycles=='all':
    return True
  elif qmr.run_in_macro_cycles=='first_only' and macro_cycle==1:
    return True
  elif qmr.run_in_macro_cycles=='test' and macro_cycle in [0, None]:
    return True
  if energy_only:
    # if macro_cycle in [0, None]: return False
    if pre_refinement:
      if qmr.calculate_starting_energy or qmr.calculate_starting_strain:
        return True
    else:
      if qmr.calculate_final_energy or qmr.calculate_final_strain:
        if macro_cycle==number_of_macro_cycles or macro_cycle==-1:
          return True
  return False

def update_bond_restraints(ligand_model,
                           buffer_model,
                           model_grm=None, # flag for update vs checking
                           ignore_x_h_distance_protein=False,
                           include_inter_residue_restraints=False,
                           log=StringIO()):
  ligand_grm = ligand_model.get_restraints_manager()
  ligand_atoms = ligand_model.get_atoms()
  ligand_i_seqs = []
  for atom in ligand_atoms:
    ligand_i_seqs.append(atom.tmp)
  buffer_grm = buffer_model.get_restraints_manager()
  atoms = buffer_model.get_atoms()
  bond_proxies_simple, asu = buffer_grm.geometry.get_all_bond_proxies(
    sites_cart=buffer_model.get_sites_cart())
  sorted_table, n_not_shown = bond_proxies_simple.get_sorted(
    'delta',
    buffer_model.get_sites_cart())
  i=0
  for info in sorted_table:
    #
    # loop over buffer restraints
    #
    (i_seq, j_seq, i_seqs, distance_ideal, distance_model, slack, delta, sigma, weight, residual, sym_op_j, rt_mx) = info
    i_atom=atoms[i_seq]
    j_atom=atoms[j_seq]
    if model_grm:
      if i_atom.element_is_hydrogen(): continue
      if j_atom.element_is_hydrogen(): continue
    i_seqs=[i_seq, j_seq]
    bond=buffer_grm.geometry.bond_params_table.lookup(*list(i_seqs))
    i+=1
    if model_grm:
      if include_inter_residue_restraints:
        if not(i_atom.tmp in ligand_i_seqs or j_atom.tmp in ligand_i_seqs):
          continue
      else:
        if not(i_atom.tmp in ligand_i_seqs and j_atom.tmp in ligand_i_seqs):
          continue
      print('    %-2d %s - %s %5.3f ~> %5.3f' % (
        i,
        i_atom.id_str().replace('pdb=',''),
        j_atom.id_str().replace('pdb=',''),
        bond.distance_ideal,
        distance_model), file=log)
      bond.distance_ideal=distance_model
      i_seqs=[i_atom.tmp, j_atom.tmp]
      bond=model_grm.geometry.bond_params_table.lookup(*list(i_seqs))
      bond.distance_ideal=distance_model
    else:
      if ( i_atom.element_is_hydrogen() or j_atom.element_is_hydrogen()):
        if distance_model>1.5:
          print('    %-2d %s - %s %5.3f ~> %5.3f' % (
            i,
            i_atom.id_str().replace('pdb=',''),
            j_atom.id_str().replace('pdb=',''),
            bond.distance_ideal,
            distance_model), file=log)
          if not ignore_x_h_distance_protein:
            raise Sorry('''
  The QM optimisation has caused a X-H covalent bond to exceed 1.5 Angstrom.
  This may be because the input geometry was not appropriate or the QM method
  is not providing the biological answer. Check the QM result for issues. This
  check can be skipped by using ignore_x_h_distance_protein.
  ''')

def update_bond_restraints_simple(model):
  """Update bond restraints in model to match the actual model values

  Args:
      model (model): model!
  """
  grm = model.get_restraints_manager()
  atoms = model.get_atoms()
  bond_proxies_simple, asu = grm.geometry.get_all_bond_proxies(
    sites_cart=model.get_sites_cart())
  sorted_table, n_not_shown = bond_proxies_simple.get_sorted(
    'delta',
    model.get_sites_cart())
  # for single ions there is no bond table
  if sorted_table is None: return
  for info in sorted_table:
    (i_seq, j_seq, i_seqs, distance_ideal, distance_model, slack, delta, sigma, weight, residual, sym_op_j, rt_mx) = info
    i_atom=atoms[i_seq]
    j_atom=atoms[j_seq]
    # exclude X-H because x-ray...
    if i_atom.element_is_hydrogen(): continue
    if j_atom.element_is_hydrogen(): continue
    i_seqs=[i_seq, j_seq]
    bond=grm.geometry.bond_params_table.lookup(*list(i_seqs))
    bond.distance_ideal=distance_model

def get_program_goal(qmr, macro_cycle=None, energy_only=False):
  program_goal='opt' # can be 'opt', 'energy', 'strain'
  if not energy_only: return program_goal
  if macro_cycle==1:
    if qmr.calculate_starting_energy:
      program_goal='energy'
    if qmr.calculate_starting_strain:
      program_goal='strain'
  else: # only called with final energy on final macro cycle
    if qmr.calculate_final_energy:
      program_goal='energy'
    if qmr.calculate_final_strain:
      program_goal='strain'
  return program_goal

def setup_qm_jobs(model,
                  params,
                  macro_cycle,
                  energy_only=False,
                  pre_refinement=True,
                  log=StringIO()):
  prefix = get_prefix(params)
  objects = []
  for i, qmr in enumerate(params.qi.qm_restraints):
    number_of_macro_cycles = 1
    if hasattr(params, 'main'):
      number_of_macro_cycles = params.main.number_of_macro_cycles
    if macro_cycle is not None and not running_this_macro_cycle(
        qmr,
        macro_cycle,
        energy_only=energy_only,
        number_of_macro_cycles=number_of_macro_cycles,
        pre_refinement=pre_refinement):
      print('    Skipping this selection in this macro_cycle : %s' % qmr.selection)
      continue
    preamble = quantum_interface.get_preamble(macro_cycle, i, qmr)
    #
    # get ligand and buffer region models
    #
    ligand_model, buffer_model = get_ligand_buffer_models(model, qmr)
    assert buffer_model.restraints_manager_available()
    if not energy_only: # only write PDB files for restraints update
      if qmr.write_pdb_core:
        write_pdb_file(ligand_model, '%s_ligand_%s.pdb' % (prefix, preamble), log)
      if qmr.write_pdb_buffer:
        write_pdb_file(buffer_model, '%s_cluster_%s.pdb' % (prefix, preamble), log)
    #
    # get appropriate QM manager
    #
    program_goal = get_program_goal(qmr, macro_cycle, energy_only=energy_only)
    if not energy_only and qmr.do_not_even_calculate_qm_restraints:
      print('    Skipping QM calculation : %s' % qmr.selection)
      continue
    qmm = get_qm_manager(ligand_model, buffer_model, qmr, program_goal, log=log)
    qmm.preamble='%s_%s' % (prefix, preamble)
    objects.append([ligand_model, buffer_model, qmm, qmr])
  print('',file=log)
  return objects

def run_jobs(objects, macro_cycle, nproc=1, log=StringIO()):
  if nproc==1:
    xyzs = []
    xyzs_buffer = []
    for i, (ligand_model, buffer_model, qmm, qmr) in enumerate(objects):
      #
      # run QM program
      #
      xyz, xyz_buffer = qm_manager.qm_runner(
        qmm,
        cleanup=qmr.cleanup,
        file_read=qmr.package.read_output_to_skip_opt_if_available,
        check_file_read_safe=not(qmr.package.ignore_input_differences),
        log=log,
        )
      print('  Time for calculation of "%s" using %s %s %s: %s' % (
        qmr.selection,
        qmr.package.method,
        qmr.package.basis_set,
        qmr.package.solvent_model,
        qmm.get_timings().split(':')[-1],
        ), file=log)
      xyzs.append(xyz)
      xyzs_buffer.append(xyz_buffer)
      # if qmm.program_goal in ['energy', 'strain']:
        # print('  Energy = %f' % xyz, file=log)
  else:
    print('nproc',nproc)
    assert 0
    xyzs, junk = get_all_xyz(objects, nproc)
  print('',file=log)
  return xyzs, xyzs_buffer

def run_energies(model,
                 params,
                 macro_cycle=None,
                 run_program=True,
                 pre_refinement=True,
                 nproc=1,
                 log=StringIO(),
                 ):
  t0 = time.time()
  energy_only=True
  if not model.restraints_manager_available():
    model.log=null_out()
    model.process(make_restraints=True)
  if macro_cycle in [None, 0]: run_program=False
  qi_array = quantum_interface.get_qi_macro_cycle_array(params)
  if quantum_interface.is_quantum_interface_active_this_macro_cycle(params,
                                                                    macro_cycle,
                                                                    energy_only=energy_only,
                                                                    # pre_refinement=pre_refinement,
                                                                    ) and run_program:
    print('  QM energy calculations for macro cycle %s' % macro_cycle,
      file=log)
  #
  # setup QM jobs
  #
  objects = setup_qm_jobs(model,
                          params,
                          macro_cycle,
                          energy_only=energy_only,
                          pre_refinement=pre_refinement,
                          log=log)
  if not run_program: return
  assert macro_cycle is not None
  #
  # run jobs
  #
  xyzs, xyzs_buffer = run_jobs(objects, macro_cycle=macro_cycle, nproc=nproc, log=log)
  print('  Total time for QM energies: %0.1fs' % (time.time()-t0), file=log)
  print('%s%s' % ('<'*40, '>'*40), file=log)

def update_restraints(model,
                      params,
                      macro_cycle=None,
                      # run_program=True,
                      # transfer_internal_coordinates=True,
                      never_write_restraints=False,
                      nproc=1,
                      log=StringIO(),
                      ):
  t0 = time.time()
  energy_only=False
  if not model.restraints_manager_available():
    model.log=null_out()
    model.process(make_restraints=True)
  if quantum_interface.is_quantum_interface_active_this_macro_cycle(params,
                                                                    macro_cycle,
                                                                    energy_only=energy_only,
                                                                    ):
    print('  QM restraints calculations for macro cycle %s' % macro_cycle,
      file=log)
  #
  # setup QM jobs
  #
  objects = setup_qm_jobs(model, params, macro_cycle, energy_only=energy_only, log=log)
  #
  # run jobs
  #
  xyzs, xyzs_buffer = run_jobs(objects, macro_cycle=macro_cycle, nproc=nproc, log=log)
  #
  # update model restraints
  #
  prefix = get_prefix(params)
  for i, ((ligand_model, buffer_model, qmm, qmr), xyz, xyz_buffer) in enumerate(
    zip(objects,
        xyzs,
        xyzs_buffer,
        )):
    if qmr.package.view_output: qmm.view(qmr.package.view_output)
    if i: print(' ',file=log)
    print('  Updating QM restraints: "%s"' % qmr.selection, file=log)
    print(show_ligand_buffer_models(ligand_model, buffer_model), file=log)
    gs = ligand_model.geometry_statistics()
    print('  Starting stats: %s' % gs.show_bond_and_angle_and_dihedral(), file=log)
    #
    # update coordinates of ligand
    #
    old = ligand_model.get_hierarchy().atoms().extract_xyz()
    if not qmr.include_nearest_neighbours_in_optimisation:
      rmsd = old.rms_difference(xyz)
      if rmsd>5:
        print('  QM minimisation has large rms difference in cartesian coordinates: %0.1f' % (rmsd),
              file=log)
        print('  Check the QM minimisation for errors or incorrect protonation.',
              file=log)
        if rmsd>20:
          print('  Movement of cartesian coordinates is very large.', file=log)
      ligand_model.get_hierarchy().atoms().set_xyz(xyz)
    old = buffer_model.get_hierarchy().atoms().extract_xyz()
    # rmsd = old.rms_difference(xyz_buffer)
    buffer_model.get_hierarchy().atoms().set_xyz(xyz_buffer)
    gs = ligand_model.geometry_statistics()
    print('  Interim stats : %s' % gs.show_bond_and_angle_and_dihedral(), file=log)
    preamble = quantum_interface.get_preamble(macro_cycle, i, qmr)
    if qmr.write_final_pdb_core:
      write_pdb_file(ligand_model, '%s_ligand_final_%s.pdb' % (prefix, preamble), log)
    if qmr.write_final_pdb_buffer:
      write_pdb_file(buffer_model, '%s_cluster_final_%s.pdb' % (prefix, preamble), log)
    if qmr.do_not_update_restraints:
      print('  Skipping updating restaints')
      continue
    #
    # transfer geometry to proxies
    #  - bonds
    #
    model_grm = model.get_restraints_manager()
    print('\n  Checking', file=log)
    update_bond_restraints(buffer_model,
                           buffer_model,
                           ignore_x_h_distance_protein=qmr.ignore_x_h_distance_protein,
                           log=log)
    print('\n  Transfer', file=log)
    update_bond_restraints(ligand_model,
                           buffer_model,
                           model_grm=model_grm,
                           include_inter_residue_restraints=qmr.include_nearest_neighbours_in_optimisation,
                           log=log)
    update_bond_restraints_simple(ligand_model)
    #
    #    - angles
    #
    ligand_grm = ligand_model.get_restraints_manager()
    atoms = ligand_model.get_atoms()
    sorted_table, n_not_shown = ligand_grm.geometry.angle_proxies.get_sorted(
      'delta',
      ligand_model.get_sites_cart())
    ligand_lookup = {}
    model_lookup = {}
    i=0
    if sorted_table is None: sorted_table=[]
    for info in sorted_table:
      (i_seqs, angle_ideal, angle_model, delta, sigma, weight, residual) = info
      i_atom=atoms[int(i_seqs[0])]
      j_atom=atoms[int(i_seqs[1])]
      k_atom=atoms[int(i_seqs[2])]
      key = (int(i_seqs[0]), int(i_seqs[1]), int(i_seqs[2]))
      ligand_lookup[key]=angle_model
      i+=1
      print('    %-2d %s - %s - %s %5.1f ~> %5.1f' % (
        i,
        atoms[key[0]].id_str().replace('pdb=',''),
        atoms[key[1]].id_str().replace('pdb=',''),
        atoms[key[2]].id_str().replace('pdb=',''),
        angle_ideal,
        angle_model), file=log)
      key = (atoms[key[0]].tmp, atoms[key[1]].tmp, atoms[key[2]].tmp)
      model_lookup[key]=angle_model
      key = (int(i_seqs[2]), int(i_seqs[1]), int(i_seqs[0]))
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
    #
    #    - torsions
    #
    sorted_table, n_not_shown = ligand_grm.geometry.dihedral_proxies.get_sorted(
      'delta',
      ligand_model.get_sites_cart())
    ligand_lookup = {}
    model_lookup = {}
    if sorted_table is None: sorted_table=[]
    i=0
    for info in sorted_table:
      (i_seqs, angle_ideal, angle_model, delta, period, sigma, weight, residual) = info
      i_atom=atoms[int(i_seqs[0])]
      j_atom=atoms[int(i_seqs[1])]
      k_atom=atoms[int(i_seqs[2])]
      l_atom=atoms[int(i_seqs[3])]
      key = (int(i_seqs[0]), int(i_seqs[1]), int(i_seqs[2]), int(i_seqs[3]))
      ligand_lookup[key]=angle_model
      i+=1
      print('    %-2d %s - %s - %s - %s %6.1f ~> %6.1f' % (
        i,
        atoms[key[0]].id_str().replace('pdb=',''),
        atoms[key[1]].id_str().replace('pdb=',''),
        atoms[key[2]].id_str().replace('pdb=',''),
        atoms[key[3]].id_str().replace('pdb=',''),
        angle_ideal,
        angle_model), file=log)
      key = (atoms[key[0]].tmp, atoms[key[1]].tmp, atoms[key[2]].tmp, atoms[key[3]].tmp)
      model_lookup[key]=angle_model
      key = (int(i_seqs[3]), int(i_seqs[2]), int(i_seqs[1]), int(i_seqs[0]))
      ligand_lookup[key]=angle_model
      key = (atoms[key[3]].tmp, atoms[key[2]].tmp, atoms[key[1]].tmp, atoms[key[0]].tmp)
      model_lookup[key]=angle_model
    for angle_proxy in ligand_grm.geometry.dihedral_proxies:
      angle = ligand_lookup.get(angle_proxy.i_seqs, None)
      if angle is None: continue
      angle_proxy.angle_ideal=angle
    for angle_proxy in model_grm.geometry.dihedral_proxies:
      angle = model_lookup.get(angle_proxy.i_seqs, None)
      if angle is None: continue
      angle_proxy.angle_ideal=angle

    print('', file=log)
    #
    # final stats
    #
    gs = ligand_model.geometry_statistics()
    print('  Finished stats : %s' % gs.show_bond_and_angle_and_dihedral(assert_zero=True),
          file=log)
    print('%s%s' % (' '*19, gs.show_planarity_details()), file=log)
    r=gs.result()
    if r.planarity.mean>0.02 or r.planarity.max>0.05:
      print('  %s\n   rmsd values for planarity restraints are high. Check QM minimisation. \n  %s' % (
            '-'*71,
            '-'*71,
            ),
            file=log)
    if qmr.write_restraints and not never_write_restraints:
      header='''
Restraints written by QMR process in phenix.refine
      ''' % ()
      if qmr.restraints_filename is not Auto:
        tmp_cif_filename = qmr.restraints_filename
      else:
        tmp_cif_filename = '%s.cif' % qmm.preamble
      if not tmp_cif_filename.endswith('.cif'):
        tmp_cif_filename = '%s.cif' % tmp_cif_filename
      write_restraints(ligand_model,
                       tmp_cif_filename,
                       header=header,
                       log=log,
                       )
  print('\n  Total time for QM restaints: %0.1fs\n' % (time.time()-t0), file=log)
  print('%s%s' % ('/'*39, '\\'*40))
  print('%s%s' % ('\\'*39, '/'*40))

if __name__ == '__main__':
  print(quantum_chemistry_scope)
