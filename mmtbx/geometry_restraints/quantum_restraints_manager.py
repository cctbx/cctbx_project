import os
from io import StringIO

from cctbx.geometry_restraints.manager import manager as standard_manager

class manager(standard_manager):
  def __init__(self,
               params,
               log=StringIO()
               ):
    pass

  def __repr__(self):
    return 'QM restraints manager'

  def energies_sites1(self,
                     sites_cart,
                     flags=None,
                     custom_nonbonded_function=None,
                     compute_gradients=False,
                     gradients=None,
                     disable_asu_cache=False,
                     normalization=False,
                     external_energy_function=None,
                     extension_objects=[],
                     site_labels=None,
                     log=None):
    assert 0
    result = standard_manager.energies_sites(
      self,
      sites_cart,
      flags=flags,
      custom_nonbonded_function=custom_nonbonded_function,
      compute_gradients=False, #compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization,
      external_energy_function=external_energy_function,
      extension_objects=extension_objects,
      site_labels=site_labels,
      )
    assert 0

  def cleanup(self):
    make_header('Cleaning up - QM restraints')
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
  # Create selection lists
  #
  from mmtbx.geometry_restraints import qm_manager
  qm_managers = []
  for i, qmr in enumerate(params.qi.qm_restraints):
    ligand_selection = model.selection(qmr.selection)
    ligand_model = model.select(ligand_selection)

    # qmm = qm_manager.orca_manager(
    #   ligand_model.get_atoms(),
    #   qmr.package.method,
    #   qmr.package.basis_set,
    #   '', #qmr.package.solvent_model,
    #   qmr.package.charge,
    #   qmr.package.multiplicity,
    #   preamble='%02d' % (i+1),
    #   )
    # qm_managers.append(qmm)
  #
  # Add to QI GRM
  #
  # qm_grm.qm_managers = qm_managers
  return qm_grm

def get_ligand_buffer_models(model, qmr):
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
        print(atom.quote(), atom.tmp)
    return mod

  ligand_selection = model.selection(qmr.selection)
  ligand_model = model.select(ligand_selection)
  rc = _reindexing(ligand_model, ligand_selection)

  buffer_selection_string = 'residues_within(%s, %s)' % (qmr.buffer,
                                                         qmr.selection)
  buffer_selection = model.selection(buffer_selection_string)
  buffer_model = model.select(buffer_selection)
  rc = _reindexing(buffer_model, buffer_selection)
  return ligand_model, buffer_model

def show_ligand_buffer_models(ligand_model, buffer_model):
  outl = '    Core atoms\n'
  for atom in ligand_model.get_atoms():
    outl += '      %s\n' % (atom.quote().replace('"', ''))
  return outl

def get_qm_manager(ligand_model, buffer_model, qmr):
  from mmtbx.geometry_restraints import qm_manager
  program = qmr.package.program
  if program=='test':
    qmm = qm_manager.base_manager
  elif program=='orca':
    qmm = qm_manager.base_manager
  else:
    assert 0
  qmm = qmm(ligand_model.get_atoms(),
            qmr.package.method,
            qmr.package.basis_set,
            '', #qmr.package.solvent_model,
            qmr.package.charge,
            qmr.package.multiplicity,
            # preamble='%02d' % (i+1),
            )
  qmm.program=program
  return qmm

def update_restraints(model,
                      params,
                      log=StringIO(),
                      ):
  # atoms = model.get_atoms()
  for i, qmr in enumerate(params.qm_restraints):
    print('  Updating QM restraints: %s' % qmr.selection, file=log)
    #
    # get ligand and buffer region models
    #
    ligand_model, buffer_model = get_ligand_buffer_models(model, qmr)
    print(show_ligand_buffer_models(ligand_model, buffer_model))
    gs = ligand_model.geometry_statistics()
    print('  Starting stats: %s' % gs.show_bond_and_angle(), file=log)
    #
    # get appropriate QM manager
    #
    qmm = get_qm_manager(ligand_model, buffer_model, qmr)
    #
    # optimise
    #
    xyz = qmm.get_opt()
    #
    # update coordinates of ligand
    #
    ligand_model.get_hierarchy().atoms().set_xyz(xyz)
    gs = ligand_model.geometry_statistics()
    print('  Interim stats : %s' % gs.show_bond_and_angle(), file=log)
    #
    # transfer geometry to proxies
    #
    for atom in ligand_model.get_atoms():
      print(atom.quote(),atom.tmp)
    model_grm = model.get_restraints_manager()
    ligand_grm = ligand_model.get_restraints_manager()
    print(dir(ligand_grm))
    print(dir(ligand_grm.geometry))
    bond_proxies_simple, asu = ligand_grm.geometry.get_all_bond_proxies(
      sites_cart=ligand_model.get_sites_cart())
    sorted_table, n_not_shown = bond_proxies_simple.get_sorted(
      'delta',
      ligand_model.get_sites_cart())
    for info in sorted_table:
      print(info)
      (i_seq, j_seq, i_seqs, distance_ideal, distance_model, slack, delta, sigma, weight, residual, sym_op_j, rt_mx) = info
      print(i_seq, j_seq, delta)
      atoms = ligand_model.get_atoms()
      i_atom=atoms[i_seq]
      j_atom=atoms[j_seq]
      i_seqs=[i_seq, j_seq]
      print(i_atom.quote(), j_atom.quote())
      bond=ligand_grm.geometry.bond_params_table.lookup(*list(i_seqs))
      bond.distance_ideal=distance_model
      i_seqs=[i_atom.tmp, j_atom.tmp]
      if 0 :
        print(i_seqs)
        atoms=model.get_atoms()
        print(atoms[i_seqs[0]].quote())
        print(atoms[i_seqs[1]].quote())
      bond=model_grm.geometry.bond_params_table.lookup(*list(i_seqs))
      bond.distance_ideal=distance_model
    #
    # final stats
    #
    gs = ligand_model.geometry_statistics()
    print('  Finished stats : %s' % gs.show_bond_and_angle(), file=log)

if __name__ == '__main__':
  print(quantum_chemistry_scope)


