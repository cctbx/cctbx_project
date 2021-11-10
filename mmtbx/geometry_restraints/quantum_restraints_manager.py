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

  def energies_sites(self,
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
  ligand_selection = model.selection(qmr.selection)
  print(list(ligand_selection))
  ligand_model = model.select(ligand_selection)
  assert 0

def get_qm_manager(qmr):
  print(dir(qmr))
  assert 0

def update_restraints(model,
                      params,
                      log=StringIO(),
                      ):
  print(model)
  atoms = model.get_atoms()
  print(dir(params))
  for i, qmr in enumerate(params.qm_restraints):
    print(i,qmr.selection)
    print(qmr)
    print('-'*80)
    print(ligand_model)
    for sel in ['within', 'residues_within']:
      buffer_selection_string = '%s(%s, %s)' % (sel, qmr.buffer,
                                                    qmr.selection)
      buffer_selection = model.selection(buffer_selection_string)
      print(list(buffer_selection))
      buffer_model = model.select(buffer_selection)
      print(buffer_model)
      for i_seq, i_bool in enumerate(buffer_selection):
        if i_bool:
          print(i_seq, atoms[i_seq].quote())
      for atom in buffer_model.get_atoms():
        print(atom.quote())

if __name__ == '__main__':
  print(quantum_chemistry_scope)


