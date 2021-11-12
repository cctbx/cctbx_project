from __future__ import absolute_import,division, print_function
from io import StringIO

from libtbx import Auto
from libtbx.str_utils import make_header

from cctbx.geometry_restraints.manager import manager as standard_manager
from mmtbx.geometry_restraints import quantum_interface

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
  # Create selection lists
  #
  from mmtbx.geometry_restraints import qm_manager
  # qm_managers = []
  make_header('QM Restraints Initialisation', out=log)
  for i, qmr in enumerate(params.qi.qm_restraints):
    print('  Selection %2d: %s' % (i+1, qmr.selection), file=log)
    ligand_model, buffer_model = get_ligand_buffer_models(model, qmr)
    qmm = get_qm_manager(ligand_model, buffer_model, qmr)
    # print(qmm)
    print(show_ligand_buffer_models(ligand_model, buffer_model), file=log)
    # qm_managers.append(qmm)
  #
  # Add to QI GRM
  #
  # qm_grm.qm_managers = qm_managers
  return qm_grm

def select_and_reindex(model, selection_str=None, selection_array=None):
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
  assert (selection_array, selection_str).count(None)==1
  if selection_str:
    selection_array = model.selection(selection_str)
  selected_model = model.select(selection_array)
  rc = _reindexing(selected_model, selection_array)
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
    outl += '      %s\n' % (atom.id_str().replace('pdb=',''))
  for atom in buffer_model.get_atoms():
    agi = atom.parent().id_str()
    if agi not in ags: ags.append(agi)
  outl += '    Buffer residues\n'
  for agi in ags:
    outl += '      %s\n' % agi
  return outl

def get_qm_manager(ligand_model, buffer_model, qmr):
  from mmtbx.geometry_restraints import qm_manager
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
      qmr.package.charge=-2
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

def get_all_xyz(inputs, nproc):
  assert 0
  from libtbx import easy_mp
  from mmtbx.geometry_restraints import qm_manager
  cmds = []
  for i, (ligand_model, buffer_model, qmm) in enumerate(inputs):
    cmd = qmm.get_cmd()
    cmds.append([cmd])
  for args, res, err_str in easy_mp.multi_core_run(qm_manager.run_orca_cmd,
                                                   tuple(cmds),
                                                   nproc,
                                                   ):
    if res:
      print (args, res)
      # results.update(res)
      # results[args[0]]=res
  xyzs = []
  return xyzs

def update_restraints(model,
                      params, # just the qi scope
                      macro_cycle=None,
                      nproc=1, # not used
                      log=StringIO(),
                      ):
  objects = []
  for i, qmr in enumerate(params.qm_restraints):
    preamble = quantum_interface.get_preamble(macro_cycle, i, qmr.selection)
    #
    # get ligand and buffer region models
    #
    ligand_model, buffer_model = get_ligand_buffer_models(model, qmr)
    if qmr.write_pdb_core:
      outl = ligand_model.model_as_pdb()
      f=open('ligand_%s.pdb' % preamble, 'w')
      f.write(outl)
      del f
    if qmr.write_pdb_buffer:
      outl = buffer_model.model_as_pdb()
      f=open('cluster_%s.pdb' % preamble, 'w')
      f.write(outl)
      del f
    #
    # get appropriate QM manager
    #
    qmm = get_qm_manager(ligand_model, buffer_model, qmr)
    qmm.preamble=preamble
    objects.append([ligand_model, buffer_model, qmm])
  #
  # optimise
  #
  # xyzs = get_all_xyz(objects, nproc)
  xyzs = []
  for i, (ligand_model, buffer_model, qmm) in enumerate(objects):
    xyz = qmm.get_opt(cleanup=qmr.cleanup)
    # qmm.print_timings(log=log)
    xyzs.append(xyz)
  #
  # update model restraints
  #
  for i, ((ligand_model, buffer_model, qmm), xyz) in enumerate(zip(objects, xyzs)):
    if i: print(' ',file=log)
    print('  Updating QM restraints: %s' % qmr.selection, file=log)
    print(show_ligand_buffer_models(ligand_model, buffer_model), file=log)
    gs = ligand_model.geometry_statistics()
    print('  Starting stats: %s' % gs.show_bond_and_angle(), file=log)
    #
    # update coordinates of ligand
    #
    ligand_model.get_hierarchy().atoms().set_xyz(xyz)
    gs = ligand_model.geometry_statistics()
    print('  Interim stats : %s' % gs.show_bond_and_angle(), file=log)
    #
    # remove H/D before transfer
    #
    if ligand_model.has_hd():
      hd_selection = ligand_model.get_hd_selection()
      ligand_model = select_and_reindex(ligand_model,
                                        selection_array=~hd_selection)
    #
    # transfer geometry to proxies
    #  - bonds
    #
    model_grm = model.get_restraints_manager()
    ligand_grm = ligand_model.get_restraints_manager()
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
      angle = ligand_lookup.get(angle_proxy.i_seqs)
      angle_proxy.angle_ideal=angle
    for angle_proxy in model_grm.geometry.angle_proxies:
      angle = model_lookup.get(angle_proxy.i_seqs)
      if angle is None: continue
      angle_proxy.angle_ideal=angle
    print('', file=log)
    #
    # final stats
    #
    gs = ligand_model.geometry_statistics()
    print('  Finished stats : %s' % gs.show_bond_and_angle(), file=log)

if __name__ == '__main__':
  print(quantum_chemistry_scope)
