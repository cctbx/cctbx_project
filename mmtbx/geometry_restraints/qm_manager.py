from __future__ import absolute_import, division, print_function
from io import StringIO

from libtbx import adopt_init_args
from libtbx.str_utils import make_header

from libtbx.utils import Sorry
from mmtbx.geometry_restraints import base_qm_manager, mopac_manager, orca_manager

harkcal = 627.50946900
bohrang = 0.52918
#
# QM runner
#
def qm_runner(qmm,
              cleanup=True,
              file_read=False,
              check_file_read_safe=True,
              log=None,
              ):
  def get_func(manager, attr):
    return getattr(manager, 'get_%s' % attr, None)
  redirect_output=True
  # fake
  coordinate_filename_ext='.nwm'
  log_filename_ext='.nwm'
  if qmm.program=='test':
    func = get_func(base_qm_manager.base_manager, qmm.program_goal)
  elif qmm.program=='orca':
    func = get_func(orca_manager.orca_manager, qmm.program_goal)
    coordinate_filename_ext='.xyz'
    log_filename_ext='.log'
    # raise Sorry('Orca temporarily unsupported. Consider using MOPAC.')
  elif qmm.program=='mopac':
    func = get_func(mopac_manager.mopac_manager, qmm.program_goal)
    coordinate_filename_ext='.arc' # maybe better from .out?
    log_filename_ext='.out'
    redirect_output=False
  else:
    raise Sorry('QM program not found or set "%s"' % qmm.program)
  if func is None:
    raise Sorry('QM manager does not have get_%s' % qmm.program_goal)
  ligand_xyz, buffer_xyz = func(qmm,
                                cleanup=cleanup,
                                file_read=file_read,
                                check_file_read_safe=check_file_read_safe,
                                coordinate_filename_ext=coordinate_filename_ext,
                                log_filename_ext=log_filename_ext,
                                redirect_output=redirect_output,
                                log=log)
  return ligand_xyz, buffer_xyz

from cctbx.geometry_restraints.manager import manager as standard_manager

class manager(standard_manager):
  def __init__(self,
               params,
               # energy_components=None,
               # gradients_factory=flex.vec3_double,
               log=StringIO()):
    # self.gradients_factory = gradients_factory
    adopt_init_args(self, locals()) #, exclude=["log"])

    self.initialisation(params, log)

  def __repr__(self):
    return 'QM Gradients manager'

  def initialisation(self, params, log=None):
    make_header("Initializing QM Gradients", out=log)

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
    result = standard_manager.energies_sites(
      self,
      sites_cart,
      flags=flags,
      custom_nonbonded_function=custom_nonbonded_function,
      compute_gradients=compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization,
      external_energy_function=external_energy_function,
      extension_objects=extension_objects,
      site_labels=site_labels,
      )
    if result.gradients is None:
      print('ENERGY NONE')
      return result
    print(len(result.gradients))
    from mmtbx.geometry_restraints.quantum_restraints_manager import run_gradients
    self.model.set_sites_cart(sites_cart)
    rc = run_gradients(
      self.model,
      self.params,
      # macro_cycle=self.macro_cycle,
      # pre_refinement=pre_refinement,
      nproc=self.params.main.nproc,
      # log=self.log,
      )
    energies=rc.energies
    print('ENERGY',energies)
    objects = rc.objects
    gradients = rc.gradients
    factor=10
    for i, (ligand_model, buffer_model, qmm, qmr) in enumerate(objects):
      print(ligand_model)
      hierarchy=buffer_model.get_hierarchy()
      for j, atom in enumerate(hierarchy.atoms()):
        if atom.tmp==-1: continue
        # print(atom.quote(), atom.tmp, result.gradients[atom.tmp], gradients[i][j])
        tmp=[]
        for k in range(3):
          tmp.append(gradients[i][j][k]*factor)
        result.gradients[atom.tmp]=tmp
    # assert 0
    return result

def digester(model,
             standard_geometry_restraints_manager,
             # geometry,
             params,
             log=StringIO(),
             ):
  sgrm = standard_geometry_restraints_manager
  agrm = manager(params, log=log)
  for attr, value in list(vars(sgrm).items()):
    if attr.startswith('__'): continue
    setattr(agrm, attr, value)
  agrm.standard_geometry_restraints_manager = sgrm
  agrm.model=model
  return agrm

def main():
  from iotbx import pdb
  pdb_lines = '''
HETATM   97  S   SO4 A  13      31.477  38.950  15.821  0.50 25.00           S
HETATM   98  O1  SO4 A  13      31.243  38.502  17.238  0.50 25.00           O
HETATM   99  O2  SO4 A  13      30.616  40.133  15.527  0.50 25.00           O
HETATM  100  O3  SO4 A  13      31.158  37.816  14.905  0.50 25.00           O
HETATM  101  O4  SO4 A  13      32.916  39.343  15.640  0.50 25.00           O
'''
  pdb_inp = pdb.input(lines=pdb_lines, source_info='lines')
  qi_grm = orca_manager(pdb_inp.atoms(),
                        'PM3',
                        '',
                        '',
                        -2,
                        1,
                        preamble='test',
                        )
  print(qi_grm)
  energy, gradients = qi_grm.get_engrad()
  print(energy, list(gradients))
  coordinates = qi_grm.get_opt()
  print(list(coordinates))

if __name__ == '__main__':
  main()
