from __future__ import absolute_import,division, print_function
import os
from libtbx import Auto

def env_exists_exists(env, var, check=True):
  if check:
    orca_env = env.get(var, False)
    if not orca_env: return orca_env
    if orca_env.find('LD_LIBRARY_PATH')>-1:
      lib, orca_env = orca_env.split()
    if os.path.exists(orca_env): return True
    # return (env.get(var, False) and os.path.exists(env[var]))
  else:
    return env.get(var, False)

def is_orca_installed(env, var):
  return env_exists_exists(env, var)

def is_mopac_installed(env, var):
  return env_exists_exists(env, var)

def is_qm_test_installed(env, var):
  return env_exists_exists(env, var, check=False)

program_options = {
  'orca' : (is_orca_installed, 'PHENIX_ORCA'),
  'mopac' : (is_mopac_installed, 'PHENIX_MOPAC'),
  'test' : (is_qm_test_installed, 'PHENIX_QM_TEST'),
  }

def get_qm_restraints_scope(validate=True, verbose=False):
  qm_package_scope = '''
  package
  {
    program = %s
      .type = choice
    charge = Auto
      .type = int
    multiplicity = Auto
      .type = int
    method = Auto
      .type = str
    basis_set = Auto
      .type = str
    solvent_model = None
      .type = str
    nproc = 1
      .type = int
    read_output_to_skip_opt_if_available = False
      .type = bool
    ignore_input_differences = False
      .type = bool
    view_output = None
      .type = str
  }
'''

  qm_restraints_scope = '''
qm_restraints
  .multiple = True
{
  selection = None
    .type = atom_selection
    .help = selection for core of atoms to calculate new restraints via a QM \
            geometry minimisation
  buffer = 3.5
    .type = float
    .help = distance to include entire residues into the enviroment of the core
  calculate_starting_energy = False
    .type = bool
  calculate_final_energy = False
    .type = bool
  calculate_starting_strain = False
    .type = bool
  calculate_final_strain = False
    .type = bool
  write_pdb_core = False
    .type = bool
  write_pdb_buffer = False
    .type = bool
  write_final_pdb_core = False
    .type = bool
  write_final_pdb_buffer = False
    .type = bool
  write_restraints = True
    .type = bool
  restraints_filename = Auto
    .type = path
    .style = new_file
  cleanup = all *most None
    .type = choice
  run_in_macro_cycles = *first_only all test
    .type = choice
  ignore_x_h_distance_protein = False
    .type = bool
  capping_groups = True
    .type = bool
  include_nearest_neighbours_in_optimisation = False
    .type = bool
  do_not_update_restraints = False
    .type = bool
    .style = hidden
    .help = For testing and maybe getting strain energy of standard restraints
  do_not_even_calculate_qm_restraints = False
    .type = bool
    .style = hidden
    .help = For testing and maybe getting strain energy of standard restraints
  %s
}
'''
  programs = ''
  for package, (func, var) in program_options.items():
    if func(os.environ, var):
      if package=='mopac':
        programs += ' *%s' % package
      else:
        programs += '   %s' % package
  if verbose: print(programs)
  if validate:
    assert programs, 'Need to set some parameters for QM programs %s' % program_options
  qm_package_scope = qm_package_scope % programs
  qm_restraints_scope = qm_restraints_scope % qm_package_scope
  return qm_restraints_scope

def electrons(model, log=None):
  from elbow.quantum import electrons
  atom_valences = electrons.electron_distribution(
    model.get_hierarchy(), # needs to be altloc free
    model.get_restraints_manager().geometry,
    log=log,
    verbose=False,
  )
  atom_valences.validate(ignore_water=True,
                         raise_if_error=False)
  charged_atoms = atom_valences.get_charged_atoms()
  return atom_valences.get_total_charge()

def get_safe_filename(s):
  while s.find('  ')>-1:
    s=s.replace('  ',' ')
  s=s.replace(' ','_')
  s=s.replace("'",'_prime_')
  s=s.replace('*','_star_')
  s=s.replace('(','_lb_')
  s=s.replace(')','_rb_')
  return s

def get_preamble(macro_cycle, i, qmr, old_style=False):
  s=''
  if macro_cycle is not None:
    s+='%02d_' % macro_cycle
  if old_style:
    s+='%02d_%s_%s' % (i+1, get_safe_filename(qmr.selection), qmr.buffer)
  else:
    s+='%s_%s' % (get_safe_filename(qmr.selection), qmr.buffer)
  if qmr.capping_groups:
    s+='_C'
  if qmr.include_nearest_neighbours_in_optimisation:
    s+='_N'
  if qmr.package.method is not Auto:
    s+='_%s' % get_safe_filename(qmr.package.method)
  if qmr.package.basis_set is not Auto and qmr.package.basis_set:
    s+='_%s' % get_safe_filename(qmr.package.basis_set)
  if qmr.package.solvent_model is not Auto and qmr.package.solvent_model:
    s+='_%s' % get_safe_filename(qmr.package.solvent_model)
  return s

def is_any_quantum_package_installed(env):
  installed = []
  actions = []
  outl = ''
  for key, (question, var) in program_options.items():
    if question(os.environ, var):
      installed.append(key)
  if installed:
    # refine_buffer_hydrogen_atoms = False
    #   .type = bool
    #   .style = hidden
    outl = '''
  qi
    .help = QM
    .expert_level = 3
  {
    %s
  }
''' % get_qm_restraints_scope()
  return outl

def validate_qm_restraints(qm_restraints, verbose=False):
  """Simple check for active QM restraints

  Args:
      qm_restraints (PHIL list): List of QM restraints PHIL scopes
      verbose (bool, optional): D'oh

  Returns:
      TYPE: Description
  """
  for i, qmr in enumerate(qm_restraints):
    if verbose: print(i, qmr.selection)
    if i==0 and qmr.selection is None:
      return False
  return True

def is_quantum_interface_active(params, verbose=False):
  """Checks whether the QI is active at all

  Args:
      params (PHIL): PHIL scope with a possible 'qi' scope
      verbose (bool, optional): D'oh

  Returns:
      TYPE: False or True and the type of QI active
  """
  if not hasattr(params, 'qi'):
    if verbose: assert 0
    return False
  if verbose: print('  len(qm_restraints)=%d' % len(params.qi.qm_restraints))
  if len(params.qi.qm_restraints):
    if validate_qm_restraints(params.qi.qm_restraints, verbose=verbose):
      return True, 'qm_restraints' # includes restraints and energy
  return False

def is_qi_energy_pre_refinement(params,
                                macro_cycle,
                                ):
  assert 0
  qi = is_quantum_interface_active(params)
  if qi:
    rc = []
    if qi[1]=='qm_restraints':
      for i, qmr in enumerate(params.qi.qm_restraints):
        if macro_cycle==1:
          if qmr.calculate_starting_energy or qmr.calculate_starting_strain:
            rc.append(True)
    return True in rc
  else:
    return False

def is_qi_energy_post_refinement(params,
                                macro_cycle,
                                ):
  assert 0
  qi = is_quantum_interface_active(params)
  if qi:
    rc = []
    if qi[1]=='qm_restraints':
      for i, qmr in enumerate(params.qi.qm_restraints):
        if macro_cycle==params.main.number_of_macro_cycles:
          if qmr.calculate_final_energy or qmr.calculate_final_strain:
            rc.append(True)
    return True in rc
  else:
    return False

def is_quantum_interface_active_this_macro_cycle(params,
                                                 macro_cycle,
                                                 energy_only=False,
                                                 verbose=False):
  from mmtbx.geometry_restraints.quantum_restraints_manager import running_this_macro_cycle
  qi = is_quantum_interface_active(params)
  if qi:
    rc = []
    if qi[1]=='qm_restraints':
      number_of_macro_cycles = 1
      if hasattr(params, 'main'):
        number_of_macro_cycles = params.main.number_of_macro_cycles
      for i, qmr in enumerate(params.qi.qm_restraints):
        rc.append(running_this_macro_cycle(qmr,
                                           macro_cycle,
                                           number_of_macro_cycles,
                                           energy_only=energy_only))
    else:
      assert 0
    return rc
  else:
    return False

class unique_item_list(list):
  def append(self, item):
    if item not in self:
      list.append(self, item)

def get_qi_macro_cycle_array(params, verbose=False, log=None):
  qi = is_quantum_interface_active(params)
  number_of_macro_cycles = params.main.number_of_macro_cycles
  tmp=[]
  for i in range(number_of_macro_cycles+1):
    tmp.append(unique_item_list())
  if qi:
    for i, qmr in enumerate(params.qi.qm_restraints):
      rc=[]
      for i in range(number_of_macro_cycles+1):
        rc.append(unique_item_list())
      if qmr.calculate_starting_strain:
        rc[1].append('strain')
      elif qmr.calculate_starting_energy:
        rc[1].append('energy')
      if not qmr.do_not_even_calculate_qm_restraints:
        if qmr.run_in_macro_cycles=='first_only':
          rc[1].append('restraints')
        elif qmr.run_in_macro_cycles=='all':
          for j in range(1,number_of_macro_cycles+1):
            rc[j].append('restraints')
        elif qmr.run_in_macro_cycles=='test':
          rc[1].append('test')
      if qmr.calculate_final_strain:
        rc[-1].append('strain')
      elif qmr.calculate_final_energy:
        rc[-1].append('energy')
    if verbose:
      print('    %s' % qmr.selection, file=log)
      for j, actions in enumerate(rc):
        if actions:
          print('      %2d : %s' % (j, ' '.join(actions)), file=log)
    for j, actions in enumerate(rc):
      for action in actions:
        tmp[j].append(action)
  return tmp

def digester(model, geometry, params, log=None):
  active, choice = is_quantum_interface_active(params)
  assert active
  if not model.has_hd():
    from libtbx.utils import Sorry
    raise Sorry('Model must has Hydrogen atoms')
  if choice=='qm_restraints':
    from mmtbx.geometry_restraints import quantum_restraints_manager
    geometry = quantum_restraints_manager.digester(model,
                                                   geometry,
                                                   params,
                                                   log=log)
  else:
    assert 0
  return geometry

def main():
  print('testing QI')
  for var, item in program_options.items():
    if item[1] in os.environ: os.environ.pop(item[1])
  assert 'PHENIX_ORCA' not in os.environ
  rc = is_any_quantum_package_installed(os.environ)
  assert not rc
  for var1, item1 in program_options.items():
    os.environ[item1[1]]=os.getcwd()
    for var2, item2 in program_options.items():
      os.environ[item2[1]]=os.getcwd()
      rc = is_any_quantum_package_installed(os.environ)
      rc = get_qm_restraints_scope(verbose=True)
      # print(rc)
      os.environ.pop(item2[1])
    if item1[1] in os.environ: os.environ.pop(item1[1])

if __name__ == '__main__':
  main()
