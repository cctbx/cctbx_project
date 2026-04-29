from __future__ import absolute_import,division, print_function
import os
from libtbx import Auto

from mmtbx.geometry_restraints import mopac_manager

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

def is_mopac_installed(env, var, verbose=False):
  if mopac_manager.get_exe(verbose=verbose):
    return True
  else:
    return env_exists_exists(env, var)

def is_qm_test_installed(env, var):
  return True #env_exists_exists(env, var, check=False)

program_options = {
  'orca' : (is_orca_installed, 'PHENIX_ORCA'),
  'mopac' : (is_mopac_installed, 'PHENIX_MOPAC'),
  'test' : (is_qm_test_installed, 'PHENIX_QM_TEST'),
  }

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

qm_selection_specs = '''
  selection = None
    .type = atom_selection
    .help = selection for core of atoms to calculate new restraints via a QM \
            geometry minimisation
  run_in_macro_cycles = %sfirst_only first_and_last %sall last_only test
    .type = choice
    .help = the steps of the refinement that the restraints generation is run
  buffer = 3.5
    .type = float
    .help = distance to include entire residues into the enviroment of the core
  specific_atom_charges
    .optional = True
    .multiple = True
    .short_caption = Specify the charge for a specific atom (mostly metal ions)
    .style = auto_align
  {
    atom_selection = None
      .type = atom_selection
      .input_size = 400
    charge = None
      .type = int
  }
  specific_atom_multiplicities
    .optional = True
    .multiple = True
    .short_caption = Specify the multiplicity for a specific atom (mostly metal ions). \
    Selection is not needed but is great for bookkeeping.
    .style = auto_align
  {
    atom_selection = None
      .type = atom_selection
      .input_size = 400
    multiplicity = None
      .type = int
  }
  ignore_lack_of_h_on_ligand = False
    .type = bool
    .help = skip check on protonation of ligand for entities such as MgF3
  capping_groups = True
    .type = bool
  cleanup = all *most None
    .type = choice
'''

def get_qm_package_scope(validate=True, verbose=False):
  programs = ''
  for package, (func, var) in program_options.items():
    if func(os.environ, var):
      if package=='mopac':
        programs += ' *%s' % package
      else:
        programs += ' %s' % package
  if verbose: print(programs)
  if validate:
    assert programs, 'Need to set some parameters for QM programs %s' % program_options
  return qm_package_scope % programs

def get_qm_restraints_scope(validate=True, verbose=False):
  qm_restraints_scope = '''
qm_restraints
  .multiple = True
{
  %s
  calculate = *in_situ_opt starting_energy final_energy \
starting_strain final_strain starting_bound final_bound \
starting_binding final_binding \
gradients
    .type = choice(multi=True)
    .help = Choose QM calculations to run
    .caption = in_situ_optimisation_of_selection \
      starting_energy_of_isolated_ligand final_energy_of_isolated_ligand \
      strain_energy_of_starting_ligand_geometry \
      strain_energy_of_final_ligand_geometry \
      starting_energy_of_bound_ligand_cluster \
      starting_binding_energy_of)ligand final_binding_energy_of_ligand \
      final_energy_of_bound_ligand_cluster gradients

  write_files = *restraints pdb_core pdb_buffer pdb_final_core pdb_final_buffer
    .type = choice(multi=True)
    .help = which ligand or cluster files to write
    .caption = restraints_file \
      input_ligand_in_PDB_format input_cluster_in_PDB_format \
      final_ligand_in_PDB_format final_cluster_in_PDB_format

  protein_optimisation_freeze = *all None main_chain main_chain_to_beta main_chain_to_delta torsions
    .type = choice(multi=True)
    .help = the parts of protein residues that are frozen when an amino \
            acid is the main selection
    .caption = all None N_CA_C_O N_CA_CB_C_O N_CA_CB_CD_C_O all_torsions

  remove_water = False
    .type = bool
    .style = hidden

  restraints_filename = Auto
    .type = path
    .style = new_file
    .help = restraints filename is based on model name if not specified
  ignore_x_h_distance_protein = False
    .type = bool
    .help = skip check on transfer of proton during QM optimisation

  freeze_specific_atoms
    .optional = True
    .multiple = True
    .short_caption = specify atoms of the main selection frozen in optimisation
    .caption = can be used to freeze a ligand from moving too far. use Auto \
               to freeze the atom closest to the centre of mass.
    .style = auto_align
  {
    atom_selection = None
      .type = atom_selection
      .input_size = 400
  }

  include_nearest_neighbours_in_optimisation = False
    .type = bool
    .short_caption = include protein side chain in ligand optimisation
    .help = include the side chains of protein in the QM optimisation

  include_inter_residue_restraints = False
    .type = bool
    .short_caption = include residue (including metal) linking in restraints \
                     update

  do_not_update_restraints = False
    .type = bool
    .style = hidden
    .help = For testing and maybe getting strain energy of standard restraints
  buffer_selection = None
    .type = atom_selection
    .help = use this instead of distance from selection
    .style = hidden
  %s
}
'''
  qm_package_scope_filled = get_qm_package_scope(validate=validate ,verbose=verbose)
  qm_restraints_scope = qm_restraints_scope % (qm_selection_specs % ('*',''),
                                               qm_package_scope_filled)
  return qm_restraints_scope

def get_qm_gradients_scope(validate=True, verbose=False):
  qm_gradients_scope = '''
qm_gradients
  .multiple = True
{
  %s
  %s
}
'''
  qm_package_scope_filled = get_qm_package_scope(validate=validate, verbose=verbose)
  qm_gradients_scope = qm_gradients_scope % (qm_selection_specs % ('','*'),
                                               qm_package_scope_filled)
  return qm_gradients_scope

master_phil_str = get_qm_restraints_scope()
# master_phil_str = get_qm_gradients_scope()

def electrons(model,
              specific_atom_charges=None,
              specific_atom_multiplicities=None,
              return_atom_valences=False,
              log=None):
  from libtbx.utils import Sorry
  from mmtbx.ligands import electrons
  atom_valences = electrons.electron_distribution(
    model.get_hierarchy(), # needs to be altloc free
    model.get_restraints_manager().geometry,
    specific_atom_charges=specific_atom_charges,
    specific_atom_multiplicities=specific_atom_multiplicities,
    log=log,
    verbose=False,
  )
  if return_atom_valences: return atom_valences
  rc = atom_valences.validate(ignore_water=True,
                              raise_if_error=False)
  # rc = atom_valences.report(ignore_water)
  if rc:
    print('''
  Unusual atom valences
%s
    ''' % str(atom_valences), file=log)
    print(atom_valences)
    for key, item in rc.items():
      print('  %s' % key, file=log)
      for i in item:
        print('    %s' % i[0], file=log)
    # raise Sorry('Unusual charges found')
  charged_atoms = atom_valences.get_charged_atoms()
  print('''
  Complete valence picture
%s
  '''% (atom_valences.show()), file=log)
  return atom_valences.get_total_charge()

def get_safe_filename(s, compact_selection_syntax=True):
  import string
  assert compact_selection_syntax
  if compact_selection_syntax:
    s=s.replace('chain', '')
    s=s.replace('resid', '')
    s=s.replace('resseq', '')
    s=s.replace('resname', '')
    s=s.replace('and', '')
  while s.find('  ')>-1:
    s=s.replace('  ',' ')
  if s[0]==' ': s=s[1:]
  s=s.replace(' ','_')
  for i in range(26):
    a=string.ascii_uppercase[i]
    s=s.replace("'%s'"%a, a)
  s=s.replace("'",'_prime_')
  s=s.replace('*','_star_')
  s=s.replace('(','_lb_')
  s=s.replace(')','_rb_')
  s=s.replace('=', '_equals_')
  s=s.replace(':', '_colon_')
  s=s.replace('"', '_quote_')
  while s.find('__')>-1:
    s=s.replace('__','_')
  return s

def populate_qmr_defaults(qmr):
  def default_defaults(qmr):
    if qmr.package.basis_set is Auto:
      qmr.package.basis_set=''
    if qmr.package.solvent_model is Auto:
      qmr.package.solvent_model=''
    # if qmr.package.multiplicity is Auto:
    #   qmr.package.multiplicity=1
    # if qmr.package.charge is Auto:
    #   qmr.package.charge=0
  program = qmr.package.program
  if program=='test':
    pass
  elif program=='orca':
    default_defaults(qmr)
    if qmr.package.method is Auto:
      qmr.package.method='AM1'
      qmr.package.method='PBEh-3c'
  elif program=='mopac':
    default_defaults(qmr)
    if qmr.package.method is Auto:
      qmr.package.method='PM7'
      qmr.package.method='PM6-D3H4'
  else:
    assert 0
  return qmr

def get_working_directory(model, params, prefix=None):
  rc = 'qm_work_dir'
  return rc
  if prefix is None:
    prefix = getattr(params.output, 'prefix', None)
  if prefix is not None:
    rc='%s_%s' % (prefix, rc)
  return rc

def get_total_multiplicity(qmr):
  tm = 0
  macs = qmr.specific_atom_multiplicities
  for mac in macs:
    tm += mac.multiplicity
  return max(tm,1)

def get_preamble(macro_cycle, i, qmr, old_style=False, compact_selection_syntax=True):
  qmr = populate_qmr_defaults(qmr)
  s=''
  if macro_cycle is not None:
    s+='%02d_' % macro_cycle
  # else:
  #   s+='00_'
  if old_style:
    s+='%02d_%s_%s' % (i+1, get_safe_filename(qmr.selection), qmr.buffer)
  else:
    s+='%s_%s' % (get_safe_filename(qmr.selection,
                                    compact_selection_syntax=compact_selection_syntax),
                  qmr.buffer)
  if qmr.capping_groups:
    s+='_C'
  if qmr.include_nearest_neighbours_in_optimisation:
    s+='_N'
  if 'main_chain_to_delta' in qmr.protein_optimisation_freeze:
    s+='_D'
  elif 'main_chain_to_beta' in qmr.protein_optimisation_freeze:
    s+='_B'
  elif 'main_chain' in qmr.protein_optimisation_freeze:
    s+='_S'
  if 'torsions' in qmr.protein_optimisation_freeze:
    s+='_T'
  if qmr.package.method is not Auto:
    s+='_%s' % get_safe_filename(qmr.package.method)
  if qmr.package.basis_set is not Auto and qmr.package.basis_set:
    s+='_%s' % get_safe_filename(qmr.package.basis_set)
  if qmr.package.solvent_model is not Auto and qmr.package.solvent_model:
    s+='_%s' % get_safe_filename(qmr.package.solvent_model)
  multiplicity=get_total_multiplicity(qmr)
  if multiplicity!=1:
    s+='_%s' % (multiplicity)
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
    working_directory = None
      .type = path
      .style = hidden
      .caption = not implemented
    %s
    %s
  }
''' % (get_qm_restraints_scope(), get_qm_gradients_scope())
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
  if len(params.qi.qm_restraints):
    if validate_qm_restraints(params.qi.qm_restraints, verbose=verbose):
      return True, 'qm_restraints' # includes restraints and energy
  elif len(params.qi.qm_gradients):
    if validate_qm_restraints(params.qi.qm_gradients, verbose=verbose):
      return True, 'qm_gradients'
  return False

def is_quantum_interface_active_this_macro_cycle(params,
                                                 macro_cycle,
                                                 energy_only=False,
                                                 gradients_only=False,
                                                 verbose=False):
  from mmtbx.geometry_restraints.quantum_restraints_manager import running_this_macro_cycle
  qi = is_quantum_interface_active(params)
  if verbose: print('qi',qi,'energy_only',energy_only,'gradients_only',gradients_only,'macro_cycle',macro_cycle)
  if qi:
    rc = []
    if gradients_only:
      rc.append(True)
    elif qi[1]=='qm_restraints':
      number_of_macro_cycles = 1
      if hasattr(params, 'main'):
        number_of_macro_cycles = params.main.number_of_macro_cycles
      for i, qmr in enumerate(params.qi.qm_restraints):
        pre_refinement=True
        if energy_only and macro_cycle==number_of_macro_cycles:
          pre_refinement=False
        tmp = running_this_macro_cycle(qmr,
                                       macro_cycle,
                                       number_of_macro_cycles,
                                       energy_only=energy_only,
                                       pre_refinement=pre_refinement,
                                       verbose=verbose)
        if verbose: print(tmp)
        if tmp: rc.append(True)
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
  from mmtbx.geometry_restraints.quantum_restraints_manager import running_this_macro_cycle
  qi = is_quantum_interface_active(params)
  if not qi: return {}
  if hasattr(params, 'main'):
    number_of_macro_cycles = params.main.number_of_macro_cycles
  else:
    number_of_macro_cycles = 1
  if qi:
    data={}
    for i, qmr in enumerate(params.qi.qm_restraints):
      data[qmr.selection]=[]
      rc=[]
      for j in range(number_of_macro_cycles+1):
        rc.append(unique_item_list())
        yn = running_this_macro_cycle(qmr, j, number_of_macro_cycles)
        if yn:
          rc[j].append(yn)
        pre_refinement=(j!=number_of_macro_cycles)
        yn = running_this_macro_cycle(qmr, j, number_of_macro_cycles,
                                           energy_only=True,
                                           pre_refinement=pre_refinement)
        if yn:
          for k in range(len(yn)):
            yn[k]=yn[k].split('_')[-1]
          rc[j]+=yn
      data[qmr.selection] = rc
      if verbose:
        print('    %s' % qmr.selection, file=log)
        for j, actions in enumerate(rc):
          if actions:
            print('      %2d : %s' % (j, ' '.join(actions)), file=log)
    if 0:
      for key, item in data.items():
        print(key, item)
  return data

def digester(model, geometry, params, log=None):
  active, choice = is_quantum_interface_active(params)
  assert active
  if not model.has_hd():
    from libtbx.utils import Sorry
    raise Sorry('Model must have Hydrogen atoms')
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
