from __future__ import absolute_import,division, print_function
import os

def env_exists_exists(env, var, check=True):
  if check:
    return (env.get(var, False) and os.path.exists(env[var]))
  else:
    return env.get(var, False)

def is_orca_installed(env, var):
  return env_exists_exists(env, var)

def is_qm_test_installed(env, var):
  return env_exists_exists(env, var, check=False)

program_options = {
  'orca' : (is_orca_installed, 'PHENIX_ORCA'),
  'test' : (is_qm_test_installed, 'PHENIX_QM_TEST'),
  }


def get_qm_restraints_scope(verbose=False):
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
  buffer = 5.
    .type = float
    .help = distance to include entire residues into the enviroment of the core
  write_pdb_core = False
    .type = bool
  write_pdb_buffer = False
    .type = bool
  cleanup = True
    .type = bool
  %s
}
'''
  programs = ''
  for package, (func, var) in program_options.items():
    if func(os.environ, var):
      if package=='orca':
        programs += ' *%s' % package
      else:
        programs += '   %s' % package
  if verbose: print(programs)
  qm_package_scope = qm_package_scope % programs
  qm_restraints_scope = qm_restraints_scope % qm_package_scope
  return qm_restraints_scope

def orca_action():
  outl = '''
    orca
      .help = Orca
    {
      include scope mmtbx.geometry_restraints.qm_manager.orca_master_phil_str
    }
  '''
  return outl

def get_safe_filename(s):
  s=s.replace(' ','_')
  s=s.replace("'",'_prime_')
  s=s.replace('*','_star_')
  return s

def get_preamble(macro_cycle, i, selection):
  s=''
  if macro_cycle is not None:
    s+='%02d_' % macro_cycle
  s+='%02d_%s' % (i+1, get_safe_filename(selection))
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

def validate_qm_restraints(qm_restraints):
  for qmr in qm_restraints:
    print ('...',qmr)
    for attr, item in qmr.__dict__.items():
      print(attr,item)

def is_quantum_interface_active(params):
  if not hasattr(params, 'qi'): return False
  if len(params.qi.qm_restraints):
    # validate_qm_restraints(params.qi.qm_restraints)
    return True, 'qm_restraints'
  return False

def digester(model, geometry, params, log=None):
  active, choice = is_quantum_interface_active(params)
  assert active
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
