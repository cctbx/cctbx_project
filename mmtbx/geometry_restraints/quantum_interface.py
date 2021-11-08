import os

def is_orca_installed(env):
  return (env.get('PHENIX_ORCA', False) and os.path.exists(env['PHENIX_ORCA']))

programs = ' *orca'

qm_package_scope = '''
  package
  {
    program = %(programs)s
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
''' % locals()

qm_restraints_scope = '''
qm_restraints
  .multiple = True
{
  selection = None
    .type = atom_selection
  buffer = 5.
    .type = float
  %(qm_package_scope)s
}
''' % locals()

def orca_action():
  outl = '''
    orca
      .help = Orca
    {
      include scope mmtbx.geometry_restraints.qm_manager.orca_master_phil_str
    }
  '''
  return outl

def is_any_quantum_package_installed(env):
  installed = []
  actions = []
  for key, (question, action) in {'orca' : (is_orca_installed, orca_action),
                                  }.items():
    if question(os.environ):
      rc = action()
      actions.append(rc)
      installed.append(key)
  if installed:
    outl = '''
  qi
    .help = QM
    .expert_level = 3
  {
    use_quantum_interface = *None qm_restraints qm_gradients
      .type = bool
    selection = None
      .type = atom_selection
    charge = 0
      .type = int
    multiplicity = 1
      .type = int
    buffer = 0.
      .type = float
      .style = hidden
    update_metal_coordination = False
      .type = bool
      .style = hidden
    refine_buffer_hydrogen_atoms = False
      .type = bool
      .style = hidden
'''
    for action in actions:
      outl += action
    outl += '}'
    outl = '''
  qi
    .help = QM
    .expert_level = 3
  {
    %s
  }
''' % qm_restraints_scope
  return outl

def validate_qm_restraints(qm_restraints):
  for qmr in qm_restraints:
    print ('...',qmr)
    for attr, item in qmr.__dict__.items():
      print(attr,item)

def is_quantum_interface_active(params):
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
  pass

if __name__ == '__main__':
  main()
