from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.geometry_restraints.quantum_interface import get_qm_restraints_scope

import iotbx.pdb
import iotbx.phil

get_class = iotbx.pdb.common_residue_names_get_class

def get_selection_from_user(hierarchy):
  j=0
  opts = []
  for residue_group in hierarchy.residue_groups():
    # if len(residue_group.atom_groups())>1:
    #   print('skip')
    #   assert 0
    atom_group = residue_group.atom_groups()[0]
    rc = get_class(atom_group.resname)
    if rc in ['common_amino_acid', 'common_water']: continue
    print(j, residue_group.id_str(), atom_group.resname, rc)
    opts.append('chain %s and resid %s and resname %s' % (
      residue_group.parent().id,
      residue_group.resid(),
      atom_group.resname,
    ))
    j+=1
  print('\n\n')
  for i, sel in enumerate(opts):
    print('    %2d : "%s"' % (i+1,sel))
  rc = input('\n  Enter selection by choosing number or typing a new one ~> ')
  try:
    rc = int(rc)
    rc = opts[rc-1]
  except ValueError:
    pass
  return rc

class Program(ProgramTemplate):

  description = '''
phenix.quantum_inteface: tool for selecting some atoms for QI

Usage examples:
  phenix.quantum_inteface model.pdb "chain A"
  '''

  datatypes = ['model', 'phil', 'restraint']

  master_phil_str = """
  qi_helper {
    selection = None
      .type = atom_selection
      .help = what to select
      .multiple = True
  }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)

  # ---------------------------------------------------------------------------
  def run(self):
    model = self.data_manager.get_model()
    #
    # get selection
    #
    cif_object = None
    if self.data_manager.has_restraints():
      cif_object = self.data_manager.get_restraint()
    if not self.params.qi_helper.selection:
      rc = get_selection_from_user(model.get_hierarchy())
      self.params.qi_helper.selection = [rc]
    #
    # validate selection
    #
    selection_array = model.selection(self.params.qi_helper.selection[0])
    selected_model = model.select(selection_array)
    print(selected_model)

    # for i, line in enumerate(get_qm_restraints_scope().splitlines()):
    #   print(i, line)

    qi_phil_string = get_qm_restraints_scope()
    qi_phil_string = qi_phil_string.replace('selection = None',
                                            'selection = "%s"' % self.params.qi_helper.selection[0])
    qi_phil_string = qi_phil_string.replace('read_output_to_skip_opt_if_available = False',
                                            'read_output_to_skip_opt_if_available = True')
    qi_phil_string = qi_phil_string.replace('capping_groups = False',
                                            'capping_groups = True')

    outl = ''
    for line in qi_phil_string.splitlines():
      if line.find(' write_')>-1 or line.find(' calculate_')>-1:
        tmp=line.split()
        line = '  %s = True' % tmp[0]
      outl += '%s\n' % line
    qi_phil_string = outl
    qi_phil = iotbx.phil.parse(qi_phil_string,
                             # process_includes=True,
                             )
    qi_phil.show()

    qi_phil_string = qi_phil_string.replace('qm_restraints',
                                            'refinement.qi.qm_restraints',
                                            1)
    def safe_filename(s):
      s=s.replace('chain ','')
      s=s.replace('resname ','')
      s=s.replace('resid ','')
      s=s.replace('and ','')
      s=s.replace('   ','_')
      s=s.replace('  ','_')
      s=s.replace(' ','_')
      return s
    pf = '%s_%s.phil' % (
      self.data_manager.get_default_model_name().replace('.pdb',''),
      safe_filename(self.params.qi_helper.selection[0]),
      )
    f=open(pf, 'w')
    for line in qi_phil_string.splitlines():
      if line.strip().startswith('.'): continue
      f.write('%s\n' % line)
    del f

  # ---------------------------------------------------------------------------
  def get_results(self):
    return None
