# LIBTBX_SET_DISPATCHER_NAME phenix.development.qi
from __future__ import absolute_import, division, print_function
import os
from libtbx.program_template import ProgramTemplate

from mmtbx.geometry_restraints.quantum_restraints_manager import update_restraints
from mmtbx.geometry_restraints.quantum_interface import get_qm_restraints_scope

import iotbx.pdb
import iotbx.phil
from libtbx.utils import Sorry

get_class = iotbx.pdb.common_residue_names_get_class

def _add_HIS_H_atom_to_atom_group(ag, name):
  # from scitbx.array_family import flex
  from mmtbx.ligands.ready_set_basics import construct_xyz
  from mmtbx.ligands.ready_set_basics import get_hierarchy_h_atom
 # move to basics
  bonded = {'HD1' : ['ND1', 'CE1', 'NE2'],
            'HE2' : ['NE2', 'CD2', 'CG'],
           }
  atoms = []
  for i in range(3):
    atoms.append(ag.get_atom(bonded[name.strip()][i]))
  ro2 = construct_xyz(atoms[0], 0.9,
                      atoms[1], 126.,
                      atoms[2], 180.,
                     )
  atom = get_hierarchy_h_atom(name, ro2[0], atoms[0])
  ag.append_atom(atom)
  return ag

def add_histidine_H_atoms(hierarchy):
  '''
  HIS      ND1    HD1       coval       0.860    0.020    1.020
  HIS      NE2    HE2       coval       0.860    0.020    1.020
  '''
  for i, ag in enumerate(hierarchy.atom_groups()): pass
  assert i==0
  ag = ag.detached_copy()
  for name in [' HD1', ' HE2']:
    atom = ag.get_atom(name.strip())
    if atom is None:
      ag = _add_HIS_H_atom_to_atom_group(ag, name)
  return ag

def assert_histidine_double_protonated(ag):
  count = 0
  for atom in ag.atoms():
    if atom.name.strip() in ['HD1', 'HE2', 'DD1', 'DE2']:
      count+=1
  if count not in [1]:
    raise Sorry('incorrect protonation of %s' % ag.id_str())

def generate_flipping_his(ag,
                          return_hierarchy=False,
                          include_unprotonated=False,
                          chain_id=None,
                          resseq=None):
  # assume double protonated HIS
  assert_histidine_double_protonated(ag)
  booleans = [[1,1], [1,0], [0,1]]
  if include_unprotonated: booleans = [[1,1], [1,0], [0,1], [0,0]]
  for flip in range(2):
    for i, (hd, he) in enumerate([[1,1], [1,0], [0,1], [0,0]]):
      if i==0 and flip:
        for n1, n2 in [[' ND1', ' CD2'],
                       [' CE1', ' NE2'],
                       [' HD1', ' HD2'],
                       [' HE1', ' HE2'],
                       [' DD1', ' DD2'],
                       [' DE1', ' DE2'],
                      ]:
          a1 = ag.get_atom(n1.strip())
          a2 = ag.get_atom(n2.strip())
          tmp = a1.xyz
          a1.xyz = a2.xyz
          a2.xyz = tmp
      rc = iotbx.pdb.hierarchy.atom_group()
      rc.resname='HIS'
      for atom in ag.atoms():
        if hd==0 and atom.name in [' HD1', ' DD1']: continue
        if he==0 and atom.name in [' HE2', ' DE2']: continue
        atom = atom.detached_copy()
        rc.append_atom(atom)
      if return_hierarchy:
        ph = iotbx.pdb.hierarchy.root()
        m = iotbx.pdb.hierarchy.model()
        c = iotbx.pdb.hierarchy.chain()
        if chain_id is None: c.id='A'
        else: c.id=chain_id
        r = iotbx.pdb.hierarchy.residue_group()
        if resseq is None: r.resseq='1'
        else: r.resseq=resseq
        r.append_atom_group(rc)
        c.append_residue_group(r)
        m.append_chain(c)
        ph.append_model(m)
        yield ph
      else:
        yield rc

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
    print(dir(residue_group))
    # opts.append('chain %s and resid %s and resname %s' % (
    #   residue_group.parent().id,
    #   residue_group.resid(),
    #   atom_group.resname,
    # ))
    for conformer in residue_group.conformers():
      print(dir(conformer))
      for residue in conformer.residues():
        print(dir(residue))
        sel_str = 'chain %s and resid %s and resname %s' % (
            residue_group.parent().id,
            residue_group.resid(),
            residue.resname,
          )
        print(sel_str)
        if residue.is_pure_main_conf:
          opts.append(sel_str)
        else:
          altlocs=[]
          print(dir(altlocs))
          for atom in residue.atoms():
            print(atom.quote())
            print(dir(atom))
            print(atom.parent().altloc)
            altloc = atom.parent().altloc
            if altloc not in altlocs: altlocs.append(altloc)
          print(altlocs)
          ts = []
          for altloc in altlocs:
            ts.append('(%s and altloc %s)' % (sel_str, altloc))
          print(ts)
          opts.append(' or '.join(ts))
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
phenix.qi: tool for selecting some atoms for QI

Usage examples:
  phenix.quantum_inteface model.pdb "chain A"
  '''

  datatypes = ['model', 'phil', 'restraint']

  master_phil_str = """
  qi {
    %s
    selection = None
      .type = atom_selection
      .help = what to select
      .multiple = True
    format = *phenix_refine quantum_interface
      .type = choice
    write_qmr_phil = False
      .type = bool
    run_qmr = False
      .type = bool
    iterate_histidine = None
      .type = atom_selection
  }
""" % (get_qm_restraints_scope())

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    model = self.data_manager.get_model()
    if not model.has_hd():
      raise Sorry('Model must has Hydrogen atoms')
    if self.params.output.prefix is None:
      prefix = os.path.splitext(self.data_manager.get_default_model_name())[0]
      print('  Setting output prefix to %s' % prefix, file=self.logger)
      self.params.output.prefix = prefix
    if self.params.qi.iterate_histidine:
      self.params.qi.selection = [self.params.qi.iterate_histidine]

  # ---------------------------------------------------------------------------
  def run(self, log=None):
    model = self.data_manager.get_model()
    #
    # get selection
    #
    # cif_object = None
    # if self.data_manager.has_restraints():
    #   cif_object = self.data_manager.get_restraint()
    if (not self.params.qi.selection and
        self.params.qi.iterate_histidine is None and
        len(self.params.qi.qm_restraints)==0):
      rc = get_selection_from_user(model.get_hierarchy())
      self.params.qi.selection = [rc]
    #
    # validate selection
    #
    if len(self.params.qi.qm_restraints)!=0:
      selection = self.params.qi.qm_restraints[0].selection
    elif self.params.qi.selection:
      selection=self.params.qi.selection[0]
    selection_array = model.selection(selection)
    selected_model = model.select(selection_array)
    print('Selected model  %s' % selected_model, file=log)
    self.data_manager.add_model('ligand', selected_model)

    if self.params.qi.write_qmr_phil:
      self.write_qmr_phil(self.params.qi.format)

    if self.params.qi.run_qmr:
      self.params.qi.qm_restraints.selection=self.params.qi.selection
      self.run_qmr(self.params.qi.format)

    if self.params.qi.iterate_histidine:
      self.iterate_histidine(self.params.qi.iterate_histidine)

  def iterate_histidine(self, selection, log=None):
    if len(self.params.qi.qm_restraints)<1:
      self.write_qmr_phil(iterate_histidine=True)
      print('Restart command with PHIL file')
      return
    model = self.data_manager.get_model()
    selection_array = model.selection(selection)
    selected_model = model.select(selection_array)
    hierarchy = selected_model.get_hierarchy()
    # add all H atoms
    his_ag = add_histidine_H_atoms(hierarchy)
    for atom in hierarchy.atoms(): break
    rg_resseq = atom.parent().parent().resseq
    chain_id = atom.parent().parent().parent().id
    for i, flipping_his in enumerate(generate_flipping_his(his_ag)):
      model = self.data_manager.get_model()
      hierarchy = model.get_hierarchy()
      for chain in hierarchy.chains():
        if chain.id!=chain_id: continue
        for rg in chain.residue_groups():
          if rg.resseq!=rg_resseq: continue
          for j, ag in enumerate(rg.atom_groups()): pass
          assert j==0
          rg.remove_atom_group(ag)
          rg.insert_atom_group(0, flipping_his)
      # hierarchy.write_pdb_file('his_%02d.pdb' % i)
      self.params.output.prefix='iterate_histidine_%02d' % (i+1)
      # self.params.output.prefix='his_%02d' % (i+1)
      update_restraints(model,
                        self.params,
                        never_write_restraints=True,
                        log=log)
      # run clashscore

  def run_qmr(self, format, log=None):
    model = self.data_manager.get_model()
    rc = update_restraints( model,
                            self.params,
                            log=log,
                            )

  def write_qmr_phil(self, iterate_histidine=False, log=None):
    qi_phil_string = get_qm_restraints_scope()
    qi_phil_string = qi_phil_string.replace('selection = None',
                                            'selection = "%s"' % self.params.qi.selection[0])
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
    if iterate_histidine:
      qi_phil_string = qi_phil_string.replace('refinement.', '')
      qi_phil_string = qi_phil_string.replace('ignore_x_h_distance_protein = False',
                                              'ignore_x_h_distance_protein = True')

    if format=='quantum_inteface':
      qi_phil_string = qi_phil_string.replace('refinement.qi', 'qi')

    def safe_filename(s):
      s=s.replace('chain ','')
      s=s.replace('resname ','')
      s=s.replace('resid ','')
      s=s.replace('and ','')
      s=s.replace('   ','_')
      s=s.replace('  ','_')
      s=s.replace(' ','_')
      s=s.replace('(','')
      s=s.replace(')','')
      return s

    pf = '%s_%s.phil' % (
      self.data_manager.get_default_model_name().replace('.pdb',''),
      safe_filename(self.params.qi.selection[0]),
      )
    print('  Writing QMR phil scope to %s' % pf, file=log)
    f=open(pf, 'w')
    for line in qi_phil_string.splitlines():
      if line.strip().startswith('.'): continue
      f.write('%s\n' % line)
    del f

  # ---------------------------------------------------------------------------
  def get_results(self):
    return None
