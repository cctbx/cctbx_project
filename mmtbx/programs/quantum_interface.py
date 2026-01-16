# LIBTBX_SET_DISPATCHER_NAME phenix.development.qi
from __future__ import absolute_import, division, print_function
import os
import sys
import time
import copy

from libtbx.program_template import ProgramTemplate

from mmtbx.monomer_library.linking_setup import ad_hoc_single_metal_residue_element_types

from mmtbx.geometry_restraints.quantum_restraints_manager import run_energies
from mmtbx.geometry_restraints.quantum_restraints_manager import run_gradients
from mmtbx.geometry_restraints.quantum_restraints_manager import update_restraints
from mmtbx.geometry_restraints.quantum_restraints_manager import min_dist2
from mmtbx.geometry_restraints.quantum_interface import get_qm_restraints_scope
from mmtbx.geometry_restraints.quantum_interface import get_qm_gradients_scope
from mmtbx.geometry_restraints.quantum_interface import get_working_directory
from mmtbx.geometry_restraints.qi_utils import classify_histidine
from mmtbx.geometry_restraints.qi_utils import run_serial_or_parallel
from mmtbx.geometry_restraints.qi_utils import get_hbonds_via_filenames
from mmtbx.geometry_restraints.qi_utils import get_rotamers_via_filenames

import iotbx.pdb
import iotbx.phil
from libtbx import Auto
from libtbx.utils import Sorry
from libtbx.utils import null_out
from libtbx.utils import multi_out

from mmtbx.monomer_library.linking_setup import ad_hoc_single_metal_residue_element_types

get_class = iotbx.pdb.common_residue_names_get_class

def safe_filename(s):
  while s.find('  ')>-1:
    s=s.replace('  ',' ')
  s=s.replace('chain ','')
  s=s.replace('resname ','')
  s=s.replace('resid ','')
  s=s.replace('and ','')
  s=s.replace('altloc ','')
  s=s.replace(' ','_')
  s=s.replace('(','')
  s=s.replace(')','')
  s=s.replace(':','_')
  s=s.replace("'",'')
  return s

def stepper(b,e,s):
  b=float(b)
  e=float(e)
  s=float(s)
  for i in range(100):
    yield b
    b+=s
    if b>e: break

def substitute_ag_into_hierarchy(hierarchy, atom_group):
  for ag in hierarchy.atom_groups():
    if ag.id_str()==atom_group.id_str():
      rg=ag.parent()
      # not alt loc aware
      rg.remove_atom_group(ag)
      rg.append_atom_group(atom_group.detached_copy())
      break

def merge_water(model, filenames, chain_id='A', carve_selection=None):
  hierarchy=model.get_hierarchy()
  selected_model=None
  if carve_selection:
    selection_array = model.selection(carve_selection)
    selected_model = model.select(selection_array)
    for ag in hierarchy.atom_groups():
      if ag.resname=='DOD': assert 0
      if ag.resname=='HOH':
        remove=[]
        for atom in ag.atoms():
          if atom.element.strip()=='H':
            remove.append(atom)
        for atom in remove:
          ag.remove_atom(atom)
  hierarchies = []
  for filename in filenames:
    ph = iotbx.pdb.input(filename).construct_hierarchy()
    hierarchies.append(ph)
  waters = {}
  for ph in hierarchies:
    # ph.show()
    if selected_model:
      smp=selected_model.get_hierarchy()
      for hoh in ph.atom_groups(): break
      for ag in smp.atom_groups():
        if ag.id_str()==hoh.id_str(): break
      else:
        continue
    for ag in ph.atom_groups():
      substitute_ag_into_hierarchy(hierarchy, ag)
  return

  outl = ''
  for id_str, ags in waters.items():
    if len(ags)>1 or 1:
      for i, ag in enumerate(ags):
        # ag.altloc='ABCDEFGHIJKLMNOPQRSTUVWXYZ'[i]
        for atom in ag.atoms():
          outl += '%s\n' % atom.format_atom_record()
  print('-'*80)
  print(outl)
  f=open('water.pdb', 'w')
  f.write(outl)
  del f

def _add_HIS_H_atom_to_atom_group(ag, name):
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

def _add_NQ_H_atom_to_atom_group(ag, name):
  from mmtbx.ligands.ready_set_basics import construct_xyz
  from mmtbx.ligands.ready_set_basics import get_hierarchy_h_atom
  bonded = {'HD21' : ['ND2', 'CG', 'CB'], #ASN
            'HD22' : ['ND2', 'CG', 'OD1'],#ASN
            'HE21' : ['NE2', 'CD', 'CG'], #GLN
            'HE22' : ['NE2', 'CD', 'OE1'],#GLN
           }
  atoms = []
  for i in range(3):
    atoms.append(ag.get_atom(bonded[name.strip()][i]))
  ro2 = construct_xyz(atoms[0], 0.9,
                      atoms[1], 120.,
                      atoms[2], 180.,
                     )
  atom = get_hierarchy_h_atom(name, ro2[0], atoms[0])
  ag.append_atom(atom)
  return ag

def _merge_atom_groups(hierarchy):
  bag=None
  for i, ag in enumerate(hierarchy.atom_groups()):
    if i:
      rg=bag.parent()
      tag=ag.detached_copy()
      for atom in tag.atoms():
        bag.append_atom(atom.detached_copy())
      rg.remove_atom_group(ag)
    else:
      bag=ag
  return hierarchy

def get_first_ag(hierarchy):
  hierarchy=_merge_atom_groups(hierarchy)
  for i, ag in enumerate(hierarchy.atom_groups()): pass
  assert i==0
  return ag.detached_copy()

def add_histidine_H_atoms(hierarchy):
  '''
  HIS      ND1    HD1       coval       0.860    0.020    1.020
  HIS      NE2    HE2       coval       0.860    0.020    1.020
  '''
  ag = get_first_ag(hierarchy)
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
  if count not in [2]:
    raise Sorry('incorrect protonation of %s' % ag.id_str())

def construct_hierarchy(ag, chain_id=None, resseq=None):
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
  return ph

def generate_flipping_his(ag,
                          return_hierarchy=False,
                          include_unprotonated=False,
                          chain_id=None,
                          resseq=None):
  assert_histidine_double_protonated(ag)
  ag=ag.detached_copy()
  booleans = [[1,1], [1,0], [0,1]]
  if include_unprotonated: booleans = [[1,1], [1,0], [0,1], [0,0]]
  for flip in range(2):
    for i, (hd, he) in enumerate(booleans):
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
          if a1 is None or a2 is None: continue
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
        yield construct_hierarchy(rc, chain_id=chain_id, resseq=resseq)
      else:
        yield rc

def generate_flipping_NQ(ag,
                         return_hierarchy=False,
                         chain_id=None,
                         resseq=None):
  ag=ag.detached_copy()
  if ag.resname=='ASN':
    hs = ['HD21', 'HD22']
    nos = [' ND2', ' OD1']
  elif ag.resname=='GLN':
    hs = ['HE21', 'HE22']
    nos = [' NE2', ' OE1']
  else:
    assert ag.resname in ['ASN', 'GLN'], 'resname %s' % ag.resname
  for atom in ag.atoms():
    if atom.name.strip() in hs:
      ag.remove_atom(atom)
  for flip in range(2):
    if flip:
      for n1, n2 in [nos]:
        a1 = ag.get_atom(n1.strip())
        a2 = ag.get_atom(n2.strip())
        if a1 is None or a2 is None: continue
        tmp = a1.xyz
        a1.xyz = a2.xyz
        a2.xyz = tmp
    rc = iotbx.pdb.hierarchy.atom_group()
    rc.resname=ag.resname
    for atom in ag.atoms():
      atom = atom.detached_copy()
      rc.append_atom(atom)
    for name in hs:
      _add_NQ_H_atom_to_atom_group(rc, name)
    if return_hierarchy:
      yield construct_hierarchy(rc, chain_id=chain_id, resseq=resseq)
    else:
      yield rc

def should_get_selection_from_user(params):
  def qm_restraints_has_selection(params):
    if len(params.qi.qm_restraints)==0:
      return True
    if not params.qi.qm_restraints[0].selection:
      return True
    return False
  bools = [
    qm_restraints_has_selection(params),
    not params.qi.each_amino_acid,
    not params.qi.each_water,
    not params.qi.merge_water,
    ]
  bc=0
  if not params.qi.selection:
    for i, boolean in enumerate(bools):
      if boolean: bc+=1
  else:
    return False
  return bc==len(bools)

def get_selection_from_user(hierarchy, include_amino_acids=None, return_list=False, log=None):
  j=0
  opts = []
  if include_amino_acids is not None:
    for a,b in enumerate(include_amino_acids):
      include_amino_acids[a]=b.upper()
  for residue_group in hierarchy.residue_groups():
    atom_group = residue_group.atom_groups()[0]
    rc = get_class(atom_group.resname)
    if include_amino_acids and atom_group.resname in include_amino_acids: pass
    elif (rc!='common_rna_dna' and
          atom_group.resname.strip() in ad_hoc_single_metal_residue_element_types):
      pass # ions
    elif rc in ['common_amino_acid', 'common_water', 'common_rna_dna']: continue
    for conformer in residue_group.conformers():
      for residue in conformer.residues():
        sel_str = 'chain %s and resid %s and resname %s' % (
            residue_group.parent().id,
            residue_group.resid(),
            residue.resname,
          )
        if residue.is_pure_main_conf:
          opts.append(sel_str)
        else:
          altlocs=[]
          for atom in residue.atoms():
            altloc = atom.parent().altloc
            if altloc not in altlocs: altlocs.append(altloc)
          ts = []
          for altloc in altlocs:
            if not altloc: altloc=' '
            ts.append("(%s and altloc '%s')" % (sel_str, altloc))
          opts.append(' or '.join(ts))
    j+=1
  print('\n\n', file=log)
  if return_list: return opts
  for i, sel in enumerate(opts):
    print('    %2d : "%s"' % (i+1,sel), file=log)
  if len(opts)==1:
    print('\n  Automatically selecting', file=log)
    rc=[opts[0]]
    return rc
  else:
    rc = input('\n  Enter selection by choosing number or typing a new one ~> ')
  try:
    rc = int(rc)
    rc = [opts[rc-1]]
  except ValueError:
    pass
  except IndexError:
    rc = ['resid 1']
  return rc

class Program(ProgramTemplate):

  description = '''
phenix.qi: tool for selecting some atoms for QI

Usage examples:
  phenix.quantum_interface model.pdb "chain A"
  '''

  datatypes = ['model', 'phil', 'restraint']

  master_phil_str = """
  qi {
    %s
    selection = None
      .type = atom_selection
      .help = what to select
      .multiple = True
    buffer_selection = None
      .type = atom_selection
      .help = what to select for buffer
      .style = hidden
    format = *phenix_refine qi
      .type = choice
    write_qmr_phil = Auto
      .type = bool
    write_qmg_phil = False
      .type = bool
    run_qmr = False
      .type = bool
    run_directory = None
      .type = str
    iterate_NQH = HIS ASN GLN
      .type = choice
    proton_energy_difference = None
      .type = float
    only_i = None
      .type = int
    step_buffer_radius = None
      .type = str
    iterate_metals = None
      .type = str
    include_amino_acids = None
      .type = str
    carve_selection = None
      .type = atom_selection
      .help = only work on this selection
    each_amino_acid = False
      .type = bool
    each_water = False
      .type = bool
    merge_water = False
      .type = bool
    nproc = 1
      .type = int
    randomise_selection = None
      .type = float
    verbose = False
      .type = bool
    debug = False
      .type = bool
      .style = hidden
      .help = more output including PDB files from each step

  }
""" % (get_qm_restraints_scope())

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    model = self.data_manager.get_model()
    if not model.has_hd():
      raise Sorry('Model must have Hydrogen atoms')
    if self.params.output.prefix is None:
      prefix = os.path.splitext(self.data_manager.get_default_model_name())[0]
      print('  Setting output prefix to %s' % prefix, file=self.logger)
      self.params.output.prefix = prefix
    if self.params.qi.randomise_selection and self.params.qi.randomise_selection>0.5:
      raise Sorry('Random select value %s is too large' % self.params.qi.randomise_selection)
    if self.params.qi.write_qmr_phil is Auto:
      if self.params.qi.run_qmr:
        self.params.qi.write_qmr_phil=False

  # ---------------------------------------------------------------------------
  def run(self, log=None):
    self.results=[]
    self.logger = multi_out()
    self.logger.register("stdout", sys.stdout)
    log_filename = '%s.log' % self.data_manager.get_default_output_filename()
    if self.params.qi.run_qmr:
      if os.path.exists(log_filename):
        for i in range(1,1000):
          tf = log_filename.replace('000', '%03d'% i)
          if not os.path.exists(tf):
            log_filename=tf
            break
      log_file = open(log_filename, 'w')
      self.logger.register('logfile', log_file)
    print('''
    Running mmtbx.quantum
      excutable: %s
      arguments: %s
          ''' % (sys.executable, ' '.join(sys.argv[1:])))
    model = self.data_manager.get_model()
    if self.params.qi.each_water:
      model.add_hydrogens(1., occupancy=1.)
    self.restraint_filenames = []
    rc = self.data_manager.get_restraint_names()
    for f in rc:
      self.restraint_filenames.append(os.path.abspath(f))
    #
    # get selection
    #
    if self.params.qi.iterate_NQH and len(self.params.qi.qm_restraints)==0:
      resname = self.params.qi.iterate_NQH
      if not self.params.qi.selection:
        rc = get_selection_from_user(model.get_hierarchy(),
                                     include_amino_acids=[resname])
        if rc[0].find('resname %s' % resname)>-1:
          self.params.qi.format='qi'
        self.params.qi.selection = rc
    #
    include_amino_acids=self.params.qi.include_amino_acids
    if include_amino_acids:
      include_amino_acids=[include_amino_acids]

    rc = should_get_selection_from_user(self.params)
    if rc:
      rc = get_selection_from_user(model.get_hierarchy(),
                                   include_amino_acids=include_amino_acids)
      self.params.qi.selection = rc
      if len(self.params.qi.qm_restraints)!=0:
        for attr, item in self.params.qi.qm_restraints[0].__dict__.items():
          print(attr, item)
    #
    # validate selection
    #
    include_nearest_neighbours=False
    selection=None
    if len(self.params.qi.qm_restraints)!=0:
      selection = self.params.qi.qm_restraints[0].selection
    elif self.params.qi.selection:
      selection=self.params.qi.selection[0]
    if selection:
      selection_array = model.selection(selection)
      selected_model = model.select(selection_array)
      print('Selected model:\n%s' % selected_model, file=self.logger)
      self.data_manager.add_model('ligand', selected_model)
      ags = selected_model.get_hierarchy().atom_groups()
      names = []
      for ag in ags: names.append(ag.resname.strip())
      if len(names)==1:
        if names[0] in ad_hoc_single_metal_residue_element_types:
          include_nearest_neighbours=True

    qm_work_dir = get_working_directory(model, self.params)

    if self.params.qi.step_buffer_radius:
      step_buffer_radius = self.params.qi.step_buffer_radius
      assert len(step_buffer_radius.split(','))==3
    #
    # options
    #
    if self.params.qi.each_amino_acid:
      hierarchy = model.get_hierarchy()
      outl = ''
      for rg in hierarchy.residue_groups():
        if len(rg.atom_groups())!=1: continue
        resname=rg.atom_groups()[0].resname
        include_amino_acids=self.params.qi.include_amino_acids
        if include_amino_acids and include_amino_acids!=resname: continue
        gc = get_class(resname)
        if gc not in ['common_amino_acid', 'modified_amino_acid']: continue
        selection = 'chain %s and resid %s' % (rg.parent().id, rg.resseq.strip())
        qi_phil_string = self.get_single_qm_restraints_scope(selection)
        qi_phil_string = qi_phil_string.replace('ignore_x_h_distance_protein = False',
                                                'ignore_x_h_distance_protein = True')
        qi_phil_string = qi_phil_string.replace(' pdb_final_buffer', ' *pdb_final_buffer')
        print('  writing phil for %s %s' % (rg.id_str(), rg.atom_groups()[0].resname), file=self.logger)
        outl += '%s' % qi_phil_string
      pf = '%s_all.phil' % (
        self.data_manager.get_default_model_name().replace('.pdb',''))
      f=open(pf, 'w')
      f.write('refinement.qi {\n')
      for line in outl.splitlines():
        if line.strip().startswith('.'): continue
        f.write('%s\n' % line)
      f.write('}\n')
      del f
      print('  phenix.refine %s %s qi.nproc=6' % (self.data_manager.get_default_model_name(),
                                                  pf), file=self.logger)
      return

    if self.params.qi.merge_water:
      pf = self.data_manager.get_default_model_name().replace('.pdb', '')
      if os.path.exists(qm_work_dir):
        os.chdir(qm_work_dir)
        filenames=[]
        for filename in os.listdir('.'):
          if not filename.startswith(pf): continue
          if not filename.endswith('pdb'): continue
          if filename.find('HOH')==-1: continue
          if filename.find('_cluster_final_')==-1: continue
          filenames.append(filename)
        merge_water(filenames)
      else:
        print('no working directory', file=self.logger)

    if self.params.qi.each_water:
      def _water_scope_generator(hierarchy, carve_selection=None):
        if carve_selection:
          selection_array = model.selection(carve_selection)
          selected_model = model.select(selection_array)
          hierarchy = selected_model.get_hierarchy()
        else:
          hierarchy = model.get_hierarchy()
        for rg in hierarchy.residue_groups():
          if len(rg.atom_groups())!=1: continue
          resname=rg.atom_groups()[0].resname
          include_amino_acids=self.params.qi.include_amino_acids
          if include_amino_acids and include_amino_acids!=resname: continue
          gc = get_class(resname)
          if gc not in ['common_water']: continue
          selection = 'chain %s and resid %s and resname HOH' % (rg.parent().id, rg.resseq.strip())
          qi_phil_string = self.get_single_qm_restraints_scope(selection)
          qi_phil_string = self.set_one_write_to_true(qi_phil_string, 'pdb_final_core')
          qi_phil_string = qi_phil_string.replace('ignore_x_h_distance_protein = False',
                                                  'ignore_x_h_distance_protein = True')
          qi_phil_string = qi_phil_string.replace('do_not_update_restraints = False',
                                                  'do_not_update_restraints = True')
          # print('  writing phil for %s %s' % (rg.id_str(), rg.atom_groups()[0].resname), file=self.logger)
          yield qi_phil_string
      if not self.params.qi.run_qmr:
        outl = ''
        for qi_phil_string in _water_scope_generator(model,
                                                     carve_selection=self.params.qi.carve_selection):
          outl += '%s' % qi_phil_string
        pf = '%s_water.phil' % (
          self.data_manager.get_default_model_name().replace('.pdb',''))
        f=open(pf, 'w')
        f.write('qi {\n')
        for line in outl.splitlines():
          if line.strip().startswith('.'): continue
          f.write('%s\n' % line)
        f.write('}\n')
        del f

        ih = 'each_water=True'
        if self.params.qi.carve_selection:
          ih += ' carve_selection="%s"' % self.params.qi.carve_selection
        print('''

        mmtbx.quantum_interface %s run_qmr=True %s %s
        ''' % (self.data_manager.get_default_model_name(),
               ih,
               pf), file=self.logger)

      else:
        nproc=self.params.qi.nproc
        argstuples=[]
        # for qmr in self.params.qi.qm_restraints:
        for i, qi_phil_string in enumerate(_water_scope_generator(model,
                                                                  carve_selection=self.params.qi.carve_selection)):
          if nproc==-1:
            print('  Running %s flip %d' % (nq_or_h, i+1), file=self.logger)
            res = update_restraints(model,
                                    self.params,
                                    never_write_restraints=True,
                                    # never_run_strain=True,
                                    log=arg_log)
            energies.append(energies)
            units=res.units
            rmsds.append(res.rmsds[0][1])
            print('    Energy : %s %s' % (energies[-1],units), file=self.logger)
            print('    Time   : %ds' % (time.time()-t0), file=self.logger)
            print('    RMSD   : %8.3f' % res.rmsds[0][1], file=self.logger)
          else:
            params = copy.deepcopy(self.params)
            for j in range(len(params.qi.qm_restraints)-1,-1,-1):
              if i==j: continue
              del params.qi.qm_restraints[j]
            if len(params.qi.qm_restraints)==0: continue
            print('  Queued selection "%s"' % params.qi.qm_restraints[0].selection, file=self.logger)
            argstuples.append(( model,
                                params,
                                None, # macro_cycle=None,
                                True, #never_write_restraints=False,
                                1, # nproc=1,
                                None, #null_out(),
                                ))
        results = run_serial_or_parallel(update_restraints, argstuples, nproc)
        args = []
        for result in results:
          for filenames in result.final_pdbs:
            args.append(filenames[-1])

        cwd=os.getcwd()
        os.chdir(qm_work_dir)
        merge_water(model, args, carve_selection=self.params.qi.carve_selection)
        os.chdir(cwd)
        model = self.data_manager.get_model()
        fn=self.data_manager.get_default_model_name()
        fnw=fn.replace('.pdb', '_water.pdb')
        self.data_manager.write_model_file(model, fnw, overwrite=True)
        cmd = '\n\tphenix.start_coot %s %s' % (fn, fnw)
        print(cmd)

      return

    if self.params.qi.randomise_selection:
      import random
      sites_cart = model.get_sites_cart()
      for i, (b, xyz) in enumerate(zip(selection_array, sites_cart)):
        if b:
          xyz=list(xyz)
          for j in range(3):
            xyz[j] += random.gauss(0, self.params.qi.randomise_selection)
          sites_cart[i]=tuple(xyz)
      model.set_sites_cart(sites_cart)

    if self.params.qi.write_qmg_phil:
      pf = self.write_qmg_phil()
      program='phenix.refine'
      ih=''
      ih2 = self.data_manager.get_default_model_name()
      if ih2.endswith('.updated.pdb'):
        ih2 = ih2.replace('.updated.pdb', '.mtz')
      elif ih2.endswith('H.pdb'):
        ih2 = ih2.replace('H.pdb', '-sf.cif')
      else:
        ih2 = ' %s' % 'test.mtz'

      print('''

      %s %s %s %s %s
      ''' % (program,
             self.data_manager.get_default_model_name(),
             ih,
             pf,
             ih2,
             ), file=self.logger)
      return

    if self.params.qi.write_qmr_phil:
      # for attr, item in self.params.qi.qm_restraints[0].__dict__.items():
      #   print(attr, item)
      pf = self.write_qmr_phil(iterate_NQH=self.params.qi.iterate_NQH,
                               iterate_metals=self.params.qi.iterate_metals,
                               step_buffer_radius=self.params.qi.step_buffer_radius,
                               output_format=self.params.qi.format,
                               include_nearest_neighbours=include_nearest_neighbours,
                               )
      ih = ''
      if self.params.qi.iterate_NQH:
        ih = 'iterate_NQH=%s' % self.params.qi.iterate_NQH
        ih = ih.replace('  ',' ')
      if self.params.qi.iterate_metals:
        ih = 'iterate_metals="%s"' % self.params.qi.iterate_metals
      if self.params.qi.step_buffer_radius:
        ih = 'step_buffer_radius="%s"' % self.params.qi.step_buffer_radius

      program = 'mmtbx.quantum_interface'
      ih2 = ' run_qmr=True'
      if self.params.qi.format=='qi':
        ih2 += ' qi.nproc=%s' % self.params.qi.nproc
      else:
        program='phenix.refine'
        ih2 = self.data_manager.get_default_model_name()
        if ih2.endswith('.updated.pdb'):
          ih2 = ih2.replace('.updated.pdb', '.mtz')
        else:
          ih2 = ' %s' % 'test.mtz'

      print('''

      %s %s %s %s %s
      ''' % (program,
             self.data_manager.get_default_model_name(),
             ih,
             pf,
             ih2,
             ), file=self.logger)
      return

    if self.params.qi.run_directory:
      if not os.path.exists(self.params.qi.run_directory):
        os.mkdir(self.params.qi.run_directory)
      os.chdir(self.params.qi.run_directory)

    if self.params.qi.run_qmr and self.params.qi.step_buffer_radius:
      ph=selected_model.get_hierarchy()
      for i, atom_group in enumerate(ph.atom_groups()):
        id_str = atom_group.id_str()
      assert i==0
      rc = self.step_thru_buffer_radii(id_str=id_str, log=log)
      return

    if ( self.params.qi.run_qmr and
         not (self.params.qi.iterate_NQH or
              self.params.qi.iterate_metals)):
      self.params.qi.qm_restraints.selection=self.params.qi.selection
      rcs = self.run_qmr(self.params.qi.format)
      self.results.append(rcs)
      dmn = self.data_manager.get_default_model_name()
      cmd = '\n\n\tphenix.start_coot %s' % dmn
      for rc in rcs:
        filenames = getattr(rc, 'final_pdbs', [])
        if filenames:
          for f in filenames:
            for g in f:
              cmd += ' %s' % os.path.join('qm_work_dir', g)
      for r in ['.pdb', '.updated.pdb']:
        mf = dmn.replace(r, '_map_coeffs.mtz')
        if os.path.exists(mf):
          cmd += ' --auto=%s' % mf
          break
      print('%s\n\n' % cmd)

    if self.params.qi.iterate_NQH:
      print('"%s"' % self.params.qi.iterate_NQH, file=self.logger)
      assert self.params.qi.randomise_selection==None
      if self.params.qi.iterate_NQH=='HIS':
        rc = self.iterate_histidine()
      elif self.params.qi.iterate_NQH=='ASN':
        rc = self.iterate_ASN()
      elif self.params.qi.iterate_NQH=='GLN':
        rc = self.iterate_GLN()
      else:
        assert 0
      if rc:
        update, filename = rc
        self.update_residue(filename)
      model = self.data_manager.get_model()
      fn=self.data_manager.get_default_model_name()
      self.data_manager.write_model_file(model,
                                         fn.replace('.pdb', '_best.pdb'),
                                         overwrite=True)

    if self.params.qi.iterate_metals:
      self.iterate_metals()

  def update_residue(self, filename):
    def dist2(r1,r2): return (r1[0]-r2[0])**2+(r1[1]-r2[1])**2+(r1[2]-r2[2])**2
    selection = self.params.qi.qm_restraints[0].selection
    model = self.data_manager.get_model()
    from iotbx import pdb
    pdb_inp = pdb.input(filename)
    hierarchy = pdb_inp.construct_hierarchy()
    asc = hierarchy.atom_selection_cache()
    sel = asc.selection(selection)
    hierarchy = hierarchy.select(sel)
    for replace in hierarchy.residue_groups(): pass
    hierarchy = model.get_hierarchy()
    for h_model in hierarchy.models():
      for chain in h_model.chains():
        for rgi, residue_group in enumerate(chain.residue_groups()):
          for atom1 in replace.atoms(): break
          for atom2 in residue_group.atoms(): break
          d2 = dist2(atom1.xyz, atom2.xyz)
          if d2<1e-3:
            chain.remove_residue_group(residue_group)
            chain.insert_residue_group(rgi, replace)

  def get_selected_hierarchy(self):
    selection = self.params.qi.qm_restraints[0].selection
    model = self.data_manager.get_model()
    selection_array = model.selection(selection)
    selected_model = model.select(selection_array)
    hierarchy = selected_model.get_hierarchy()
    return hierarchy

  def iterate_metals(self, log=None):
    from mmtbx.geometry_restraints.quantum_interface import get_preamble
    if len(self.params.qi.qm_restraints)<1:
      self.write_qmr_phil(iterate_metals=True)
      print('Restart command with PHIL file', file=self.logger)
      return
    def generate_metals(s):
      if s.lower()=='all': s='LI,NA,MG,K,CA,CU,ZN'
      s=s.replace(' ',',')
      for t in s.split(','):
        yield t.upper()
    selection = self.params.qi.qm_restraints[0].selection
    nproc = self.params.qi.nproc
    preamble = get_preamble(None, 0, self.params.qi.qm_restraints[0])
    model = self.data_manager.get_model()
    selection_array = model.selection(selection)
    selected_model = model.select(selection_array)
    hierarchy = selected_model.get_hierarchy()
    for atom in hierarchy.atoms(): break
    rg_resseq = atom.parent().parent().resseq
    chain_id = atom.parent().parent().parent().id
    energies = []
    rmsds = []
    argstuples = []
    filenames = []
    t0=time.time()
    for element in generate_metals(self.params.qi.iterate_metals):
      print('Substituting element : %s' % element, file=self.logger)
      model = self.data_manager.get_model()
      if nproc>1: model=model.deep_copy()
      hierarchy = model.get_hierarchy()
      for chain in hierarchy.chains():
        if chain.id!=chain_id: continue
        for rg in chain.residue_groups():
          if rg.resseq!=rg_resseq: continue
          for j, ag in enumerate(rg.atom_groups()): pass
          assert j==0
          eag = ag.detached_copy()
          eag.resname = element
          eag.atoms()[0].element=element
          eag.atoms()[0].name=element
          tselection = selection.replace(' %s' % ag.resname.upper().strip(),
                                        ' %s' % element)
          self.params.qi.qm_restraints[0].selection = tselection
          rg.remove_atom_group(ag)
          rg.insert_atom_group(0, eag)
      self.params.output.prefix='iterate_metals'
      arg_log=null_out()
      if self.params.qi.verbose:
        arg_log=log
      filenames.append('%s_cluster_final_%s.pdb' % (self.params.output.prefix,
                                                   preamble))
      if nproc==-1:
        print('  Running metal swap %s' % (element), file=self.logger)
        res = update_restraints(model,
                                self.params,
                                never_write_restraints=True,
                                # never_run_strain=True,
                                log=arg_log)
        energies.append(energies)
        units=res.units
        rmsds.append(res.rmsds[0][1])
        print('    Energy : %s %s' % (energies[-1],units), file=self.logger)
        print('    Time   : %ds' % (time.time()-t0), file=self.logger)
        print('    RMSD   : %8.3f' % res.rmsds[0][1], file=self.logger)
      else:
        params = copy.deepcopy(self.params)
        argstuples.append(( model,
                            params,
                            None, # macro_cycle=None,
                            True, #never_write_restraints=False,
                            1, # nproc=1,
                            null_out(),
                            ))
    results = run_serial_or_parallel(update_restraints, argstuples, nproc)
    for i, filename in enumerate(filenames):
      assert os.path.exists(filename), ' Output %s missing' % filename
      results[i].filename = filename
    return results

  def iterate_NQH(self, nq_or_h, classify_nqh, add_nqh_H_atoms, generate_flipping, log=None):
    from mmtbx.geometry_restraints.quantum_interface import get_preamble
    if len(self.params.qi.qm_restraints)<1:
      pf = self.write_qmr_phil(iterate_NQH=True)
      print('Restart command with PHIL file : %s' % pf, file=self.logger)
      return
    qm_work_dir = get_working_directory(self.data_manager.get_model(), self.params)
    nproc = self.params.qi.nproc
    preamble = get_preamble(None, 0, self.params.qi.qm_restraints[0])
    hierarchy = self.get_selected_hierarchy()
    original_ch = classify_nqh(hierarchy)
    # add all H atoms
    his_ag = add_nqh_H_atoms(hierarchy)
    for atom in hierarchy.atoms(): break
    rg_resseq = atom.parent().parent().resseq
    chain_id = atom.parent().parent().parent().id
    buffer_selection = ''
    buffer = self.params.qi.qm_restraints[0].buffer
    buffer *= buffer
    t0=time.time()
    argstuples = []
    filenames = []
    for i, flipping_his in enumerate(generate_flipping(his_ag)):
      model = self.data_manager.get_model()
      model=model.deep_copy()
      hierarchy = model.get_hierarchy()
      for chain in hierarchy.chains():
        if chain.id!=chain_id: continue
        for rg in chain.residue_groups():
          if rg.resseq!=rg_resseq: continue
          for j, ag in enumerate(rg.atom_groups()): pass
          assert j==0, 'alt. loc. not supported'
          rg.remove_atom_group(ag)
          rg.insert_atom_group(0, flipping_his)
      self.params.output.prefix='iterate_%s_%02d' % (nq_or_h, i+1)
      # self.params.output.prefix='iterate_histidine_%02d' % (i+1)
      arg_log=null_out()
      if self.params.qi.verbose:
        arg_log=log
      #
      # need the same buffer for energy
      #
      if not buffer_selection:
        for rg in hierarchy.residue_groups():
          min_d2 = min_dist2(rg,ag)
          if min_d2[0]>=buffer: continue
          buffer_selection += ' (chain %s and resname %s and resid %s) or' % (
            rg.parent().id,
            rg.atom_groups()[0].resname,
            rg.resseq,
            )
        assert buffer_selection
        buffer_selection = buffer_selection[:-3]
      self.params.qi.qm_restraints[0].buffer_selection=buffer_selection
      #
      if self.params.qi.only_i is not None and self.params.qi.only_i!=i+1:
        continue
      filenames.append('%s_cluster_final_%s.pdb' % (self.params.output.prefix,
                                                   preamble))
      #
      if nproc==-1:
        print('  Running %s flip %d' % (nq_or_h, i+1), file=self.logger)
        res = update_restraints(model,
                                self.params,
                                never_write_restraints=True,
                                # never_run_strain=True,
                                log=arg_log)
        energies.append(energies)
        units=res.units
        rmsds.append(res.rmsds[0][1])
        print('    Energy : %s %s' % (energies[-1],units), file=self.logger)
        print('    Time   : %ds' % (time.time()-t0), file=self.logger)
        print('    RMSD   : %8.3f' % res.rmsds[0][1], file=self.logger)
      else:
        params = copy.deepcopy(self.params)
        argstuples.append(( model,
                            params,
                            None, # macro_cycle=None,
                            True, #never_write_restraints=False,
                            1, # nproc=1,
                            None, #null_out(),
                            ))
    results = run_serial_or_parallel(update_restraints, argstuples, nproc)
    for i, filename in enumerate(filenames):
      filename=os.path.join(qm_work_dir, filename)
      assert os.path.exists(filename), ' Output %s missing' % filename
      results[i].filename = filename
    return results

  def _process_energies(self, energies, units, resname='HIS', energy_adjustment=None, log=None):
    te=[]
    adjust = []
    if resname=='HIS': adjust=[0,3]
    for i, energy in enumerate(energies):
      energy=energy[1]
      if i in adjust:
        if units.lower() in ['kcal/mol']:
          if energy_adjustment is None:
            energy_adjustment=94.51
        elif units.lower() in ['hartree']:
          energy+=0.5
          assert 0
        elif units.lower() in ['ev']:
          if energy_adjustment is None:
            # energy_adjustment=13.61 # ???
            energy_adjustment=4.098 # 94.51 kcal/mol JM
        else:
          assert 0
        #
        # ADJUST!!!
        #
        energy-=energy_adjustment
      te.append(energy)
    outl=''
    if resname=='HIS':
      outl = '  Proton energy used : %0.4f %s' % (energy_adjustment, units)
    return te, outl

  def iterate_ASN(self, log=None):
    def classify_NQ(args): pass
    rc=self.iterate_NQH('ASN',
                        classify_NQ,
                        get_first_ag,
                        generate_flipping_NQ,
                        log=log)
    if rc is None: return None
    protonation = ['original', 'flipped']
    nproc = self.params.qi.nproc
    rc = self.process_flipped_jobs('ASN', rc, protonation=protonation, nproc=nproc, log=log)
    return rc

  def iterate_GLN(self, log=None):
    def classify_NQ(args): pass
    rc=self.iterate_NQH('GLN',
                        classify_NQ,
                        get_first_ag,
                        generate_flipping_NQ,
                        log=log)
    if rc is None: return None
    protonation = ['original', 'flipped']
    nproc = self.params.qi.nproc
    rc = self.process_flipped_jobs('GLN', rc, protonation=protonation, nproc=nproc, log=log)
    return rc

  def iterate_histidine(self, log=None):
    rc=self.iterate_NQH('HIS',
                        classify_histidine,
                        add_histidine_H_atoms,
                        generate_flipping_his,
                        log=log)
    if rc is None: return None
    protonation = [ 'HD1, HE2',
                    'HD1 only',
                    'HE2 only',
                    'HD1, HE2 flipped',
                    'HD1 only flipped',
                    'HE2 only flipped',
    ]
    nproc = self.params.qi.nproc
    energy_adjustment=self.params.qi.proton_energy_difference
    rc = self.process_flipped_jobs('HIS', rc, protonation=protonation, nproc=nproc, energy_adjustment=energy_adjustment, log=log)
    return rc

  def process_flipped_jobs(self, resname, rc, protonation=None, id_str=None, nproc=-1, energy_adjustment=None, log=None):
    energies = []
    units = None
    rmsds = []
    rotamers = []
    filenames = []
    for i, res in enumerate(rc):
      for selection, te in res.energies.items(): pass
      te=te[0]
      print('  Energy %d %s : %9.1f %s # ligand atoms : %d # cluster atoms : %d cluster charge :%3d' % (
        i+1,
        te[0],
        te[1],
        res.units,
        te[2],
        te[3],
        te[4],
        ), file=self.logger)
      energies.append(te)
      units=res.units
      rmsds.append(res.rmsds[0][1])
      filenames.append(res.filename)

    if protonation is None: protonation=filenames
    rc=get_hbonds_via_filenames(filenames,
                                resname,
                                nproc=nproc,
                                restraint_filenames=self.restraint_filenames)
    hbondss, pymols = rc
    selection = self.params.qi.qm_restraints[0].selection
    if resname not in ['radius']:
      rotamers=get_rotamers_via_filenames(filenames, selection, resname)
    else:
      rotamers=['None']*len(filenames)

    energies, outl = self._process_energies(energies, units, resname=resname, energy_adjustment=energy_adjustment, log=log)
    me=min(energies)
    cmd = '\n\n  phenix.start_coot'
    #
    if resname not in ['radius']:
      hierarchy = self.get_selected_hierarchy()
      original_ch = classify_histidine(hierarchy, resname=resname)
      print('\n\nEnergies in units of %s' % units, file=self.logger)
      print('%s\n' % outl, file=self.logger)
      print('  %i. %-20s : rotamer "%s"' % (
        0,
        original_ch[1],
        original_ch[0]),
      file=self.logger)
    #
    final_result=None
    close_result=[]
    outl = '  %i. %-20s : %12.3f %s ~> %10.2f kcal/mol. H-Bonds : %2d rmsd : %7.2f rotamer "%s"'
    # outl = '%i|%-20s|%7.5f|%s|%10.2f|%2d|%7.2f|%s'
    min_i=None
    for i, filename in enumerate(filenames):
      assert os.path.exists(filename), '"%s"' % filename
      if self.params.qi.run_directory:
        cmd += ' %s' % os.path.join(self.params.qi.run_directory, filename)
      else:
        cmd += ' %s' % filename
      #
      nci = hbondss[i].get_counts(filter_id_str=id_str, min_data_size=1)
      n=-1
      if nci:
        n=nci.n_filter
      energy = energies[i]
      #
      # convert to kcal/mol
      #
      if units.lower() in ['hartree']:
        de = (energy-me)*627.503
      elif units.lower() in ['kcal/mol', 'dirac']:
        de = (energy-me)
      elif units.lower() in ['ev']:
        de = (energy-me)*23.0605
      args = (
        i+1,
        protonation[i],
        energy,
        units,
        de,
        n,
        rmsds[i],
        rotamers[i],
        )
      print(outl % args, file=self.logger)
      update=False
      if de<1e-3:
        final_result = outl % args
        min_i=i
        tmp = protonation[i]
        tmp=tmp.replace('only','')
        # tmp=tmp.replace('flipped', '')
        if original_ch[1]:
          if tmp.strip()==original_ch[1].strip():
            final_result += ' SAME'
          elif tmp.replace('flipped', '').strip()==original_ch[1].strip():
            final_result += ' FLIPPED'
          else:
            final_result += ' DIFFERENT'
            update=True
        else:
          update=i
        j=i
      elif de<5.:
        close_result.append(outl % args)

    if units.lower()=='dirac':
      raise Sorry('MOPAC not installed! Please install or update to Python3.')

    import shutil
    nf = filenames[min_i].replace('%02d_cluster' % (min_i+1), 'best_cluster')
    shutil.copy(filenames[min_i], nf)

    cmd += '\n\n'
    print(cmd, file=self.logger)
    print(pymols, file=self.logger)
    print('!!! %s' % final_result, file=self.logger)
    if close_result:
      print('\nClose', file=self.logger)
      for cl in close_result:
        print('  >< %s' % cl, file=self.logger)
    return update, filenames[j]

  def step_thru_buffer_radii(self, id_str=None, log=None):
    from mmtbx.geometry_restraints.quantum_interface import get_preamble
    if len(self.params.qi.qm_restraints)<1:
      self.write_qmr_phil(step_buffer_radius=True)
      print('Restart command with PHIL file', file=self.logger)
      return
    nproc = self.params.qi.nproc
    #
    start, end, step = self.params.qi.step_buffer_radius.split(',')
    steps=[]
    for r in stepper(start, end, step): steps.append(r)
    steps.reverse()
    buffer = self.params.qi.qm_restraints[0].buffer
    t0=time.time()
    argstuples = []
    filenames = []
    for i, r in enumerate(steps):
      model = self.data_manager.get_model()
      model=model.deep_copy()
      self.params.output.prefix='iterate_radii'
      arg_log=null_out()
      if self.params.qi.verbose:
        arg_log=log
      #
      if self.params.qi.only_i is not None and self.params.qi.only_i!=i+1:
        continue
      params = copy.deepcopy(self.params)
      params.qi.qm_restraints[0].buffer=r
      preamble = get_preamble(None, 0, params.qi.qm_restraints[0])
      filenames.append('%s_cluster_final_%s.pdb' % (self.params.output.prefix,
                                                   preamble))
      #
      if nproc==-1:
        print('  Running radius %d' % (r), file=self.logger)
        res = update_restraints(model,
                                params,
                                never_write_restraints=True,
                                # never_run_strain=True,
                                log=arg_log)
        energies.append(res.energies)
        units=res.units
        rmsds.append(res.rmsds[0][1])
        print('    Energy : %s %s' % (energies[-1],units), file=self.logger)
        print('    Time   : %ds' % (time.time()-t0), file=self.logger)
        print('    RMSD   : %8.3f' % res.rmsds[0][1], file=self.logger)
      else:
        argstuples.append(( model,
                            params,
                            None, # macro_cycle=None,
                            False, #never_write_restraints=False,
                            1, # nproc=1,
                            None, #null_out(),
                            ))
    results = run_serial_or_parallel(update_restraints, argstuples, nproc)
    for i, filename in enumerate(filenames):
      assert os.path.exists(filename), ' Output %s missing' % filename
      results[i].filename = filename
    if results is None: return
    self.process_flipped_jobs('radius', results, id_str=id_str, nproc=nproc, log=log)
    return results

  def run_qmr(self, format, return_minimised=False, log=None):
    from mmtbx.refinement.energy_monitor import digest_return_energy_object
    model = self.data_manager.get_model()
    qmr = self.params.qi.qm_restraints[0]
    checks = 'starting_strain starting_energy starting_bound starting_binding'
    energies = None
    rcs=[]
    if any(item in checks for item in qmr.calculate):
      rc = run_energies(
        model,
        self.params,
        macro_cycle=1,
        pre_refinement=True,
        nproc=self.params.qi.nproc,
        log=log,
        )
      energies = digest_return_energy_object(rc, 1, energy_only=True)
      outl = energies.as_string()
      print(outl, file=self.logger)
      rcs.append(rc)
    if any(item in 'gradients' for item in qmr.calculate):
      rc = run_gradients(
        model,
        self.params,
        macro_cycle=1,
        pre_refinement=True,
        nproc=self.params.qi.nproc,
        log=log,
        )
      assert 0
    #
    # minimise ligands geometry
    #
    if return_minimised:
      for qmr in self.params.qi.qm_restraints:
        if 'pdb_final_buffer' in qmr.write_files:
          qmr.write_files.remove('pdb_final_buffer')
        qmr.write_files.append('pdb_final_core')
    rc = update_restraints( model,
                            self.params,
                            log=log,
                            )
    rcs.append(rc)
    if energies is None:
      energies = digest_return_energy_object(rc, 1, energy_only=False)
    else:
      digest_return_energy_object(rc, 1, False, energies)
    outl = energies.as_string()
    print(outl, file=self.logger)
    #
    cannot_run_final_energies=False
    if len(rc.final_pdbs)>1:
      cannot_run_final_energies=True
    checks = 'final_strain final_energy final_bound'
    if rc.final_pdbs[0]:
      fn = os.path.join('qm_work_dir', rc.final_pdbs[0][0])
      self.data_manager.process_model_file(fn)
      model = self.data_manager.get_model(fn)
    if any(item in checks for item in qmr.calculate):
      if cannot_run_final_energies:
        print('Cannot run final energies. Restart with only one ligand selection.',
              file=self.logger)
      else:
        rc = run_energies(
          model,
          self.params,
          macro_cycle=99,
          pre_refinement=False,
          nproc=self.params.qi.nproc,
          log=log,
          )
        energies = digest_return_energy_object(rc, 99, True, energies)
        outl = energies.as_string()
        print(outl, file=self.logger)
        rcs.append(rc)
    return rcs

  def get_single_qm_restraints_scope(self, selection):
    qi_phil_string = get_qm_restraints_scope()
    qi_phil_string = qi_phil_string.replace(' selection = None',
                                            ' selection = "%s"' % selection)
    qi_phil_string = qi_phil_string.replace('read_output_to_skip_opt_if_available = False',
                                            'read_output_to_skip_opt_if_available = True')
    qi_phil_string = qi_phil_string.replace('capping_groups = False',
                                            'capping_groups = True')
    return qi_phil_string

  def get_single_qm_gradients_scope(self, selection):
    qi_phil_string = get_qm_gradients_scope()
    qi_phil_string = qi_phil_string.replace(' selection = None',
                                            ' selection = "%s"' % selection)
    qi_phil_string = qi_phil_string.replace('read_output_to_skip_opt_if_available = False',
                                            'read_output_to_skip_opt_if_available = True')
    qi_phil_string = qi_phil_string.replace('capping_groups = False',
                                            'capping_groups = True')
    return qi_phil_string


  def set_all_calculate_to_true(self, qi_phil_string):
    outl = ''
    for line in qi_phil_string.splitlines():
      if line.find(' calculate =')>-1:
        tmp=line.split()
        line=''
        for i, t in enumerate(tmp):
          if i>1 and t.find('*')==-1:
            line += ' *%s' % tmp[i]
          else:
            line += ' %s' % tmp[i]
      outl += '%s\n' % line
    return outl

  def set_all_write_to_true(self, qi_phil_string):
    outl = ''
    # write_files = *restraints pdb_core pdb_buffer pdb_final_core pdb_final_buffer
    for line in qi_phil_string.splitlines():
      if line.find(' write_')>-1:
        tmp=line.split()
        line = '  '
        for t in tmp:
          if t.startswith('pdb'): t='*%s'%t
          line+='%s ' % t
      outl += '%s\n' % line
    return outl

  def set_one_write_to_true(self, qi_phil_string, attr):
    return qi_phil_string.replace(' %s' % attr, ' *%s' % attr)

  def write_qmg_phil(self):
    qi_phil_string = ''
    for i in range(len(self.params.qi.selection)):
      qi_phil_string += self.get_single_qm_gradients_scope(self.params.qi.selection[i])
    qi_phil_string = qi_phil_string.replace('qm_gradients',
                                            'refinement.qi.qm_gradients',
                                            1)
    qi_phil_string = qi_phil_string.replace(' *mopac',
                                            ' ')
    qi_phil_string = qi_phil_string.replace(' mopac',
                                            ' ')
    qi_phil_string = qi_phil_string.replace(' orca',
                                            ' *orca')
    qi_phil_string = qi_phil_string.replace(' method = Auto',
                                            ' method = PM3')
    qi_phil_string = qi_phil_string.replace(' read_output_to_skip_opt_if_available = True',
                                            ' read_output_to_skip_opt_if_available = False')

    qi_phil = iotbx.phil.parse(qi_phil_string,
                             # process_includes=True,
                             )

    pf = '%s_%s_grad.phil' % (
      self.data_manager.get_default_model_name().replace('.pdb',''),
      safe_filename(self.params.qi.selection[0]),
      )
    print('  Writing QMG phil scope to %s' % pf, file=self.logger)
    f=open(pf, 'w')
    for line in qi_phil_string.splitlines():
      if line.strip().startswith('.'): continue
      print('%s' % line, file=self.logger)
      f.write('%s\n' % line)
    del f
    return pf

  def write_qmr_phil(self,
                     iterate_NQH=False,
                     iterate_metals=False,
                     step_buffer_radius=False,
                     output_format=None,
                     include_nearest_neighbours=False,
                     log=None):
    qi_phil_string = ''
    for i in range(len(self.params.qi.selection)):
      qi_phil_string += self.get_single_qm_restraints_scope(self.params.qi.selection[i])
    # qi_phil_string = self.get_single_qm_restraints_scope(self.params.qi.selection[0])
    # qi_phil_string = self.set_all_calculate_to_true(qi_phil_string)
    # qi_phil_string = self.set_all_write_to_true(qi_phil_string)
    qi_phil = iotbx.phil.parse(qi_phil_string,
                             # process_includes=True,
                             )
    qi_phil_string = qi_phil_string.replace(' pdb_final_buffer',
                                            ' *pdb_final_buffer')
    # qi_phil_string = qi_phil_string.replace(' pdb_buffer',
    #                                         ' *pdb_buffer')
    # qi_phil.show()

    qi_phil_string = qi_phil_string.replace('qm_restraints',
                                            'refinement.qi.qm_restraints',
                                            1)
    if step_buffer_radius:
      # start, end, step = step_buffer_radius.split(',')
      qi_phil_string=qi_phil_string.replace('refinement.qi.qm_restraints',
                                            'qm_restraints')
      qi_phil_string=qi_phil_string.replace('ignore_x_h_distance_protein = False',
                                            'ignore_x_h_distance_protein = True')
      start=3.5
      end=3.5
      step=1
      tmp = 'qi {\n'
      for r in stepper(start, end, step):
        tmp+=qi_phil_string.replace('buffer = 3.5', 'buffer = %s' % r)
      tmp += '\n}\n'
      qi_phil_string = tmp

    if iterate_NQH:
      qi_phil_string = qi_phil_string.replace('refinement.', '')
      qi_phil_string = qi_phil_string.replace('ignore_x_h_distance_protein = False',
                                              'ignore_x_h_distance_protein = True')
      qi_phil_string = qi_phil_string.replace(
        'protein_optimisation_freeze = *all None main_chain main_chain_to_beta main_chain_to_delta torsions',
        'protein_optimisation_freeze = all None main_chain main_chain_to_beta *main_chain_to_delta *torsions')
      qi_phil_string = qi_phil_string.replace(
        'solvent_model = None', 'solvent_model = EPS=78.4',
        )
      # qi_phil_string = qi_phil_string.replace('buffer = 3.5', 'buffer = 4.')

    if iterate_metals:
      qi_phil_string = qi_phil_string.replace('refinement.', '')
      qi_phil_string = qi_phil_string.replace('include_nearest_neighbours_in_optimisation = False',
                                              'include_nearest_neighbours_in_optimisation = True')

    if include_nearest_neighbours:
      qi_phil_string = qi_phil_string.replace('include_nearest_neighbours_in_optimisation = False',
                                              'include_nearest_neighbours_in_optimisation = True')
      qi_phil_string = qi_phil_string.replace('include_inter_residue_restraints = False',
                                              'include_inter_residue_restraints = True')

    if output_format=='qi':
      qi_phil_string = qi_phil_string.replace('refinement.qi', 'qi')

    # qi_phil_string = qi_phil_string.replace('nproc = 1', 'nproc = 6')


    print(self.params.qi.selection)
    pf = '%s_%s.phil' % (
      self.data_manager.get_default_model_name().replace('.pdb',''),
      safe_filename(self.params.qi.selection[0]),
      )
    print('  Writing QMR phil scope to %s' % pf, file=self.logger)
    f=open(pf, 'w')
    for line in qi_phil_string.splitlines():
      if line.strip().startswith('.'): continue
      print('%s' % line, file=self.logger)
      f.write('%s\n' % line)
    del f
    return pf

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results
