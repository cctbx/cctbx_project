# LIBTBX_SET_DISPATCHER_NAME phenix.development.qi
from __future__ import absolute_import, division, print_function
import os
import time
from libtbx.program_template import ProgramTemplate

from mmtbx.monomer_library.linking_setup import ad_hoc_single_metal_residue_element_types

from mmtbx.geometry_restraints.quantum_restraints_manager import run_energies
from mmtbx.geometry_restraints.quantum_restraints_manager import update_restraints
from mmtbx.geometry_restraints.quantum_restraints_manager import min_dist2
from mmtbx.geometry_restraints.quantum_interface import get_qm_restraints_scope
from mmtbx.geometry_restraints.quantum_interface import classify_histidine

import iotbx.pdb
import iotbx.phil
from libtbx.utils import Sorry
from libtbx.utils import null_out

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
  if count not in [2]:
    raise Sorry('incorrect protonation of %s' % ag.id_str())

def generate_flipping_his(ag,
                          return_hierarchy=False,
                          include_unprotonated=False,
                          chain_id=None,
                          resseq=None):
  assert_histidine_double_protonated(ag)
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
    if atom_group.resname in ['HIS']: pass
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
    format = *phenix_refine quantum_interface
      .type = choice
    write_qmr_phil = False
      .type = bool
    run_qmr = False
      .type = bool
    randomise_selection = None
      .type = float
    run_directory = None
      .type = str
    iterate_histidine = False
      .type = bool
    only_i = None
      .type = int
    iterate_metals = None
      .type = str
    each_amino_acid = False
      .type = bool
    nproc = 1
      .type = int
    verbose = False
      .type = bool
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
    if self.params.qi.randomise_selection and self.params.qi.randomise_selection>0.5:
      raise Sorry('Random select value %s is too large' % self.params.qi.randomise_selection)

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
        len(self.params.qi.qm_restraints)==0 and
        not self.params.qi.each_amino_acid):
      rc = get_selection_from_user(model.get_hierarchy())
      if rc.find('resname HIS')>-1:
        self.params.qi.iterate_histidine=rc
        self.params.qi.format='quantum_interface'
      self.params.qi.selection = [rc]
    #
    # validate selection
    #
    selection=None
    if len(self.params.qi.qm_restraints)!=0:
      selection = self.params.qi.qm_restraints[0].selection
    elif self.params.qi.selection:
      selection=self.params.qi.selection[0]
    if selection:
      selection_array = model.selection(selection)
      selected_model = model.select(selection_array)
      print('Selected model  %s' % selected_model, file=log)
      self.data_manager.add_model('ligand', selected_model)

    if self.params.qi.each_amino_acid:
      hierarchy = model.get_hierarchy()
      outl = ''
      for rg in hierarchy.residue_groups():
        if len(rg.atom_groups())!=1: continue
        gc = get_class(rg.atom_groups()[0].resname)
        if gc not in ['common_amino_acid']: continue
        selection = 'chain %s and resid %s' % (rg.parent().id, rg.resseq.strip())
        qi_phil_string = self.get_single_qm_restraints_scope(selection)
        qi_phil_string = self.set_all_write_to_true(qi_phil_string)
        print('  writing phil for %s %s' % (rg.id_str(), rg.atom_groups()[0].resname))
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

    if self.params.qi.write_qmr_phil:
      pf = self.write_qmr_phil(iterate_histidine=self.params.qi.iterate_histidine,
                               iterate_metals=self.params.qi.iterate_metals,
                               output_format=self.params.qi.format,
                               )
      ih = ''
      # if self.params.qi.iterate_histidine:
      #   ih = 'iterate_histidine="%s"' % self.params.qi.iterate_histidine
      #   ih = ih.replace('  ',' ')
      print('''

      mmtbx.quantum_interface %s run_qmr=True %s %s
      ''' % (self.data_manager.get_default_model_name(),
             ih,
             pf))
      return

    if self.params.qi.run_directory:
      if not os.path.exists(self.params.qi.run_directory):
        os.mkdir(self.params.qi.run_directory)
      os.chdir(self.params.qi.run_directory)

    if ( self.params.qi.run_qmr and
         not (self.params.qi.iterate_histidine or
              self.params.qi.iterate_metals)):
      self.params.qi.qm_restraints.selection=self.params.qi.selection
      self.run_qmr(self.params.qi.format)

    if self.params.qi.iterate_histidine:
      assert self.params.qi.randomise_selection==None
      self.iterate_histidine()

    if self.params.qi.iterate_metals:
      self.iterate_metals()

  def iterate_metals(self, log=None):
    if len(self.params.qi.qm_restraints)<1:
      self.write_qmr_phil(iterate_metals=True)
      print('Restart command with PHIL file')
      return
    def generate_metals(s):
      s=s.replace(' ',',')
      for t in s.split(','):
        yield t.upper()
    selection = self.params.qi.qm_restraints[0].selection
    nproc = self.params.qi.nproc
    model = self.data_manager.get_model()
    selection_array = model.selection(selection)
    selected_model = model.select(selection_array)
    hierarchy = selected_model.get_hierarchy()
    for atom in hierarchy.atoms(): break
    rg_resseq = atom.parent().parent().resseq
    chain_id = atom.parent().parent().parent().id
    energies = []
    argstuples = []
    # logs = []
    # buffer_selection = ''
    # buffer = self.params.qi.qm_restraints[0].buffer
    # buffer *= buffer
    # for i, flipping_his in enumerate(generate_flipping_his(his_ag)):
    for element in generate_metals(self.params.qi.iterate_metals):
      print('Substituting element : %s' % element)
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
          selection = selection.replace(' %s' % ag.resname.upper().strip(),
                                        ' %s' % element)
          self.params.qi.qm_restraints[0].selection = selection
          rg.remove_atom_group(ag)
          rg.insert_atom_group(0, eag)
      self.params.output.prefix='iterate_metals_%s' % (element)
      arg_log=null_out()
      if self.params.qi.verbose:
        arg_log=log
      if nproc==1:
        print('  Running metal swap %s' % (element), file=log)
        res = update_restraints(model,
                                self.params,
                                never_write_restraints=True,
                                # never_run_strain=True,
                                log=arg_log)
        energies.append(res[2][0])
        units=res[3]
        # print('    Energy %s %s' % (energies[-1],units), file=log)
      else:
        import copy
        params = copy.deepcopy(self.params)
        argstuples.append(( model,
                            params,
                            None, # macro_cycle=None,
                            True, #never_write_restraints=False,
                            1, # nproc=1,
                            null_out(),
                            ))
    from libtbx import easy_mp
    from mmtbx.geometry_restraints.quantum_interface import get_preamble
    preamble = get_preamble(None, 0, self.params.qi.qm_restraints[0])
    if argstuples:
      print('  Running %d jobs in %d procs' % (len(argstuples), nproc), file=log)
    i=0
    for args, res, err_str in easy_mp.multi_core_run( update_restraints,
                                                      argstuples,
                                                      nproc,
                                                      keep_input_order=True):
      assert not err_str, '\n\nDebug in serial :\n%s' % err_str
      print('  Running metal substitution %d' % (i+1), file=log)
      energies.append(res[0][0])
      units = res[1]
      print('    Energy %s %s' % (energies[-1],units), file=log)
      i+=1

  def iterate_histidine(self, log=None):
    if len(self.params.qi.qm_restraints)<1:
      self.write_qmr_phil(iterate_histidine=True)
      print('Restart command with PHIL file')
      return
    selection = self.params.qi.qm_restraints[0].selection
    nproc = self.params.qi.nproc
    model = self.data_manager.get_model()
    selection_array = model.selection(selection)
    selected_model = model.select(selection_array)
    hierarchy = selected_model.get_hierarchy()
    original_ch = classify_histidine(hierarchy)
    # add all H atoms
    his_ag = add_histidine_H_atoms(hierarchy)
    for atom in hierarchy.atoms(): break
    rg_resseq = atom.parent().parent().resseq
    chain_id = atom.parent().parent().parent().id
    energies = []
    rmsds = []
    argstuples = []
    # logs = []
    buffer_selection = ''
    buffer = self.params.qi.qm_restraints[0].buffer
    buffer *= buffer
    t0=time.time()
    for i, flipping_his in enumerate(generate_flipping_his(his_ag)):
      model = self.data_manager.get_model()
      if nproc>1: model=model.deep_copy()
      hierarchy = model.get_hierarchy()
      for chain in hierarchy.chains():
        if chain.id!=chain_id: continue
        for rg in chain.residue_groups():
          if rg.resseq!=rg_resseq: continue
          for j, ag in enumerate(rg.atom_groups()): pass
          assert j==0
          rg.remove_atom_group(ag)
          rg.insert_atom_group(0, flipping_his)
      self.params.output.prefix='iterate_histidine_%02d' % (i+1)
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
      #
      if nproc==1:
        print('  Running HIS flip %d' % (i+1), file=log)
        res = update_restraints(model,
                                self.params,
                                never_write_restraints=True,
                                # never_run_strain=True,
                                log=arg_log)
        energies.append(res.energies[0])
        units=res.units
        rmsds.append(res.rmsds[0][1])
        if 0:
          print('    Energy : %s %s' % (energies[-1],units), file=log)
          print('    Time   : %ds' % (time.time()-t0), file=log)
          print('    RMSD   : %8.3f' % res.rmsds[0][1], file=log)
      else:
        import copy
        params = copy.deepcopy(self.params)
        argstuples.append(( model,
                            params,
                            None, # macro_cycle=None,
                            True, #never_write_restraints=False,
                            1, # nproc=1,
                            null_out(),
                            ))

    from libtbx import easy_mp
    from mmtbx.geometry_restraints.quantum_interface import get_preamble
    preamble = get_preamble(None, 0, self.params.qi.qm_restraints[0])
    if nproc>1:
      print('  Running %d jobs in %d procs' % (len(argstuples), nproc), file=log)
    i=0
    for args, res, err_str in easy_mp.multi_core_run( update_restraints,
                                                      argstuples,
                                                      nproc,
                                                      keep_input_order=True):
      assert not err_str, '\n\nDebug in serial :\n%s' % err_str
      print('  Running HIS flip %d' % (i+1), file=log)
      # energies.append(res[0][0])
      # units = res[1]
      # print('    Energy %s %s' % (energies[-1],units), file=log)
      assert 0
      i+=1

    protonation = [ 'ND1,NE2',
                    'ND1 only',
                    'NE2 only',
                    'ND1,NE2 flipped',
                    'ND1 only flipped',
                    'NE2 only flipped',
    ]
    cmd = '\n\n  phenix.start_coot'
    te=[]
    for i, (pro, energy) in enumerate(zip(protonation, energies)):
      if i in [0,3]:
        if units.lower() in ['kcal/mol']:
          # energy-=247.80642
          energy-=156.9
        elif units.lower() in ['hartree']:
          energy+=0.5
      te.append(energy)
    me=min(te)
    #
    # parallel HBond
    #
    argstuples = []
    for i, (pro, energy) in enumerate(zip(protonation, energies)):
      prefix='iterate_histidine_%02d' % (i+1)
      filename = '%s_cluster_final_%s.pdb' % (prefix, preamble)
      assert os.path.exists(filename), '"%s"' % filename
      argstuples.append([[filename,
                         '--quiet',
                         'output_pymol_file=True',
                         'output_restraint_file=False',
                         'output_skew_kurtosis_plot=False',
                         'prefix=%s_%s' % (prefix, preamble),
                         ]])
    def run_hbond(args):
      from iotbx.cli_parser import run_program
      from mmtbx.programs.hbond import Program
      hbonds = run_program(program_class=Program,
                           args=tuple(args),
                           )
      return hbonds
    # print('  Running %d jobs in %d procs' % (len(argstuples), nproc), file=log)
    i=0
    hbondss=[]
    rotamers=[]
    pymols = ''
    for args, res, err_str in easy_mp.multi_core_run( run_hbond,
                                                      argstuples,
                                                      max(nproc,6),
                                                      keep_input_order=True):
      assert not err_str, '\n\nDebug in serial :\n%s' % err_str
      hbondss.append(res)
      prefix='iterate_histidine_%02d' % (i+1)
      pf = '%s_%s.pml' % (prefix, preamble)
      f=open(pf, 'a')
      f.write('\n')
      # f.write('zoom resn HIS\n')
      f.write('show sticks, resn HIS\n')
      del f
      pymols += '  phenix.pymol %s &\n' % pf
      from iotbx import pdb
      hierarchy = pdb.input(pf.replace('.pml', '.pdb')).construct_hierarchy()
      rc = classify_histidine(hierarchy)
      rotamers.append(rc[0])
      i+=1
    #
    results = {}
    print('\n\nEnergies in units of %s\n' % units, file=log)
    print('  %i. %-20s : rotamer "%s"' % (
      0,
      original_ch[1],
      original_ch[0])
    )
    for i, (pro, energy) in enumerate(zip(protonation, energies)):
      energy=te[i]
      prefix='iterate_histidine_%02d' % (i+1)
      filename = '%s_cluster_final_%s.pdb' % (prefix, preamble)
      assert os.path.exists(filename), '"%s"' % filename
      if self.params.qi.run_directory:
        cmd += ' %s' % os.path.join(self.params.qi.run_directory, filename)
      else:
        cmd += ' %s' % filename
      #
      n = hbondss[i].get_counts(min_data_size=1).n
      #
      if units.lower() in ['hartree']:
        de = (energy-me)*627.503
        print('  %i. %-20s : %7.5f %s ~> %10.2f kcal/mol. H-Bonds : %2d rmsd : %7.2f rotamer "%s"' % (
          i+1,
          pro,
          energy,
          units,
          de,
          n,
          rmsds[i],
          rotamers[i],
          ), file=log)
      elif units.lower() in ['kcal/mol']:
        de = (energy-me)
        print('  %i. %-20s : %7.2f %s ~> %7.2f kcal/mol. H-Bonds : %2d rmsd : %7.2f rotamer "%s"' % (
          i+1,
          pro,
          energy,
          units,
          de,
          n,
          rmsds[i],
          rotamers[i],
          ), file=log)
      results.setdefault(pro, {})
      results[pro]['delta E'] = de
      results[pro]['H bonds'] = n
      results[pro]['rmsd'] = rmsds[i]

    # d = dict(sorted(results.items(), key=lambda item: item[1]['delta E']))
    # for i, (key, item) in enumerate(d.items()):
    #   print('%-20s %s' % (key,item))

    cmd += '\n\n'
    print(cmd)
    print(pymols)

  def run_qmr(self, format, log=None):
    model = self.data_manager.get_model()
    qmr = self.params.qi.qm_restraints[0]
    if qmr.calculate_starting_strain:
      rc = run_energies(
        model,
        self.params,
        # macro_cycle=self.macro_cycle,
        pre_refinement=True,
        # nproc=self.params.main.nproc,
        log=log,
        )
      print('starting strain',rc)
    #
    # minimise ligands geometry
    #
    rc = update_restraints( model,
                            self.params,
                            log=log,
                            )

  def get_single_qm_restraints_scope(self, selection):
    qi_phil_string = get_qm_restraints_scope()
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
      if line.find(' calculate_')>-1:
        tmp=line.split()
        line = '  %s = True' % tmp[0]
      outl += '%s\n' % line
    return outl

  def set_all_write_to_true(self, qi_phil_string):
    outl = ''
    for line in qi_phil_string.splitlines():
      if line.find(' write_')>-1:
        tmp=line.split()
        line = '  %s = True' % tmp[0]
      outl += '%s\n' % line
    return outl

  def write_qmr_phil(self,
                     iterate_histidine=False,
                     iterate_metals=False,
                     output_format=None,
                     log=None):
    qi_phil_string = self.get_single_qm_restraints_scope(self.params.qi.selection[0])
    qi_phil_string = self.set_all_calculate_to_true(qi_phil_string)
    # qi_phil_string = self.set_all_write_to_true(qi_phil_string)
    qi_phil = iotbx.phil.parse(qi_phil_string,
                             # process_includes=True,
                             )
    qi_phil_string = qi_phil_string.replace('write_final_pdb_buffer = False',
                                            'write_final_pdb_buffer = True')
    # qi_phil.show()

    qi_phil_string = qi_phil_string.replace('qm_restraints',
                                            'refinement.qi.qm_restraints',
                                            1)
    if iterate_histidine:
      qi_phil_string = qi_phil_string.replace('refinement.', '')
      qi_phil_string = qi_phil_string.replace('ignore_x_h_distance_protein = False',
                                              'ignore_x_h_distance_protein = True')
      # qi_phil_string = qi_phil_string.replace('write_restraints = True',
      #                                         'write_restraints = False')
      qi_phil_string = qi_phil_string.replace('exclude_protein_main_chain_to_delta_from_optimisation = False',
                                              'exclude_protein_main_chain_to_delta_from_optimisation = True')
      # qi_phil_string = qi_phil_string.replace('exclude_protein_main_chain_from_optimisation = False',
      #                                         'exclude_protein_main_chain_from_optimisation = True')

    if iterate_metals:
      qi_phil_string = qi_phil_string.replace('include_nearest_neighbours_in_optimisation = False',
                                              'include_nearest_neighbours_in_optimisation = True')

    if output_format=='quantum_interface':
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
      print('%s' % line)
      f.write('%s\n' % line)
    del f
    return pf

  # ---------------------------------------------------------------------------
  def get_results(self):
    return None
