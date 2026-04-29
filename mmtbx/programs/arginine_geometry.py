"""Check geometry of arginines"""
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate
from cctbx import geometry_restraints

from cctbx.array_family import flex
from libtbx.utils import null_out

limits = {
  ('CD', 'NE', 'CZ', 'NH1')  : 10,
  ('CD', 'NE', 'CZ', 'NH2')  : 10,
  ('NE', 'CZ', 'NH1', 'NH2') : 2,
  ('NE', 'CZ', 'NH1', 'HH11'): 1,
  ('NE', 'CZ', 'NH1', 'HH12'): 1,
  ('NE', 'CZ', 'NH2', 'HH21'): 1,
  ('NE', 'CZ', 'NH2', 'HH22'): 1,
  #
  ('CD', 'NE', 'CZ', 'NH?')  : 10,
  ('NE', 'CZ', 'NHn', 'HHn?'): 1,
  }

def get_torsion_deviations(v):
  if v<-90: v+=360
  if v>90: v-=180
  return v

def get_atoms_key(atoms):
  if atoms[3].find('HH')>-1: return ('NE', 'CZ', 'NHn', 'HHn?')
  if atoms[0]=='CD': return ('CD', 'NE', 'CZ', 'NH?')
  return atoms

def format_atoms(t):
  outl = ''
  for atom in t:
    outl += '%4s -' % atom
  return outl[:-1]

def format_values(v):
  outl = ''
  for attr, value in zip(['mean', 'min', 'max'], v):
    outl += ' %s: %4.1f' % (attr, value)
  return outl

def format_annotation(values):
  limit = limits.get(atoms, 1e9)
  ann = ''
  if abs(values[1])>limit: ann+='<'
  if values[2]>limit: ann+='>'
  return ann

class arginine_deviations(dict):
  def __repr__(self):
    outl = 'ARG devs\n'
    outl += '  Number of ARG : %s\n' % len(self)
    current = getattr(self, 'atoms',{})
    if current: outl += '  Stats\n'
    for atoms, values in current.items():
      outl += '  %s ~> %s %s\n' % (
        format_atoms(atoms),
        format_values(values),
        format_annotation(values),
        )
    return outl

  def process(self):
    self.data = {}
    for residue, chis in self.items():
      for atoms, chi in chis.items():
        key = get_atoms_key(atoms)
        self.data.setdefault(key, flex.float())
        self.data[key].append(get_torsion_deviations(chi))
    self.atoms = {}
    for atoms, chis in self.data.items():
      self.atoms[atoms]=[flex.mean(chis), flex.min(chis), flex.max(chis)]

class plane_torsion_deviations(arginine_deviations):
  def __init__(self, print_torsion_limit=1., print_torsion_number=10):
    self.print_torsion_limit=print_torsion_limit
    self.print_torsion_number=print_torsion_number

  def __repr__(self):
    tmp = dict(sorted(self.items(), reverse=True, key=lambda item: abs(item[1][0])))
    outl = 'Plane Torsion Deviations\n'
    outl += '  Number : %s\n' % len(self)
    i=0
    for key, item in tmp.items():
      if abs(item[0])<self.print_torsion_limit: continue
      i+=1
      outl += '    %s : %5.1f %s\n' % (key, item[0], '*' if item[1] else '')
      if i>= self.print_torsion_number: break
    return outl

  def is_possible_reduce(self):
    rc=0
    for key, item in self.items():
      if not item[1]: continue
      if abs(item[0])>1.:
        rc+=1
    return rc

class Program(ProgramTemplate):

  description = '''
mmtbx.arginine_geometry:

Usage examples:
  mmtbx.arginine_geometry model.pdb
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """
  arginine {
    selection = None
      .type = atom_selection
      .help = what to select
      .multiple = True
    exclude_hydrogen = False
      .type = bool
    print_torsion_limit = 1.
      .type = float
    print_torsion_number = 10
      .type = int
  }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    # print('Validating inputs', file=self.logger)
    pass

  def arginine_simple(self):
    model = self.data_manager.get_model()
    hierarchy = model.get_hierarchy()
    for residue_group in hierarchy.residue_groups():
      if len(residue_group.atom_groups())>1: continue
      atom_group = residue_group.atom_groups()[0]
      if atom_group.resname!='ARG': continue
      current = self.results.setdefault(atom_group.id_str(), {})
      for torsion_atom_names in limits:
        torsion_xyzs = []
        for atom in torsion_atom_names:
          ta = atom_group.get_atom(atom)
          if ta is None:
            torsion_xyzs=[]
            break
          torsion_xyzs.append(ta.xyz)
        if not torsion_xyzs: continue
        v = geometry_restraints.dihedral(sites=torsion_xyzs, angle_ideal=0, weight=1).angle_model
        current[tuple(torsion_atom_names)]=v
    self.results.process()

  def torsions_in_planes(self, exclude_hydrogen=False):
    model = self.data_manager.get_model()
    model.set_log(null_out())
    model.process(make_restraints=True)
    grm = model.get_restraints_manager()
    atoms = model.get_hierarchy().atoms()
    plane_i_seqs = []
    remove = []
    for i, plane in enumerate(grm.geometry.planarity_proxies):
      plane_i_seqs.append(set(plane.i_seqs))
      torsion_xyzs = []
      torsion_quotes = []
      for i_seq in plane.i_seqs:
        torsion_xyzs.append(atoms[i_seq].xyz)
        torsion_quotes.append(atoms[i_seq].quote()[8:])
        if len(torsion_xyzs)==4:
          v = geometry_restraints.dihedral(sites=torsion_xyzs, angle_ideal=0, weight=1).angle_model
          v=get_torsion_deviations(v)
          if abs(v)>1: break
          del torsion_xyzs[0]
          del torsion_quotes[0]
      else:
        remove.append(i)
    if remove:
      remove.reverse()
      for r in remove:
        del plane_i_seqs[r]
    data = {}
    for torsion in grm.geometry.dihedral_proxies:
      for plane in plane_i_seqs:
        if len(plane.intersection(set(torsion.i_seqs)))==4:
          torsion_xyzs = []
          torsion_strs = ''
          h_atom = False
          for i_seq in torsion.i_seqs:
            atom = atoms[i_seq]
            if atom.element_is_hydrogen(): h_atom=True
            torsion_xyzs.append(atom.xyz)
            torsion_strs += '%s - ' % (atom.id_str()[4:])
          if exclude_hydrogen and h_atom: continue
          v = geometry_restraints.dihedral(sites=torsion_xyzs, angle_ideal=0, weight=1).angle_model
          v = get_torsion_deviations(v)
          self.results[torsion_strs[:-2]]=[v, h_atom]

  # ---------------------------------------------------------------------------
  def run(self, log=None):
    # self.results = arginine_deviations()
    # self.arginine_simple()
    self.results = plane_torsion_deviations(self.params.arginine.print_torsion_limit,
                                            self.params.arginine.print_torsion_number)
    self.torsions_in_planes(self.params.arginine.exclude_hydrogen)
    print('\n%s' % self.results, file=log)
    reduce_H_atoms = self.results.is_possible_reduce()
    print('  Has reduce-like H atoms : %s\n' % reduce_H_atoms if reduce_H_atoms else '', file=log)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results
