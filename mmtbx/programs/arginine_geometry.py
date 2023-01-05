from __future__ import absolute_import, division, print_function
import os
from libtbx.program_template import ProgramTemplate
from cctbx import geometry_restraints

from cctbx.array_family import flex

limits = {
  ('CD', 'NE', 'CZ', 'NH1') : 5,
  ('CD', 'NE', 'CZ', 'NH2') : 5,
  ('NE', 'CZ', 'NH1', 'NH2'): 2,
  ('NE', 'CZ', 'NH1', 'HH11'): 1,
  ('NE', 'CZ', 'NH1', 'HH12'): 1,
  ('NE', 'CZ', 'NH2', 'HH21'): 1,
  ('NE', 'CZ', 'NH2', 'HH22'): 1,
  }

def get_atoms_key(atoms):
  if atoms[3].find('HH')>-1: return ('NE', 'CZ', 'NHn', 'HHn?')
  if atoms[0]=='CD': return ('CD', 'NE', 'CZ', 'NH?')
  return atoms

class arginine_deviations(dict):
  def __repr__(self):
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
    def positive(v):
      if v<-90: v+=360
      if v>90: v-=180
      return v
    self.data = {}
    for residue, chis in self.items():
      for atoms, chi in chis.items():
        key = get_atoms_key(atoms)
        self.data.setdefault(key, flex.float())
        self.data[key].append(positive(chi))
    self.atoms = {}
    for atoms, chis in self.data.items():
      self.atoms[atoms]=[flex.mean(chis), flex.min(chis), flex.max(chis)]

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
  }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    # print('Validating inputs', file=self.logger)
    pass

  # ---------------------------------------------------------------------------
  def run(self, log=None):
    model = self.data_manager.get_model()
    hierarchy = model.get_hierarchy()
    self.results = arginine_deviations()
    for residue_group in hierarchy.residue_groups():
      if len(residue_group.atom_groups())>1: continue
      atom_group = residue_group.atom_groups()[0]
      if atom_group.resname!='ARG': continue
      current = self.results.setdefault(atom_group.id_str(), {})
      for torsion_atom_names in limits:
        torsion_xyzs = []
        for atom in torsion_atom_names:
          ta = atom_group.get_atom(atom)
          torsion_xyzs.append(ta.xyz)
        v = geometry_restraints.dihedral(sites=torsion_xyzs, angle_ideal=0, weight=1).angle_model
        current[tuple(torsion_atom_names)]=v
    self.results.process()
    print('\n%s' % self.results, file=log)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results
