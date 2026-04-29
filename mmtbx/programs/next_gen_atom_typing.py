"""Report on Ramachandran values for a model"""
from __future__ import absolute_import, division, print_function

import os
import iotbx.phil
from libtbx.program_template import ProgramTemplate
# from six.moves import range
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass
from libtbx.utils import Sorry
# from datetime import datetime

master_phil_str = """
verbose = False
  .type = bool
  .help = Create graphics of plots (if Matplotlib is installed)
"""

def master_params():
  return iotbx.phil.parse(master_phil_str)

class running_proxy_table(dict):
  def __init__(self, atoms=None, plane=False):
    self.atoms=atoms
    self.plane=plane

  def __repr__(self):
    if self.plane:
      return self._plane_repr()
    return self._bond_repr()

  def _bond_repr(self):
    outl='Bonding\n'
    for k1, item in self.items():
      for k2 in item:
        outl+='  %3d - %3d ' % (k1,k2)
        if self.atoms:
          atom1=self.atoms[k1]
          atom2=self.atoms[k2]
          outl+=' %s - %s' % (atom1.quote(), atom2.quote())
          altloc1=atom1.parent().altloc
          altloc2=atom2.parent().altloc
          if altloc1!=altloc2:
            outl+=' * "%s" "%s" *' % (altloc1, altloc2)
        outl+='\n'
    return outl

  def _plane_repr(self):
    outl='Planes\n'
    for k1, item in self.items():
      outl+='  %3d ' % (k1)
      if self.atoms:
        atom1=self.atoms[k1]
        outl+='%s' % atom1.quote()
      outl+=' : ['
      for k2 in item:
        outl+=' %d' % k2
      outl+=']\n'
    return outl

  def add(self, i, j):
    self.setdefault(i, [])
    if j not in self[i]: self[i].append(j)
    self.setdefault(j, [])
    if i not in self[j]: self[j].append(i)
    if self.plane:
      if i not in self[i]: self[i].append(i)
      if j not in self[j]: self[j].append(j)

  def get_bonded(self, atom, as_atoms=False):
    if atom.i_seq not in self: return None
    if as_atoms:
      assert self.atoms
      rc=[]
      for i_seq in self[atom.i_seq]:
        rc.append(self.atoms[i_seq])
      return rc
    return self[atom.i_seq]

  def get_bonded_altlocs(self, atom):
    if atom.i_seq not in self: return None
    assert self.atoms
    rc=[]
    for i_seq in self[atom.i_seq]:
      rc.append(self.atoms[i_seq].parent().altloc)
    return rc

  def get_plane(self, atom):
    if atom.i_seq not in self: return None
    assert self.plane
    return self[atom.i_seq]

def compute(hierarchy, grm, verbose=False):
  if verbose: hierarchy.show()
  atoms=hierarchy.atoms()
  geometry=grm.geometry
  rpt=running_proxy_table(atoms=atoms, plane=True)
  for plane in geometry.planarity_proxies:
    for i, i_seq in enumerate(plane.i_seqs):
      for j, j_seq in enumerate(plane.i_seqs):
        if i==j: continue
        rpt.add(i_seq, j_seq)
  rbt=running_proxy_table(atoms=atoms)
  for i, bd in enumerate(geometry.bond_params_table):
    for key, item in bd.items():
      rbt.add(i,key)
  bond_simple, bond_asu=geometry.get_all_bond_proxies()
  sbt=running_proxy_table()
  for bond in bond_asu:
    sbt.add(bond.i_seq, bond.j_seq)

  for chain in hierarchy.chains():
    if verbose: print('  Chain %s' % chain.id)
    for conformer in chain.conformers():
      if verbose: print('    Conformer "%s" (protein %s, NA %s)' % (conformer.altloc, conformer.is_protein(), conformer.is_na()))
      for atom in conformer.atoms():
        bonding=rbt.get_bonded(atom)
        altlocs=rbt.get_bonded_altlocs(atom)
        plane=rpt.get_plane(atom)
        syms=sbt.get_bonded(atom)
        print('      %s bonds: %s altlocs: %s plane: %s syms: %s' % (atom.quote(), bonding, altlocs, plane, syms))
        if 0:
          bonding=rbt.get_bonded(atom, as_atoms=True)
          if bonding is not None:
            for atom in bonding:
              print(atom.quote())

  # think of an object to send back

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB file

Example:

  %(prog)s model=1ubq.pdb
""" % locals()

  master_phil_str = master_phil_str

  datatypes = ['model','phil']
  # data_manager_options = ['model_skip_expand_with_mtrix']
  # known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)
    model=self.data_manager.get_model()
    if not model.has_hd():
      raise Sorry('Model must have Hydrogen atoms')


  def run(self):
    model = self.data_manager.get_model()
    # p = m.get_default_pdb_interpretation_params()
    # p.pdb_interpretation.nonbonded_distance_cutoff = nonbonded_distance_cutoff
    model.process(make_restraints=True) #pdb_interpretation_params=p
    hierarchy = model.get_hierarchy()
    grm = model.get_restraints_manager()
    self.results = compute(hierarchy, grm)

  def get_results(self):
    return self.results
