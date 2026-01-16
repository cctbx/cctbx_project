from __future__ import absolute_import, division, print_function
class LinkedResidues(list):
  def __init__(self,
               geometry,
               length=3, # CDL & other psi/phi apps
               allow_poly_ca=False,
               registry=None,
               include_non_linked=False,
              ):
    assert registry is not None
    self.length = length
    self.geometry = geometry
    self.registry = registry
    if geometry is None:
      self.bond_params_table = None
    else:
      self.bond_params_table = geometry.bond_params_table
    self.errors = []
    self.start = None
    self.end = None
    self.include_non_linked = include_non_linked
    self.allow_poly_ca = allow_poly_ca

  def __repr__(self):
    if 1: return self.show()
    outl = ''
    for residue in self:
      outl += '%s ' % residue.resname
    return '"%s"\n' % outl

  def show(self): assert 0

  def show_detailed(self): assert 0

  def atoms(self):
    for residue in self:
      for atom in residue.atoms():
        yield atom

  def is_pure_main_conf(self):
    tmp = [rg.is_pure_main_conf for rg in self]
    return len(list(filter(None, tmp)))==self.length

  def are_linked(self, *args, **kwds): assert 0

  def append(self, residue):
    list.append(self, residue)
    while len(self)>self.length:
      del self[0]
    if self.include_non_linked: return
    if len(self)>=self.length-1:
      while not self.are_linked():
        del self[0]
        if len(self)==0: break

  def get_i_seqs(self):
    rc=[]
    for residue in self:
      for atom in residue.atoms():
        rc.append(atom.i_seq)
    return rc

  def get_resnames(self):
    rc = []
    for residue in self: rc.append(residue.resname)
    return rc

  def is_pure_main_conf(self):
    for one in self:
      if not one.is_pure_main_conf: return False
    return True

  def altloc(self):
    if self.is_pure_main_conf(): return ' '
    rc=[]
    for one in self:
      rc.append(self[0].parent().altloc)
    rc = list(filter(None,rc))
    assert rc
    return rc[0]
