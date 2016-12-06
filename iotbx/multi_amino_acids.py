from __future__ import division
#import copy

#from scitbx.math import dihedral_angle
#from libtbx.utils import Sorry

class multi_residue_class(list):
  def __init__(self,
               residues=None,
               length=3, # CDL & other psi/phi apps
               ):
    if residues:
      for residue in residues: self.append(residue)
    self.length = length
    #self.errors = []
    #self.start = None
    #self.end = None

  def __repr__(self):
    return self.show()

  def show(self):
    outl = "multi_residue_class"
    for residue in self:
      if residue is not None: outl += " %s(%s)" % (residue.resname, residue.resseq)
      else: outl += ' "%s"' % residue
    outl += " %s" % self.are_linked(return_value=True)
    if self.start is not None: outl += " start=T"
    if self.end is not None: outl += " end=T"
    return outl

  def show_detailed(self):
    outl = "multi_residue_class"
    for residue in self:
      outl += "\nREMARK"
      for atom in residue.atoms():
        outl += "\n%s" % atom.format_atom_record()
    return outl

  def atoms(self):
    for residue in self:
      for atom in residue.atoms():
        yield atom

  def get_omega_value(self,
                      omega_cdl=False,
                     ):
    #
    # this is very poor! there needs to be a better way to check for cis-
    #
    for i, residue in enumerate(self):
      if i==0: continue
      if omega_cdl:
        if len(self)==3:
          if i==1: continue
      else:
        if i==2: continue
      ccn1, outl1 = get_c_ca_n(residue, return_subset=True)
      ccn2, outl2 = get_c_ca_n(self[i-1], return_subset=True)
      ca1 = ccn1[1]
      n = ccn1[2]
      c = ccn2[0]
      ca2 = ccn2[1]
      omega_atoms = [ca1, n, c, ca2]
      if None in omega_atoms: return None
      omega = dihedral_angle(sites=[atom.xyz for atom in omega_atoms], deg=True)
      return omega

  def cis_group(self,
                limit=45., # based on the default in pdb_interpretation
                verbose=False):
    cis_peptide_bond = False
    omega = self.get_omega_values()
    if omega is None: return None
    if (180.-abs(omega))>limit:
      cis_peptide_bond = True
    if verbose:
      if cis_peptide_bond:
        print 'cis peptide bond', cis_peptide_bond, omega
        print self
    return cis_peptide_bond

  def is_pure_main_conf(self):
    tmp = [rg.is_pure_main_conf for rg in self]
    return len(filter(None, tmp))==self.length

  def are_linked(self,
                 return_value=False,
                 verbose=True):
    d2 = None
    for i, residue in enumerate(self):
      if i==0: continue
      ccn1, outl1 = get_c_ca_n(residue, return_subset=True)
      if self[i-1] is None: # place holder for omega CDL
        return False
      ccn2, outl2 = get_c_ca_n(self[i-1], return_subset=True)
      if ccn1 is None:
        for line in outl1:
          if line not in self.errors:
            self.errors.append(line)
        break
      if ccn2 is None:
        for line in outl2:
          if line not in self.errors:
            self.errors.append(line)
        break
      n = ccn1[2]
      c = ccn2[0]
      if n is None or c is None: return False
      if self.bond_params_table is None:
        d2 = distance2(n,c)
        if d2<4: bond=True
        else: bond=False
      else:
        bond=self.bond_params_table.lookup(c.i_seq, n.i_seq)
      if not bond:
        #assert c.i_seq
        #assert n.i_seq
        break
    else:
      return True
    #assert d2
    if return_value: return d2
    return False

  def append(self, residue):
    list.append(self, residue)
    while len(self)>3:
      del self[0]
    if len(self)>=2:
      while not self.are_linked():
        del self[0]
        if len(self)==0: break

  def get_resnames(self):
    return self[0].resname, self[1].resname, self[2].resname

  def get_phi_psi_atoms(self,
                        only_psi_phi_pairs=True,
                        force_plus_one=False,
                        omega_cdl=False,
                        ):
    if omega_cdl:
      if len(self) not in [self.length, self.length-1]:
        return None, None
      if len(self)==2:
        self.insert(0, None)
    else:
      if len(self)!=self.length: return None, None
    if force_plus_one: only_psi_phi_pairs=False
    if self[0] is None:
      backbone_i_minus_1 = None
    else:
      backbone_i_minus_1, junk = get_c_ca_n(self[0], return_subset=True)
      assert len(backbone_i_minus_1)==3
    backbone_i, junk = get_c_ca_n(self[1], return_subset=True)
    if None in backbone_i: return None
    backbone_i_plus_1, junk = get_c_ca_n(self[2], return_subset=True)
    if None in backbone_i_plus_1: return None
    assert len(backbone_i)==3
    assert len(backbone_i_plus_1)==3
    if omega_cdl: # phi(+1)
      phi_atoms = [
        backbone_i[0],
        backbone_i_plus_1[2],
        backbone_i_plus_1[1],
        backbone_i_plus_1[0],
        ]
    else:
      phi_atoms = [
        backbone_i_minus_1[0],
        backbone_i[2],
        backbone_i[1],
        backbone_i[0],
        ]
    psi_atoms = [
      backbone_i[2],
      backbone_i[1],
      backbone_i[0],
      backbone_i_plus_1[2],
      ]
    atoms = [phi_atoms, psi_atoms]
    if not only_psi_phi_pairs:
      if self.start:
        psi_atoms = [
          backbone_i_minus_1[2],
          backbone_i_minus_1[1],
          backbone_i_minus_1[0],
          backbone_i[2],
          ]
        atoms.insert(0, psi_atoms)
      if self.end or force_plus_one:
        phi_atoms = [
          backbone_i[0],
          backbone_i_plus_1[2],
          backbone_i_plus_1[1],
          backbone_i_plus_1[0],
          ]
        atoms.append(phi_atoms)
    if 0:
      for dihedral in atoms:
        print '-'*80
        for atom in dihedral:
          print atom.quote()
    return atoms

  def get_phi_psi_angles(self,
                         only_psi_phi_pairs=True,
                         force_plus_one=False,
                         omega_cdl=False,
                         verbose=False,
                         ):
    atoms = self.get_phi_psi_atoms(only_psi_phi_pairs=only_psi_phi_pairs,
                                   force_plus_one=force_plus_one,
                                   omega_cdl=omega_cdl,
                                  )
    if atoms is None: return None
    dihedrals = []
    for dihedral in atoms:
      phi_or_psi=dihedral_angle(sites=[atom.xyz for atom in dihedral], deg=True)
      dihedrals.append(phi_or_psi)
    if verbose:
      for phi_or_psi in dihedrals:
        print 'phi_or_psi',phi_or_psi
    return dihedrals

  def get_ramalyze_key(self,
                       limit=30.,
                       verbose=False,
                       ):
    from mmtbx.validation import ramalyze
    # defined in mmtbx.validation.ramalyze:
    # res_types = ["general", "glycine", "cis-proline", "trans-proline",
    #              "pre-proline", "isoleucine or valine"]
    #
    # This should be consistent with mmtbx/validation/ramalyze.py,
    # lines 219-240. Particularly, prepro comes before ile/val
    if self[1].resname == "PRO":
      if self.cis_group(limit=limit): return ramalyze.RAMA_CISPRO
      else: return ramalyze.RAMA_TRANSPRO
    elif self[2].resname == "PRO": return ramalyze.RAMA_PREPRO
    elif self[1].resname in ["ILE", "VAL"]: return ramalyze.RAMA_ILE_VAL
    elif self[1].resname == "GLY": return ramalyze.RAMA_GLYCINE
    else: return ramalyze.RAMA_GENERAL

  def get_dummy_dihedral_proxies(self, only_psi_phi_pairs=True):
    from cctbx.geometry_restraints import dihedral_proxy
    atoms = self.get_phi_psi_atoms(only_psi_phi_pairs=only_psi_phi_pairs)
    proxies = []
    if atoms is None: return proxies
    for dihedral in atoms:
      proxy = dihedral_proxy(
          i_seqs=[atom.i_seq for atom in dihedral],
          angle_ideal=0,
          weight=1)
      proxies.append(proxy)
    return proxies

