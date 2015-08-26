from __future__ import division
import copy

from scitbx.math import dihedral_angle
from libtbx.utils import Sorry

from mmtbx.conformation_dependent_library.cdl_utils import \
  get_c_ca_n, distance2, round_to_ten
from mmtbx.conformation_dependent_library.cdl_setup import columns

class RestraintsRegistry(dict):
  def __init__(self):
    self.n = {}

  def __repr__(self):
    outl = "RestraintsRegistry"
    outl += "\n  %s(%d)" % (self.keys(), len(self))
    outl += "\n  %s" % self.n
    return outl

  def __setitem__(self, key, item):
    if key in self:
      if self[key]!=item:
        self.n.setdefault(key,1)
        self.n[key]+=1
        dict.__setitem__(self, key, (self[key]+item))
    else:
      dict.__setitem__(self, key, item)

class ThreeProteinResidues(list):
  def __init__(self,
               geometry,
               length=3, # CDL & other psi/phi apps
               registry=None,
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

  def __repr__(self):
    return self.show()

  def show(self):
    outl = "ThreeProteinResidues"
    for residue in self:
      outl += " %s(%s)" % (residue.resname, residue.resseq)
    outl += " %s" % self.are_linked(return_value=True)
    if self.start is not None: outl += " start=T"
    if self.end is not None: outl += " end=T"
    return outl

  def show_detailed(self):
    outl = "ThreeProteinResidues"
    for residue in self:
      outl += "\nREMARK"
      for atom in residue.atoms():
        outl += "\n%s" % atom.format_atom_record()
    return outl

  def cis_group(self, limit=45., verbose=False):
    cis_peptide_bond = False
    for i, residue in enumerate(self):
      if i==0: continue
      ccn1, outl1 = get_c_ca_n(residue)
      ccn2, outl2 = get_c_ca_n(self[i-1])
      ca1 = ccn1[1]
      n = ccn1[2]
      c = ccn2[0]
      ca2 = ccn2[1]
      omega_atoms = [ca1, n, c, ca2]
      omega = dihedral_angle(sites=[atom.xyz for atom in omega_atoms], deg=True)
      if (180.-abs(omega))>limit:
        cis_peptide_bond = True
        break
    if verbose:
      if cis_peptide_bond:
        print 'cis peptide bond', cis_peptide_bond, omega
        print self
    return cis_peptide_bond

  def are_linked(self, return_value=False, verbose=True):
    d2 = None
    for i, residue in enumerate(self):
      if i==0: continue
      ccn1, outl1 = get_c_ca_n(residue)
      ccn2, outl2 = get_c_ca_n(self[i-1])
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

  def get_i_seqs(self):
    atoms = {}
    # i-1
    for name in ["C", "CA"]: # need CA_minus_1 for omega-CDL
      for atom in self[0].atoms():
        if atom.name.strip()==name:
          atoms["%s_minus_1" % name] = atom
          break
    # i
    for name in ["N", "CA", "CB", "C", "O"]:
      for atom in self[1].atoms():
        if atom.name.strip()==name:
          atoms["%s_i" % name] = atom
          break
      #else:
      #  if name not in ["CB", "O"]:
      #    print self
      #    for atom in self[1].atoms():
      #      print atom.name, atom.xyz
      #    assert 0
    # i+1
    for name in ["N", "CA"]: # need CA_plus_1 for omega-CDL
      for atom in self[2].atoms():
        if atom.name.strip()==name:
          atoms["%s_plus_1" % name] = atom
          break
    return atoms

  def get_resnames(self):
    return self[0].resname, self[1].resname, self[2].resname

  def get_phi_psi_atoms(self, only_psi_phi_pairs=True, force_plus_one=False):
    if len(self)!=self.length: return None, None
    if force_plus_one: only_psi_phi_pairs=False
    backbone_i_minus_1, junk = get_c_ca_n(self[0])
    backbone_i, junk = get_c_ca_n(self[1])
    backbone_i_plus_1, junk = get_c_ca_n(self[2])
    assert len(backbone_i_minus_1)==3
    assert len(backbone_i)==3
    assert len(backbone_i_plus_1)==3
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

  def get_cdl_key(self,
                  exact=False,
                  only_psi_phi_pairs=True,
                  force_plus_one=False,
                  verbose=False):
    atoms = self.get_phi_psi_atoms(only_psi_phi_pairs=only_psi_phi_pairs,
                                   force_plus_one=force_plus_one)
    dihedrals = []
    for dihedral in atoms:
      phi_or_psi=dihedral_angle(sites=[atom.xyz for atom in dihedral], deg=True)
      dihedrals.append(phi_or_psi)
    if verbose:
      for phi_or_psi in dihedrals:
        print 'phi_or_psi',phi_or_psi
    key = []
    for phi_or_psi in dihedrals:
      if exact:
        key.append(phi_or_psi)
      else:
        key.append(round_to_ten(phi_or_psi))
    return tuple(key)

  def get_dummy_dihedral_proxies(self, only_psi_phi_pairs=True):
    from cctbx.geometry_restraints import dihedral_proxy
    atoms = self.get_phi_psi_atoms(only_psi_phi_pairs=only_psi_phi_pairs)
    proxies = []
    for dihedral in atoms:
      proxy = dihedral_proxy(
          i_seqs=[atom.i_seq for atom in dihedral],
          angle_ideal=0,
          weight=1)
      proxies.append(proxy)
    return proxies

  def apply_updates(self,
                    restraint_values,
                    cdl_proxies,
                    ideal=True,
                    esd=True,
                    esd_factor=1.,
                    average=True,
                    verbose=False,
                    ):
    if not average:
      if restraint_values[0]=="I":
        print restraint_values
        assert 0
        return
    atoms = self.get_i_seqs()
    for i, value in enumerate(restraint_values):
      if i<2: continue
      if columns[i][0]=="s": continue
      code = columns[i][1:]
      names = []
      if code=="CNA":   names = ["C_minus_1", "N_i",  "CA_i"      ]
      elif code=="NAB": names = ["N_i",       "CA_i", "CB_i"      ]
      elif code=="NAC": names = ["N_i",       "CA_i", "C_i"       ]
      elif code=="BAC": names = ["CB_i",      "CA_i", "C_i"       ]
      elif code=="ACO": names = ["CA_i",      "C_i",  "O_i"       ]
      elif code=="ACN": names = ["CA_i",      "C_i",  "N_plus_1"  ]
      elif code=="OCN": names = ["O_i",       "C_i",  "N_plus_1"  ]
      elif code=="CN":  names = ["C_minus_1",  "N_i" ]
      elif code=="NA":  names = ["N_i",  "CA_i" ]
      elif code=="AB":  names = ["CA_i", "CB_i" ]
      elif code=="AC":  names = ["CA_i", "C_i" ]
      elif code=="CO":  names = ["C_i",  "O_i" ]
      # not all amino acids have a CB
      if "CB_i" in names and not "CB_i" in atoms: continue
      # sometimes the O is not in the model
      if "O_i" in names and not "O_i" in atoms: continue
      for j in range(len(names)):
        names[j] = atoms[names[j]].i_seq
      if len(names)==3:
        rnames = copy.deepcopy(names)
        rnames.reverse()
        angle_proxy = cdl_proxies.get(tuple(names), None)
        if angle_proxy is None:
          angle_proxy = cdl_proxies.get(tuple(rnames), None)
        if angle_proxy is None: continue
        if 0:
          outl=""
          for key in atoms:
            outl += "\n    %-10s %s" % ( key, atoms[key].quote())
          raise Sorry("""CDL angle to be changed not set in model.
  Possible problems:
    Residue on special positions.

  Check:%s""" % outl)
        if verbose:
          print " i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            angle_proxy.i_seqs,
            angle_proxy.angle_ideal,
            angle_proxy.weight,
            restraint_values[i],
            1/restraint_values[i+1]**2,
            )
        names.sort()
        self.registry[tuple(names)] = restraint_values[i]
        if ideal: angle_proxy.angle_ideal = restraint_values[i]
        if esd: angle_proxy.weight = esd_factor * 1/restraint_values[i+1]**2
      elif len(names)==2:
        bond=self.bond_params_table.lookup(*names)
        assert bond
        if verbose:
          print " i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            names,
            bond.distance_ideal,
            bond.weight,
            restraint_values[i],
            1/restraint_values[i+1]**2,
            )
        names.sort()
        self.registry[tuple(names)] = restraint_values[i]
        #print "BOND", 1/restraint_values[i+1]**2/bond.weight,1/restraint_values[i+1]**2, bond.weight
        if ideal: bond.distance_ideal = restraint_values[i]
        if esd: bond.weight = esd_factor * 1/restraint_values[i+1]**2
      else:
        assert 0

  def apply_average_updates(self, averages, verbose=False):
    if verbose:
      print averages
      print averages.n
    for key in averages.n:
      if len(key)==2:
        bond=self.bond_params_table.lookup(*key)
        bond.distance_ideal = averages[key]/averages.n[key]
      elif len(key)==3:
        rkey = list(copy.deepcopy(key))
        rkey.reverse()
        rkey=tuple(rkey)
#        for angle in self.restraints_manager.geometry.angle_proxies:
        for angle in self.geometry.angle_proxies:
          # could be better!
          akey = list(copy.deepcopy(angle.i_seqs))
          akey.sort()
          akey=tuple(akey)
          if akey==key or akey==rkey:
            angle.angle_ideal = averages[key]/averages.n[key]
            break
        else:
          print key,rkey
          print averages[key]
          print averages.n[key]
          assert 0
