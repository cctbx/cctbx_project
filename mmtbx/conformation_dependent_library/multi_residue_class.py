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

class ProteinResidues(list):
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
    if 0: return self.show()
    outl = ''
    for residue in self:
      outl += '%s ' % residue.resname
    return '"%s"\n' % outl

  def show(self):
    outl = "%sProteinResidues" % self.length
    for residue in self:
      if residue is not None: outl += " %s(%s)" % (residue.resname, residue.resseq)
      else: outl += ' "%s"' % residue
    outl += " %s" % self.are_linked(return_value=True)
    if self.start is not None: outl += " start=T"
    if self.end is not None: outl += " end=T"
    return outl

  def show_detailed(self):
    outl = "%sProteinResidues" % self.length
    for residue in self:
      outl += "\nREMARK"
      for atom in residue.atoms():
        outl += "\n%s" % atom.format_atom_record()
    return outl

  def atoms(self):
    for residue in self:
      for atom in residue.atoms():
        yield atom

  def get_omega_value(self): assert 0

  def _define_omega_a_la_duke_using_limit(self,
                                          omega,
                                          limit=45,
                                          ):
    if abs(omega)<limit: return 'cis'
    elif 180-abs(omega)<limit: return 'trans'
    else: return 'twisted'

  def cis_group(self,
                limit=45.,
                omega_cdl=False, # need last not middle
                verbose=False):
    cis_peptide_bond = False
    omega = self.get_omega_value(omega_cdl=omega_cdl)
    if omega is None: return None
    if self._define_omega_a_la_duke_using_limit(omega, limit=limit)=='cis':
      cis_peptide_bond = True
    if verbose:
      if cis_peptide_bond:
        print 'cis peptide bond', cis_peptide_bond, omega
        print self
    return cis_peptide_bond

  def trans_group(self, limit=30.):
    return not self.cis_group(limit=limit)

  def is_pure_main_conf(self):
    tmp = [rg.is_pure_main_conf for rg in self]
    return len(filter(None, tmp))==self.length

  def are_linked(self,
                 return_value=False,
                 use_distance_always=False,
                 bond_cut_off=2.,
                 verbose=True):
    d2 = None
    bond_cut_off *= bond_cut_off
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
        if d2<bond_cut_off: bond=True
        else: bond=False
      else:
        bond=self.bond_params_table.lookup(c.i_seq, n.i_seq)
        if not bond and use_distance_always:
          # needed for situations where atoms are added and the i_seq is updated
          if distance2(n,c)<bond_cut_off: bond=True
      if not bond:
        break
    else:
      return True
    if return_value: return d2
    return False

  def append(self, residue):
    list.append(self, residue)
    while len(self)>self.length:
      del self[0]
    if len(self)>=self.length-1:
      while not self.are_linked():
        del self[0]
        if len(self)==0: break

  def get_i_seqs(self): assert 0

  def get_resnames(self):
    rc = []
    for residue in self: rc.append(residue.resname)
    return rc

  def get_phi_psi_angles(self): assert 0

  def is_pure_main_conf(self):
    for one in self:
      if not one.is_pure_main_conf: return False
    return True

  def altloc(self):
    if self.is_pure_main_conf(): return ' '
    rc=[]
    for one in self:
      rc.append(self[0].parent().altloc)
    rc = filter(None,rc)
    assert rc
    return rc[0]

class TwoProteinResidues(ProteinResidues):
  def get_omega_value(self,
                      omega_cdl=False,
                      ):
    ccn1, outl1 = get_c_ca_n(self[0], return_subset=True)
    ccn2, outl2 = get_c_ca_n(self[1], return_subset=True)
    ca1 = ccn1[1]
    n = ccn1[2]
    c = ccn2[0]
    ca2 = ccn2[1]
    omega_atoms = [ca1, n, c, ca2]
    if None in omega_atoms: return None
    omega = dihedral_angle(sites=[atom.xyz for atom in omega_atoms], deg=True)
    return omega

class ThreeProteinResidues(ProteinResidues):
  def get_omega_value(self,
                      omega_cdl=False,
                     ):
    #
    # this is very poor! there needs to be a better way to check for cis-
    #
    for i, residue in enumerate(self):
      if i==0: continue
      if omega_cdl:
        if len(self)==self.length:
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

  def trans_group(self,
                  limit=45.,
                  verbose=False,
                  ):
    omega = self.get_omega_value() #omega_cdl=omega_cdl)
    if omega is None: return None
    if self._define_omega_a_la_duke_using_limit(omega, limit=limit)=='trans':
      return True
    return False

  def provide_second_sub_unit_if_unlinked(self):
    # used if residue is appended using superclass method
    if not self.are_linked():
      sub_unit = copy.copy(self) # calls append to delete first sub unit
      while not self.are_linked():
        del self[-1]
      return sub_unit
    return None

  def get_i_seqs(self):
    atoms = {}
    # i-1
    if self[0]:
      for name in [" C  ", " CA "]: # need CA_minus_1 for omega-CDL
        atom = self[0].find_atom_by(name=name)
        if atom: atoms["%s_minus_1" % name.strip()] = atom
    # i
    for name in [" N  ", " CA ", " CB ", " C  ", " O  "]:
      atom = self[1].find_atom_by(name=name)
      if atom: atoms["%s_i" % name.strip()] = atom
    # i+1
    for name in [" N  ", " CA "]: # need CA_plus_1 for omega-CDL
      atom = self[2].find_atom_by(name=name)
      if atom: atoms["%s_plus_1" % name.strip()] = atom
    return atoms

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
    # lines 179-195. Particularly, prepro comes before ile/val
    if self[1].resname == "PRO":
      if self.cis_group(limit=limit): return ramalyze.RAMA_CISPRO
      else: return ramalyze.RAMA_TRANSPRO
    elif self[1].resname == "GLY": return ramalyze.RAMA_GLYCINE
    elif self[2].resname == "PRO": return ramalyze.RAMA_PREPRO
    elif self[1].resname in ["ILE", "VAL"]: return ramalyze.RAMA_ILE_VAL
    else: return ramalyze.RAMA_GENERAL

  def get_phi_psi_atoms(self,
                        only_psi_phi_pairs=True,
                        force_plus_one=False,
                        omega_cdl=False,
                        verbose=False,
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
    if verbose: print backbone_i
    if None in backbone_i: return None
    backbone_i_plus_1, junk = get_c_ca_n(self[2], return_subset=True)
    if verbose: print backbone_i_plus_1, junk
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
    if verbose:
      print atoms
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
                                   verbose=verbose,
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

  def get_dummy_dihedral_proxies(self, only_psi_phi_pairs=True):
    #
    # Needs testing. One of the candidates is 3j0d, chain I, the first
    # residue is missing CA atom.
    #
    from cctbx.geometry_restraints import dihedral_proxy
    atoms = self.get_phi_psi_atoms(only_psi_phi_pairs=only_psi_phi_pairs)
    proxies = []
    if atoms is None: return proxies
    for dihedral in atoms:
      if None not in dihedral:
        proxy = dihedral_proxy(
            i_seqs=[atom.i_seq for atom in dihedral],
            angle_ideal=0,
            weight=1)
        proxies.append(proxy)
    return proxies
  #
  # CDL specific methods
  #
  def get_cdl_key(self,
                  exact=False,
                  only_psi_phi_pairs=True,
                  force_plus_one=False,
                  omega_cdl=False,
                  verbose=False):
    dihedrals=self.get_phi_psi_angles(only_psi_phi_pairs=only_psi_phi_pairs,
                                      omega_cdl=omega_cdl,
                                      verbose=verbose,
                                      )
    if dihedrals is None: return None
    key = []
    for phi_or_psi in dihedrals:
      if exact:
        key.append(phi_or_psi)
      else:
        key.append(round_to_ten(phi_or_psi))
    return tuple(key)

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
        angle_proxy = cdl_proxies.get(tuple(names), None)
        if angle_proxy is None:
          rnames = copy.deepcopy(names)
          rnames.reverse()
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
        if not bond:
          atoms = []
          for atom in self.atoms():
            if atom.i_seq in names: atoms.append(atom)
          outl = 'CDL error: bond not found between %s - %s' % (
            atoms[0].quote(),
            atoms[1].quote(),
            )
          raise Sorry(outl)
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
        assert restraint_values[i+1]<.1, 'CDL bond restraint larger than 0.1'
      else:
        assert 0

  def apply_average_updates(self, averages, verbose=False):
    if verbose:
      print averages
      print averages.n
    if not averages.n: return
    keys = averages.n.keys()
    for key in keys:
      if len(key)==2:
        bond=self.bond_params_table.lookup(*key)
        bond.distance_ideal = averages[key]/averages.n[key]
      elif len(key)==3:
        rkey = (key[2],key[1],key[0])
        averages.n[rkey]=averages.n[key]
    for angle in self.geometry.angle_proxies:
      if angle.i_seqs in averages.n:
        key = angle.i_seqs
        if key not in averages:
          assert 0
        angle.angle_ideal = averages[key]/averages.n[key]

if __name__=="__main__":
  import sys
  from iotbx import pdb
  from test_rdl import get_geometry_restraints_manager
  filename=sys.argv[1]
  pdb_inp = pdb.input(filename)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  geometry_restraints_manager = get_geometry_restraints_manager(filename)
  pdb_hierarchy.reset_i_seq_if_necessary()
  from mmtbx.conformation_dependent_library import generate_protein_twos
  from mmtbx.conformation_dependent_library import generate_protein_threes
  if 1:
    generate_protein_tuples = generate_protein_twos
  else:
    generate_protein_tuples = generate_protein_threes
  for threes in generate_protein_tuples(pdb_hierarchy,
                                        geometry_restraints_manager,
                                        #verbose=verbose,
                                        ):
    print threes
    print '  omega   %5.1f' % threes.get_omega_value()
    print "  cis?    %-5s %s" % (threes.cis_group(), threes.cis_group(limit=30))
    print "  trans?  %-5s %s" % (threes.trans_group(), threes.trans_group(limit=30))
    try: print "  rama    %s" % threes.get_ramalyze_key()
    except: print '  rama not specified'
    print '  conf    %s' % threes.is_pure_main_conf()
  print "OK"
