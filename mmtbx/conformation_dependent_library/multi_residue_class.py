from __future__ import division
import copy

from mmtbx.conformation_dependent_library.cdl_utils import \
  get_c_ca_n
from mmtbx.conformation_dependent_library.cdl_utils import \
  distance2, get_omega_value, get_phi_psi_angles
from mmtbx.conformation_dependent_library.cdl_utils import \
  get_ca_dihedrals

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

class LinkedResidues(list):
  def __init__(self,
               geometry,
               length=3, # CDL & other psi/phi apps
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
    return len(filter(None, tmp))==self.length

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

  def get_i_seqs(self): assert 0

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
    rc = filter(None,rc)
    assert rc
    return rc[0]

class ProteinResidues(LinkedResidues):
  def __init__(self,
               geometry,
               length=3, # CDL & other psi/phi apps
               registry=None,
               include_non_linked=False,
              ):
    LinkedResidues.__init__(self,
                            geometry,
                            length=3,
                            registry=registry,
                            include_non_linked=include_non_linked,
                            )

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
    outl += "\nREMARK"
    for residue in self:
      for atom in residue.atoms():
        outl += "\n%s" % atom.format_atom_record()
    return outl

  def get_omega_value(self): assert 0

  def _define_omega_a_la_duke_using_limit(self,
                                          omega,
                                          limit=45.,
                                          ):
    if omega is None: return None
    if abs(omega)<limit: return 'cis'
    elif 180-abs(omega)<limit: return 'trans'
    else: return 'twisted'

  def cis_group(self,
                limit=45.,
#                omega_cdl=False, # need last not middle
                verbose=False):
    # is any omega a cis angle?
    # assert not omega_cdl
    #cis_peptide_bond = False
    #omega = self.get_omega_value(omega_cdl=omega_cdl)
    #if omega is None: return None
    omegas = self.get_omega_values()
    assert omegas
    def _is_cis(angle):
      return self._define_omega_a_la_duke_using_limit(angle, limit=limit)=='cis'
    if filter(_is_cis, omegas): return True
    return False
    #if self._define_omega_a_la_duke_using_limit(omega, limit=limit)=='cis':
    #  cis_peptide_bond = True

  def trans_group(self, limit=45.):
    return not self.cis_group(limit=limit)

  def cis_trans_twisted_list(self, limit=45.):
    omegas = self.get_omega_values()
    def _is_cis_trans_twisted(angle):
      return self._define_omega_a_la_duke_using_limit(angle, limit=limit)
    return map(_is_cis_trans_twisted, omegas)

  def are_linked(self,
                 return_value=False,
                 use_distance_always=False,
                 bond_cut_off=3., # Same as link_distance_cutoff of pdb_interpretation
                 allow_poly_ca=False,
                 poly_ca_cut_off=4.,
                 verbose=True):
    '''
    Need to add poly-Calpha chains
      CA-CA 4.5 is use in CaBLAM, maybe shorter
    '''
    if allow_poly_ca:
      assert 0
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

  def get_phi_psi_angles(self): assert 0

  def get_omega_values(self, verbose=False):
    rc=[]
    for i, residue in enumerate(self):
      if i==0: continue
      omega = get_omega_value(residue, self[i-1], verbose=verbose)
      rc.append(omega)
    return rc

class TwoProteinResidues(ProteinResidues):
  def get_omega_value(self):
    return get_omega_value(self[1], self[0])

class ThreeProteinResidues(ProteinResidues):
  def get_omega_values(self,
                       #omega_cdl=None,
                       verbose=False,
                       ):
    #assert omega_cdl is None, 'can not use omega_cdl for %sProteinResidues' % self.length
    return ProteinResidues.get_omega_values(self, verbose=verbose)

  def get_phi_psi_atoms(self,
                        only_psi_phi_pairs=True,
                        force_plus_one=False,
                        verbose=False,
                        ):
    if len(self)!=self.length: return None, None
    if force_plus_one: only_psi_phi_pairs=False
    if self[0] is None:
      backbone_i_minus_1 = None
    else:
      backbone_i_minus_1, junk = get_c_ca_n(self[0], return_subset=True)
      assert len(backbone_i_minus_1)==self.length
    backbone_i, junk = get_c_ca_n(self[1], return_subset=True)
    if verbose: print backbone_i
    if None in backbone_i: return None
    backbone_i_plus_1, junk = get_c_ca_n(self[2], return_subset=True)
    if verbose: print backbone_i_plus_1, junk
    if None in backbone_i_plus_1: return None
    assert len(backbone_i)==self.length
    assert len(backbone_i_plus_1)==self.length
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
    if verbose: print atoms
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
    if verbose:
      for dihedral in atoms:
        print '-'*80
        for atom in dihedral:
          print atom.quote()
    return atoms

  def get_phi_psi_angles(self, verbose=False):
    if verbose:
      for residue in self:
        print residue.id_str()
    return get_phi_psi_angles(self, verbose=verbose)

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

  def provide_second_sub_unit_if_unlinked(self):
    # used if residue is appended using superclass method
    if not self.are_linked():
      sub_unit = copy.copy(self) # calls append to delete first sub unit
      while not self.are_linked():
        del self[-1]
      return sub_unit
    return None

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

class FourProteinResidues(ThreeProteinResidues):
  def get_ca_dihedrals(self, verbose=False):
    if verbose:
      for residue in self:
        print residue.id_str()
    return get_ca_dihedrals(self)

class FiveProteinResidues(FourProteinResidues):
  def get_cablam_info(self):
    assert 0

def _get_atoms(atom_group, atom_names):
  atoms, outl = get_c_ca_n(atom_group, atom_names)
  if atoms is None:
    for i in range(len(atom_names)):
      atom_names[i] = atom_names[i].replace("'", '*')
    atoms, outl = get_c_ca_n(atom_group, atom_names)
  return atoms

def calc_pseudorotation(t0,t1,t2,t3,t4):
  import math
  if t0 > 180.0: t0 = t0 - 360.0
  if t1 > 180.0: t1 = t1 - 360.0
  if t2 > 180.0: t2 = t2 - 360.0
#JC hack
  if t2 == 0.0: t2 = 0.1
#/JC
  if t3 > 180.0: t3 = t3 - 360.0
  if t4 > 180.0: t4 = t4 - 360.0

  taus = [t0, t1, t2, t3, t4]

  tanP = ((taus[4] + taus[1]) - (taus[3] + taus[0]))/(2 * taus[2] * (math.sin(36.0*math.pi/180.0) + math.sin(72.0*math.pi/180.0)))

  P = math.atan(tanP)*180.0/math.pi
  if taus[2] < 0: P = P + 180.0
  elif tanP < 0: P = P + 360.0
  #P = "%.1f" % P
  return P

def get_distance(ag1, ag2, an1, an2):
  atoms = _get_atoms(ag1, an1) + _get_atoms(ag2, an2)
  # for atom in atoms: print atom.quote()
  return atoms[0].distance(atoms[1])

def get_torsion(ag1, ag2, an1, an2, limits='-180-180'):
  from scitbx.math import dihedral_angle
  atoms = _get_atoms(ag1, an1) + _get_atoms(ag2, an2)
  omega = dihedral_angle(sites=[atom.xyz for atom in atoms], deg=True)
  if limits=='-180-180':
    if omega>180:
      print omega, limits
      assert 0
  elif limits=='0-360':
    if omega<0:
      omega+=360
  # for atom in atoms: print atom.quote()
  return omega

class TwoNucleicResidues(LinkedResidues):
  def show(self):
    outl = "%sNucleicResidues" % self.length
    for residue in self:
      if residue is not None: outl += " %s(%s)" % (residue.resname, residue.resseq)
      else: outl += ' "%s"' % residue
    outl += " %s" % self.are_linked(return_value=True)
    if self.start is not None: outl += " start=T"
    if self.end is not None: outl += " end=T"
    return outl

  @staticmethod
  def get_o3prime_p(residue, return_subset=False):
    rc = get_c_ca_n(residue, atom_name_list=[' O3', ' P  '], return_subset=return_subset)
    if rc[0] is None:
      rc = get_c_ca_n(residue, atom_name_list=[' O3*', ' P  '], return_subset=return_subset)
    return rc

  def are_linked(self,
                 return_value=False,
                 use_distance_always=False,
                 bond_cut_off=3.5, # Same as link_distance_cutoff of pdb_interpretation
                 verbose=True,
                 ):
    bond_cut_off *= bond_cut_off
    for i, residue in enumerate(self):
      if i==0: continue
      op1, outl1 = self.get_o3prime_p(residue, return_subset=False)
      # if self[i-1] is None: # place holder for omega CDL
      #   return False
      op2, outl2 = self.get_o3prime_p(self[i-1], return_subset=False)
      # if ccn1 is None:
      #   for line in outl1:
      #     if line not in self.errors:
      #       self.errors.append(line)
      #   break
      # if ccn2 is None:
      #   for line in outl2:
      #     if line not in self.errors:
      #       self.errors.append(line)
      #   break
      p = op1[1]
      o3prime = op2[0]
      if p is None or o3prime is None: return False
      if self.bond_params_table is None:
        d2 = distance2(p,o3prime)
        if d2<bond_cut_off: bond=True
        else: bond=False
      else:
        bond=self.bond_params_table.lookup(p.i_seq, o3prime.i_seq)
        if not bond and use_distance_always:
          # needed for situations where atoms are added and the i_seq is updated
          if distance2(p,o3prime)<bond_cut_off: bond=True
      if not bond:
        break
    else:
      return True
    if return_value: return d2
    return False

  def get_base_types(self):
    rc = []
    for base in self:
      for atom in base.atoms():
        if atom.name==' N9 ':
          rc.append('R')
          break
      else:
        rc.append('Y')
    return rc

  def get_id(self):
    outl = []
    outl.append(self[0].parent().parent().id)
    outl.append(self[0].resname.strip())
    outl.append(self[0].resseq.strip())
    assert not self[0].parent().altloc
    outl.append(self[1].resname.strip())
    outl.append(self[1].resseq.strip())
    assert not self[1].parent().altloc
    return '_'.join(outl)

  def get_ntc_angles(self):
    angles = {
      'd' :[[" C5'", " C4'", " C3'", " O3'"],[]], # delta0
      'e' :[[" C4'", " C3'", " O3'" ],       [" P  "]], # epsilon
      'z' :[[" C3'", " O3'"],                [" P  ", " O5'"]], # zeta
      'a1':[[" O3'"],                        [" P  ", " O5'", " C5'"]], # alpha
      'b1':[[],                              [" P  ", " O5'", " C5'", " C4'"]], # beta
      'g1':[[],                              [" O5'", " C5'", " C4'", " C3'"]], # gamma
      'd1':[[],                              [" C5'", " C4'", " C3'", " O3'"]], # delta1
    }
    types = self.get_base_types()
    if types[0]=='R':
      angles['ch'] = [[" O4'", " C1'", " N9 ", " C4 "],[]] # chi0
      N0 = ' N9 '
    else:
      angles['ch'] = [[" O4'", " C1'", " N1 ", " C2 "],[]] # chi0
      N0 = ' N1 '
    if types[1]=='R':
      angles['ch1'] = [[], [" O4'", " C1'", " N9 ", " C4 "]] # chi1
      N1 = ' N9 '
    else:
      angles['ch1'] = [[], [" O4'", " C1'", " N1 ", " C2 "]] # chi1
      N1 = ' N1 '
    angles['NCCN'] = [[N0, " C1'"], [" C1'", N1]]
    rc = {}
    for angle, atom_names in angles.items():
      rc[angle] = get_torsion(self[0], self[1], atom_names[0], atom_names[1], limits='0-360')
    rc['NN'] = get_distance(self[0], self[1], [N0], [N1])
    rc['CC'] = get_distance(self[0], self[1], [" C1'"], [" C1'"])
    # tau
    args1 = []
    args2 = []
    for atom_names in [
      [" C4'", " O4'", " C1'", " C2'"],
      [" O4'", " C1'", " C2'", " C3'"],
      [" C1'", " C2'", " C3'", " C4'"],
      [" C2'", " C3'", " C4'", " O4'"],
      [" C3'", " C4'", " O4'", " C1'"],
      ]:
      args1.append(get_torsion(self[0], self[1], atom_names, []))
      args2.append(get_torsion(self[0], self[1], [], atom_names))
    rc['P'] = calc_pseudorotation(*tuple(args1))
    rc['P1'] = calc_pseudorotation(*tuple(args2))
    for label, item in rc.items():
      # print '  %s : %0.2f' % (label, item)
      rc[label] = '%0.1f' % item
    rc['step_id'] = self.get_id()
    return rc

  def get_ntc_coordinates(self):
    query = {}
    for atom_key in ['C5pa',
                     'C4pa',
                     'O4pa',
                     'C3pa',
                     'O3pa',
                     'C2pa',
                     'C1pa',
                     'N19a',
                     'C24a',
                     'Pb',
                     'O5pb',
                     'C5pb',
                     'C4pb',
                     'O4pb',
                     'C3pb',
                     'O3pb',
                     'C2pb',
                     'C1pb',
                     'N19b',
                     'C24b',
      ]:
      if atom_key[-1]=='a': atom_group = self[0]
      elif atom_key[-1]=='b': atom_group = self[1]
      else: assert 0
      if atom_key.find('P')>-1: names = [' P  ']
      elif atom_key.find('N19')>-1: names = [' N1 ', ' N9 ']
      elif atom_key.find('C24')>-1: names = [' C2 ', ' C4 ']
      else: names = ['%4s' % atom_key[:-1].replace('p',"'")]
      for name in names:
        atom = atom_group.find_atom_by(name=name)
        if atom is None:
          atom = atom_group.find_atom_by(name=name.replace("'", '*'))
        if atom: break
      else:
        assert atom
      query[atom_key]= ['%s'%atom.xyz[0], '%s'%atom.xyz[1], '%s'%atom.xyz[2]]
    query['step_id'] = self.get_id()
    return query

if __name__=="__main__":
  import sys
  from iotbx import pdb
  from test_rdl import get_geometry_restraints_manager
  filename=sys.argv[1]
  pdb_inp = pdb.input(filename)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  geometry_restraints_manager = get_geometry_restraints_manager(filename)
  pdb_hierarchy.reset_i_seq_if_necessary()
  from mmtbx.conformation_dependent_library import generate_protein_fragments
  for i in range(2,6):
    for threes in generate_protein_fragments(pdb_hierarchy,
                                             geometry_restraints_manager,
                                             length=i,
                                             #verbose=verbose,
                                             ):
      print threes
      try: print '  omega   %5.1f' % threes.get_omega_value()
      except: print '  omega is not valid' # intentional
      print '  omegas  %s' % threes.get_omega_values()
      try: print "  cis?    %-5s %s" % (threes.cis_group(), threes.cis_group(limit=30))
      except: print '  cis? is not valid' # intentional
      try: print "  trans?  %-5s %s" % (threes.trans_group(), threes.trans_group(limit=30))
      except: print '  tran? is not valid' # intentional
      print '  cis/trans/twisted? %s' % ' '.join(threes.cis_trans_twisted_list())
      try: print "  rama    %s" % threes.get_ramalyze_key()
      except: print '  rama not specified' # intentional
      print '  conf    %s' % threes.is_pure_main_conf()
      try: print '  phi/psi %s' % threes.get_phi_psi_angles()
      except: print '  phi/psi not specified' # intentional
      try: print '  CA dihedrals %s' % threes.get_ca_dihedrals()
      except: print '  CA dihedrals not specified' # intentional
    print "OK",i+2
