from __future__ import absolute_import, division, print_function
import copy

from scitbx.math import dihedral_angle
from libtbx.utils import Sorry
from mmtbx.conformation_dependent_library.multi_residue_class import \
  ThreeProteinResidues
from mmtbx.conformation_dependent_library.cdl_utils import \
  get_c_ca_n, round_to_ten, get_omega_value
from mmtbx.conformation_dependent_library.cdl_setup import columns

from cctbx.geometry_restraints.linking_class import linking_class
from six.moves import range
origin_ids = linking_class()

class ThreeProteinResiduesWithCDL(ThreeProteinResidues):
  #
  # CDL specific methods
  #
  def get_i_seqs(self):
    atoms = {}
    # i-1
    if self[0]:
      for name in [" C  ", " CA "]: # need CA_minus_1 for omega-CDL
        atom = self[0].find_atom_by(name=name)
        if atom: atoms["%s_minus_1" % name.strip()] = atom
    # i
    for name in [" N  ", " CA ", " CB ", " C  ", " O  ", ' H  ', ' CD ', ' CG ']:
      atom = self[1].find_atom_by(name=name)
      if atom: atoms["%s_i" % name.strip()] = atom
    # i+1
    for name in [" N  ", " CA "]: # need CA_plus_1 for omega-CDL
      atom = self[2].find_atom_by(name=name)
      if atom: atoms["%s_plus_1" % name.strip()] = atom
    return atoms

  def cis_group(self, limit=45, omega_cdl=False):
    omegas = self.get_omega_values()
    assert omegas
    if omega_cdl: del omegas[0]
    else: del omegas[1]
    def _is_cis(angle):
      return self._define_omega_a_la_duke_using_limit(angle, limit=limit)=='cis'
    if list(filter(_is_cis, omegas)): return True
    return False

  def get_omega_value(self, omega_cdl=False):
    #
    # this is very poor! there needs to be a better way to check for cis-
    #
    assert not omega_cdl
    for i, residue in enumerate(self):
      if i==0: continue
      if omega_cdl:
        if len(self)==self.length:
          if i==1: continue
      else:
        if i==2: continue
      omega = get_omega_value(residue, self[i-1])
      return omega

  def get_phi_psi_atoms(self,
                        only_psi_phi_pairs=True,
                        force_plus_one=False,
                        omega_cdl=False,
                        verbose=False,
                        ):
    if omega_cdl:
      if len(self) not in [self.length, self.length-1]:
        return None, None
      if len(self)==self.length-1:
        self.insert(0, None)
    else:
      if len(self)!=self.length: return None, None
    if force_plus_one: only_psi_phi_pairs=False
    if self[0] is None:
      backbone_i_minus_1 = None
    else:
      backbone_i_minus_1, junk = get_c_ca_n(self[0], return_subset=True)
      assert len(backbone_i_minus_1)==self.length
    backbone_i, junk = get_c_ca_n(self[1], return_subset=True)
    if verbose: print(backbone_i)
    if None in backbone_i: return None
    backbone_i_plus_1, junk = get_c_ca_n(self[2], return_subset=True)
    if verbose: print(backbone_i_plus_1, junk)
    if None in backbone_i_plus_1: return None
    assert len(backbone_i)==self.length
    assert len(backbone_i_plus_1)==self.length
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
      print(atoms)
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
        print('-'*80)
        for atom in dihedral:
          print(atom.quote())
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
        print('phi_or_psi',phi_or_psi)
    return dihedrals

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
    if None in dihedrals: return None
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
    def _get_angle_proxy(names):
      angle_proxy = cdl_proxies.get(tuple(names), None)
      if angle_proxy is None:
        rnames = copy.deepcopy(names)
        rnames.reverse()
        angle_proxy = cdl_proxies.get(tuple(rnames), None)
      return angle_proxy
    def _get_i_seqs(names):
      for j in range(len(names)):
        if names[j] not in atoms:
          # print('names not found',names)
          return None
        names[j] = atoms[names[j]].i_seq
      return names
    ####################
    if not average:
      if restraint_values[0]=="I":
        print(restraint_values)
        assert 0
        return
    atoms = self.get_i_seqs()
    for i, value in enumerate(restraint_values):
      if i<2: continue
      if columns[i][0]=="s": continue
      if restraint_values[i] is None: continue
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
      # needed for cis_127
      elif code=='CND': names = ['C_minus_1', 'N_i', 'CD_i']
      elif code=='AND': names = ['CA_i', 'N_i', 'CD_i']
      elif code=='NDG': names = ['N_i', 'CD_i', 'CG_i']
      elif code=='ABG': names = ['CA_i','CB_i', 'CG_i']
      elif code=='BGD': names = ['CB_i','CG_i', 'CD_i']
      elif code=='BG':  names = ['CB_i','CG_i']
      elif code=='GD':  names = ['CG_i','CD_i']
      elif code=='ND':  names = ['N_i', 'CD_i']
      # not all amino acids have a CB
      if "CB_i" in names and not "CB_i" in atoms: continue
      # sometimes the O is not in the model
      if "O_i" in names and not "O_i" in atoms: continue
      names = _get_i_seqs(names)
      if names is None: continue
      if len(names)==3:
        angle_proxy = _get_angle_proxy(names)
        if angle_proxy is None: continue
        if angle_proxy.origin_id==origin_ids.get_origin_id('edits'): continue
        if 0:
          outl=""
          for key in atoms:
            outl += "\n    %-10s %s" % ( key, atoms[key].quote())
          raise Sorry("""CDL angle to be changed not set in model.
  Possible problems:
    Residue on special positions.

  Check:%s""" % outl)
        if verbose:
          print(" %s i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            code,
            angle_proxy.i_seqs,
            angle_proxy.angle_ideal,
            angle_proxy.weight,
            restraint_values[i],
            1/restraint_values[i+1]**2,
            ))
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
          print(" i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
            names,
            bond.distance_ideal,
            bond.weight,
            restraint_values[i],
            1/restraint_values[i+1]**2,
            ))
        names.sort()
        self.registry[tuple(names)] = restraint_values[i]
        #print "BOND", 1/restraint_values[i+1]**2/bond.weight,1/restraint_values[i+1]**2, bond.weight
        if ideal: bond.distance_ideal = restraint_values[i]
        if esd: bond.weight = esd_factor * 1/restraint_values[i+1]**2
        assert restraint_values[i+1]<.1, 'CDL bond restraint larger than 0.1'
      else:
        assert 0, 'names %s not found' % names
    # adjust X-N-H angles to obtain planar
    nh_atoms = 0
    if self[0].find_atom_by(name=" C  "):
      nh_atoms += 1
    for name in [" N  ", " CA ", ' H  ']:
      if self[1].find_atom_by(name=name):
        nh_atoms += 1
    if nh_atoms == 4:
      CNCA = _get_angle_proxy(_get_i_seqs(["C_minus_1", "N_i", "CA_i"]))
      CNH = _get_angle_proxy(_get_i_seqs(["C_minus_1", "N_i", "H_i"]))
      CANH = _get_angle_proxy(_get_i_seqs(["CA_i", "N_i", "H_i"]))
      if not (CNH and CANH):
        error_atoms = []
        for atom in self[0].atoms():
          error_atoms.append('%s\n' %atom.quote())
        for atom in self[1].atoms():
          error_atoms.append('%s\n' %atom.quote())
        raise Sorry('''
  Certain angles in the protein chain (C-N-H or CA-N-H) have not been found
  in the restraints by the Conformational Dependent Library. This usually
  means that the protein backbone is traversing a special position.

  This is unlikely.

  However, to proceed, set cdl=False.

%s
                    ''' % ''.join(error_atoms))
      total = CNCA.angle_ideal + CNH.angle_ideal + CANH.angle_ideal
      diff = (total-360)/2
      CNH.angle_ideal-=diff
      CANH.angle_ideal-=diff
      assert 360 - (CNCA.angle_ideal + CNH.angle_ideal + CANH.angle_ideal)<0.1

  def apply_average_updates(self, averages, verbose=False):
    if verbose:
      print(averages)
      print(averages.n)
    if not averages.n: return
    for key in list(averages.n.keys()):
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
