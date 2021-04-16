from __future__ import absolute_import, division, print_function
from mmtbx.conformation_dependent_library.LinkedResidues import LinkedResidues
from mmtbx.conformation_dependent_library.cdl_utils import \
  get_c_ca_n, distance2
from six.moves import range

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

def _get_atoms(atom_group, atom_names):
  atoms, outl = get_c_ca_n(atom_group, atom_names)
  if atoms is None:
    for i in range(len(atom_names)):
      atom_names[i] = atom_names[i].replace("'", '*')
    atoms, outl = get_c_ca_n(atom_group, atom_names)
  return atoms

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
      print(omega, limits)
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
    rc = get_c_ca_n(residue, atom_name_list=[" O3'", ' P  '], return_subset=return_subset)
    if rc[0] is None:
      rc = get_c_ca_n(residue, atom_name_list=[" O3'", ' P  '], return_subset=return_subset)
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
      op1, outl1 = self.get_o3prime_p(residue, return_subset=True)
      # if self[i-1] is None: # place holder for omega CDL
      #   return False
      op2, outl2 = self.get_o3prime_p(self[i-1], return_subset=True)
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
