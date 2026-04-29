from __future__ import absolute_import, division, print_function
import math
import copy
import mmtbx.model
from libtbx.utils import null_out
from collections import OrderedDict
from scitbx.array_family import flex
from mmtbx.secondary_structure.build import ss_idealization as ssb

NaN = float('nan')

def format_HELIX_records_from_AEV(aev_values_dict):
  threshold1 = [2.48, 0.9, 4]
  threshold2 = [2.38, 0.9, 4]
  threshold3 = [2.28, 0.85, 2]
  helix_num = 1
  helix_list = [[], [], [], []]
  result = []
  cut_num = threshold1
  #This is select rules.This function judges if this atom is a B,M or E?
  def select(value, cut_num, residue, judge_list):
    #Is this atom a B?
    if value['B']>=value['M'] and value['B']>=value['E'] and value['B']>=cut_num:
      if judge_list[0]==[]:
        judge_list[0] = residue
      elif judge_list[1]!=[] and judge_list[2]!=[] and judge_list[3]==[]:
        judge_list[3] = residue
        judge_list[1] = residue
    # Is this atom a M
    elif value['M']>=value['E'] and value['M']>=cut_num:
      if judge_list[0]!=[]:
        judge_list[1] = residue
      elif judge_list[1]!=[] and judge_list[2]!=[] and judge_list[3]==[]:
        judge_list[3] = residue
    # Is this atom a E?
    elif value['E']>=value['B'] and value['E']>=value['M'] and value['E']>=cut_num:
      if judge_list[0]!=[] and judge_list[1]!=[] :
        judge_list[2] = residue
    else:
      if judge_list[0]!=[] and judge_list[1]!=[] and judge_list[2]!=[]:
        judge_list[3] = residue
      else:
        #reset
        judge_list = [[], [], [], []]
    return judge_list
  #Judge all atoms
  for key, value in aev_values_dict.items():
    if value['all'] >= cut_num[0]:
      helix_list = select(value, cut_num[1], key, helix_list)
    elif value['all'] < cut_num[0]:
      if NaN in value:
        helix_list = select(value, cut_num[1], key, helix_list)
      else:
        helix_list = select(value, cut_num[1], key, helix_list)
    # helices records input
    if not [] in helix_list:
      length = int(helix_list[2].split()[-1]) - int(helix_list[0].split()[-1]) + 1
      if length > cut_num[2]:
        if helix_list[1]==helix_list[3]:
          fmt = "HELIX   {0:>2}  {0:>2} {1:>}  {2:>}   {3:>36}"
          result.append(fmt.format(helix_num, helix_list[0], helix_list[2], length))
          start = copy.copy(helix_list[3])
          helix_list = [start, [], [], []]
        else:
          fmt = "HELIX   {0:>2}  {0:>2} {1:>}  {2:>}   {3:>36}"
          result.append(fmt.format(helix_num, helix_list[0], helix_list[2], length))
          helix_list = [[], [], [], []]
        helix_num += 1
  return result

def generate_perfect_helix(rs_values,
                           ts_values,
                           angular_rs_values,
                           radial_eta,
                           angular_eta,
                           angular_zeta,
                           crystal_symmetry,
                           n_residues=10,
                           residue_code="G",
                           ):
  """
  Compute AEV values for the perfect helix.
  """
  perfect_helix_ph = ssb.secondary_structure_from_sequence(ssb.alpha_helix_str,
    residue_code*n_residues)
  model = mmtbx.model.manager(
    model_input   = None,
    crystal_symmetry = crystal_symmetry,
    pdb_hierarchy = perfect_helix_ph,
    log           = null_out())
  model.crystal_symmetry()
  return AEV(model = model,
             rs_values=rs_values,
             ts_values=ts_values,
             angular_rs_values=angular_rs_values,
             radial_eta=radial_eta,
             angular_eta=angular_eta,
             angular_zeta=angular_zeta,
             ).get_values()

def compare(aev_values_dict):
  """
  Compare perfect helix with a target structure and get correlation coefficient
  values. The result include 3 direction AEVs CC values. If the c-alpha doesn't
  have BAEVs or EAEVs the CC values are None.
  """
  result = diff_class()
  perfect_helix = generate_perfect_helix(
                    crystal_symmetry=aev_values_dict.crystal_symmetry,
                    rs_values=aev_values_dict.rs_values,
                    ts_values=aev_values_dict.ts_values,
                    angular_rs_values=aev_values_dict.angular_rs_values,
                    radial_eta=aev_values_dict.radial_eta,
                    angular_eta=aev_values_dict.angular_eta,
                    angular_zeta=aev_values_dict.angular_zeta,
                    )

  def pretty_aev(v):
    outl = 'AEV'
    for i in v:
      outl += ' %0.3f' % i
    return outl


  def set_vals(result, d, verbose=False):
    for key1, value1 in d.items():
      result.setdefault(key1, OrderedDict())
      result[key1].setdefault(key, NaN)
      if value1 != [] and value != []:
        cc = flex.linear_correlation(
        x=flex.double(value), y=flex.double(value1)).coefficient()
        result[key1][key] = cc
      if verbose:
        print('comparing\n%s\n%s' % (pretty_aev(value), pretty_aev(value1)))
        print('  CC = %0.3f' % cc)
  for key,value in perfect_helix.items():
    if   key == 'B': set_vals(result=result, d=aev_values_dict.BAEVs)
    elif key == 'M': set_vals(result=result, d=aev_values_dict.MAEVs)
    elif key == 'E': set_vals(result=result, d=aev_values_dict.EAEVs)
  for key1, value1 in result.items():
    sum = 0
    for num in value1.values():
      if num != NaN:
        sum += num
    result[key1]['all'] = sum
  return result

# It is format calss of corelation coefficient values.
class diff_class(OrderedDict):
  def __repr__(self):
    outl = '...\n'
    for key, item in self.items():
      outl += '  %s :' % (key)
      for key1,value in item.items():
        outl += ' %s: '%key1
        if value is NaN:
          outl += 'NaN, '
        else:
          outl += '%0.2f, ' % value
      outl += '\n'
    return outl

# This is format class of AEV. It makes print of AEV more clearly.
class format_class(OrderedDict):
  def __init__(self, length_of_radial=None):
    OrderedDict.__init__(self)
    self.length_of_radial=length_of_radial


  def __repr__(self):
    outl = '...\n'
    for key, item in sorted(self.items()):
      outl += '  %s :' % (key)
      for i, v in enumerate(item):
        if i==self.length_of_radial: outl+='|'
        outl += '%0.4f, ' % v
      outl += '\n'
    return outl

class AEV(object):
  """
  Smith J S, Isayev O, Roitberg A E. ANI-1: an extensible neural network potential with DFT
  accuracy at force field computational cost[J]. Chemical science, 2017, 8(4): 3192-3203.
  """
  def __init__( self,
                model,
                rs_values = [2.0, 3.8, 5.2, 5.5, 6.2, 7.0, 8.6, 10.0],
                # probe distances (A) for radial
                radial_eta = 4,
                cutoff = 8.1,
                # radial cutoff distance
                ts_values = [0.392699, 1.178097, 1.963495, 2.748894],
                # probe angles (rad) for angular
                angular_rs_values = [3.8, 5.2, 5.5, 6.2],
                # probe distances (A) for angular
                angular_eta = 4,
                angular_zeta = 8,
                # parameters for probe angles
                ):
    self.hierarchy = model.get_hierarchy()
    self.rs_values = rs_values
    self.crystal_symmetry = model.crystal_symmetry()
    self.radial_eta = radial_eta
    self.cutoff = cutoff
    self.ts_values = ts_values
    self.angular_rs_values = angular_rs_values
    self.angular_eta = angular_eta
    self.angular_zeta = angular_zeta
    self.EAEVs = format_class(length_of_radial=len(self.rs_values))
    self.MAEVs = format_class(length_of_radial=len(self.rs_values))
    self.BAEVs = format_class(length_of_radial=len(self.rs_values))
    self.center_atom     = None
    self.chain_hierarchy = None
    self.generate_AEV()

  def get_values(self):
    result = OrderedDict()
    result['B'] = list(self.BAEVs.values())[0]
    result['M'] = list(self.MAEVs.values())[5]
    result['E'] = list(self.EAEVs.values())[-1]
    return result

  def empty_dict(self):
    results = OrderedDict()
    for atom in self.hierarchy.atoms():
      if atom.name == ' CA ':
        res_name = atom.format_atom_record()[17:20] + '  ' + atom.format_atom_record()[21:26]
        results.setdefault(res_name, [])
      self.BAEVs.update(results)
      self.MAEVs.update(results)
      self.EAEVs.update(results)
    return 0

  def generate_ca(self, length = 5):
    hierarchy = self.chain_hierarchy
    for chain in hierarchy.chains():
      con = chain.conformers()[0]
      b = int(length)
      for a in range(len(con.atoms())):
        rc = []
        if b<=len(con.atoms()):
          b = a + int(length)
          for atom in con.atoms()[a:b]:
            if atom.name == ' CA ':
              rc.append(atom)
        if len(rc) == int(b-a):
          yield  rc

  def generate_AEV(self):
    """
    Traversal the C-alpha atoms list and generate AEV values of every chain's
    every atom. Every C-alpha atom has 3 direction AEVs: forward(BAEVs),
    backward(EAEVs) and all direction(MAEVS).
    """
    chain_list = []
    self.empty_dict()
    for chain in self.hierarchy.chains():
      chain_list.append(chain.id)
    chain_list = list(set(chain_list))
    for chain in chain_list:
      chain_hierarchy = self.hierarchy.deep_copy()
      asc = chain_hierarchy.atom_selection_cache()
      sel = asc.selection("protein and name CA and chain " + chain )
      self.chain_hierarchy = chain_hierarchy.select(sel)
      for atomlist in self.generate_ca():
        self.center_atom = atomlist[0]
        self.BAEVs.update(self.calculate(atomlist))
        self.center_atom = atomlist[-1]
        self.EAEVs.update(self.calculate(atomlist))
        self.center_atom = atomlist[2]
        self.MAEVs.update(self.calculate(atomlist))

  def cutf(self, distance):
    """
    Formula (2), page 3194
    """
    if distance <= self.cutoff:
      Fc = 0.5 * math.cos(math.pi * distance / self.cutoff) + 0.5
    else:
      Fc = 0
    return Fc

  def calculate(self, atom_list):
    """
    Formula (3) and (4), page 3194
    """
    n = self.radial_eta
    assert  self.radial_eta==self.angular_eta
    l = self.angular_zeta
    AEVs = format_class()
    res_name = self.center_atom.format_atom_record()[17:20]+'  '+\
               self.center_atom.format_atom_record()[21:26]
    AEVs.setdefault(res_name, [])
    if atom_list != []:
      atom1 = self.center_atom
      atomlist = copy.copy(atom_list)
      if atom1 in atomlist:
        atomlist.remove(atom1)
      for Rs in self.rs_values:
        GmR = 0
        for atom2 in atomlist:
          R = atom1.distance(atom2)
          f = self.cutf(R)
          if f != 0:
            mR = math.exp(- n * ((R - Rs) ** 2)) * f
            GmR += mR
        AEVs[res_name].append(GmR)
      for Rs in self.angular_rs_values:
        for zetas in self.ts_values:
          i = 0
          GmA = 0
          zeta_list = []
          atomlist = copy.copy(atom_list)
          if atom1 in atomlist:
            atomlist.remove(atom1)
          for atom2 in atomlist[::-1]:
            if atom2 in atomlist:
              atomlist.remove(atom2)
            for atom3 in atomlist:
              Rij = atom1.distance(atom2)
              Rik = atom1.distance(atom3)
              ZETAijk = atom1.angle(atom2, atom3)
              if ZETAijk != 0:
                i += 1
                fk = self.cutf(Rik)
                fj = self.cutf(Rij)
                if fk != 0 and fj != 0:
                  zeta_list.append(ZETAijk)
                  mA = (((1 + math.cos(ZETAijk - zetas))) ** l) * \
                       math.exp(- n * ((((Rij + Rik) / 2) - Rs) ** 2)) * fj * fk
                  GmA += mA
          GmA = GmA * (2**(1-l))
          AEVs[res_name].append(GmA)
    return AEVs
