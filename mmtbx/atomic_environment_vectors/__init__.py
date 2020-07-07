from __future__ import absolute_import, division, print_function
import math
import copy
import mmtbx.model
from libtbx.utils import null_out
from collections import OrderedDict
from scitbx.array_family import flex
from mmtbx.conformation_dependent_library import generate_protein_fragments
from mmtbx.secondary_structure.build import ss_idealization as ssb

def format_HELIX_records_from_AEV(aev_values_dict, cc_cutoff):
  """
  Geting HELIX records for target model using AEV values and specified cutoff.
  """
  start = []
  end = []
  M = 0
  i = 0
  result = []
  for key,value in aev_values_dict.items():
    if start==[] and value['B'] > cc_cutoff and value['M'] > cc_cutoff - 0.1:
      start = key
      length = 1
    elif start and value['M'] > cc_cutoff:
      M = 1
      length += 1
      if value['E'] > cc_cutoff:
        end = key
    elif start and M == 1 and value['E'] + value['M'] > 2 * cc_cutoff:
      end = key
      length += 1
    else:
      if start and end and M ==1 and length > 3:
        i += 1
        fmt = "HELIX   {0:>2}  {0:>2} {1:>}  {2:>}   {3:>36}"
        result.append(fmt.format(i, start, end, length))
      start = []
      end = []
      M = 0
  return result

def generate_perfect_helix(n_residues=10, residue_code="G"):
  """
  Compute AEV values for the perfect helix.
  """
  perfect_helix_ph = ssb.secondary_structure_from_sequence(ssb.alpha_helix_str,
    residue_code*n_residues)
  model = mmtbx.model.manager(
    model_input   = None,
    pdb_hierarchy = perfect_helix_ph,
    build_grm     = True,
    log           = null_out())
  return AEV(model = model).get_values()

def compare(aev_values_dict):
  """
  Compare perfect helix with a target structure and get correlation coefficient
  values. The result include 3 direction AEVs CC values. If the c-alpha doesn't
  have BAEVs or EAEVs the CC values are 0.
  """
  result = diff_class()
  perfect_helix = generate_perfect_helix()
  def set_vals(result, d):
    for key1, value1 in d.items():
      if value1 != []:
        cc = flex.linear_correlation(
          x=flex.double(value), y=flex.double(value1)).coefficient()
        result.setdefault(key1, OrderedDict())
        result[key1].setdefault(key, cc)
      else:
        result.setdefault(key1, OrderedDict())
        result[key1].setdefault(key, 1)
  for key,value in perfect_helix.items():
    if   key == 'B': set_vals(result=result, d=aev_values_dict.BAEVs)
    elif key == 'M': set_vals(result=result, d=aev_values_dict.MAEVs)
    elif key == 'E': set_vals(result=result, d=aev_values_dict.EAEVs)
  return result

# It is format calss of corelation coefficient values.
class diff_class(OrderedDict):
  def __repr__(self):
    outl = '...\n'
    for key, item in self.items():
      outl += '  %s :' % (key)
      for key1,value in item.items():
        outl += ' %s: '%key1
        outl += '%0.2f, ' % value
      outl += '\n'
    return outl

# This is format class of AEV. It makes print of AEV more clearly.
class format_class(OrderedDict):
  def __repr__(self):
    outl = '...\n'
    for key, item in self.items():
      outl += '  %s :' % (key)
      for v in item:
        print (v)
        outl += '%0.4f, ' % v
      outl += '\n'
    return outl

class AEV(object):
  """
  Cite paper...
  """

  def __init__(self, model):
    self.hierarchy = model.get_hierarchy()
    self.geometry_restraints_manager = model.get_restraints_manager().geometry
    self.rs_values = [2.0, 3.8, 5.2, 5.5, 6.2, 7.0, 8.6, 10.0]
    self.Rj = [2.1, 2.2, 2.5]
    self.cutoff = 8.1
    self.ts_values = [0.392699, 1.178097, 1.963495, 2.748894]
    self.angular_rs_values = [3.8, 5.2, 5.5, 6.2]
    self.angular_zeta = 8
    self.EAEVs = format_class()
    self.MAEVs = format_class()
    self.BAEVs = format_class()
    self.center_atom     = None
    self.chain_hierarchy = None
    self.generate_AEV()

  def get_values(self):
    result = OrderedDict()
    result['B'] = self.BAEVs.values()[0]
    result['M'] = self.MAEVs.values()[5]
    result['E'] = self.EAEVs.values()[-1]
    return result

  def generate_ca(self):
    """
    ???
    """
    protain_fragments = generate_protein_fragments(
      hierarchy = self.chain_hierarchy,
      geometry = self.geometry_restraints_manager,
      include_non_standard_peptides=True,
      length=5)
    for five in protain_fragments:
      rc = []
      for residue in five:
        for atom in residue.atoms():
          if atom.name == ' CA ':
            rc.append(atom)
      if len(rc) == 5:
        yield rc

  def generate_AEV(self):
    """
    Traversal the C-alpha atoms list and generate AEV values of every chain's
    every atom. Every C-alpha atom has 3 direction AEVs: forward(BAEVs),
    backward(EAEVs) and all direction(MAEVS).
    """
    chain_list = []
    for chain in self.hierarchy.chains():
      chain_list.append(chain.id)
    chain_list = list(set(chain_list))
    for chain in chain_list:
      chain_hierarchy = self.hierarchy.deep_copy()
      asc = chain_hierarchy.atom_selection_cache()
      sel = asc.selection("chain " + chain)
      self.chain_hierarchy = chain_hierarchy.select(sel)
      begin = next(self.generate_ca())
      self.center_atom = begin[0]
      self.EAEVs.update(self.calculate([]))
      self.MAEVs.update(self.calculate(begin))
      self.center_atom = begin[1]
      self.EAEVs.update(self.calculate(begin[0:1]))
      self.MAEVs.update(self.calculate(begin))
      self.center_atom = begin[2]
      self.EAEVs.update(self.calculate(begin[0:2]))
      self.center_atom = begin[3]
      self.EAEVs.update(self.calculate(begin[0:3]))
      for atomlist in self.generate_ca():
        end = atomlist
        self.center_atom = atomlist[0]
        self.BAEVs.update(self.calculate(atomlist))
        self.center_atom = atomlist[-1]
        self.EAEVs.update(self.calculate(atomlist))
        self.center_atom = atomlist[2]
        self.MAEVs.update(self.calculate(atomlist))
      self.center_atom = end[1]
      self.BAEVs.update(self.calculate(end[1:]))
      self.center_atom = end[2]
      self.BAEVs.update(self.calculate(end[2:]))
      self.center_atom = end[3]
      self.BAEVs.update(self.calculate(end[3:]))
      self.MAEVs.update(self.calculate(end))
      self.center_atom = end[4]
      self.BAEVs.update(self.calculate([]))
      self.MAEVs.update(self.calculate(end))

  def cutf(self, distance):
    """
    Formula number ???, page number ???
    """
    if distance <= self.cutoff:
      Fc = 0.5 * math.cos(math.pi * distance / self.cutoff) + 0.5
    else:
      Fc = 0
    return Fc

  def calculate(self, atom_list):
    """
    Formula number ???, page number ???
    """
    n = 4.0
    l = 8.0
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
