from __future__ import absolute_import, division, print_function
import math
import copy
import mmtbx.model
from libtbx.utils import null_out
from collections import OrderedDict
from scitbx.array_family import flex
from mmtbx.conformation_dependent_library import generate_protein_fragments
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import uctbx

def format_HELIX_records_from_AEV(aev_values_dict, cc_cutoff):
  """
  Geting HELIX records for target model using AEV values and specified cutoff.
  """
  cc_cutoff = float(cc_cutoff)
  start = []
  end = []
  M = 0
  i = 0
  result = []
  length = 0
  for key,value in aev_values_dict.items():
    if start==[] and value['B'] > cc_cutoff and value['M'] > cc_cutoff - 0.1:
      start = key
    elif start and value['M'] > cc_cutoff:
      M = 1
      if value['E'] > cc_cutoff:
        end = key
    elif start and M == 1 and value['E'] > cc_cutoff:
      end = key
      length = int(end.split()[-1]) - int(start.split()[-1]) + 1
    if value['B'] < 0.6 and end and M == 1 and length > 3:
      i += 1
      fmt = "HELIX   {0:>2}  {0:>2} {1:>}  {2:>}   {3:>36}"
      length = int(end.split()[-1]) - int(start.split()[-1]) + 1
      result.append(fmt.format(i, start, end, length))
      start = []
      end = []
      M = 0
      length = 0
  return result

def generate_perfect_helix(rs_values,
                           ts_values,
                           angular_rs_values,
                           radial_eta,
                           angular_eta,
                           angular_zeta,
                           n_residues=10,
                           residue_code="G",
                           ):
  """
  Compute AEV values for the perfect helix.
  """
  perfect_helix_ph = ssb.secondary_structure_from_sequence(ssb.alpha_helix_str,
    residue_code*n_residues)
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = perfect_helix_ph.atoms().extract_xyz(),
    buffer_layer = 5)
  model = mmtbx.model.manager(
    model_input      = None,
    crystal_symmetry = box.crystal_symmetry(),
    pdb_hierarchy    = perfect_helix_ph,
    build_grm        = True,
    log              = null_out())
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
      result[key1].setdefault(key, None)
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
  return result

# It is format calss of corelation coefficient values.
class diff_class(OrderedDict):
  def __repr__(self):
    outl = '...\n'
    for key, item in self.items():
      outl += '  %s :' % (key)
      for key1,value in item.items():
        outl += ' %s: '%key1
        if value is None:
          outl += 'None, '
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
        # print (v)
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
    self.geometry_restraints_manager = model.get_restraints_manager().geometry
    self.rs_values = rs_values
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

  def generate_ca(self, length=5):
    """
    Geting all C-alpha atoms from an whole pdb file.
    """
    # faster with atom selection to CA re CJW update
    protain_fragments = generate_protein_fragments(
      hierarchy = self.chain_hierarchy,
      geometry = self.geometry_restraints_manager,
      include_non_linked=True,
      include_non_standard_peptides=True,
      length=length)
    for five in protain_fragments:
      rc = []
      for residue in five:
        for atom in residue.atoms():
          if atom.name == ' CA ':
            rc.append(atom)
      if len(rc) == length:
        yield rc

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
      sel = asc.selection("chain " + chain)
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
          GmA = 10 * GmA * (2**(1-l))
          AEVs[res_name].append(GmA)
    return AEVs
