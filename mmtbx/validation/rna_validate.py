
# TODO reduce to one outlier per residue

from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import rna_sugar_pucker_analysis
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.validation import utils
from mmtbx import monomer_library
from mmtbx.validation import rna_geometry
from mmtbx.validation import residue
from mmtbx.validation import atoms
from mmtbx.validation import get_atoms_info
from iotbx.pdb import common_residue_names_get_class as get_res_class
from cctbx import geometry_restraints
from scitbx.array_family import flex
from libtbx.str_utils import make_sub_header, format_value
from libtbx import slots_getstate_setstate
from math import sqrt
from mmtbx.suitename import suitealyze
import sys
import json

rna_backbone_atoms = set([
  "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C1'",
  "C3'", "O3'", "C2'", "O2'", "N1", "N9" ]) #version 3.x naming

# individual validation results
class rna_bond(atoms):
  __slots__ = atoms.__slots__ + ["sigma", "delta"]

  @staticmethod
  def header():
    return "%-20s  %6s  %6s  %6s" % ("residue", "atom 1", "atom 2", "sigmas")

  def id_str(self, spacer=" "):
    return "%s%s%s" % (self.atoms_info[0].id_str(), spacer,
      self.atoms_info[1].id_str())

  def format_values(self):
    return "%-20s  %6s  %6s  %6.2f" % (self.id_str(), self.atoms_info[0].name,
      self.atoms_info[1].name, self.score)

  def as_string(self, prefix=""):
    return prefix + self.format_values()

  def as_JSON(self):
    atoms_dict = {}
    for s in self.atoms_info[0].__slots__:
      atoms_dict["atoms_"+s] = [getattr(self.atoms_info[0], s), getattr(self.atoms_info[1], s)]
    serializable_slots = [s for s in self.__slots__ if s != 'atoms_info' and hasattr(self, s)]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    return json.dumps(self.merge_two_dicts(slots_as_dict, atoms_dict), indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  def as_table_row_phenix(self):
    return [ self.id_str(), self.atoms_info[0].name, self.atoms_info[1].name,
             self.score ]

class rna_angle(atoms):
  __slots__ = atoms.__slots__ + ["sigma", "delta"]

  @staticmethod
  def header():
    return "%-20s  %6s  %6s  %6s  %6s" % ("residue", "atom 1", "atom 2",
      "atom 3", "sigmas")

  def id_str(self, spacer=" "):
    return "%s%s%s%s%s" % (self.atoms_info[0].id_str(), spacer,
      self.atoms_info[1].id_str(), spacer,
      self.atoms_info[2].id_str())

  def format_values(self):
    return "%-20s  %6s  %6s  %6s  %6.2f" % (self.id_str(),
      self.atoms_info[0].name, self.atoms_info[1].name,
      self.atoms_info[2].name, self.score)

  def as_string(self, prefix=""):
    return prefix + self.format_values()

  def as_JSON(self):
    atoms_dict = {}
    for s in self.atoms_info[0].__slots__:
      atoms_dict["atoms_"+s] = [getattr(self.atoms_info[0], s), getattr(self.atoms_info[1], s), getattr(self.atoms_info[2], s)]
    serializable_slots = [s for s in self.__slots__ if s != 'atoms_info' and hasattr(self, s)]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    return json.dumps(self.merge_two_dicts(slots_as_dict, atoms_dict), indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  def as_table_row_phenix(self):
    return [ self.id_str(), self.atoms_info[0].name, self.atoms_info[1].name,
             self.atoms_info[2].name, self.score ]

class rna_pucker(residue):
  """
  Validation using pucker-specific restraints library.
  """
  __slots__ = residue.__slots__ + [
    "delta_angle",
    "is_delta_outlier",
    "epsilon_angle",
    "is_epsilon_outlier",
    "Pperp_distance",
    "probable_pucker",
    "model_id",
  ]

  @staticmethod
  def header():
    return "%-20s  %8s  %8s  %8s  %8s" % ("residue", "delta", "outlier",
      "epsilon", "outlier")

  def format_values(self):
    def format_outlier_flag(flag):
      if (flag) : return "yes"
      else : return "no"
    def format_angle(val):
      return format_value("%8.1f", val, replace_none_with="---")
    return "%-20s  %8s  %8s  %8s  %8s" % (self.id_str(),
      format_angle(self.delta_angle),
      format_outlier_flag(self.is_delta_outlier),
      format_angle(self.epsilon_angle),
      format_outlier_flag(self.is_epsilon_outlier))

  def as_string(self, prefix=""):
    return prefix + self.format_values()

  def as_JSON(self):
    serializable_slots = [s for s in self.__slots__ if hasattr(self, s)]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    return json.dumps(slots_as_dict, indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  def as_table_row_phenix(self):
    return [ self.id_str(), self.delta_angle, self.epsilon_angle ]

# analysis objects
class rna_bonds(rna_geometry):
  output_header = "#residue:atom_1:atom_2:num_sigmas"
  label = "Backbone bond lenths"
  gui_list_headers = ["Residue", "Atom 1", "Atom 2", "Sigmas"]
  gui_formats = ["%s", "%s", "%s", "%.2f"]
  wx_column_widths = [160] * 4
  def __init__(self, pdb_hierarchy, pdb_atoms, geometry_restraints_manager,
                outliers_only=True):
    self.n_outliers_large_by_model = {}
    self.n_outliers_small_by_model = {}
    rna_geometry.__init__(self)
    cutoff = 4
    sites_cart = pdb_atoms.extract_xyz()
    flags = geometry_restraints.flags.flags(default=True)
    pair_proxies = geometry_restraints_manager.pair_proxies(
      flags=flags,
      sites_cart=sites_cart)
    bond_proxies = pair_proxies.bond_proxies
    for proxy in bond_proxies.simple:
      restraint = geometry_restraints.bond(
        sites_cart=sites_cart,
        proxy=proxy)
      atom1 = pdb_atoms[proxy.i_seqs[0]].name
      atom2 = pdb_atoms[proxy.i_seqs[1]].name
      labels = pdb_atoms[proxy.i_seqs[0]].fetch_labels()
      model_id = labels.model_id
      if model_id not in self.n_total_by_model.keys():
        self.n_total_by_model[model_id] = 0
        self.n_outliers_by_model[model_id] = 0
        self.n_outliers_small_by_model[model_id] = 0
        self.n_outliers_large_by_model[model_id] = 0
      if (atom1.strip() not in rna_backbone_atoms or
          atom2.strip() not in rna_backbone_atoms):
        continue
      self.n_total += 1
      self.n_total_by_model[model_id] += 1
      sigma = sqrt(1 / restraint.weight)
      num_sigmas = - restraint.delta / sigma
      is_outlier = (abs(num_sigmas) >= cutoff)
      if is_outlier:
        self.n_outliers += 1
        self.n_outliers_by_model[model_id] += 1
        if num_sigmas < 0:
          self.n_outliers_small_by_model[model_id] += 1
        else:
          self.n_outliers_large_by_model[model_id] += 1
      if (is_outlier or not outliers_only):
        self.results.append(rna_bond(
          atoms_info=get_atoms_info(pdb_atoms, proxy.i_seqs),
          sigma=sigma,
          score=num_sigmas,
          delta=restraint.delta,
          xyz=flex.vec3_double([pdb_atoms[proxy.i_seqs[0]].xyz, pdb_atoms[proxy.i_seqs[1]].xyz]).mean(),
          outlier=is_outlier))

  def get_result_class(self) : return rna_bond

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "rna_bonds"
    data = addon_json
    flat_results = []
    hierarchical_results = {}
    summary_results = {}
    for result in self.results:
      flat_results.append(json.loads(result.as_JSON()))
      hier_result = json.loads(result.as_hierarchical_JSON())
      hierarchical_results = self.merge_dict(hierarchical_results, hier_result)

    data['flat_results'] = flat_results
    data['hierarchical_results'] = hierarchical_results
    for mod_id in self.n_total_by_model.keys():
      summary_results[mod_id] = {"num_outliers": self.n_outliers_by_model[mod_id],
                                 "num_total": self.n_total_by_model[mod_id],
                                 "num_outliers_too_small": self.n_outliers_small_by_model[mod_id],
                                 "num_outliers_too_large": self.n_outliers_large_by_model[mod_id]}
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def show_summary(self, out=sys.stdout, prefix=""):
    if (self.n_total == 0):
      print(prefix + "No RNA backbone atoms found.", file=out)
    elif (self.n_outliers == 0):
      print(prefix + "All bonds within 4.0 sigma of ideal values.", file=out)
    else :
      print(prefix + "%d/%d bond outliers present" % (self.n_outliers,
        self.n_total), file=out)

class rna_angles(rna_geometry):
  output_header = "#residue:atom_1:atom_2:atom_3:num_sigmas"
  label = "Backbone bond angles"
  gui_list_headers = ["Residue", "Atom 1", "Atom 2", "Atom 3", "Sigmas"]
  gui_formats = ["%s", "%s", "%s", "%s", "%.2f"]
  wx_column_widths = [160] * 5
  def __init__(self, pdb_hierarchy, pdb_atoms, geometry_restraints_manager,
                outliers_only=True):
    self.n_outliers_large_by_model = {}
    self.n_outliers_small_by_model = {}
    rna_geometry.__init__(self)
    cutoff = 4
    sites_cart = pdb_atoms.extract_xyz()
    flags = geometry_restraints.flags.flags(default=True)
    i_seq_name_hash = utils.build_name_hash(pdb_hierarchy=pdb_hierarchy)
    for proxy in geometry_restraints_manager.angle_proxies:
      restraint = geometry_restraints.angle(
        sites_cart=sites_cart,
        proxy=proxy)
      atom1 = pdb_atoms[proxy.i_seqs[0]].name
      atom2 = pdb_atoms[proxy.i_seqs[1]].name
      atom3 = pdb_atoms[proxy.i_seqs[2]].name
      labels = pdb_atoms[proxy.i_seqs[0]].fetch_labels()
      model_id = labels.model_id
      if model_id not in self.n_total_by_model.keys():
        self.n_total_by_model[model_id] = 0
        self.n_outliers_by_model[model_id] = 0
        self.n_outliers_small_by_model[model_id] = 0
        self.n_outliers_large_by_model[model_id] = 0
      if (atom1.strip() not in rna_backbone_atoms or
          atom2.strip() not in rna_backbone_atoms or
          atom3.strip() not in rna_backbone_atoms):
        continue
      self.n_total += 1
      self.n_total_by_model[model_id] += 1
      sigma = sqrt(1 / restraint.weight)
      num_sigmas = - restraint.delta / sigma
      is_outlier = (abs(num_sigmas) >= cutoff)
      if is_outlier:
        self.n_outliers += 1
        self.n_outliers_by_model[model_id] += 1
        if num_sigmas < 0:
          self.n_outliers_small_by_model[model_id] += 1
        else:
          self.n_outliers_large_by_model[model_id] += 1
      if (is_outlier or not outliers_only):
        self.results.append(rna_angle(
          atoms_info=get_atoms_info(pdb_atoms, proxy.i_seqs),
          sigma=sigma,
          score=num_sigmas,
          delta=restraint.delta,
          xyz=pdb_atoms[proxy.i_seqs[1]].xyz,
          outlier=is_outlier))

  def get_result_class(self) : return rna_angle

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "rna_angles"
    data = addon_json
    flat_results = []
    hierarchical_results = {}
    summary_results = {}
    for result in self.results:
      flat_results.append(json.loads(result.as_JSON()))
      hier_result = json.loads(result.as_hierarchical_JSON())
      hierarchical_results = self.merge_dict(hierarchical_results, hier_result)

    data['flat_results'] = flat_results
    data['hierarchical_results'] = hierarchical_results
    for mod_id in self.n_total_by_model.keys():
      summary_results[mod_id] = {"num_outliers": self.n_outliers_by_model[mod_id],
                                 "num_total": self.n_total_by_model[mod_id],
                                 "num_outliers_too_small": self.n_outliers_small_by_model[mod_id],
                                 "num_outliers_too_large": self.n_outliers_large_by_model[mod_id]}
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def show_summary(self, out=sys.stdout, prefix=""):
    if (self.n_total == 0):
      print(prefix + "No RNA backbone atoms found.", file=out)
    elif (self.n_outliers == 0):
      print(prefix + "All angles within 4.0 sigma of ideal values.", file=out)
    else :
      print(prefix + "%d/%d angle outliers present" % (self.n_outliers,
        self.n_total), file=out)

class rna_puckers(rna_geometry):
  __slots__ = rna_geometry.__slots__ + [
    "pucker_states",
    "pucker_perp_xyz",
    "pucker_dist",
  ]
  output_header = "#residue:delta_angle:is_delta_outlier:epsilon_angle:is_epsilon_outler"
  label = "Sugar pucker"
  gui_list_headers = ["Residue", "Delta", "Epsilon"]
  gui_formats = ["%s", "%.2f", "%.2f"]
  wx_column_widths = [200]*3
  def __init__(self, pdb_hierarchy, params=None, outliers_only=True):
    if (params is None):
      params = rna_sugar_pucker_analysis.master_phil.extract()
    self.pucker_states = []
    self.pucker_perp_xyz = {}
    self.pucker_dist = {}
    rna_geometry.__init__(self)
    from iotbx.pdb.rna_dna_detection import residue_analysis
    for model in pdb_hierarchy.models():
      self.n_outliers_by_model[model.id] = 0
      self.n_total_by_model[model.id] = 0
      for chain in model.chains():
        first_altloc = [conformer.altloc for conformer in chain.conformers()][0]
        #can skip some calculations on later altlocs
        for conformer in chain.conformers():
          residues = conformer.residues()
          for i_residue,residue in enumerate(residues):
            def _get_next_residue():
              j = i_residue + 1
              if (j == len(residues)): return None
              return residues[j]
            ra1 = residue_analysis(
              residue_atoms=residue.atoms(),
              distance_tolerance=params.bond_detection_distance_tolerance)
            if (ra1.problems is not None): continue
            if (not ra1.is_rna): continue
            residue_2_p_atom = None
            next_pdb_residue = _get_next_residue()
            if (next_pdb_residue is not None):
              residue_2_p_atom = next_pdb_residue.find_atom_by(name=" P  ")
            if (get_res_class(residue.resname) != "common_rna_dna"):
              continue
            if conformer.altloc.strip() == '':
              local_altloc = ''
            else:
              local_altloc = self.local_altloc_from_atoms(ra1.deoxy_ribo_atom_dict, ra1.c1p_outbound_atom, residue_2_p_atom)
            if local_altloc == '' and local_altloc != conformer.altloc and conformer.altloc != first_altloc:
              #if the pucker atoms contain no alternates, then the calculation can be skipped if this isn't the first
              #  conformer to be encountered
              continue
            ana = rna_sugar_pucker_analysis.evaluate(
              params=params,
              residue_1_deoxy_ribo_atom_dict=ra1.deoxy_ribo_atom_dict,
              residue_1_c1p_outbound_atom=ra1.c1p_outbound_atom,
              residue_2_p_atom=residue_2_p_atom)
            self.pucker_states.append(ana)
            self.n_total += 1
            self.n_total_by_model[model.id] += 1
            is_outlier = ana.is_delta_outlier or ana.is_epsilon_outlier
            if (is_outlier):
              self.n_outliers += 1
              self.n_outliers_by_model[model.id] += 1
            if (is_outlier or not outliers_only):
              pucker = rna_pucker(
                model_id=model.id,
                chain_id=chain.id,
                resseq=residue.resseq,
                icode=residue.icode,
                altloc=local_altloc, #'' if none
                resname=residue.resname,
                delta_angle=ana.delta,
                is_delta_outlier=ana.is_delta_outlier,
                epsilon_angle=ana.epsilon,
                is_epsilon_outlier=ana.is_epsilon_outlier,
                xyz=self.get_sugar_xyz_mean(residue),
                outlier=is_outlier)
              self.results.append(pucker)
              key = pucker.id_str() #[8:-1]
              self.pucker_perp_xyz[key] = [ana.p_perp_xyz, ana.o3p_perp_xyz]
              self.pucker_dist[key] = [ana.p_distance_c1p_outbound_line,
                                       ana.o3p_distance_c1p_outbound_line]
              pucker.Pperp_distance=ana.p_distance_c1p_outbound_line
              if not pucker.Pperp_distance:
                pucker.probable_pucker=None
              elif pucker.Pperp_distance and pucker.Pperp_distance < 2.9:
                pucker.probable_pucker="C2'-endo"
              else:
                pucker.probable_pucker="C3'-endo"

  def get_sugar_xyz_mean(self, residue):
    sugar_atoms = [" C1'", " C2'", " C3'", " C4'", " O4'"]
    atom_xyzs = []
    sums = [0]*3
    for at_name in sugar_atoms:
      sug_atom = residue.find_atom_by(name = at_name)
      if sug_atom is not None:
        atom_xyzs.append(sug_atom.xyz)
    assert len(atom_xyzs) > 0, "RNA sugar pucker validation found zero sugar atoms, probably due to non-standard atom names"
    mean = flex.vec3_double(atom_xyzs).mean()
    return mean

  def get_result_class(self) : return rna_pucker

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "rna_puckers"
    data = addon_json
    flat_results = []
    hierarchical_results = {}
    summary_results = {}
    for result in self.results:
      flat_results.append(json.loads(result.as_JSON()))
      hier_result = json.loads(result.as_hierarchical_JSON())
      hierarchical_results = self.merge_dict(hierarchical_results, hier_result)

    data['flat_results'] = flat_results
    data['hierarchical_results'] = hierarchical_results
    for mod_id in self.n_total_by_model.keys():
      summary_results[mod_id] = {"num_outliers": self.n_outliers_by_model[mod_id],
                                 "num_residues": self.n_total_by_model[mod_id]}
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def show_summary(self, out=sys.stdout, prefix=""):
    if (self.n_total == 0):
      print(prefix + "No RNA sugar groups found.", file=out)
    elif (self.n_outliers == 0):
      print(prefix + "All puckers have reasonable geometry.", file=out)
    else :
      print(prefix + "%d/%d pucker outliers present" % (self.n_outliers,
        self.n_total), file=out)

  def local_altloc_from_atoms(self, residue_1_deoxy_ribo_atom_dict, residue_1_c1p_outbound_atom, residue_2_p_atom):
    #conformer.altloc masks whether a residue has true alternate conformations
    #check atom.parent().altloc for whether any atoms have alt positions
    #only run this if conformer.altloc != ''
    for atom in [residue_1_c1p_outbound_atom, residue_2_p_atom]:
      if atom is None: continue
      altloc = atom.parent().altloc
      if altloc != "":
        return altloc
    for atom in residue_1_deoxy_ribo_atom_dict.values():
      if atom is None: continue
      altloc = atom.parent().altloc
      if altloc != "":
        return altloc
    return ""

class rna_validation(slots_getstate_setstate):
  __slots__ = ["bonds", "angles", "puckers", "suites"]

  def __init__(self,
      pdb_hierarchy,
      geometry_restraints_manager=None,
      params=None,
      outliers_only=True):
    if (geometry_restraints_manager is None):
      mon_lib_srv = monomer_library.server.server()
      ener_lib = monomer_library.server.ener_lib()
      processed_pdb_file = pdb_interpretation.process(
        mon_lib_srv=mon_lib_srv,
        ener_lib=ener_lib,
        pdb_hierarchy=pdb_hierarchy,
        substitute_non_crystallographic_unit_cell_if_necessary=True)
      geometry_restraints_manager = \
        processed_pdb_file.geometry_restraints_manager()
      pdb_hierarchy = \
        processed_pdb_file.all_chain_proxies.pdb_hierarchy
    pdb_atoms = pdb_hierarchy.atoms()
    self.bonds = rna_bonds(
      pdb_hierarchy=pdb_hierarchy,
      pdb_atoms=pdb_atoms,
      geometry_restraints_manager=geometry_restraints_manager,
      outliers_only=outliers_only)
    self.angles = rna_angles(
      pdb_hierarchy=pdb_hierarchy,
      pdb_atoms=pdb_atoms,
      geometry_restraints_manager=geometry_restraints_manager,
      outliers_only=outliers_only)
    self.puckers = rna_puckers(
      pdb_hierarchy=pdb_hierarchy,
      params=getattr(params, "rna_sugar_pucker_analysis", None),
      outliers_only=outliers_only)
    self.suites = suitealyze.suitealyze(
      pdb_hierarchy=pdb_hierarchy,
      outliers_only=outliers_only)
    #self.suites = rna_suites(
    #  pdb_hierarchy=pdb_hierarchy,
    #  geometry_restraints_manager=geometry_restraints_manager,
    #  outliers_only=outliers_only)

  def show_summary(self, out=sys.stdout, prefix=""):
    pass

  def show(self, out=sys.stdout, prefix="", outliers_only=None,
      verbose=True):
    for geo_type in self.__slots__ :
      rv = getattr(self, geo_type)
      if (rv.n_outliers > 0) or (not outliers_only):
        make_sub_header(rv.label, out=out)
        rv.show(out=out)

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    rna_json = addon_json
    for slot in self.__slots__:
      slot_json = json.loads(getattr(self, slot).as_JSON())
      rna_json["rna_"+slot] = slot_json
    return json.dumps(rna_json, indent=2)

