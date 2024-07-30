# TODO reduce to one outlier per residue
# CDL on by default?

from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.validation import utils
from mmtbx import monomer_library
from mmtbx.validation import validation
from mmtbx.validation import residue
from mmtbx.validation import atoms
from mmtbx.validation import get_atoms_info
from scitbx.array_family import flex
from cctbx import geometry_restraints
from libtbx.str_utils import make_sub_header
from libtbx import slots_getstate_setstate
from math import sqrt
import sys
import json

# individual validation results
class mp_bond(atoms):
  __slots__ = atoms.__slots__ + ["sigma", "delta", "target", "distance_value", "macromolecule_type"]

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

class mp_angle(atoms):
  __slots__ = atoms.__slots__ + ["sigma", "delta", "target", "angle_value", "macromolecule_type"]

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

# analysis objects
class mp_bonds(validation):
  output_header = "#residue:atom_1:atom_2:num_sigmas"
  label = "Bond lengths"
  gui_list_headers = ["Residue", "Atom 1", "Atom 2", "Sigmas"]
  gui_formats = ["%s", "%s", "%s", "%.2f"]
  wx_column_widths = [160] * 4
  def __init__(self, pdb_hierarchy, pdb_atoms, geometry_restraints_manager,
                outliers_only=True):
    validation.__init__(self)
    self.n_outliers_large_by_model = {}
    self.n_outliers_small_by_model = {}
    self.n_outliers_protein_by_model = {}
    self.n_outliers_na_by_model = {}
    self.n_outliers_other_by_model = {}
    self.n_total_protein_by_model = {}
    self.n_total_na_by_model = {}
    self.n_total_other_by_model = {}
    for m in pdb_hierarchy.models():
      self.n_total_by_model[m.id] = 0
      self.n_outliers_by_model[m.id] = 0
      self.n_outliers_small_by_model[m.id] = 0
      self.n_outliers_large_by_model[m.id] = 0
      self.n_total_protein_by_model[m.id] = 0
      self.n_outliers_protein_by_model[m.id] = 0
      self.n_total_na_by_model[m.id] = 0
      self.n_outliers_na_by_model[m.id] = 0
      self.n_total_other_by_model[m.id] = 0
      self.n_outliers_other_by_model[m.id] = 0
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
      self.n_total += 1
      #iotbx.pdb.common_residue_names_get_class
      self.n_total_by_model[model_id] += 1
      mm_type = utils.get_mmtype_from_resname(pdb_atoms[proxy.i_seqs[0]].parent().resname)
      if mm_type=="PROTEIN":
        self.n_total_protein_by_model[model_id] += 1
      elif mm_type=="NA":
        self.n_total_na_by_model[model_id] += 1
      else:
        self.n_total_other_by_model[model_id] += 1
      sigma = sqrt(1 / restraint.weight)
      num_sigmas = - restraint.delta / sigma
      is_outlier = (abs(num_sigmas) >= cutoff)
      if is_outlier:
        self.n_outliers += 1
        self.n_outliers_by_model[model_id] += 1
        if mm_type == "PROTEIN":
          self.n_outliers_protein_by_model[model_id] += 1
        elif mm_type == "NA":
          self.n_outliers_na_by_model[model_id] += 1
        else:
          self.n_outliers_other_by_model[model_id] += 1
        if num_sigmas < 0:
          self.n_outliers_small_by_model[model_id] += 1
        else:
          self.n_outliers_large_by_model[model_id] += 1
      if (is_outlier or not outliers_only):
        self.results.append(mp_bond(
          atoms_info=get_atoms_info(pdb_atoms, proxy.i_seqs),
          target=restraint.distance_ideal,
          distance_value=restraint.distance_model,
          sigma=sigma,
          score=num_sigmas,
          delta=restraint.delta,
          xyz=flex.vec3_double([pdb_atoms[proxy.i_seqs[0]].xyz, pdb_atoms[proxy.i_seqs[1]].xyz]).mean(),
          outlier=is_outlier,
          macromolecule_type=mm_type))

  def get_result_class(self) : return mp_bond

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "mp_bonds"
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
                                 "num_outliers_too_large": self.n_outliers_large_by_model[mod_id],
                                 "num_total_protein": self.n_total_protein_by_model[mod_id],
                                 "num_outliers_protein": self.n_outliers_protein_by_model[mod_id],
                                 "num_total_na": self.n_total_na_by_model[mod_id],
                                 "num_outliers_na": self.n_outliers_na_by_model[mod_id],
                                 "num_total_other": self.n_total_other_by_model[mod_id],
                                 "num_outliers_other": self.n_outliers_other_by_model[mod_id]}
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def show_summary(self, out=sys.stdout, prefix=""):
    if (self.n_total == 0):
      print(prefix + "No bond lengths found.", file=out)
    elif (self.n_outliers == 0):
      print(prefix + "All bonds within 4.0 sigma of ideal values.", file=out)
    else :
      print(prefix + "%d/%d bond outliers present" % (self.n_outliers,
        self.n_total), file=out)

class mp_angles(validation):
  output_header = "#residue:atom_1:atom_2:atom_3:num_sigmas"
  label = "Bond angles"
  gui_list_headers = ["Residue", "Atom 1", "Atom 2", "Atom 3", "Sigmas"]
  gui_formats = ["%s", "%s", "%s", "%s", "%.2f"]
  wx_column_widths = [160] * 5
  def __init__(self, pdb_hierarchy, pdb_atoms, geometry_restraints_manager,
                outliers_only=True):
    validation.__init__(self)
    self.n_outliers_large_by_model = {}
    self.n_outliers_small_by_model = {}
    self.n_outliers_protein_by_model = {}
    self.n_outliers_na_by_model = {}
    self.n_outliers_other_by_model = {}
    self.n_total_protein_by_model = {}
    self.n_total_na_by_model = {}
    self.n_total_other_by_model = {}
    for m in pdb_hierarchy.models():
      self.n_total_by_model[m.id] = 0
      self.n_outliers_by_model[m.id] = 0
      self.n_outliers_small_by_model[m.id] = 0
      self.n_outliers_large_by_model[m.id] = 0
      self.n_total_protein_by_model[m.id] = 0
      self.n_outliers_protein_by_model[m.id] = 0
      self.n_total_na_by_model[m.id] = 0
      self.n_outliers_na_by_model[m.id] = 0
      self.n_total_other_by_model[m.id] = 0
      self.n_outliers_other_by_model[m.id] = 0
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
      self.n_total += 1
      self.n_total_by_model[model_id] += 1
      mm_type = utils.get_mmtype_from_resname(pdb_atoms[proxy.i_seqs[0]].parent().resname)
      if mm_type=="PROTEIN":
        self.n_total_protein_by_model[model_id] += 1
      elif mm_type=="NA":
        self.n_total_na_by_model[model_id] += 1
      else:
        self.n_total_other_by_model[model_id] += 1
      sigma = sqrt(1 / restraint.weight)
      num_sigmas = - restraint.delta / sigma
      is_outlier = (abs(num_sigmas) >= cutoff)
      if is_outlier:
        self.n_outliers += 1
        self.n_outliers_by_model[model_id] += 1
        if mm_type == "PROTEIN":
          self.n_outliers_protein_by_model[model_id] += 1
        elif mm_type == "NA":
          self.n_outliers_na_by_model[model_id] += 1
        else:
          self.n_outliers_other_by_model[model_id] += 1
        if num_sigmas < 0:
          self.n_outliers_small_by_model[model_id] += 1
        else:
          self.n_outliers_large_by_model[model_id] += 1

      if (is_outlier or not outliers_only):
        self.results.append(mp_angle(
          atoms_info=get_atoms_info(pdb_atoms, proxy.i_seqs),
          target=restraint.angle_ideal,
          angle_value=restraint.angle_model,
          sigma=sigma,
          score=num_sigmas,
          delta=restraint.delta,
          xyz=pdb_atoms[proxy.i_seqs[1]].xyz,
          outlier=is_outlier,
          macromolecule_type=mm_type))

  def get_result_class(self) : return mp_angle

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "mp_angles"
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
                                 "num_outliers_too_large": self.n_outliers_large_by_model[mod_id],
                                 "num_total_protein": self.n_total_protein_by_model[mod_id],
                                 "num_outliers_protein": self.n_outliers_protein_by_model[mod_id],
                                 "num_total_na": self.n_total_na_by_model[mod_id],
                                 "num_outliers_na": self.n_outliers_na_by_model[mod_id],
                                 "num_total_other": self.n_total_other_by_model[mod_id],
                                 "num_outliers_other": self.n_outliers_other_by_model[mod_id]}
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def show_summary(self, out=sys.stdout, prefix=""):
    if (self.n_total == 0):
      print(prefix + "No bond angles found.", file=out)
    elif (self.n_outliers == 0):
      print(prefix + "All angles within 4.0 sigma of ideal values.", file=out)
    else :
      print(prefix + "%d/%d angle outliers present" % (self.n_outliers,
        self.n_total), file=out)

class mp_validate_bonds(slots_getstate_setstate):
  __slots__ = ["bonds", "angles"]

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
    self.bonds = mp_bonds(
      pdb_hierarchy=pdb_hierarchy,
      pdb_atoms=pdb_atoms,
      geometry_restraints_manager=geometry_restraints_manager,
      outliers_only=outliers_only)
    self.angles = mp_angles(
      pdb_hierarchy=pdb_hierarchy,
      pdb_atoms=pdb_atoms,
      geometry_restraints_manager=geometry_restraints_manager,
      outliers_only=outliers_only)

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
    mp_json = addon_json
    for slot in self.__slots__:
      slot_json = json.loads(getattr(self, slot).as_JSON())
      mp_json["mp_"+slot] = slot_json
    return json.dumps(mp_json, indent=2)
