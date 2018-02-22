
"""
Module for managing and manipulating all properties (and underlying objects)
associated with an atomic (or pseudo-atomic) model, including both basic
attributes stored in a PDB file, scattering information, and geometry
restraints.
"""

from __future__ import division

from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry, user_plus_sys_time, null_out
from libtbx import group_args, str_utils

import iotbx.pdb
import iotbx.cif.model
from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter
from iotbx.pdb.atom_selection import AtomSelectionError
from iotbx.pdb.misc_records_output import link_record_output
from iotbx.cif import category_sort_function
from iotbx.pdb.mmcif import tls_as_cif_block

from cctbx.array_family import flex
from cctbx import xray
from cctbx import adptbx
from cctbx import geometry_restraints
from cctbx import adp_restraints
from cctbx import crystal

import mmtbx.restraints
import mmtbx.hydrogens
from mmtbx.hydrogens import riding
import mmtbx.model.statistics
import mmtbx.monomer_library.server
import mmtbx.geometry_restraints.torsion_restraints.utils as torsion_utils
import mmtbx.tls.tools as tls_tools
from mmtbx import ias
from mmtbx import utils
from mmtbx import ncs
from mmtbx.command_line import find_tls_groups
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from mmtbx.geometry_restraints.torsion_restraints.reference_model import \
    add_reference_dihedral_restraints_if_requested, reference_model_str
from mmtbx.geometry_restraints.torsion_restraints.torsion_ncs import torsion_ncs
from mmtbx.refinement import print_statistics
from mmtbx.refinement import anomalous_scatterer_groups
from mmtbx.refinement import geometry_minimization

import boost.python
ext = boost.python.import_ext("mmtbx_validation_ramachandran_ext")
from mmtbx_validation_ramachandran_ext import rama_eval
from mmtbx.rotamer.rotamer_eval import RotamerEval
from mmtbx.rotamer.rotamer_eval import RotamerID

ext2 = boost.python.import_ext("iotbx_pdb_hierarchy_ext")
from iotbx_pdb_hierarchy_ext import *

from cStringIO import StringIO
from copy import deepcopy
import sys
import math

time_model_show = 0.0

def find_common_water_resseq_max(pdb_hierarchy):
  get_class = iotbx.pdb.common_residue_names_get_class
  hy36decode = iotbx.pdb.hy36decode
  result = None
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if (get_class(name=ag.resname) == "common_water"):
            try: i = hy36decode(width=4, s=rg.resseq)
            except (RuntimeError, ValueError): pass
            else:
              if (result is None or result < i):
                result = i
            break
  return result

class xh_connectivity_table(object):
  # XXX need angle information as well
  def __init__(self, geometry, xray_structure):
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
        sites_cart=xray_structure.sites_cart())
    scatterers = xray_structure.scatterers()
    self.table = []
    for proxy in bond_proxies_simple:
      i_seq, j_seq = proxy.i_seqs
      i_x, i_h = None, None
      if(scatterers[i_seq].element_symbol() in ["H", "D"]):
        i_h = i_seq
        i_x = j_seq
        site_x = scatterers[i_x].site
        site_h = scatterers[i_h].site
        const_vect = flex.double(site_h)-flex.double(site_x)
        distance_model = xray_structure.unit_cell().distance(site_x, site_h)
        self.table.append([i_x, i_h, const_vect, proxy.distance_ideal,
                           distance_model])
      if(scatterers[j_seq].element_symbol() in ["H", "D"]):
        i_h = j_seq
        i_x = i_seq
        site_x = scatterers[i_x].site
        site_h = scatterers[i_h].site
        const_vect = flex.double(site_h)-flex.double(site_x)
        distance_model = xray_structure.unit_cell().distance(site_x, site_h)
        self.table.append([i_x, i_h, const_vect, proxy.distance_ideal,
                           distance_model])

class xh_connectivity_table2(object):
  def __init__(self, geometry, xray_structure):
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
        sites_cart=xray_structure.sites_cart())
    scatterers = xray_structure.scatterers()
    self.table = {}
    for proxy in bond_proxies_simple:
      i_seq, j_seq = proxy.i_seqs
      i_x, i_h = None, None
      if(scatterers[i_seq].element_symbol().upper() in ["H", "D"]):
        i_h = i_seq
        i_x = j_seq
      if(scatterers[j_seq].element_symbol().upper() in ["H", "D"]):
        i_h = j_seq
        i_x = i_seq
      if([i_x, i_h].count(None)==0):
        site_x = scatterers[i_x].site
        site_h = scatterers[i_h].site
        const_vect = flex.double(site_h)-flex.double(site_x)
        distance_model = xray_structure.unit_cell().distance(site_x, site_h)
        self.table.setdefault(i_h, []).append([i_x, i_h, const_vect,
          proxy.distance_ideal, distance_model])
    for p in geometry.geometry.angle_proxies:
      k,l,m = p.i_seqs
      els = [scatterers[k].element_symbol().upper(),
             scatterers[l].element_symbol().upper(),
             scatterers[m].element_symbol().upper()]
      o = flex.double()
      h = []
      ih=None
      if(els.count("H")<2 and els.count("D")<2):
        for i in p.i_seqs:
          s = scatterers[i]
          o.append(s.occupancy)
          sct = s.scattering_type.strip().upper()
          h.append(sct)
          if(sct in ["H","D"]): ih = i
        if("H" in h and not o.all_eq(o[0])):
          self.table.setdefault(ih).append(p.i_seqs)

class manager(object):
  """
  Wrapper class for storing and manipulating an iotbx.pdb.hierarchy object and
  a cctbx.xray.structure object, plus optional restraints-related objects.
  Being refactored to handle all model needs.
  Refactoring roadmap:
  1. Be able to create it in a new way (from model_input object) and keep
    everything else working. This constraints changing names of variables.
  2. Use it in some small places.
  3. Use it in phenix.refine and make sure it is really working
  4. Check all places where it is used and convert them to new way of creating it.
  5. Start refactoring with phenix.refine to remove unnecessary things.
  """

  def __init__(self,
      model_input,
      crystal_symmetry = None,
      restraint_objects = None, # ligand restraints in cif format
      monomer_parameters = None, # mmtbx.utils.cif_params scope # temporarily for phenix.refine
      pdb_interpretation_params = None,
      process_input = False, # obtain processed_pdb_file straight away
      build_grm = False,  # build GRM straight away, without waiting for get_restraints_manager() call
      stop_for_unknowns = True,
      log = None,
      expand_with_mtrix = True,
      # for GRM, selections etc. Do not use when creating an object.
                     processed_pdb_file = None, # Temporary, for refactoring phenix.refine
                     xray_structure = None, # remove later
                     pdb_hierarchy = None,  # remove later
                     restraints_manager = None, # remove later
                     tls_groups = None,         # remove later? ! used in def select()
                     ):

    self._xray_structure = xray_structure
    self._pdb_hierarchy = pdb_hierarchy
    self._model_input = model_input
    self._restraint_objects = restraint_objects
    self._monomer_parameters = monomer_parameters
    self._pdb_interpretation_params = None
    self.set_pdb_interpretation_params(pdb_interpretation_params)
    # Important! if shift_manager is not None, model_input - in original coords,
    # self.get_hierarchy(), self.get_xray_structure, self.get_sites_cart - in shifted coords.
    # self.crystal_symmetry() - shifted (boxed) one
    # If shift_manager is None - everything is consistent.
    self._shift_manager = None # mmtbx.utils.extract_box_around_model_and_map

    self._stop_for_unknowns = stop_for_unknowns
    self._processed_pdb_file = processed_pdb_file
    self._processed_pdb_files_srv = None
    self._crystal_symmetry = crystal_symmetry

    self.restraints_manager = restraints_manager
    self.refinement_flags = None
    # IAS related, need some cleaning!
    self.ias_manager = None
    self.use_ias = False    # remove later, use presence of ias_manager
                          # somewhat challenging because of ias.build_only parameter in phenix.refine
    self.tls_groups = tls_groups         # remove later? No action performed on them here
    self._ncs_groups = None
    self._anomalous_scatterer_groups = []
    self.log = log
    self.exchangable_hd_groups = []
    self.original_xh_lengths = None
    self.riding_h_manager = None
    self.header_tls_groups = None
    # for reprocessing. Probably select() will also use...
    self.scattering_dict_info = None

    self._ss_annotation = None
    self.link_records_in_pdb_format = None # Nigel's stuff up to refactoring
    self.all_chain_proxies = None # probably will keep it. Need to investigate
          # 'smart' selections provided by it. If they useless, remove.

    # These are new
    self.xray_scattering_dict = None
    self.neutron_scattering_dict = None
    self.has_hd = None
    self.model_statistics_info = None
    self._master_sel = flex.size_t([]) # selection of master part if NCS constraints are in use
    self._atom_selection_cache = None
    self._ncs_obj = None
    self._mon_lib_srv = None
    self._ener_lib = None
    self._rotamer_eval = None
    self._rotamer_id = None
    self._rama_eval = None
    self._original_model_format = None
    self._ss_manager = None
    self._site_symmetry_table = None

    # here we start to extract and fill appropriate field one by one
    # depending on what's available.
    if self._model_input is not None:
      s = str(type(model_input))
      if s.find("cif") > 0:
        self._original_model_format = "mmcif"
      elif s.find("pdb") > 0:
        self._original_model_format = "pdb"
      # input xray_structure most likely don't have proper crystal symmetry
      if self.crystal_symmetry() is None:
        self._crystal_symmetry = self._model_input.crystal_symmetry()
      if expand_with_mtrix:
        self.expand_with_MTRIX_records()

    if process_input or build_grm:
      assert self._processed_pdb_file is None
      assert self.all_chain_proxies is None
      self._process_input_model()

    # do pdb_hierarchy
    if self._pdb_hierarchy is None: # got nothing in parameters
      if self._processed_pdb_file is not None:
        self.all_chain_proxies = self._processed_pdb_file.all_chain_proxies
        self._pdb_hierarchy = self.all_chain_proxies.pdb_hierarchy
      elif self._model_input is not None:
        # self._pdb_hierarchy = deepcopy(self._model_input).construct_hierarchy()
        self._pdb_hierarchy = deepcopy(self._model_input).construct_hierarchy(
            self._pdb_interpretation_params.pdb_interpretation.sort_atoms)
    self._update_atom_selection_cache()
    self._update_pdb_atoms()

    if not self.all_chain_proxies and self._processed_pdb_file:
      self.all_chain_proxies = self._processed_pdb_file.all_chain_proxies

    # do GRM
    if self.restraints_manager is None and build_grm:
      self.setup_restraints_manager()

    # do xray_structure

    # if self._xray_structure is None:
    #   self._create_xray_structure()
    if self._xray_structure is not None:
      if self.crystal_symmetry() is None:
        self._crystal_symmetry = self._xray_structure.crystal_symmetry()

    if self._xray_structure is not None:
      assert self._xray_structure.scatterers().size() == self._pdb_hierarchy.atoms_size()

    # if self.crystal_symmetry() is None:
    #   self._crystal_symmetry = self._xray_structure.crystal_symmetry()
    # print "self._crystal_symmetry3", self._crystal_symmetry

    if self.has_hd is None:
      self._update_has_hd()

    if self.has_hd:
      self.exchangable_hd_groups = utils.combine_hd_exchangable(
        hierarchy = self._pdb_hierarchy)

  @staticmethod
  def get_default_pdb_interpretation_scope():
    """
    Get parsed parameters (libtbx.phil.scope). Use this function to
    avoid importing pdb_interpretation phil strings and think about how to
    parse it. Does not need the instance of class (staticmethod).
    Then modify what needed to be modified and init this class normally.
    """
    return iotbx.phil.parse(
          input_string=grand_master_phil_str+reference_model_str,
          process_includes=True)

  @staticmethod
  def get_default_pdb_interpretation_params():
    """
    Get the default extract object (libtbx.phil.scope_extract)
    """
    return manager.get_default_pdb_interpretation_scope().extract()

  def set_log(self, log):
    self.log = log

  def set_shift_manager(self, shift_manager):
    self._shift_manager = shift_manager
    if shift_manager is not None:
      self.set_xray_structure(self._shift_manager.xray_structure_box)
      self.set_crystal_symmetry_if_undefined(self._shift_manager.get_shifted_cs())
      self.unset_restraints_manager()

  def get_shift_manager(self):
    return self._shift_manager

  def set_pdb_interpretation_params(self, params):
    #
    # Consider invalidating self.restraints_manager here, because it could be already
    # constructed with different params. Done.
    #
    # check if we got only inside of pdb_interpretation scope.
    # For mmtbx.command_line.load_model_and_data
    if params is None:
      self._pdb_interpretation_params = manager.get_default_pdb_interpretation_params()
    else:
      # if getattr(params, "sort_atoms", None) is not None:
      full_params = manager.get_default_pdb_interpretation_params()
      if getattr(params, "sort_atoms", None) is not None:
        full_params.pdb_interpretation = params
      if hasattr(params, "pdb_interpretation"):
        full_params.pdb_interpretation = params.pdb_interpretation
      if hasattr(params, "geometry_restraints"):
        full_params.geometry_restraints = params.geometry_restraints
      if hasattr(params, "reference_model"):
        full_params.reference_model = params.reference_model
      self._pdb_interpretation_params = full_params
    self.unset_restraints_manager()

  def check_consistency(self):
    """
    Primarilly for debugging
    """
    # print "NUMBER OF ATOMS", self.get_number_of_atoms()
    assert self.pdb_atoms.size() == self.get_xray_structure().scatterers().size()
    assert self.pdb_atoms.size() == self.get_hierarchy().atoms().size()
    # finally, check all_chain_proxies
    if self.all_chain_proxies is not None:
      self.all_chain_proxies.extract_xray_structure().scatterers().size() == self.pdb_atoms.size()
    # checking label consistency
    # for ha, sc in zip(self.get_hierarchy().atoms_with_labels(),
    #     self.get_xray_structure().scatterers()):
    #   if ha.id_str() != sc.label and sc.scattering_type[:2] != "IS":
    #     print "XXXXXXX: '%s' != '%s'" % (ha.id_str(), sc.label)
    #     assert 0


  def get_model_input(self):
    return self._model_input

  def crystal_symmetry(self):
    return self._crystal_symmetry

  def get_restraint_objects(self):
    return self._restraint_objects

  def setup_ss_annotation(self, log=null_out()):
    self._ss_annotation = self._model_input.extract_secondary_structure()

  def get_ss_annotation(self, log=null_out()):
    if self._ss_annotation is None:
      self.setup_ss_annotation(log=log)
    return self._ss_annotation

  def set_ss_annotation(self, ann):
    self._ss_annotation = ann

  def set_crystal_symmetry_if_undefined(self, cs):
    """
    Function to set crystal symmetry if it is not defined yet.
    Special case when incoming cs is the same to self._crystal_symmetry,
    then just do nothing for compatibility with existing code.
    It would be better to remove this function at all.
    """
    is_empty_cs = self._crystal_symmetry is None or (
        self._crystal_symmetry.unit_cell() is None and
        self._crystal_symmetry.space_group() is None)
    if not is_empty_cs:
      if self._crystal_symmetry.is_similar_symmetry(cs):
        return
    assert is_empty_cs, "%s,%s" % (
          self._crystal_symmetry.show_summary(), cs.show_summary())
    self._crystal_symmetry = cs

  def set_refinement_flags(self, flags):
    self.refinement_flags = flags

  def get_number_of_atoms(self):
    return self.pdb_atoms.size()

  def get_atoms(self):
    return self.pdb_atoms

  def get_sites_cart(self):
    if self.get_xray_structure() is None:
      self._create_xray_structure()
    return self._xray_structure.sites_cart()

  def get_atom_selection_cache(self):
    if self._atom_selection_cache is not None:
      return self._atom_selection_cache
    return None

  def get_number_of_models(self):
    if self._model_input is not None:
      return self._model_input.model_ids().size()
    return 0

  def get_site_symmetry_table(self):
    if self._site_symmetry_table is not None:
      return self._site_symmetry_table
    elif self.all_chain_proxies is not None:
      self._site_symmetry_table = self.all_chain_proxies.site_symmetry_table()
    return self._site_symmetry_table

  def initialize_anomalous_scatterer_groups(
      self,
      find_automatically=True,
      groups_from_params=None):
    self.n_anomalous_total = 0
    result = []
    if find_automatically:
      result = anomalous_scatterer_groups.find_anomalous_scatterer_groups(
        pdb_atoms=self.get_atoms(),
        xray_structure=self.get_xray_structure(),
        group_same_element=False,
        out=self.log)
      for group in result:
        self.n_anomalous_total += group.iselection.size()
    else:
      if len(groups_from_params) != 0:
        chain_proxies = self.all_chain_proxies
        sel_cache = self.get_atom_selection_cache()
        for group in groups_from_params:
          if (group.f_prime is None): group.f_prime = 0
          if (group.f_double_prime is None): group.f_double_prime = 0
          aag = xray.anomalous_scatterer_group(
            iselection=chain_proxies.phil_atom_selection(
              cache=sel_cache,
              scope_extract=group,
              attr="selection",
              raise_if_empty_selection=True).iselection(),
            f_prime=group.f_prime,
            f_double_prime=group.f_double_prime,
            refine=group.refine,
            selection_string=group.selection)
          # aag.show_summary(out=log)
          result.append(aag)
          self.n_anomalous_total += aag.iselection.size()
    self._anomalous_scatterer_groups = result
    for group in self._anomalous_scatterer_groups:
      group.copy_to_scatterers_in_place(scatterers=self._xray_structure.scatterers())
    return self._anomalous_scatterer_groups, self.n_anomalous_total

  def set_anomalous_scatterer_groups(self, groups):
    self._anomalous_scatterer_groups = groups
    new_n_anom_total = 0
    for g in self._anomalous_scatterer_groups:
      new_n_anom_total += g.iselection.size()
    self.n_anomalous_total = new_n_anom_total

  def get_anomalous_scatterer_groups(self):
    return self._anomalous_scatterer_groups

  def have_anomalous_scatterer_groups(self):
    return (self._anomalous_scatterer_groups is not None and
        len(self._anomalous_scatterer_groups) > 0)

  def update_anomalous_groups (self, out=sys.stdout) :
    if self.have_anomalous_scatterer_groups():
      modified = False
      sel_cache = self.get_atom_selection_cache()
      for i_group, group in enumerate(self._anomalous_scatterer_groups) :
        if (group.update_from_selection) :
          isel = sel_cache.selection(group.selection_string).iselection()
          assert (len(isel) == len(group.iselection))
          if (not isel.all_eq(group.iselection)) :
            print >> out, "Updating %d atom(s) in anomalous group %d" % \
              (len(isel), i_group+1)
            print >> out, "  selection string: %s" % group.selection_string
            group.iselection = isel
            modified = True
      return modified

  def anomalous_scatterer_groups_as_pdb(self):
    out = StringIO()
    if self.have_anomalous_scatterer_groups():
      pr = "REMARK   3  "
      print >> out, pr+"ANOMALOUS SCATTERER GROUPS DETAILS."
      print >> out, pr+" NUMBER OF ANOMALOUS SCATTERER GROUPS : %-6d"%\
        len(self._anomalous_scatterer_groups)
      counter = 0
      for group in self._anomalous_scatterer_groups:
        counter += 1
        print >>out,pr+" ANOMALOUS SCATTERER GROUP : %-6d"%counter
        lines = str_utils.line_breaker(group.selection_string, width=45)
        for i_line, line in enumerate(lines):
          if(i_line == 0):
            print >> out, pr+"  SELECTION: %s"%line
          else:
            print >> out, pr+"           : %s"%line
        print >>out,pr+"  fp  : %-15.4f"%group.f_prime
        print >>out,pr+"  fdp : %-15.4f"%group.f_double_prime
    return out.getvalue()

  def get_restraints_manager(self):
    if self.restraints_manager is not None:
      return self.restraints_manager
    else:
      self.setup_restraints_manager()
      return self.restraints_manager

  def set_non_unit_occupancy_implies_min_distance_sym_equiv_zero(self,value):
    if self._xray_structure is not None:
      self._xray_structure.set_non_unit_occupancy_implies_min_distance_sym_equiv_zero(value)
      self.set_xray_structure(self._xray_structure.customized_copy(
          non_unit_occupancy_implies_min_distance_sym_equiv_zero=value))

  def get_hd_selection(self):
    if self._xray_structure is not None:
      return self._xray_structure.hd_selection()
    return None

  def get_ias_selection(self):
    if self.ias_manager is None:
      return None
    else:
      return self.ias_manager.get_ias_selection()

  def selection(self, selstr, optional=True):
    if self.all_chain_proxies is None:
      return self.get_atom_selection_cache().selection(selstr, optional=optional)
    else:
      return self.all_chain_proxies.selection(
          selstr,
          cache=self._atom_selection_cache,
          optional=optional)

  def iselection(self, selstr):
    result = self.selection(selstr)
    if result is None:
      return None
    return result.iselection()

  def sel_backbone(self):
    assert self.all_chain_proxies is not None
    return self.all_chain_proxies.sel_backbone_or_sidechain(True, False)

  def sel_sidechain(self):
    assert self.all_chain_proxies is not None
    return self.all_chain_proxies.sel_backbone_or_sidechain(False, True)

  def get_xray_structure(self):
    if self._xray_structure is None:
      self._create_xray_structure()
    return self._xray_structure

  def set_xray_structure(self, xray_structure, update_hierarchy=True):
    same_symmetry = True
    if self._xray_structure is not None:
      same_symmetry = self._xray_structure.crystal_symmetry().is_similar_symmetry(
          xray_structure.crystal_symmetry())
    # This is happening e.g. in iotbx/map_and_model.py where xrs and map are
    # being boxed and model is updated with this function.
    if not same_symmetry:
      self.unset_restraints_manager()
    self._xray_structure = xray_structure
    if update_hierarchy:
      self.set_sites_cart_from_xrs()
    self._update_has_hd()
    self._update_pdb_atoms()
    # Not using self.set_crystal_symmetry() because it is different.
    self._crystal_symmetry = xray_structure.crystal_symmetry()
    if(not self._xray_structure.crystal_symmetry().is_similar_symmetry(
       xray_structure.crystal_symmetry())):
      self.unset_restraints_manager()
    # XXX what else needs to be done here?

  def _create_xray_structure(self):
    if self._xray_structure is not None:
      return
    if self.all_chain_proxies is not None:
      self._xray_structure = self.all_chain_proxies.extract_xray_structure()
    else:
      self._xray_structure = self._pdb_hierarchy.extract_xray_structure(
          crystal_symmetry=self.crystal_symmetry())

  def get_mon_lib_srv(self):
    if self._mon_lib_srv is None:
      self._mon_lib_srv = mmtbx.monomer_library.server.server()
    return self._mon_lib_srv

  def get_ener_lib(self):
    if self._ener_lib is None:
      self._ener_lib = mmtbx.monomer_library.server.ener_lib()
    return self._ener_lib

  def get_rotamer_manager(self):
    if self._rotamer_eval is None:
      self._rotamer_eval = RotamerEval(mon_lib_srv=self.get_mon_lib_srv())
    return self._rotamer_eval

  def get_rotamer_id(self):
    if self._rotamer_id is None:
      self._rotamer_id = RotamerID()
    return self._rotamer_id

  def get_ramachandran_manager(self):
    if self._rama_eval is None:
      self._rama_eval = rama_eval()
    return self._rama_eval

  def get_apply_cif_links(self):
    if self.all_chain_proxies is not None:
      return self.all_chain_proxies.apply_cif_links
    return []

  def update_xrs(self, hierarchy=None):
    """
    Updates xray structure using self._pdb_hierarchy for cases when it
    was modified outside. E.g. refinement, minimization, etc.
    """
    if hierarchy is not None:
      # !!! This could fail even for very similar hierarchies.
      # see example in
      # cctbx_project/mmtbx/secondary_structure/build/tst_1.py
      assert hierarchy.is_similar_hierarchy(other=self._pdb_hierarchy)
      self._pdb_hierarchy = hierarchy
    self._xray_structure = self._pdb_hierarchy.extract_xray_structure(
        crystal_symmetry=self.crystal_symmetry())

  def model_as_pdb(self,
      output_cs = True,
      atoms_reset_serial_first_value=None,
      do_not_shift_back = False):
    """
    move all the writing here later.
    """
    if do_not_shift_back and self._shift_manager is None:
      do_not_shift_back = False
    cs_to_output = None
    if output_cs:
      if do_not_shift_back:
        cs_to_output = self._shift_manager.get_shifted_cs()
      else:
        cs_to_output = self.crystal_symmetry()

    result = StringIO()
    # outputting HELIX/SHEET records
    ss_records = ""
    ss_ann = None
    if self._ss_manager is not None:
      ss_ann = self._ss_manager.actual_sec_str
    elif self.get_ss_annotation() is not None:
      ss_ann = self.get_ss_annotation()
    if ss_ann is not None:
      ss_records = ss_ann.as_pdb_str()
    if ss_records != "":
      if ss_records[-1] != "\n":
        ss_records += "\n"
      result.write(ss_records)

    #
    # Here should be NCS output somehow
    #
    if (self.link_records_in_pdb_format is not None
        and len(self.link_records_in_pdb_format)>0):
      result.write("%s\n" % self.link_records_in_pdb_format)

    hierarchy_to_output = self.get_hierarchy()
    if hierarchy_to_output is not None:
      hierarchy_to_output = hierarchy_to_output.deep_copy()
    if self._shift_manager is not None and not do_not_shift_back:
      self._shift_manager.shift_back(hierarchy_to_output)

    if hierarchy_to_output is not None:
      result.write(hierarchy_to_output.as_pdb_string(
          crystal_symmetry=cs_to_output,
          atoms_reset_serial_first_value=atoms_reset_serial_first_value,
          append_end=True))
    return result.getvalue()

  def model_as_mmcif(self,
      cif_block_name = "default",
      additional_blocks = None,
      align_columns = False):
    out = StringIO()
    cif = iotbx.cif.model.cif()
    cif_block = None
    if self.crystal_symmetry() is not None:
      cif_block = self.crystal_symmetry().as_cif_block()
    if self._pdb_hierarchy is not None:
      if self.crystal_symmetry() is not None:
        cif_block.update(self._pdb_hierarchy.as_cif_block())
      else:
        cif_block = self._pdb_hierarchy.as_cif_block()
    # outputting HELIX/SHEET records
    ss_cif_loops = []
    ss_ann = None
    if self._ss_manager is not None:
      ss_ann = self._ss_manager.actual_sec_str
    elif self.get_ss_annotation() is not None:
      ss_ann = self.get_ss_annotation()
    if ss_ann is not None:
      ss_cif_loops = ss_ann.as_cif_loops()
    for loop in ss_cif_loops:
      cif_block.add_loop(loop)

    # self.get_model_statistics_info()
    if self.model_statistics_info is not None:
      cif_block.update(self.model_statistics_info.as_cif_block())
    if additional_blocks is not None:
      for ab in additional_blocks:
        cif_block.update(ab)
    cif_block.sort(key=category_sort_function)
    cif[cif_block_name] = cif_block
    # here we are outputting a huge number of restraints
    #
    if (self.all_chain_proxies is not None
        and self.all_chain_proxies.cif is not None):
      # should writing ALL restraints be optional?
      skip_residues = one_letter_given_three_letter.keys() + ['HOH']
      restraints = self.all_chain_proxies.cif
      keys = restraints.keys()
      for key in keys:
        # need more control
        assert key.find('UNK')==-1
        if key.replace('comp_', '') in skip_residues:
          del restraints[key]
      cif.update(self.all_chain_proxies.cif)
    cif.show(out=out, align_columns=align_columns)
    return out.getvalue()

  def restraints_as_geo(self,
      header = "# Geometry restraints\n",
      excessive_distance_limit = 1.5,
      force=False):
    """
    get geo file as string.
    force = True actually will try to build GRM. Not advised, because if
    it is not already build, most likely a program never needed it, so
    no reason to output.
    """
    result = StringIO()
    if force:
      self.get_restraints_manager()
    self.restraints_manager.write_geo_file(
        hierarchy = self.get_hierarchy(),
        sites_cart=self.get_sites_cart(),
        site_labels=self.get_xray_structure().scatterers().extract_labels(),
        header=header,
        # Stuff for outputting ncs_groups
        excessive_distance_limit = excessive_distance_limit,
        xray_structure=self.get_xray_structure(),
        file_descriptor=result)
    return result.getvalue()

  def input_format_was_cif(self):
    return self._original_model_format == "mmcif"

  def _process_input_model(self):
    # Not clear if we can handle this correctly for self._xray_structure
    # assert self.get_number_of_models() < 2
    # assert 0
    if self._processed_pdb_files_srv is None:
      self._processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
          crystal_symmetry          = self.crystal_symmetry(),
          pdb_interpretation_params = self._pdb_interpretation_params.pdb_interpretation,
          stop_for_unknowns         = self._stop_for_unknowns,
          log                       = self.log,
          cif_objects               = self._restraint_objects,
          cif_parameters            = self._monomer_parameters, # mmtbx.utils.cif_params scope - should be refactored to remove
          mon_lib_srv               = None,
          ener_lib                  = None,
          use_neutron_distances     = self._pdb_interpretation_params.pdb_interpretation.use_neutron_distances)
    if self._processed_pdb_file is None:
      self._processed_pdb_file, junk = self._processed_pdb_files_srv.process_pdb_files(
          pdb_inp = self._model_input,
          hierarchy = self._pdb_hierarchy,
          # because hierarchy already extracted
          # raw_records = flex.split_lines(self._pdb_hierarchy.as_pdb_string()),
          # stop_if_duplicate_labels = True,
          allow_missing_symmetry=True)
    if self.all_chain_proxies is None:
      self.all_chain_proxies = self._processed_pdb_file.all_chain_proxies
    self._atom_selection_cache = self._processed_pdb_file.all_chain_proxies.pdb_hierarchy.atom_selection_cache()
    if self._pdb_hierarchy is None:
      self._pdb_hierarchy = self._processed_pdb_file.all_chain_proxies.pdb_hierarchy
    if self._xray_structure is None:
      xray_structure_all = \
          self._processed_pdb_file.xray_structure(show_summary = False)
      # XXX ad hoc manipulation
      for sc in xray_structure_all.scatterers():
        lbl=sc.label.split()
        if("IAS" in lbl and sc.scattering_type=="?" and lbl[1].startswith("IS")):
          sc.scattering_type = lbl[1]
      #
      if(xray_structure_all is None):
        raise Sorry("Cannot extract xray_structure.")
      if(xray_structure_all.scatterers().size()==0):
        raise Sorry("Empty xray_structure.")
      if self.all_chain_proxies is not None:
        self.all_chain_proxies = self._processed_pdb_file.all_chain_proxies
      self._xray_structure = xray_structure_all

    self._mon_lib_srv = self._processed_pdb_files_srv.mon_lib_srv
    self._ener_lib = self._processed_pdb_files_srv.ener_lib
    self._ncs_obj = self._processed_pdb_file.ncs_obj

    # copy-paste from command_line/geometry_minimization.py: get_geometry_restraints_manager
    # should live only here and be removed there.
    self._update_has_hd()

  def _update_has_hd(self):
    if(self._xray_structure is not None):
      sctr_keys = self._xray_structure.scattering_type_registry().type_count_dict().keys()
      self.has_hd = "H" in sctr_keys or "D" in sctr_keys
    if not self.has_hd:
      self.unset_riding_h_manager()

  def _update_atom_selection_cache(self):
    if self.all_chain_proxies is not None:
      self._atom_selection_cache = self.all_chain_proxies.pdb_hierarchy.atom_selection_cache()
    elif self.crystal_symmetry() is not None:
      self._atom_selection_cache = self.get_hierarchy().atom_selection_cache(
          special_position_settings=crystal.special_position_settings(
              crystal_symmetry = self.crystal_symmetry() ))
    else:
      self._atom_selection_cache = self.get_hierarchy().atom_selection_cache()

  def unset_restraints_manager(self):
    self.restraints_manager = None

  def setup_restraints_manager(
      self,
      grm_normalization = True,
      external_energy_function = None,
      plain_pairs_radius=5.0,
      custom_nb_excl=None,
      file_descriptor_for_geo_in_case_of_failure=None,
      ):
    if(self.restraints_manager is not None): return
    if self._processed_pdb_file is None:
      self._process_input_model()

    # assert self.get_number_of_models() < 2 # one molprobity test triggered this.
    # not clear how they work with multi-model stuff...

    # if self._pdb_interpretation_params is None:
    #   self._pdb_interpretation_params = master_params().fetch().extract()
    # disabled temporarily due to architecture changes
    geometry = self._processed_pdb_file.geometry_restraints_manager(
      show_energies      = False,
      plain_pairs_radius = plain_pairs_radius,
      params_edits       = self._pdb_interpretation_params.geometry_restraints.edits,
      params_remove      = self._pdb_interpretation_params.geometry_restraints.remove,
      custom_nonbonded_exclusions  = custom_nb_excl,
      external_energy_function=external_energy_function,
      assume_hydrogens_all_missing = not self.has_hd,
      file_descriptor_for_geo_in_case_of_failure=file_descriptor_for_geo_in_case_of_failure)

    self._ss_manager = self._processed_pdb_file.ss_manager

    # For test GRM pickling
    # from cctbx.regression.tst_grm_pickling import make_geo_pickle_unpickle
    # geometry = make_geo_pickle_unpickle(
    #     geometry=geometry,
    #     xrs=xray_structure,
    #     prefix=None)

    # Link treating should be rewritten. They should not be saved in
    # all_chain_proxies and they should support mmcif.
    self.link_records_in_pdb_format = link_record_output(self._processed_pdb_file.all_chain_proxies)
    # For test GRM pickling
    # from cctbx.regression.tst_grm_pickling import make_geo_pickle_unpickle
    # geometry = make_geo_pickle_unpickle(
    #     geometry=geometry,
    #     xrs=xray_structure,
    #     prefix=None)

    restraints_manager = mmtbx.restraints.manager(
      geometry      = geometry,
      normalization = grm_normalization)
    # Torsion restraints from reference model
    if hasattr(self._pdb_interpretation_params, "reference_model") and restraints_manager is not None:
      add_reference_dihedral_restraints_if_requested(
          geometry=restraints_manager.geometry,
          processed_pdb_file=self._processed_pdb_file,
          mon_lib_srv=self.get_mon_lib_srv(),
          ener_lib=self.get_ener_lib(),
          has_hd=self.has_hd,
          params=self._pdb_interpretation_params.reference_model,
          selection=None,
          log=self.log)
    if(self._xray_structure is not None):
      restraints_manager.crystal_symmetry = self._xray_structure.crystal_symmetry()
    self.restraints_manager = restraints_manager

    # Here we do all what is necessary when GRM and all related become available
    #
    self.extract_tls_selections_from_input()

  #
  # =======================================================================
  # NCS-related features.
  # =======================================================================
  #

  def setup_torsion_ncs_restraints(self,
      fmodel,
      ncs_torsion_params,
      sites_individual,
      log):
    torsion_utils.check_for_internal_chain_ter_records(
        pdb_hierarchy = self.get_hierarchy(),
        ter_indices   = self._model_input.ter_indices())
    ncs_obj = self.get_ncs_obj()
    if ncs_obj is None: return
    geometry = self.get_restraints_manager().geometry
    if ncs_obj.number_of_ncs_groups > 0:
      print >> log, "\n"
      geometry.ncs_dihedral_manager = torsion_ncs(
          model              = self,
          fmodel             = fmodel,
          params             = ncs_torsion_params,
          selection          = sites_individual,
          log                = log)
      if geometry.ncs_dihedral_manager.get_n_proxies() == 0:
        geometry.ncs_dihedral_manager = None
    geometry.sync_reference_dihedral_with_ncs(log=log)

  def torsion_NCS_present(self):
    rm = self.get_restraints_manager()
    if rm is None:
      return False
    if rm.geometry.ncs_dihedral_manager is None:
      return False
    return True

  def torsion_NCS_as_pdb(self):
    result = StringIO()
    if self.torsion_NCS_present():
      self.get_restraints_manager().geometry.ncs_dihedral_manager.as_pdb(
          out=result)
    return result.getvalue()

  def torsion_NCS_as_cif_block(self, cif_block):
    if not self.torsion_NCS_present():
      return cif_block
    loops = manager._get_NCS_cif_loops()
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    cif_block = self.get_restraints_manager().geometry.ncs_dihedral_manager.as_cif_block(
        loops=loops, cif_block=cif_block)
    return cif_block

  def torsion_ncs_restraints_update(self, log=None):
    if self.get_restraints_manager() is not None:
      self.get_restraints_manager().geometry.update_dihedral_ncs_restraints(
          model=self,
          log=log)

  def ncs_constraints_present(self):
    g = self.get_ncs_groups()
    return g is not None and len(g)>0

  def setup_ncs_constraints_groups(self, filter_groups=False):
    """
    This will be used directly (via get_ncs_groups) in
    mmtbx/refinement/minimization.py, mmtbx/refinement/adp_refinement.py

    supposedly NCS constraints
    """
    if self.get_ncs_obj() is not None:
      self._ncs_groups = self.get_ncs_obj().get_ncs_restraints_group_list()
    if filter_groups:
      ncs_group_list = ncs_obj.get_ncs_restraints_group_list()
      # set up new groups for refinements
      self._ncs_groups = self.ncs_restr_group_list.filter_ncs_restraints_group_list(
          self.get_hierarchy())
    if len(self._ncs_groups) > 0:
      # determine master selections
      self._master_sel = flex.bool(self.get_number_of_atoms(), True)
      for ncs_gr in self._ncs_groups:
        for copy in ncs_gr.copies:
          self._master_sel.set_selected(copy.iselection, False)

  def get_master_hierarchy(self):
    if self._master_sel.size() > 0:
      return self.get_hierarchy().select(self._master_sel)
    else:
      return self.get_hierarchy()

  def get_master_selection(self):
    return self._master_sel

  def get_ncs_obj(self):
    return self._ncs_obj

  def get_ncs_groups(self):
    """
    This returns object mmtbx.ncs.ncs.ncs_group, used primarily in
    Tom's tools for historic reasons
    """
    return self._ncs_groups

  def setup_cartesian_ncs_groups(self, ncs_params=None, log=null_out()):
    import mmtbx.ncs.cartesian_restraints
    cartesian_ncs = mmtbx.ncs.cartesian_restraints.cartesian_ncs_manager(
        model=self,
        ncs_params=ncs_params)
    rm = self.get_restraints_manager()
    if cartesian_ncs.get_n_groups() > 0:
      assert rm is not None
      rm.cartesian_ncs_manager = cartesian_ncs
    else:
      print >> self.log, "No NCS restraint groups specified."
      print >> self.log

  def cartesian_NCS_as_pdb(self):
    result = StringIO()
    if (self.restraints_manager is not None and
        self.restraints_manager.cartesian_ncs_manager is not None):
      self.restraints_manager.cartesian_ncs_manager.as_pdb(
          sites_cart=self.get_sites_cart(),
          out=result)
    return result.getvalue()

  @staticmethod
  def _get_NCS_cif_loops():
    ncs_ens_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_ens.id",
      "_struct_ncs_ens.details"))
    ncs_dom_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_dom.id",
      "_struct_ncs_dom.pdbx_ens_id",
      "_struct_ncs_dom.details"))
    ncs_dom_lim_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_dom_lim.pdbx_ens_id",
      "_struct_ncs_dom_lim.dom_id",
      #"_struct_ncs_dom_lim.pdbx_component_id",
      #"_struct_ncs_dom_lim.pdbx_refine_code",
      "_struct_ncs_dom_lim.beg_auth_asym_id",
      "_struct_ncs_dom_lim.beg_auth_seq_id",
      "_struct_ncs_dom_lim.end_auth_asym_id",
      "_struct_ncs_dom_lim.end_auth_seq_id",
      "_struct_ncs_dom_lim.selection_details"))

    ncs_oper_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_oper.id",
      "_struct_ncs_oper.code",
      "_struct_ncs_oper.matrix[1][1]",
      "_struct_ncs_oper.matrix[1][2]",
      "_struct_ncs_oper.matrix[1][3]",
      "_struct_ncs_oper.matrix[2][1]",
      "_struct_ncs_oper.matrix[2][2]",
      "_struct_ncs_oper.matrix[2][3]",
      "_struct_ncs_oper.matrix[3][1]",
      "_struct_ncs_oper.matrix[3][2]",
      "_struct_ncs_oper.matrix[3][3]",
      "_struct_ncs_oper.vector[1]",
      "_struct_ncs_oper.vector[2]",
      "_struct_ncs_oper.vector[3]",
      "_struct_ncs_oper.details"))

    ncs_ens_gen_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_ens_gen.dom_id_1",
      "_struct_ncs_ens_gen.dom_id_2",
      "_struct_ncs_ens_gen.ens_id",
      "_struct_ncs_ens_gen.oper_id"))

    return ncs_ens_loop, ncs_dom_loop, ncs_dom_lim_loop, ncs_oper_loop, ncs_ens_gen_loop

  def cartesian_NCS_present(self):
    c_ncs = self.get_cartesian_NCS_manager()
    if c_ncs is not None:
      return c_ncs.get_n_groups() > 0
    return False

  def get_cartesian_NCS_manager(self):
    if self.get_restraints_manager() is not None:
      return self.get_restraints_manager().cartesian_ncs_manager
    return None

  def cartesian_NCS_as_cif_block(self, cif_block=None):
    if not self.cartesian_NCS_present():
      return cif_block
    loops = manager._get_NCS_cif_loops()
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    if (self.get_restraints_manager() is not None and
        self.get_restraints_manager().cartesian_ncs_manager is not None):
      cif_block = self.get_restraints_manager().cartesian_ncs_manager.as_cif_block(
          loops=loops, cif_block=cif_block, sites_cart=self.get_sites_cart())
    return cif_block


  def setup_scattering_dictionaries(self,
      scattering_table,
      d_min=None,
      log = None,
      set_inelastic_form_factors=None,
      iff_wavelength=None):
    self._create_xray_structure()
    self.scattering_dict_info = group_args(
        scattering_table=scattering_table,
        d_min = d_min,
        set_inelastic_form_factors=set_inelastic_form_factors,
        iff_wavelength=iff_wavelength)
    if(log is not None):
      str_utils.make_header("Scattering factors", out = log)
    known_scattering_tables = [
      "n_gaussian", "wk1995", "it1992", "electron", "neutron"]
    if(not (scattering_table in known_scattering_tables)):
      raise Sorry("Unknown scattering_table: %s\n%s"%
        (show_string(scattering_table),
        "Possible choices are: %s"%" ".join(known_scattering_tables)))
    if(scattering_table in ["n_gaussian", "wk1995", "it1992", "electron"]):
      self._xray_structure.scattering_type_registry(
        table = scattering_table,
        d_min = d_min,
        types_without_a_scattering_contribution=["?"])
      self._xray_structure.scattering_type_registry(
        custom_dict = ias.ias_scattering_dict)
      self.xray_scattering_dict = \
        self._xray_structure.scattering_type_registry().as_type_gaussian_dict()
      if(log is not None):
        print_statistics.make_sub_header("X-ray scattering dictionary",out=log)
        self._xray_structure.scattering_type_registry().show(out = log)
    if(scattering_table == "neutron"):
      try :
        self.neutron_scattering_dict = \
          self._xray_structure.switch_to_neutron_scattering_dictionary()
      except ValueError, e :
        raise Sorry("Error setting up neutron scattering dictionary: %s"%str(e))
      if(log is not None):
        print_statistics.make_sub_header(
          "Neutron scattering dictionary", out = log)
        self._xray_structure.scattering_type_registry().show(out = log)
      self._xray_structure.scattering_type_registry_params.table = "neutron"
    if self.all_chain_proxies is not None:
      scattering_type_registry = self.all_chain_proxies.scattering_type_registry
      if(scattering_type_registry.n_unknown_type_symbols() > 0):
        scattering_type_registry.report(
          pdb_atoms = self.pdb_atoms,
          log = log,
          prefix = "",
          max_lines = None)
        raise Sorry("Unknown scattering type symbols.\n"
          "  Possible ways of resolving this error:\n"
          "    - Edit columns 77-78 in the PDB file to define"
            " the scattering type.\n"
          "    - Provide custom monomer definitions for the affected residues.")
      if(log is not None):
        print >> log
    if set_inelastic_form_factors is not None and iff_wavelength is not None:
      self._xray_structure.set_inelastic_form_factors(
          photon=iff_wavelength,
          table=set_inelastic_form_factors)
    return self.xray_scattering_dict, self.neutron_scattering_dict

  def get_header_tls_selections(self):
    if "input_tls_selections" not in self.__dict__.keys():
      self.extract_tls_selections_from_input()
    return self.input_tls_selections

  def get_searched_tls_selections(self, nproc):
    if "searched_tls_selections" not in self.__dict__.keys():
      tls_params = find_tls_groups.master_phil.fetch().extract()
      tls_params.nproc = nproc
      self.searched_tls_selections = find_tls_groups.find_tls(
        params=tls_params,
        pdb_inp=self._model_input,
        pdb_hierarchy=self._pdb_hierarchy,
        xray_structure=deepcopy(self._xray_structure),
        return_as_list=True,
        ignore_pdb_header_groups=True,
        out=StringIO())
    return self.searched_tls_selections

  def determine_tls_groups(self, selection_strings, generate_tlsos):
    self.tls_groups = tls_tools.tls_groups(selection_strings = selection_strings)
    if generate_tlsos is not None:

      tlsos = tls_tools.generate_tlsos(
        selections     = generate_tlsos,
        xray_structure = self._xray_structure,
        value          = 0.0)
      self.tls_groups.tlsos = tlsos

  def tls_groups_as_pdb(self):
    out = StringIO()
    if self.tls_groups is not None:
      tls_tools.remark_3_tls(
          tlsos             = self.tls_groups.tlsos,
          selection_strings = self.tls_groups.selection_strings,
          out               = out)
    return out.getvalue()

  def tls_groups_as_cif_block(self, cif_block=None):
    if self.tls_groups is not None:
      cif_block = tls_as_cif_block(
        tlsos=self.tls_groups.tlsos,
        selection_strings=self.tls_groups.selection_strings,
        cif_block=cif_block)
    return cif_block

  def extract_tls_selections_from_input(self):
    self.input_tls_selections = []
    acp = self.all_chain_proxies
    pdb_inp_tls = None

    if self._model_input is not None:
      pdb_inp_tls = self._model_input.extract_tls_params(self._pdb_hierarchy)
    elif acp is not None and acp.pdb_inp is not None:
      pdb_inp_tls = acp.pdb_inp.extract_tls_params(acp.pdb_hierarchy)

    if(pdb_inp_tls and pdb_inp_tls.tls_present):
      print_statistics.make_header(
        "TLS group selections from PDB file header", out=self.log)
      print >> self.log, "TLS group selections:"
      atom_counts = []
      for t in pdb_inp_tls.tls_params:
        try :
          n_atoms = self.iselection(t.selection_string).size()
        except AtomSelectionError, e :
          print >> self.log, "AtomSelectionError:"
          print >> self.log, str(e)
          print >> self.log, "Ignoring PDB header TLS groups"
          self.input_tls_selections = []
          return
        print >> self.log, "  selection string:"
        print >> self.log, "    %s"%t.selection_string
        print >> self.log, "    selects %d atoms"%n_atoms
        self.input_tls_selections.append(t.selection_string)
        atom_counts.append(n_atoms)
      if(pdb_inp_tls.tls_present):
        if(pdb_inp_tls.error_string is not None):
          print >> self.log, "  %s"%pdb_inp_tls.error_string
          self.input_tls_selections = []
      if(0 in atom_counts):
        msg="""
  One of TLS selections is an empty selection: skipping TLS infromation found in
  PDB file header.
"""
        print >> self.log, msg

  def get_riding_h_manager(self, idealize=True, force=False):
    """
    Force=True: Force creating manager if hydrogens are available
    """
    if self.riding_h_manager is None and force:
      self.setup_riding_h_manager(idealize=True)
    return self.riding_h_manager

  def unset_riding_h_manager(self):
    self.riding_h_manager = None

  def setup_riding_h_manager(self, idealize=True):
    assert self.riding_h_manager is None
    if not self.has_hd: return
    if(self.restraints_manager is None): return
    sites_cart = self._xray_structure.sites_cart()
    self.riding_h_manager = riding.manager(
      pdb_hierarchy       = self.get_hierarchy(),
      geometry_restraints = self.get_restraints_manager().geometry)
    if(idealize):
      self.riding_h_manager.idealize(sites_cart = sites_cart)
      self.set_sites_cart(sites_cart)

  def get_hierarchy(self, sync_with_xray_structure=False):
    """
    Accessor for the underlying PDB hierarchy, incorporating optional update of
    properties from the xray_structure attribute.

    :param sync_with_xray_structure: update all properties of the atoms in the
      hierarchy to agree with the xray_structure scatterers.  Will fail if the
      scatterer labels do not match the atom labels.
    """
    if(sync_with_xray_structure and self._xray_structure is not None):
      self._pdb_hierarchy.adopt_xray_structure(
        xray_structure = self._xray_structure)
    return self._pdb_hierarchy

  def set_sites_cart_from_hierarchy(self, multiply_ncs=False):
    if (multiply_ncs and self.ncs_constraints_present()):
      ncs_groups = self.get_ncs_groups()
      for ncs_gr in ncs_groups:
        h = self.get_hierarchy()
        new_sites = h.select(ncs_gr.master_iselection).atoms().extract_xyz()
        for c in ncs_gr.copies:
          new_c_sites = c.r.elems * new_sites + c.t
          h.select(c.iselection).atoms().set_xyz(new_c_sites)
    self._xray_structure.set_sites_cart(self._pdb_hierarchy.atoms().extract_xyz())
    self._update_pdb_atoms()
    self.model_statistics_info = None

  def _update_pdb_atoms(self):
    self.pdb_atoms = self._pdb_hierarchy.atoms()
    self.pdb_atoms.reset_i_seq()

  def set_sites_cart_from_xrs(self):
    self._pdb_hierarchy.adopt_xray_structure(self._xray_structure)
    self._update_pdb_atoms()
    self.model_statistics_info = None

  def set_sites_cart(self, sites_cart, update_grm=False):
    assert sites_cart.size() == self._pdb_hierarchy.atoms_size() == \
        self._xray_structure.scatterers().size()
    self._xray_structure.set_sites_cart(sites_cart)
    self._pdb_hierarchy.adopt_xray_structure(self._xray_structure)
    if update_grm:
      self.restraints_manager.geometry.pair_proxies(sites_cart)
    self.model_statistics_info = None
    self._update_pdb_atoms()

  def sync_pdb_hierarchy_with_xray_structure(self):
    # to be deleted.
    self._pdb_hierarchy.adopt_xray_structure(xray_structure=self._xray_structure)
    self._update_pdb_atoms()
    self.model_statistics_info = None

  def normalize_adjacent_adp(self):
    bond_proxies_simple, asu = \
      self.restraints_manager.geometry.get_all_bond_proxies(
        sites_cart = self._xray_structure.sites_cart())
    scatterers = self._xray_structure.scatterers()
    for proxy in bond_proxies_simple:
      i_seq, j_seq = proxy.i_seqs
      bi = adptbx.u_as_b(scatterers[i_seq].u_iso)
      bj = adptbx.u_as_b(scatterers[j_seq].u_iso)
      if(abs(bi-bj)>20):
        if(bi>bj): bi=bj+20
        else:      bj=bi+20
        scatterers[i_seq].u_iso = adptbx.b_as_u(bi)
        scatterers[j_seq].u_iso = adptbx.b_as_u(bj)
    self.set_sites_cart_from_xrs()

  def xh_connectivity_table(self):
    result = None
    if(self.restraints_manager is not None):
      if self.has_hd:
        xray_structure = self._xray_structure
        ias_selection = self.get_ias_selection()
        if ias_selection and ias_selection.count(True) > 0:
          xray_structure = self._xray_structure.select(~ias_selection)
        result = xh_connectivity_table(
          geometry       = self.restraints_manager,
          xray_structure = xray_structure).table
    return result

  def xh_connectivity_table2(self):
    result = None
    if(self.restraints_manager is not None):
      if self.has_hd:
        xray_structure = self._xray_structure
        ias_selection = self.get_ias_selection()
        if ias_selection and ias_selection.count(True) > 0:
          xray_structure = self._xray_structure.select(~ias_selection)
        result = xh_connectivity_table2(
          geometry       = self.restraints_manager,
          xray_structure = xray_structure).table
    return result

  def extend_xh_bonds(self, value=1.1):
    if(self.restraints_manager is None): return
    if not self.has_hd: return
    assert self.original_xh_lengths is None
    h_i_seqs = []
    xhct = self.xh_connectivity_table()
    if(xhct is None): return
    self.original_xh_lengths = flex.double()
    for xhcti in xhct:
      h_i_seqs.append(xhcti[1])
    for bp in self.restraints_manager.geometry.bond_params_table:
      for i, k in enumerate(bp.keys()):
        if(k in h_i_seqs):
          self.original_xh_lengths.append(bp.values()[i].distance_ideal)
          bp.values()[i].distance_ideal = value

  def restore_xh_bonds(self):
    if(self.restraints_manager is None): return
    if not self.has_hd: return
    assert self.original_xh_lengths is not None
    xhct = self.xh_connectivity_table()
    if(xhct is None): return
    h_i_seqs = []
    for xhcti in xhct:
      h_i_seqs.append(xhcti[1])
    counter = 0
    for bp in self.restraints_manager.geometry.bond_params_table:
      for i, k in enumerate(bp.keys()):
        if(k in h_i_seqs):
          bp.values()[i].distance_ideal = self.original_xh_lengths[counter]
          counter += 1
    self.original_xh_lengths = None
    self.idealize_h(show=False)

  def isolated_atoms_selection(self):
    if(self.restraints_manager is None):
      raise Sorry("Geometry restraints manager must be defined.")
    selection = flex.bool(self._xray_structure.scatterers().size(), True)
    bond_proxies_simple, asu = self.restraints_manager.geometry.\
        get_all_bond_proxies(sites_cart=self._xray_structure.sites_cart())
    for proxy in bond_proxies_simple:
      i_seq, j_seq = proxy.i_seqs
      selection[i_seq] = False
      selection[j_seq] = False
    return selection

  def reset_adp_for_hydrogens(self):
    """
    Set the isotropic B-factor for all hydrogens to those of the associated
    heavy atoms (using the total isotropic equivalent) times a scale factor of
    1.2.
    """
    if(self.restraints_manager is None): return
    hd_sel = self.get_hd_selection()
    if(hd_sel.count(True) > 0):
      xh_conn_table = self.xh_connectivity_table()
      bfi = self._xray_structure.extract_u_iso_or_u_equiv()
      for t in self.xh_connectivity_table():
        i_x, i_h = t[0], t[1]
        bfi[i_h] = adptbx.u_as_b(bfi[i_x])*1.2
      self._xray_structure.set_b_iso(values = bfi, selection = hd_sel)
    self.set_sites_cart_from_xrs()

  def reset_occupancies_for_hydrogens(self):
    """
    Set hydrogen occupancies to those of the associated heavy atoms.
    """
    occupancy_refinement_selections_1d = flex.size_t()
    if(self.refinement_flags.s_occupancies is not None):
      for occsel in self.refinement_flags.s_occupancies:
        for occsel_ in occsel:
          occupancy_refinement_selections_1d.extend(occsel_)
    if(self.restraints_manager is None): return
    hd_sel = self.get_hd_selection()
    scatterers = self._xray_structure.scatterers()
    if(hd_sel.count(True) > 0):
      xh_conn_table = self.xh_connectivity_table()
      qi = self._xray_structure.scatterers().extract_occupancies()
      ct = self.xh_connectivity_table2()
      for t_ in ct.values():
        i_x, i_h = t_[0][0],t_[0][1]
        assert scatterers[i_h].element_symbol() in ["H", "D"]
        if(scatterers[i_x].element_symbol() == "N" and
           i_h in occupancy_refinement_selections_1d):
          occ = flex.double()
          for t in t_:
            if(len(t) != 5):
              for i in t:
                if(i != i_h):
                  occ.append(qi[i])
              qi[i_h] = flex.min(occ)
            else:
              qi[i_h] = qi[i_x]
        else:
          qi[i_h] = qi[i_x]
      if(self.refinement_flags.s_occupancies is not None):
        for rf1 in self.refinement_flags.s_occupancies:
          o=None
          for rf2 in rf1:
            for rf_ in rf2:
              if(not hd_sel[rf_]):
                o = qi[rf_]
            for rf_ in rf2:
              if(o is not None):
                qi[rf_] = o
      self._xray_structure.scatterers().set_occupancies(qi, hd_sel)
      self.set_sites_cart_from_xrs()

  def reset_coordinates_for_exchangable_hd(self):
    if(len(self.exchangable_hd_groups) > 0):
      scatterers =  self._xray_structure.scatterers()
      occ = scatterers.extract_occupancies()
      for g in self.exchangable_hd_groups:
        i, j = g[0][0], g[1][0]
        if(occ[i] > occ[j]):
          scatterers[j].site = scatterers[i].site
        else:
          scatterers[i].site = scatterers[j].site
    self.set_sites_cart_from_xrs()

  def rotatable_hd_selection(self):
    rmh_sel = mmtbx.hydrogens.rotatable(
      pdb_hierarchy      = self.get_hierarchy(),
      mon_lib_srv        = self.get_mon_lib_srv(),
      restraints_manager = self.get_restraints_manager(),
      log                = self.log)
    sel_i = []
    for s in rmh_sel: sel_i.extend(s[1])
    return flex.size_t(sel_i)

  def h_counts(self):
    occupancies = self._xray_structure.scatterers().extract_occupancies()
    occ_sum = flex.sum(occupancies)
    hd_selection = self.get_hd_selection()
    h_occ_sum = flex.sum(occupancies.select(hd_selection))
    sel_rot = self.rotatable_hd_selection()
    hrot_occ_sum = flex.sum(occupancies.select(sel_rot))
    return group_args(
      h_count             = hd_selection.count(True),
      h_occ_sum           = h_occ_sum,
      h_fraction_of_total = h_occ_sum/occ_sum*100.,
      hrot_count             = sel_rot.size(),
      hrot_occ_sum           = hrot_occ_sum,
      hrot_fraction_of_total = hrot_occ_sum/occ_sum*100.)

  def show_h_counts(self, prefix=""):
    hc = self.h_counts()
    print >> self.log, "%sTotal:"%prefix
    print >> self.log, "%s  count: %d"%(prefix, hc.h_count)
    print >> self.log, "%s  occupancy sum: %6.2f (%s of total atoms %6.2f)"%(
      prefix, hc.h_occ_sum, "%", hc.h_fraction_of_total)
    print >> self.log, "%sRotatable:"%prefix
    print >> self.log, "%s  count: %d"%(prefix, hc.hrot_count)
    print >> self.log, "%s  occupancy sum: %6.2f (%s of total atoms %6.2f)"%(
      prefix, hc.hrot_occ_sum, "%", hc.hrot_fraction_of_total)

  def scattering_types_counts_and_occupancy_sums(self, prefix=""):
    out = StringIO()
    st_counts_and_occupancy_sums = \
        self.get_xray_structure().scattering_types_counts_and_occupancy_sums()
    atoms_occupancy_sum = \
        flex.sum(self.get_xray_structure().scatterers().extract_occupancies())
    fmt = "   %5s               %10d        %8.2f"
    print >> out, prefix+"MODEL CONTENT."
    print >> out, prefix+" ELEMENT        ATOM RECORD COUNT   OCCUPANCY SUM"
    for item in st_counts_and_occupancy_sums:
      print >> out, prefix+fmt % (item.scattering_type, item.count,
        item.occupancy_sum)
    print >> out,prefix+fmt%("TOTAL",self.get_number_of_atoms(),atoms_occupancy_sum)
    return out.getvalue()

  def idealize_h(self, correct_special_position_tolerance=1.0,
                   selection=None, show=True, nuclear=False):
    """
    Perform geometry regularization on hydrogen atoms only.
    """
    if(self.restraints_manager is None): return
    if self.has_hd:
      hd_selection = self.get_hd_selection()
      if(selection is not None): hd_selection = selection
      if(hd_selection.count(True)==0): return
      not_hd_selection = ~hd_selection
      sol_hd = self.solvent_selection().set_selected(~hd_selection, False)
      mac_hd = hd_selection.deep_copy().set_selected(self.solvent_selection(), False)
      ias_selection = self.get_ias_selection()
      if (ias_selection is not None):
        not_hd_selection.set_selected(ias_selection, False)
      sites_cart_mac_before = \
        self._xray_structure.sites_cart().select(not_hd_selection)
      xhd = flex.double()
      if(hd_selection.count(True)==0): return
      for t in self.xh_connectivity_table():
        if(hd_selection[t[1]]):
          xhd.append(abs(t[-1]-t[-2]))
      if(show):
        print >> self.log, \
        "X-H deviation from ideal before regularization (bond): mean=%6.3f max=%6.3f"%\
        (flex.mean(xhd), flex.max(xhd))
      for sel_pair in [(mac_hd, False), (sol_hd, True)]*2:
        if(sel_pair[0].count(True) > 0):
          sel = sel_pair[0]
          if(ias_selection is not None and ias_selection.count(True) > 0):
            sel = sel.select(~ias_selection)
          minimized = geometry_minimization.run2(
              restraints_manager = self.get_restraints_manager(),
              pdb_hierarchy = self.get_hierarchy(),
              correct_special_position_tolerance = correct_special_position_tolerance,
              riding_h_manager               = None, # didn't go in original implementation
              ncs_restraints_group_list      = [], # didn't go in original implementation
              max_number_of_iterations       = 500,
              number_of_macro_cycles         = 5,
              selection                      = sel,
              bond                           = True,
              nonbonded                      = sel_pair[1],
              angle                          = True,
              dihedral                       = True,
              chirality                      = True,
              planarity                      = True,
              parallelity                    = True,
              # rmsd_bonds_termination_cutoff  = rmsd_bonds_termination_cutoff,
              # rmsd_angles_termination_cutoff = rmsd_angles_termination_cutoff,
              # alternate_nonbonded_off_on     = False, # default
              # cdl                            = False,
              # rdl                            = False,
              # correct_hydrogens              = False,
              # fix_rotamer_outliers           = True,
              # allow_allowed_rotamers         = True,
              # states_collector               = None,
              log                            = StringIO(),
              mon_lib_srv                    = self.get_mon_lib_srv())
          self.set_sites_cart_from_hierarchy()
      sites_cart_mac_after = \
        self._xray_structure.sites_cart().select(not_hd_selection)
      assert approx_equal(flex.max(sites_cart_mac_before.as_double() -
        sites_cart_mac_after.as_double()), 0)
      xhd = flex.double()
      for t in self.xh_connectivity_table():
        if(hd_selection[t[1]]):
          xhd.append(abs(t[-1]-t[-2]))
      if(show):
        print >> self.log,\
        "X-H deviation from ideal after  regularization (bond): mean=%6.3f max=%6.3f"%\
        (flex.mean(xhd), flex.max(xhd))

  #def idealize_h(self):
  #   if(self.riding_h_manager is None):
  #     self.setup_riding_h_manager()
  #   self.riding_h_manager.idealize(
  #     sites_cart = self._xray_structure.sites_cart())

  def extract_water_residue_groups(self):
    result = []
    solvent_sel = self.solvent_selection()
    get_class = iotbx.pdb.common_residue_names_get_class
    for m in self._pdb_hierarchy.models():
      for chain in m.chains():
        for rg in chain.residue_groups():
          first_water = None
          first_other = None
          for ag in rg.atom_groups():
            residue_id = "%3s%4s%1s" % (ag.resname, rg.resseq, rg.icode)
            if (get_class(name=ag.resname) == "common_water"):
              for atom in ag.atoms():
                i_seq = atom.i_seq
                assert solvent_sel[i_seq]
                if (first_water is None):
                  first_water = atom
            else:
              for atom in ag.atoms():
                assert not solvent_sel[atom.i_seq]
                if (first_other is None):
                  first_other = atom
          if (first_water is not None):
            if (first_other is not None):
              raise RuntimeError(
                "residue_group with mix of water and non-water:\n"
                + "  %s\n" % first_water.quote()
                + "  %s" % first_other.quote())
            result.append(rg)
            for r in result:
              elements = r.atoms().extract_element()
              o_found = 0
              for e in elements:
                if(e.strip().upper() == 'O'): o_found += 1
              if(o_found == 0):
                print >> self.log
                for a in r.atoms():
                  print >> self.log, a.format_atom_record()
                raise Sorry(
                  "The above waters in input PDB file do not have O atom.")
    return result

  def renumber_water(self):
    for i,rg in enumerate(self.extract_water_residue_groups()):
      rg.resseq = iotbx.pdb.resseq_encode(value=i+1)
      rg.icode = " "
    self._update_pdb_atoms()
    self._sync_xrs_labels()


  def add_hydrogens(self, correct_special_position_tolerance,
        element = "H", neutron = False, occupancy=0.01):
    result = []
    xs = self._xray_structure
    if(neutron): element = "D"
    frac = xs.unit_cell().fractionalize
    sites_cart = xs.sites_cart()
    u_isos = xs.extract_u_iso_or_u_equiv()
    next_to_i_seqs = []
    last_insert_i_seq = [sites_cart.size()]
    def insert_atoms(atom, atom_names, element):
      i_seq = atom.i_seq
      assert i_seq < last_insert_i_seq[0]
      last_insert_i_seq[0] = i_seq
      xyz = sites_cart[i_seq]
      sign = True
      for i,atom_name in enumerate(atom_names):
        h = atom.detached_copy()
        h.name = atom_name
        if(sign):
          h.xyz = [a+b for a,b in zip(xyz, (1,0.001,0))]
          sign = False
        else:
          h.xyz = [a+b for a,b in zip(xyz, (-1,0,0))]
          sign = True
        h.sigxyz = (0,0,0)
        h.occ = occupancy
        h.sigocc = 0
        h.b = adptbx.u_as_b(u_isos[i_seq])
        h.sigb = 0
        h.uij = (-1,-1,-1,-1,-1,-1)
        if (iotbx.pdb.hierarchy.atom.has_siguij()):
          h.siguij = (-1,-1,-1,-1,-1,-1)
        h.element = "%2s" % element.strip()
        ag.append_atom(atom=h)
        scatterer = xray.scatterer(
          label           = h.id_str(),
          scattering_type = h.element.strip(),
          site            = frac(h.xyz),
          u               = adptbx.b_as_u(h.b),
          occupancy       = h.occ)
        xs.add_scatterer(
          scatterer = scatterer,
          insert_at_index = i_seq+i+1)
        next_to_i_seqs.append(i_seq) # not i_seq+i because refinement_flags.add
                                     # sorts next_to_i_seqs internally :-(
    water_rgs = self.extract_water_residue_groups()
    water_rgs.reverse()
    for rg in water_rgs:
      if (rg.atom_groups_size() != 1):
        raise Sorry("Not implemented: cannot add hydrogens to water "+
                    "molecules with alternate conformations")
      ag = rg.only_atom_group()
      atoms = ag.atoms()
      # do not add H or D to O at or close to special position
      skip = False
      sps = self._xray_structure.special_position_settings(
        min_distance_sym_equiv=3.0)
      for atom in atoms:
        if (atom.element.strip() == "O"):
          sps_r = sps.site_symmetry(site_cart=atom.xyz).is_point_group_1()
          if(not sps_r):
            skip = True
            break
      #
      if(not skip):
        if (atoms.size() == 2):
          o_atom = None
          h_atom = None
          for atom in atoms:
            if (atom.element.strip() == "O"): o_atom = atom
            else:                             h_atom = atom
          assert [o_atom, h_atom].count(None) == 0
          h_name = h_atom.name.strip()
          if(len(h_name) == 1):
            atom_name = " %s1 " % h_name
          elif(len(h_name) == 2):
            if(h_name[0].isdigit()):
              if(int(h_name[0]) == 1): atom_name = " %s2 " % h_name[1]
              elif(int(h_name[0]) == 2): atom_name = " %s1 " % h_name[1]
              else: raise RuntimeError
            elif(h_name[1].isdigit()):
              if(int(h_name[1]) == 1): atom_name = " %s2 " % h_name[0]
              elif(int(h_name[1]) == 2): atom_name = " %s1 " % h_name[0]
              else: raise RuntimeError
            else: raise RuntimeError
          else: raise RuntimeError
          insert_atoms(
            atom=o_atom,
            atom_names=[atom_name],
            element=h_atom.element)
        elif (atoms.size() == 1):
          atom = atoms[0]
          assert atom.element.strip() == "O"
          insert_atoms(
            atom=atom,
            atom_names=[" "+element+n+" " for n in ["1","2"]],
            element=element)
    if(neutron):
      xs.switch_to_neutron_scattering_dictionary()
    print >> self.log, "Number of H added:", len(next_to_i_seqs)
    if (len(next_to_i_seqs) == 0): return
    if (self.refinement_flags is not None):
      self.refinement_flags.add(
        next_to_i_seqs=next_to_i_seqs,
        sites_individual = True,
        s_occupancies    = neutron)
    self.reprocess_pdb_hierarchy_inefficient()
    self.idealize_h()

  def reprocess_pdb_hierarchy_inefficient(self):
    # XXX very inefficient
    """
    Re-process PDB from scratch and create restraints.  Not recommended for
    general use.
    """
    raw_records = self.model_as_pdb()
    pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records)
    pip = self._pdb_interpretation_params
    pip.pdb_interpretation.clash_guard.nonbonded_distance_threshold = -1.0
    pip.pdb_interpretation.clash_guard.max_number_of_distances_below_threshold = 100000000
    pip.pdb_interpretation.clash_guard.max_fraction_of_distances_below_threshold = 1.0
    pip.pdb_interpretation.proceed_with_excessive_length_bonds=True
    pip.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
    flags = self.refinement_flags
    ppr = self.restraints_manager.geometry.plain_pairs_radius
    cs = self.crystal_symmetry()
    norm = self.restraints_manager.normalization
    log = self.log
    scattering_dict_info = self.scattering_dict_info
    self.all_chain_proxies = None
    self.__init__(
        model_input = pdb_inp,
        crystal_symmetry = cs,
        restraint_objects = self._restraint_objects,
        pdb_interpretation_params = pip,
        process_input = True,
        build_grm = False,
        log = StringIO()
        )
    self.setup_restraints_manager(
        plain_pairs_radius = ppr,
        grm_normalization = norm)
    self.set_refinement_flags(flags)
    if scattering_dict_info is not None:
      self.setup_scattering_dictionaries(
          scattering_table=scattering_dict_info.scattering_table,
          d_min = scattering_dict_info.d_min,
          set_inelastic_form_factors = scattering_dict_info.set_inelastic_form_factors,
          iff_wavelength = scattering_dict_info.iff_wavelength)

  def hd_group_selections(self):
    return utils.combine_hd_exchangable(hierarchy = self._pdb_hierarchy)

  def reset_adp_of_hd_sites_to_be_equal(self):
    scatterers = self._xray_structure.scatterers()
    adp_fl = self.refinement_flags.adp_individual_iso
    adp_fl_a = self.refinement_flags.adp_individual_aniso
    if(adp_fl is not None):
      for gsel in self.hd_group_selections():
        i,j = gsel[0][0], gsel[1][0]
        element_symbols = \
          [scatterers[i].element_symbol(), scatterers[j].element_symbol()]
        assert element_symbols.count('H') > 0 and element_symbols.count('D')>0
        i_seq_max_q = None
        i_seq_min_q = None
        if(scatterers[i].occupancy < scatterers[j].occupancy):
          i_seq_max_q = j
          i_seq_min_q = i
        else:
          i_seq_max_q = i
          i_seq_min_q = j
        if([adp_fl[i_seq_max_q], adp_fl[i_seq_min_q]].count(True) > 0):
          if([adp_fl_a[i_seq_max_q], adp_fl_a[i_seq_min_q]].count(True) > 0):
            continue
          adp_fl[i_seq_max_q] = True
          adp_fl[i_seq_min_q] = False
          assert [adp_fl[i_seq_max_q], adp_fl[i_seq_min_q]].count(True) > 0
          scatterers[i_seq_min_q].u_iso = scatterers[i_seq_max_q].u_iso
    self.set_sites_cart_from_xrs()

  def rms_b_iso_or_b_equiv_bonded(self):
    return utils.rms_b_iso_or_b_equiv_bonded(
      restraints_manager = self.restraints_manager,
      xray_structure     = self._xray_structure,
      ias_manager        = self.ias_manager)

  def deep_copy(self):
    return self.select(selection = flex.bool(
      self._xray_structure.scatterers().size(), True))

  def add_ias(self, fmodel=None, ias_params=None, file_name=None,
                                                             build_only=False):
    """
    Generate interatomic scatterer pseudo-atoms.
    """
    if(self.ias_manager is not None):
       self.remove_ias()
       fmodel.update_xray_structure(xray_structure = self._xray_structure,
                                    update_f_calc = True)
    print >> self.log, ">>> Adding IAS.........."
    self.old_refinement_flags = None
    if not build_only: self.use_ias = True
    self.ias_manager = ias.manager(
                    geometry             = self.restraints_manager.geometry,
                    pdb_atoms            = self.pdb_atoms,
                    xray_structure       = self._xray_structure,
                    fmodel               = fmodel,
                    params               = ias_params,
                    file_name            = file_name,
                    log                  = self.log)
    size_all = self._xray_structure.scatterers().size()
    if(not build_only):
      ias_xray_structure = self.ias_manager.ias_xray_structure
      ias_selection = self.get_ias_selection()
      self._xray_structure.concatenate_inplace(other = ias_xray_structure)
      print >> self.log, "Scattering dictionary for combined xray_structure:"
      self._xray_structure.scattering_type_registry().show(out=self.log)
      if(self.refinement_flags is not None):
         self.old_refinement_flags = self.refinement_flags.deep_copy()
         # define flags
         ssites = flex.bool(ias_xray_structure.scatterers().size(), False)
         sadp = flex.bool(ias_xray_structure.scatterers().size(), False)
         # XXX set occ refinement ONLY for involved atoms
         # XXX now it refines only occupancies of IAS !!!
         occupancy_flags = []
         ms = ias_selection.count(False)
         for i in range(1, ias_selection.count(True)+1):
           occupancy_flags.append([flex.size_t([ms+i-1])])
         # set flags
         self.refinement_flags.inflate(
           sites_individual     = ssites,
           s_occupancies        = occupancy_flags,
           adp_individual_iso   = sadp,
           adp_individual_aniso = sadp,
           size_all             = size_all)
         # adjust flags
         if(self.refinement_flags.sites_individual is not None):
           self.refinement_flags.sites_individual.set_selected(ias_selection, False)
           self.refinement_flags.sites_individual.set_selected(~ias_selection, True)
         if(self.refinement_flags.adp_individual_aniso is not None):
           self.refinement_flags.adp_individual_aniso.set_selected(ias_selection, False)
         if(self.refinement_flags.adp_individual_iso is not None):
           self.refinement_flags.adp_individual_iso.set_selected(ias_selection, True)
         #occs = flex.double(self._xray_structure.scatterers().size(), 0.9)
         #self._xray_structure.scatterers().set_occupancies(occs, ~self.ias_selection)
         # D9
         sel = self._xray_structure.scatterers().extract_scattering_types() == "IS9"
         self._xray_structure.convert_to_anisotropic(selection = sel)
         if(self.refinement_flags.adp_individual_aniso is not None):
           self.refinement_flags.adp_individual_aniso.set_selected(sel, True)
         if(self.refinement_flags.adp_individual_iso is not None):
           self.refinement_flags.adp_individual_iso.set_selected(sel, False)
    n_sites = ias_xray_structure.scatterers().size()
    atom_names = []
    for sc in ias_xray_structure.scatterers():
      new_atom_name = sc.label.strip()
      if(len(new_atom_name) < 4): new_atom_name = " " + new_atom_name
      while(len(new_atom_name) < 4): new_atom_name = new_atom_name+" "
      atom_names.append(new_atom_name)
    residue_names = [ "IAS" ] * n_sites
    self._append_pdb_atoms(
      new_xray_structure=ias_xray_structure,
      atom_names=atom_names,
      residue_names=residue_names,
      chain_id=" ",
      i_seq_start=0,
      reset_labels=True, # maybe this will allow to drop IS handling in
      # adopt_xray_structure()
      )

  def remove_ias(self):
    print >> self.log, ">>> Removing IAS..............."
    self.use_ias = False
    ias_selection = self.get_ias_selection()
    if(self.ias_manager is not None):
      self.ias_manager = None
    if(self.old_refinement_flags is not None):
      self.refinement_flags = self.old_refinement_flags.deep_copy()
      self.old_refinement_flags = None
    if(ias_selection is not None):
      self._xray_structure.select_inplace(
        selection = ~ias_selection)
      self._xray_structure.scattering_type_registry().show(out=self.log)
      self._pdb_hierarchy = self._pdb_hierarchy.select(
        atom_selection = ~ias_selection)
      self._update_pdb_atoms()

  def show_rigid_bond_test(self, out=None, use_id_str=False, prefix=""):
    if (out is None): out = sys.stdout
    scatterers = self._xray_structure.scatterers()
    unit_cell = self._xray_structure.unit_cell()
    rbt_array = flex.double()
    sites_cart = self._xray_structure.sites_cart()
    ias_selection = self.get_ias_selection()
    if ias_selection is not None:#
      sites_cart = sites_cart.select(~ias_selection)
    bond_proxies_simple, asu = self.restraints_manager.geometry.\
        get_all_bond_proxies(sites_cart=sites_cart)
    for proxy in bond_proxies_simple:
      i_seqs = proxy.i_seqs
      i,j = proxy.i_seqs
      atom_i = self.pdb_atoms[i]
      atom_j = self.pdb_atoms[j]
      if (    atom_i.element.strip() not in ["H","D"]
          and atom_j.element.strip() not in ["H","D"]):
        sc_i = scatterers[i]
        sc_j = scatterers[j]
        if (sc_i.flags.use_u_aniso() and sc_j.flags.use_u_aniso()):
          p = adp_restraints.rigid_bond_pair(
            sc_i.site, sc_j.site, sc_i.u_star, sc_j.u_star, unit_cell)
          rbt_value = p.delta_z()*10000.
          rbt_array.append(rbt_value)
          name_i, name_j = atom_i.name, atom_j.name
          if (use_id_str) :
            name_i = atom_i.id_str()
            name_j = atom_j.id_str()
          print >> out, "%s%s %s %10.3f"%(prefix, name_i, name_j, rbt_value)
    if (rbt_array.size() != 0):
      print >> out, "%sRBT values (*10000):" % prefix
      print >> out, "%s  mean = %.3f" % (prefix, flex.mean(rbt_array))
      print >> out, "%s  max  = %.3f" % (prefix, flex.max(rbt_array))
      print >> out, "%s  min  = %.3f" % (prefix, flex.min(rbt_array))

  def restraints_manager_energies_sites(self,
        geometry_flags=None,
        custom_nonbonded_function=None,
        compute_gradients=False,
        gradients=None,
        disable_asu_cache=False):
    if(self.restraints_manager is None): return None
    sites_cart = self._xray_structure.sites_cart()
    ias_selection = self.get_ias_selection()
    if(self.use_ias and ias_selection is not None and
       ias_selection.count(True) > 0):
      sites_cart = sites_cart.select(~ias_selection)
    result = self.restraints_manager.energies_sites(
      sites_cart=sites_cart,
      geometry_flags=geometry_flags,
      external_energy_function=None,
      custom_nonbonded_function=custom_nonbonded_function,
      compute_gradients=compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      hd_selection=self.get_hd_selection(),
    )
    return result

  def solvent_selection(self, include_ions=False):
    result = flex.bool()
    get_class = iotbx.pdb.common_residue_names_get_class
    for a in self.get_hierarchy().atoms():
      resname = (a.parent().resname).strip()
      if(get_class(name = resname) == "common_water"):
        result.append(True)
      elif (a.segid.strip() == "ION") and (include_ions) :
        result.append(True)
      else: result.append(False)
    return result

  def xray_structure_macromolecule(self):
    sel = self.solvent_selection(include_ions=True)
    if(self.use_ias): sel = sel | self.get_ias_selection()
    result = self._xray_structure.select(~sel)
    return result

  def select(self, selection):
    # what about 3 types of NCS and self._master_sel?
    # XXX ignores IAS
    if isinstance(selection, flex.size_t):
      selection = flex.bool(self.get_number_of_atoms(), selection)
    new_pdb_hierarchy = self._pdb_hierarchy.select(selection, copy_atoms=True)
    new_refinement_flags = None
    sdi = self.scattering_dict_info
    if(self.refinement_flags is not None):
      # XXX Tom
      try:
        new_refinement_flags = self.refinement_flags.select(selection)
      except Exception: # This is potential bug. What kind of exceptions are possible here?!
        new_refinement_flags = self.refinement_flags
#      new_refinement_flags = self.refinement_flags.select(selection)
    new_restraints_manager = None
    if(self.restraints_manager is not None):
      new_restraints_manager = self.restraints_manager.select(selection = selection)
      # XXX is it necessary ?
      # YYY yes, this keeps pair_proxies initialized and available, e.g. for
      # extracting info used in .geo files.
      new_restraints_manager.geometry.pair_proxies(
          sites_cart = self._xray_structure.sites_cart().select(selection))
    new_shift_manager = self._shift_manager
    new_riding_h_manager = None
    if self.riding_h_manager is not None:
      new_riding_h_manager = self.riding_h_manager.select(selection)
    new = manager(
      model_input                = self._model_input, # any selection here?
      crystal_symmetry           = self._crystal_symmetry,
      processed_pdb_file         = self._processed_pdb_file,
      restraints_manager         = new_restraints_manager,
      expand_with_mtrix          = False,
      xray_structure             = self.get_xray_structure().select(selection),
      pdb_hierarchy              = new_pdb_hierarchy,
      pdb_interpretation_params  = self._pdb_interpretation_params,
      tls_groups                 = self.tls_groups, # XXX not selected, potential bug
      log                        = self.log)
    if new_riding_h_manager is not None:
      new.riding_h_manager = new_riding_h_manager
    new.get_xray_structure().scattering_type_registry()
    new.set_refinement_flags(new_refinement_flags)
    new.scattering_dict_info = sdi
    new._update_has_hd()
    # selecting anomalous_scatterer_groups one by one because they are simple list!
    new_anom_groups = []
    for an_gr in self._anomalous_scatterer_groups:
      new_an_gr = an_gr.select(selection)
      if new_an_gr.iselection.size() > 0:
        new_anom_groups.append(new_an_gr)
    new.set_anomalous_scatterer_groups(new_anom_groups)
    # do not use set_shift_manager for this because shifted xray_str and hierarchy
    # are used here and being cutted. We need only shift_back and get_..._cs
    # functionality form shift_manager at the moment
    new._shift_manager = new_shift_manager
    new._processed_pdb_files_srv = self._processed_pdb_files_srv
    if self.ncs_constraints_present():
      new_ncs_groups = self._ncs_groups.select(selection)
      new._ncs_groups = new_ncs_groups
    if self._master_sel is not None:
      if self._master_sel.size() == selection.size():
        new._master_sel = self._master_sel.select(selection)
      else:
        if isinstance(self._master_sel, flex.bool):
          x = flex.bool(selection.size(), self._master_sel.iselection())
          new._master_sel = x.select(selection)
        else:
          x = flex.bool(selection.size(), self._master_sel)
          new._master_sel = x.select(selection)
    return new

  def number_of_ordered_solvent_molecules(self):
    return self.solvent_selection().count(True)

  def show_groups(self, rigid_body = None, tls = None,
                        out = None, text="Information about rigid groups"):
    global time_model_show
    timer = user_plus_sys_time()
    selections = None
    if(rigid_body is not None):
       selections = self.refinement_flags.sites_rigid_body
    if(tls is not None): selections = self.refinement_flags.adp_tls
    if(self.refinement_flags.sites_rigid_body is None and
                                 self.refinement_flags.adp_tls is None): return
    assert selections is not None
    if (out is None): out = sys.stdout
    print >> out
    line_len = len("| "+text+"|")
    fill_len = 80 - line_len-1
    upper_line = "|-"+text+"-"*(fill_len)+"|"
    print >> out, upper_line
    next = "| Total number of atoms = %-6d  Number of rigid groups = %-3d                |"
    natoms_total = self._xray_structure.scatterers().size()
    print >> out, next % (natoms_total, len(selections))
    print >> out, "| group: start point:                        end point:                       |"
    print >> out, "|               x      B  atom  residue  <>        x      B  atom  residue    |"
    next = "| %5d: %8.3f %6.2f %5s %3s %5s <> %8.3f %6.2f %5s %3s %5s   |"
    sites = self._xray_structure.sites_cart()
    b_isos = self._xray_structure.extract_u_iso_or_u_equiv() * math.pi**2*8
    for i_seq, selection in enumerate(selections):
      if (isinstance(selection, flex.bool)):
        i_selection = selection.iselection()
      else:
        i_selection = selection
      start = i_selection[0]
      final = i_selection[i_selection.size()-1]
      first = self.pdb_atoms[start]
      last  = self.pdb_atoms[final]
      first_ag = first.parent()
      first_rg = first_ag.parent()
      last_ag = last.parent()
      last_rg = last_ag.parent()
      print >> out, next % (i_seq+1,
        sites[start][0], b_isos[start],
          first.name, first_ag.resname, first_rg.resid(),
        sites[final][0], b_isos[final],
          last.name, last_ag.resname, last_rg.resid())
    print >> out, "|"+"-"*77+"|"
    print >> out
    out.flush()
    time_model_show += timer.elapsed()

  def remove_alternative_conformations(self, always_keep_one_conformer):
    # XXX This is not working correctly when something was deleted.
    # Need to figure out a way to update everything so GRM
    # construction will not fail.
    self.geometry_restraints = None
    self._pdb_hierarchy.remove_alt_confs(
        always_keep_one_conformer=always_keep_one_conformer)
    self._xray_structure = self._pdb_hierarchy.extract_xray_structure(crystal_symmetry=self.crystal_symmetry())
    self._atom_selection_cache = None
    n_old_atoms = self.get_number_of_atoms()
    self.pdb_atoms = self._pdb_hierarchy.atoms()
    n_new_atoms = self.get_number_of_atoms()
    return n_old_atoms - n_new_atoms

  def remove_solvent(self):
    result = self.select(selection = ~self.solvent_selection())
    return result

  def remove_hydrogens(self):
    if self.has_hd:
      noh_selection = self.selection("not (element H or element D)")
      return self.select(noh_selection)
    else:
      return self

  def show_occupancy_statistics(self, out=None, text=""):
    global time_model_show
    timer = user_plus_sys_time()
    # XXX make this more complete and smart
    if(out is None): out = sys.stdout
    print >> out, "|-"+text+"-"*(80 - len("| "+text+"|") - 1)+"|"
    occ = self._xray_structure.scatterers().extract_occupancies()
    # this needs to stay the same size - mask out HD selection
    less_than_zero = (occ < 0.0)
    occ_min = flex.min(occ)
    occ_max = flex.max(occ)
    n_zeros = (occ < 0.1).count(True)
    percent_small = n_zeros * 100. / occ.size()
    n_large = (occ > 2.0).count(True)
    if(percent_small > 30.0):
       print >> out, "| *** WARNING: more than 30 % of atoms with small occupancy (< 0.1)       *** |"
    if(n_large > 0):
       print >> out, "| *** WARNING: there are some atoms with large occupancy (> 2.0) ***          |"
    if(abs(occ_max-occ_min) >= 0.01):
       print >> out, "| occupancies: max = %-6.2f min = %-6.2f number of "\
                     "occupancies < 0.1 = %-6d |"%(occ_max,occ_min,n_zeros)
    else:
       print >> out, "| occupancies: max = %-6.2f min = %-6.2f number of "\
                     "occupancies < 0.1 = %-6d |"%(occ_max,occ_min,n_zeros)
    print >> out, "|"+"-"*77+"|"
    out.flush()
    time_model_show += timer.elapsed()

  def add_solvent(self, solvent_xray_structure,
                        atom_name    = "O",
                        residue_name = "HOH",
                        chain_id     = " ",
                        refine_occupancies = False,
                        refine_adp = None):
    assert refine_adp is not None
    n_atoms = solvent_xray_structure.scatterers().size()
    new_atom_name = atom_name.strip()
    if(len(new_atom_name) < 4): new_atom_name = " " + new_atom_name
    while(len(new_atom_name) < 4): new_atom_name = new_atom_name+" "
    atom_names = [ new_atom_name ] * n_atoms
    residue_names = [ residue_name ] * n_atoms
    nonbonded_types = flex.std_string([ "OH2" ] * n_atoms)
    i_seq = find_common_water_resseq_max(pdb_hierarchy=self._pdb_hierarchy)
    if (i_seq is None or i_seq < 0): i_seq = 0
    self.append_single_atoms(
      new_xray_structure=solvent_xray_structure,
      atom_names=atom_names,
      residue_names=residue_names,
      nonbonded_types=nonbonded_types,
      i_seq_start=i_seq,
      chain_id=chain_id,
      refine_adp=refine_adp,
      refine_occupancies=refine_occupancies,
      reset_labels=True, # Not clear if one forgot to do this or this
      # was not available at the time. Need investigation. Probably will
      # help to eliminate special water treatment in adopt_xray_structure()
      )
    self._update_has_hd()

  def append_single_atoms(self,
      new_xray_structure,
      atom_names,
      residue_names,
      nonbonded_types,
      refine_adp,
      refine_occupancies=None,
      nonbonded_charges=None,
      segids=None,
      i_seq_start = 0,
      chain_id     = " ",
      reset_labels=False) :
    assert refine_adp in ["isotropic", "anisotropic"]
    assert new_xray_structure.scatterers().size() == len(atom_names) == \
        len(residue_names) == len(nonbonded_types)
    if segids is not None:
      assert len(atom_names) == len(segids)
    ms = self._xray_structure.scatterers().size() #
    number_of_new_atoms = new_xray_structure.scatterers().size()
    self._xray_structure = \
      self._xray_structure.concatenate(new_xray_structure)
    occupancy_flags = None
    if(refine_occupancies):
      occupancy_flags = []
      for i in range(1, new_xray_structure.scatterers().size()+1):
        occupancy_flags.append([flex.size_t([ms+i-1])])
    if(self.refinement_flags is not None and
       self.refinement_flags.individual_sites):
      ssites = flex.bool(new_xray_structure.scatterers().size(), True)
    else: ssites = None
    # add flags
    if(self.refinement_flags is not None and
       self.refinement_flags.torsion_angles):
      ssites_tors = flex.bool(new_xray_structure.scatterers().size(), True)
    else: ssites_tors = None
    #
    sadp_iso, sadp_aniso = None, None
    if (refine_adp=="isotropic"):
      nxrs_ui = new_xray_structure.use_u_iso()
      if((self.refinement_flags is not None and
          self.refinement_flags.adp_individual_iso) or nxrs_ui.count(True)>0):
        sadp_iso = nxrs_ui
        sadp_aniso = flex.bool(sadp_iso.size(), False)
      else: sadp_iso = None
    if(refine_adp=="anisotropic"):
      nxrs_ua = new_xray_structure.use_u_aniso()
      if((self.refinement_flags is not None and
          self.refinement_flags.adp_individual_aniso) or nxrs_ua.count(True)>0):
        sadp_aniso = nxrs_ua
        sadp_iso = flex.bool(sadp_aniso.size(), False)
      else: sadp_aniso = None
    if(self.refinement_flags is not None):
      self.refinement_flags.inflate(
        sites_individual       = ssites,
        sites_torsion_angles   = ssites_tors,
        adp_individual_iso     = sadp_iso,
        adp_individual_aniso   = sadp_aniso,
        s_occupancies          = occupancy_flags,
        size_all               = ms)#torsion_angles
    #
    self._append_pdb_atoms(
      new_xray_structure=new_xray_structure,
      atom_names=atom_names,
      residue_names=residue_names,
      chain_id=chain_id,
      segids=segids,
      i_seq_start=i_seq_start,
      reset_labels=reset_labels)
   #
    if(self.restraints_manager is not None):
      geometry = self.restraints_manager.geometry
      if (geometry.model_indices is None):
        model_indices = None
      else:
        model_indices = flex.size_t(number_of_new_atoms, 0)
      if (geometry.conformer_indices is None):
        conformer_indices = None
      else:
        conformer_indices = flex.size_t(number_of_new_atoms, 0)
      if (geometry.sym_excl_indices is None):
        sym_excl_indices = None
      else:
        sym_excl_indices = flex.size_t(number_of_new_atoms, 0)
      if (geometry.donor_acceptor_excl_groups is None):
        donor_acceptor_excl_groups = None
      else:
        donor_acceptor_excl_groups = flex.size_t(number_of_new_atoms, 0)
      if (nonbonded_charges is None) :
        nonbonded_charges = flex.int(number_of_new_atoms, 0)
      geometry = geometry.new_including_isolated_sites(
        n_additional_sites =number_of_new_atoms,
        model_indices=model_indices,
        conformer_indices=conformer_indices,
        sym_excl_indices=sym_excl_indices,
        donor_acceptor_excl_groups=donor_acceptor_excl_groups,
        site_symmetry_table=new_xray_structure.site_symmetry_table(),
        nonbonded_types=nonbonded_types,
        nonbonded_charges=nonbonded_charges)
      self.restraints_manager = mmtbx.restraints.manager(
        geometry      = geometry,
        cartesian_ncs_manager    = self.restraints_manager.cartesian_ncs_manager,
        normalization = self.restraints_manager.normalization)
      c_ncs_m = self.get_cartesian_NCS_manager()
      if (c_ncs_m is not None):
        c_ncs_m.register_additional_isolated_sites(
          number=number_of_new_atoms)
      self.restraints_manager.geometry.update_plain_pair_sym_table(
        sites_frac = self._xray_structure.sites_frac())
    assert self.pdb_atoms.size() == self._xray_structure.scatterers().size()
    if self.riding_h_manager is not None:
      new_riding_h_manager = self.riding_h_manager.update(
        pdb_hierarchy       = self._pdb_hierarchy,
        geometry_restraints = geometry,
        n_new_atoms         = number_of_new_atoms)
      self.riding_h_manager = new_riding_h_manager

  def _append_pdb_atoms (self,
      new_xray_structure,
      atom_names,
      residue_names,
      chain_id,
      segids=None,
      i_seq_start=0,
      reset_labels=False) :
    """ Add atoms from new_xray_structure to the model in place."""
    assert (len(atom_names) == len(residue_names) ==
            len(new_xray_structure.scatterers()))
    assert (segids is None) or (len(segids) == len(atom_names))
    pdb_model = self._pdb_hierarchy.only_model()
    new_chain = iotbx.pdb.hierarchy.chain(id=chain_id)
    orth = new_xray_structure.unit_cell().orthogonalize
    n_seq = self.pdb_atoms.size()
    i_seq = i_seq_start
    for j_seq, sc in enumerate(new_xray_structure.scatterers()) :
      i_seq += 1
      element, charge = sc.element_and_charge_symbols()
      new_atom = (iotbx.pdb.hierarchy.atom()
        .set_serial(new_serial=iotbx.pdb.hy36encode(width=5, value=n_seq+i_seq))
        .set_name(new_name=atom_names[j_seq])
        .set_xyz(new_xyz=orth(sc.site))
        .set_occ(new_occ=sc.occupancy)
        .set_b(new_b=adptbx.u_as_b(sc.u_iso))
        .set_element(element)
        .set_charge(charge)
        .set_hetero(new_hetero=True))
      if (segids is not None) :
        new_atom.segid = segids[j_seq]
      new_atom_group = iotbx.pdb.hierarchy.atom_group(altloc="",
        resname=residue_names[j_seq])
      new_atom_group.append_atom(atom=new_atom)
      new_residue_group = iotbx.pdb.hierarchy.residue_group(
        resseq=iotbx.pdb.resseq_encode(value=i_seq), icode=" ")
      new_residue_group.append_atom_group(atom_group=new_atom_group)
      new_chain.append_residue_group(residue_group=new_residue_group)
    if (new_chain.residue_groups_size() != 0):
      pdb_model.append_chain(chain=new_chain)
    self._update_atom_selection_cache()
    # This deep_copy here for the following reason:
    # sometimes (particualrly in IAS refinement), after update_atom_selection_cache()
    # part of the model atoms loose their parent(). Not clear why.
    # deep_copy seems to help with it.
    self._pdb_hierarchy = self._pdb_hierarchy.deep_copy()
    self._update_pdb_atoms()
    if (reset_labels) :
      self._sync_xrs_labels()
    self.all_chain_proxies = None
    self._processed_pdb_file = None

  def _sync_xrs_labels(self):
    for sc, atom in zip(self.get_xray_structure().scatterers(), self.get_hierarchy().atoms()) :
      sc.label = atom.id_str()

  def convert_atom (self,
      i_seq,
      scattering_type,
      atom_name,
      element,
      charge,
      residue_name,
      initial_occupancy=None,
      initial_b_iso=None,
      chain_id=None,
      segid=None,
      refine_occupancies=True,
      refine_adp = None) :
    """
    Convert a single atom (usually water) to a different type, including
    adjustment of the xray structure and geometry restraints.
    """
    atom = self.pdb_atoms[i_seq]
    atom.name = atom_name
    atom.element = "%2s" % element.strip()
    assert (atom.element.strip() == element)
    if (charge != 0) :
      symbol = "+"
      if (charge < 0) : symbol = "-"
      atom.charge = str(abs(charge)) + symbol
    else :
      atom.charge = ""
    atom.parent().resname = residue_name
    if (chain_id is not None) :
      assert (len(chain_id) <= 2)
      atom.parent().parent().parent().id = chain_id
    if (segid is not None) :
      assert (len(segid) <= 4)
      atom.segid = segid
    scatterer = self._xray_structure.scatterers()[i_seq]
    scatterer.scattering_type = scattering_type
    label = atom.id_str()
    all_labels = [ s.label for s in self._xray_structure.scatterers() ]
    while (label in all_labels) :
      rg = atom.parent().parent()
      resseq = rg.resseq_as_int()
      rg.resseq = "%4d" % (resseq + 1)
      label = atom.id_str()
    # scatterer.label = atom.id_str() # will be done below for whole xrs
    if (initial_occupancy is not None) :
      # XXX preserve partial occupancies on special positions
      if (scatterer.occupancy != 1.0) :
        initial_occupancy = scatterer.occupancy
      scatterer.occupancy = initial_occupancy
      atom.occ = initial_occupancy
    if (initial_b_iso is not None) :
      atom.b = initial_b_iso
      scatterer.u_iso = adptbx.b_as_u(initial_b_iso)
    atom_selection = flex.size_t([i_seq])
    if(refine_adp == "isotropic"):
      scatterer.convert_to_isotropic(unit_cell=self._xray_structure.unit_cell())
      if ((self.refinement_flags is not None) and
          (self.refinement_flags.adp_individual_iso is not None)) :
        self.refinement_flags.adp_individual_iso.set_selected(atom_selection,
          True)
        if (self.refinement_flags.adp_individual_aniso is not None) :
          self.refinement_flags.adp_individual_aniso.set_selected(
            atom_selection, False)
    elif(refine_adp == "anisotropic"):
      scatterer.convert_to_anisotropic(
        unit_cell=self._xray_structure.unit_cell())
      if ((self.refinement_flags is not None) and
          (self.refinement_flags.adp_individual_aniso is not None)) :
        self.refinement_flags.adp_individual_aniso.set_selected(atom_selection,
          True)
        if (self.refinement_flags.adp_individual_iso is not None) :
          self.refinement_flags.adp_individual_iso.set_selected(atom_selection,
            False)
    if ((self.refinement_flags is not None) and
        (self.refinement_flags.sites_individual is not None)) :
      self.refinement_flags.sites_individual.set_selected(atom_selection, True)
    if ((self.refinement_flags is not None) and
        (self.refinement_flags.s_occupancies is not None)) :
      flagged = False
      for occgroup in self.refinement_flags.s_occupancies:
        for occsel in occgroup :
          if (i_seq in occsel) :
            flagged = True
            break
      if (not flagged) :
        self.refinement_flags.s_occupancies.append([atom_selection])
    self.restraints_manager.geometry.update_atom_nonbonded_type(
      i_seq=i_seq,
      nonbonded_type=element,
      charge=charge)
    self._xray_structure.discard_scattering_type_registry()
    assert self.pdb_atoms.size() == self._xray_structure.scatterers().size()
    self._sync_xrs_labels()
    self.set_xray_structure(self._xray_structure)
    self._update_pdb_atoms()
    return atom

  def scale_adp(self, scale_max, scale_min):
    b_isos = self._xray_structure.extract_u_iso_or_u_equiv() * math.pi**2*8
    b_isos_mean = flex.mean(b_isos)
    max_b_iso = b_isos_mean * scale_max
    min_b_iso = b_isos_mean / scale_min
    sel_outliers_max = b_isos > max_b_iso
    sel_outliers_min = b_isos < min_b_iso
    b_isos.set_selected(sel_outliers_max, max_b_iso)
    b_isos.set_selected(sel_outliers_min, min_b_iso)
    self._xray_structure.set_b_iso(values = b_isos)
    self.set_sites_cart_from_xrs()

  def get_model_statistics_info(self,
      fmodel_x          = None,
      fmodel_n          = None,
      refinement_params = None,
      general_selection = None,
      use_molprobity    = True):
    if self.model_statistics_info is None:
      self.model_statistics_info = mmtbx.model.statistics.info(
          model             = self,
          fmodel_x          = fmodel_x,
          fmodel_n          = fmodel_n,
          refinement_params = refinement_params,
          general_selection = general_selection,
          use_molprobity    = use_molprobity)
    return self.model_statistics_info

  def geometry_statistics(self,
                          general_selection = None):
    scattering_table = \
        self.get_xray_structure().scattering_type_registry().last_table()
    if(self.restraints_manager is None): return None
    ph = self.get_hierarchy(sync_with_xray_structure=True)
    hd_selection = self.get_hd_selection()
    if(self.use_ias):
      ias_selection = self.get_ias_selection()
      hd_selection = hd_selection.select(~ias_selection)
      ph = ph.select(~ias_selection)
    rm = self.restraints_manager
    if(self.riding_h_manager is not None or
       scattering_table in ["n_gaussian","wk1995", "it1992"]):
      not_hd_sel = ~hd_selection
      ph = ph.select(not_hd_sel)
      rm = rm.select(not_hd_sel)
    return mmtbx.model.statistics.geometry(
      pdb_hierarchy               = ph,
      geometry_restraints_manager = rm.geometry)

  def occupancy_statistics(self):
    return mmtbx.model.statistics.occupancy(
      hierarchy = self.get_hierarchy(sync_with_xray_structure=True)).result

  def adp_statistics(self):
    return mmtbx.model.statistics.adp(model = self)

  def show_adp_statistics(self,
                          prefix         = "",
                          padded         = False,
                          pdb_deposition = False,
                          out            = None):
    global time_model_show
    if(out is None): out = self.log
    timer = user_plus_sys_time()
    result = self.adp_statistics()
    result.show(out = out, prefix = prefix, padded = padded,
      pdb_deposition = pdb_deposition)
    time_model_show += timer.elapsed()
    return result

  def energies_adp(self, iso_restraints, compute_gradients, use_hd):
    assert self.refinement_flags is not None
    xrs = self._xray_structure
    sel_ = xrs.use_u_iso() | xrs.use_u_aniso()
    selection = sel_
    ias_selection = self.get_ias_selection()
    if(ias_selection is not None and ias_selection.count(True) > 0):
      selection = sel_.set_selected(ias_selection, False)
    n_aniso = 0
    if(self.refinement_flags.adp_individual_aniso is not None):
      n_aniso = self.refinement_flags.adp_individual_aniso.count(True)
    if(n_aniso == 0):
      energies_adp_iso = self.restraints_manager.energies_adp_iso(
        xray_structure    = xrs,
        parameters        = iso_restraints,
        use_u_local_only  = iso_restraints.use_u_local_only,
        use_hd            = use_hd,
        compute_gradients = compute_gradients)
      target = energies_adp_iso.target
    else:
      energies_adp_aniso = self.restraints_manager.energies_adp_aniso(
        xray_structure    = xrs,
        compute_gradients = compute_gradients,
        selection         = selection,
        use_hd            = use_hd)
      target = energies_adp_aniso.target
    u_iso_gradients = None
    u_aniso_gradients = None
    if(compute_gradients):
      if(n_aniso == 0):
        u_iso_gradients = energies_adp_iso.gradients
      else:
        u_aniso_gradients = energies_adp_aniso.gradients_aniso_star
        u_iso_gradients = energies_adp_aniso.gradients_iso
    class result(object):
      def __init__(self):
        self.target = target
        self.u_iso_gradients = u_iso_gradients
        self.u_aniso_gradients = u_aniso_gradients
    return result()

  def set_refine_individual_sites(self, selection = None):
    self._xray_structure.scatterers().flags_set_grads(state=False)
    if(selection is None):
      if(self.refinement_flags is not None):
        selection = self.refinement_flags.sites_individual
    if(selection is not None):
      self._xray_structure.scatterers().flags_set_grad_site(
        iselection = selection.iselection())

  def set_refine_individual_adp(self, selection_iso = None,
                                      selection_aniso = None):
    self._xray_structure.scatterers().flags_set_grads(state=False)
    if(selection_iso is None):
      selection_iso = self.refinement_flags.adp_individual_iso
    if(selection_iso is not None):
      self._xray_structure.scatterers().flags_set_grad_u_iso(
        iselection = selection_iso.iselection())
    if(selection_aniso is None):
      selection_aniso = self.refinement_flags.adp_individual_aniso
    if(selection_aniso is not None):
      self._xray_structure.scatterers().flags_set_grad_u_aniso(
        iselection = selection_aniso.iselection())

  def _expand_symm_helper(self, records_container):
    """
    This will expand hierarchy and ss annotations. In future anything else that
    should be expanded have to be added here. e.g. TLS.
    LIMITATION: ANISOU records in resulting hierarchy will be invalid!!!
    """
    from iotbx.pdb.utils import all_chain_ids
    roots=[]
    all_cids = all_chain_ids()
    duplicate_prevention = {}
    chain_ids_match_dict = {} # {'old chain id': [new ids]}
    for m in self.get_hierarchy().models():
      for c in m.chains():
        chain_id_key = "%s%s" % (m.id, c.id)
        if chain_id_key in duplicate_prevention:
          continue
        duplicate_prevention[chain_id_key] = False
        chain_ids_match_dict[c.id] = []
        cid = c.id if len(c.id) == 2 else " "+c.id
        try:
          ind = all_cids.index(cid)
        except ValueError:
          ind = -1
        if ind >= 0:
          del all_cids[ind]
    cid_counter = 0
    for r,t in zip(records_container.r, records_container.t):
      leave_chain_ids = False
      if r.is_r3_identity_matrix() and t.is_col_zero():
        leave_chain_ids = True
      for mm in self.get_hierarchy().models():
        root = iotbx.pdb.hierarchy.root()
        m = iotbx.pdb.hierarchy.model()
        for k in duplicate_prevention.keys():
          duplicate_prevention[k] = False
        for c in mm.chains():
          c = c.detached_copy()
          if not leave_chain_ids and not duplicate_prevention["%s%s" % (mm.id, c.id)]:
            new_cid = all_cids[cid_counter]
            chain_ids_match_dict[c.id].append(new_cid)
            duplicate_prevention["%s%s" % (mm.id, c.id)] = True
            c.id = new_cid
            cid_counter += 1
          xyz = c.atoms().extract_xyz()
          new_xyz = r.elems*xyz+t
          c.atoms().set_xyz(new_xyz)
          m.append_chain(c)
        root.append_model(m)
        roots.append(root)
    result = iotbx.pdb.hierarchy.root()
    for rt in roots:
      result.transfer_chains_from_other(other=rt)
    #validation
    vals = chain_ids_match_dict.values()
    for v in vals:
      assert len(vals[0]) == len(v), chain_ids_match_dict
    result.reset_i_seq_if_necessary()
    self._pdb_hierarchy = result

    # Now deal with SS annotations
    ssa = self.get_ss_annotation()
    if ssa is not None:
      ssa.multiply_to_asu_2(chain_ids_match_dict)
      self.set_ss_annotation(ssa)
    self._update_pdb_atoms()
    self._xray_structure = None
    self._all_chain_proxies = None
    self._update_atom_selection_cache()

  def _biomt_mtrix_container_is_good(self, records_container):
    if records_container is None:
      return False
    if len(records_container.r)==1:
      r, t = records_container.r[0], records_container.t[0]
      if r.is_r3_identity_matrix() and t.is_col_zero():
        return False
    if (records_container.is_empty() or
        len(records_container.r) == 0 or
        records_container.validate()):
      return False
    return True

  def expand_with_MTRIX_records(self):
    mtrix_records_container = self._model_input.process_MTRIX_records()
    if not self._biomt_mtrix_container_is_good(mtrix_records_container):
      return
    if self.get_hierarchy() is None:
      self._pdb_hierarchy = deepcopy(self._model_input).construct_hierarchy(
          self._pdb_interpretation_params.pdb_interpretation.sort_atoms)
    self._expand_symm_helper(mtrix_records_container)

  def expand_with_BIOMT_records(self):
    """
    expanding current hierarchy and ss_annotations with BIOMT matrices.
    Known limitations: will expand everything, regardless of what selections
    were setted in BIOMT header.
    """
    biomt_records_container = self._model_input.process_BIOMT_records()
    if not self._biomt_mtrix_container_is_good(biomt_records_container):
      return
    self._expand_symm_helper(biomt_records_container)
