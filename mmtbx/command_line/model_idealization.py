from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.model_idealization

import sys, os
from iotbx.pdb import secondary_structure as ioss
from mmtbx.secondary_structure import build as ssb
from mmtbx.secondary_structure import manager, sec_str_master_phil
import mmtbx.utils
from scitbx.array_family import flex
from iotbx.pdb import write_whole_pdb_file
import iotbx.phil
from libtbx.utils import Sorry
from mmtbx.building.loop_idealization import loop_idealization
import mmtbx.building.loop_closure.utils
from mmtbx.refinement.geometry_minimization import minimize_wrapper_for_ramachandran
import iotbx.ncs
from copy import deepcopy
from mmtbx.utils import fix_rotamer_outliers
from mmtbx.command_line.geometry_minimization import get_geometry_restraints_manager
from mmtbx.model_statistics import geometry_no_grm
from time import time
import datetime
from libtbx.utils import multi_out
from cctbx import maptbx, miller
from mmtbx.refinement.real_space.individual_sites import minimize_wrapper_with_map
from mmtbx.pdbtools import truncate_to_poly_gly
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from mmtbx.rotamer.rotamer_eval import RotamerEval
from mmtbx_validation_ramachandran_ext import rama_eval
from iotbx.file_reader import any_file

turned_on_ss = ssb.ss_idealization_master_phil_str
turned_on_ss = turned_on_ss.replace("enabled = False", "enabled = True")
master_params_str = """
file_name = None
  .type = path
  .multiple = True
  .short_caption = Model file
  .style = file_type:pdb bold input_file
map_file_name = None
  .type = path
  .help = User-provided map that will be used as reference
trim_alternative_conformations = False
  .type = bool
  .help = Leave only atoms with empty altloc
additionally_fix_rotamer_outliers = True
  .type = bool
  .help = At the late stage if rotamer is still outlier choose another one \
    with minimal clash with surrounding atoms
use_ss_restraints = True
  .type = bool
  .help = Use Secondary Structure restraints
use_starting_model_for_final_gm = False
  .type = bool
  .help = Use supplied model for final geometry minimization. Otherwise just \
    use self.
output_prefix = None
  .type = str
use_map_for_reference = True
  .type = bool
run_minimization_first = False
  .type = bool
reference_map_resolution = 5
  .type = float
data_for_map = None
  .type = path
number_of_refinement_cycles = 3
  .type = int
ignore_ncs = False
  .type = bool
  .help = Don't use NCS even if it is present in model.
debug = False
  .type = bool
  .help = Output all intermediate files
%s
include scope mmtbx.secondary_structure.sec_str_master_phil_str
include scope mmtbx.building.loop_idealization.loop_idealization_master_phil_str
""" % turned_on_ss

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def format_usage_message(log):
  print >> log, "-"*79
  msg = """\
phenix.model_idealization: Idealize model geometry.
(Only protein secondary structure elements are supported currently).

Usage examples:
 phenix.model_idealization model.pdb
"""
  print >> log, msg
  print >> log, "-"*79
  print >> log, master_params().show()

class model_idealization():
  def __init__(self,
               pdb_input,
               cif_objects=None,
               map_data = None,
               crystal_symmetry = None,
               params=None,
               log=sys.stdout,
               verbose=True):
    t_0 = time()
    self.pdb_input = pdb_input
    self.cif_objects = cif_objects
    self.params = params
    self.log = log
    self.verbose = verbose

    self.rmsd_from_start = None
    self.init_model_statistics = None
    self.init_gm_model_statistics = None
    self.after_ss_idealization = None
    self.after_loop_idealization = None
    self.after_rotamer_fixing = None
    self.final_model_statistics = None
    self.reference_map = map_data
    self.master_map = None

    self.whole_grm = None
    self.master_grm = None
    self.working_grm = None

    self.mon_lib_srv = None
    self.ener_lib = None
    self.rotamer_manager = None
    self.rama_manager = rama_eval()

    self.original_hierarchy = None # original pdb_h, without any processing
    self.original_boxed_hierarchy = None # original and boxed (if needed)
    self.whole_pdb_h = None # boxed with processing (AC trimming, H trimming,...)
    self.master_pdb_h = None # master copy in case of NCS
    self.working_pdb_h = None # one to use for fixing (master_pdb_h or working_pdb_h)

    self.using_ncs = False
    self.ncs_restr_group_list = []

    # various checks, shifts, trims
    self.cs = crystal_symmetry
    # check self.cs (copy-paste from secondary_sturcure_restraints)
    corrupted_cs = False
    if self.cs is not None:
      if [self.cs.unit_cell(), self.cs.space_group()].count(None) > 0:
        corrupted_cs = True
        self.cs = None
      elif self.cs.unit_cell().volume() < 10:
        corrupted_cs = True
        self.cs = None

    self.original_hierarchy = self.pdb_input.construct_hierarchy()
    # couple checks if pdb_h is ok
    o_c = self.original_hierarchy.overall_counts()
    o_c.raise_duplicate_atom_labels_if_necessary()
    o_c.raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary()
    o_c.raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary()
    o_c.raise_improper_alt_conf_if_necessary()
    if len(self.original_hierarchy.models()) > 1:
      raise Sorry("Multi model files are not supported")
    ca_only_present = False
    for c in self.original_hierarchy.only_model().chains():
      if c.is_ca_only():
        ca_only_present = True
    if ca_only_present:
      raise Sorry("Don't support models with chains containing only CA atoms.")

    self.original_boxed_hierarchy = self.original_hierarchy.deep_copy()
    self.original_boxed_hierarchy.reset_atom_i_seqs()
    self.shift_vector = None
    if self.cs is None:
      if corrupted_cs:
        print >> self.log, "Symmetry information is corrupted, "
      else:
        print >> self.log, "Symmetry information was not found, "
      print >> self.log, "putting molecule in P1 box."
      self.log.flush()
      from cctbx import uctbx
      atoms = self.original_boxed_hierarchy.atoms()
      box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
        sites_cart=atoms.extract_xyz(),
        buffer_layer=3)
      atoms.set_xyz(new_xyz=box.sites_cart)
      self.cs = box.crystal_symmetry()
      self.shift_vector = box.shift_vector

    # self.original_boxed_hierarchy.write_pdb_file(file_name="original_boxed_h.pdb")
    if self.shift_vector is not None and self.params.debug:
      write_whole_pdb_file(
          file_name="%s_boxed.pdb" % self.params.output_prefix,
          pdb_hierarchy=self.original_boxed_hierarchy,
          crystal_symmetry=self.cs,
          ss_annotation=self.pdb_input.extract_secondary_structure())

    asc = self.original_boxed_hierarchy.atom_selection_cache()
    if self.params.trim_alternative_conformations:
      sel = asc.selection("altloc ' '")
      self.whole_pdb_h = self.original_boxed_hierarchy.select(sel).deep_copy()
      print >> self.log, "Atoms in original/working model: %d/%d" % (
          self.original_boxed_hierarchy.atoms_size(), self.whole_pdb_h.atoms_size())
    else:
      self.whole_pdb_h = self.original_boxed_hierarchy.deep_copy()
      # self.whole_pdb_h.reset_atom_i_seqs()
    # Trimming hydrogens
    # Many intermediate variables are needed due to strange behavior of
    # selections described in
    # iotbx/pdb/tst_hierarchy.py:exercise_selection_and_deep_copy()
    asc2 = self.whole_pdb_h.atom_selection_cache()
    h_sel = asc2.selection("not (element H or element D)")
    temp_h = self.whole_pdb_h.select(h_sel)
    self.whole_pdb_h = temp_h.deep_copy()
    self.whole_pdb_h.reset_atom_i_seqs()
    # self.init_model_statistics = geometry_no_grm(
    #     pdb_hierarchy=iotbx.pdb.input(
    #       source_info=None,
    #       lines=self.whole_pdb_h.as_pdb_string()).construct_hierarchy(),
    #     molprobity_scores=True)
    for_stat_h = self.get_intermediate_result_hierarchy()
    self.init_model_statistics = geometry_no_grm(
        pdb_hierarchy=for_stat_h,
        molprobity_scores=True)
    self.time_for_init = time()-t_0

  def get_intermediate_result_hierarchy(self):
    result_h = self.whole_pdb_h.deep_copy()
    result_h.reset_atom_i_seqs()
    if self.using_ncs:
      # multiply back and do geometry_minimization for the whole molecule
      for ncs_gr in self.ncs_restr_group_list:
        ssb.set_xyz_smart(result_h, self.working_pdb_h)
        new_sites = result_h.select(ncs_gr.master_iselection).atoms().extract_xyz()
        # print len(ncs_gr.master_iselection), result_h.atoms_size(), len(new_sites)
        # result_h.select(ncs_gr.master_iselection).atoms().set_xyz(new_sites)
        for c in ncs_gr.copies:
          new_c_sites = c.r.elems * new_sites + c.t
          result_h.select(c.iselection).atoms().set_xyz(new_c_sites)
    elif self.working_pdb_h is not None:
      result_h = self.working_pdb_h.deep_copy()
      result_h.reset_atom_i_seqs()
    return result_h

  def prepare_init_reference_map(self, xrs, pdb_h):
    print >> self.log, "Preparing reference map"
    # new_h = pdb_h.deep_copy()
    # truncate_to_poly_gly(new_h)
    # xrs = new_h.extract_xray_structure(crystal_symmetry=xrs.crystal_symmetry())
    xrs=xrs.set_b_iso(value=10)
    crystal_gridding = maptbx.crystal_gridding(
        unit_cell        = xrs.unit_cell(),
        space_group_info = xrs.space_group_info(),
        symmetry_flags   = maptbx.use_space_group_symmetry,
        d_min             = self.params.reference_map_resolution)
    fc = xrs.structure_factors(d_min = 2, algorithm = "fft").f_calc()
    fft_map = miller.fft_map(
        crystal_gridding=crystal_gridding,
        fourier_coefficients=fc)
    fft_map.apply_sigma_scaling()
    init_reference_map = fft_map.real_map_unpadded(in_place=False)
    if self.params.debug:
      fft_map.as_xplor_map(file_name="%s_init.map" % self.params.output_prefix)
    return init_reference_map

  def prepare_reference_map(self, xrs, pdb_h):
    print >> self.log, "Preparing reference map"
    # new_h = pdb_h.deep_copy()
    # truncate_to_poly_gly(new_h)
    # xrs = new_h.extract_xray_structure(crystal_symmetry=xrs.crystal_symmetry())
    xrs=xrs.set_b_iso(value=50)
    crystal_gridding = maptbx.crystal_gridding(
        unit_cell        = xrs.unit_cell(),
        space_group_info = xrs.space_group_info(),
        symmetry_flags   = maptbx.use_space_group_symmetry,
        d_min             = self.params.reference_map_resolution)
    fc = xrs.structure_factors(d_min = self.params.reference_map_resolution, algorithm = "fft").f_calc()
    fft_map = miller.fft_map(
        crystal_gridding=crystal_gridding,
        fourier_coefficients=fc)
    fft_map.apply_sigma_scaling()
    self.reference_map = fft_map.real_map_unpadded(in_place=False)
    if self.params.debug:
      fft_map.as_xplor_map(file_name="%s.map" % self.params.output_prefix)

  def prepare_reference_map_2(self, xrs, pdb_h):
    print >> self.log, "Preparing reference map, method 2"
    # new_h = pdb_h.deep_copy()
    # truncate_to_poly_gly(new_h)
    # xrs = new_h.extract_xray_structure(crystal_symmetry=xrs.crystal_symmetry())
    xrs=xrs.set_b_iso(value=50)

    # side_chain_no_cb_selection = ~ xrs.main_chain_selection()
    side_chain_no_cb_selection = ~ xrs.backbone_selection()
    xrs = xrs.set_b_iso(value=200, selection=side_chain_no_cb_selection)

    crystal_gridding = maptbx.crystal_gridding(
        unit_cell        = xrs.unit_cell(),
        space_group_info = xrs.space_group_info(),
        symmetry_flags   = maptbx.use_space_group_symmetry,
        d_min             = self.params.reference_map_resolution)
    fc = xrs.structure_factors(d_min = self.params.reference_map_resolution, algorithm = "direct").f_calc()
    fft_map = miller.fft_map(
        crystal_gridding=crystal_gridding,
        fourier_coefficients=fc)
    fft_map.apply_sigma_scaling()
    self.reference_map = fft_map.real_map_unpadded(in_place=False)
    if self.params.debug:
      fft_map.as_xplor_map(file_name="%s_2.map" % self.params.output_prefix)

  def prepare_reference_map_3(self, xrs, pdb_h):
    """ with ramachandran outliers """
    print >> self.log, "Preparing reference map, method 3"
    outlier_selection_txt = mmtbx.building.loop_closure.utils. \
          rama_outliers_selection(pdb_h, self.rama_manager, 1)
    asc = pdb_h.atom_selection_cache()
    print >> self.log, "rama outlier selection:", outlier_selection_txt
    rama_out_sel = asc.selection(outlier_selection_txt)
    xrs=xrs.set_b_iso(value=50)

    # side_chain_no_cb_selection = ~ xrs.main_chain_selection()
    side_chain_no_cb_selection = ~ xrs.backbone_selection()
    xrs = xrs.set_b_iso(value=200, selection=side_chain_no_cb_selection)
    xrs = xrs.set_b_iso(value=150, selection=rama_out_sel)
    # xrs = xrs.set_occupancies(value=0.3, selection=rama_out_sel)

    crystal_gridding = maptbx.crystal_gridding(
        unit_cell        = xrs.unit_cell(),
        space_group_info = xrs.space_group_info(),
        symmetry_flags   = maptbx.use_space_group_symmetry,
        d_min             = self.params.reference_map_resolution)
    fc = xrs.structure_factors(
        d_min = self.params.reference_map_resolution,
        algorithm = "direct").f_calc()
    fft_map = miller.fft_map(
        crystal_gridding=crystal_gridding,
        fourier_coefficients=fc)
    fft_map.apply_sigma_scaling()
    self.reference_map = fft_map.real_map_unpadded(in_place=False)
    if self.params.debug:
      fft_map.as_xplor_map(file_name="%s_3.map" % self.params.output_prefix)
    self.master_map = self.reference_map.deep_copy()
    if self.using_ncs and self.master_pdb_h is not None:
      # here we are negating non-master part of the model
      # self.master_sel=master_sel
      # self.master_map = self.reference_map.deep_copy()
      print type(self.reference_map), dir(self.reference_map)
      print type(self.master_map), dir(self.master_map)
      mask = maptbx.mask(
              xray_structure=xrs.select(self.master_sel),
              n_real=self.master_map.focus(),
              mask_value_inside_molecule=1,
              mask_value_outside_molecule=-1,
              solvent_radius=0,
              atom_radius=1.)
      self.master_map = self.reference_map * mask
      if self.params.debug:
        iotbx.ccp4_map.write_ccp4_map(
            file_name="%s_3_master.map" % self.params.output_prefix,
            unit_cell=xrs.unit_cell(),
            space_group=xrs.space_group(),
            map_data=self.master_map,
            labels=flex.std_string([""]))
        # fft_map.as_xplor_map(file_name="%s_3.map" % self.params.output_prefix)
    # STOP()

  def get_grm(self):
    # first make whole grm using self.whole_pdb_h
    params_line = grand_master_phil_str
    params = iotbx.phil.parse(
        input_string=params_line, process_includes=True).extract()
    params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
    params.pdb_interpretation.peptide_link.ramachandran_restraints = True
    params.pdb_interpretation.peptide_link.oldfield.weight_scale=3
    params.pdb_interpretation.peptide_link.oldfield.plot_cutoff=0.03
    params.pdb_interpretation.nonbonded_weight = 500
    params.pdb_interpretation.c_beta_restraints=True
    params.pdb_interpretation.max_reasonable_bond_distance = None
    params.pdb_interpretation.peptide_link.apply_peptide_plane = True
    params.pdb_interpretation.ncs_search.enabled = True
    params.pdb_interpretation.restraints_library.rdl = True
    if self.params.ignore_ncs:
      params.pdb_interpretation.ncs_search.enabled = False
    processed_pdb_files_srv = mmtbx.utils.\
        process_pdb_file_srv(
            crystal_symmetry= self.whole_xrs.crystal_symmetry(),
            pdb_interpretation_params = params.pdb_interpretation,
            stop_for_unknowns         = False,
            log=self.log,
            cif_objects=None)
    processed_pdb_file, junk = processed_pdb_files_srv.\
        process_pdb_files(raw_records=flex.split_lines(self.whole_pdb_h.as_pdb_string()))

    self.mon_lib_srv = processed_pdb_files_srv.mon_lib_srv
    self.ener_lib = processed_pdb_files_srv.ener_lib
    self.rotamer_manager = RotamerEval(mon_lib_srv=self.mon_lib_srv)

    self.whole_grm = get_geometry_restraints_manager(
        processed_pdb_file, self.whole_xrs, params=params)

    # set SS restratins
    if self.params.use_ss_restraints:
      ss_manager = manager(
          pdb_hierarchy=self.whole_pdb_h,
          geometry_restraints_manager=self.whole_grm.geometry,
          sec_str_from_pdb_file=self.filtered_whole_ann,
          params=None,
          mon_lib_srv=self.mon_lib_srv,
          verbose=-1,
          log=self.log)
      # self.whole_pdb_h.write_pdb_file(file_name="for_ss.pdb")
      self.whole_pdb_h.reset_atom_i_seqs()
      self.whole_grm.geometry.set_secondary_structure_restraints(
          ss_manager=ss_manager,
          hierarchy=self.whole_pdb_h,
          log=self.log)

    # now select part of it for working with master hierarchy
    if self.using_ncs:
      self.master_grm = self.whole_grm.select(self.master_sel)
      self.working_grm = self.master_grm
    else:
      self.working_grm = self.whole_grm

  def run(self):
    t_0 = time()
    self.ann = ioss.annotation.from_phil(
        phil_helices=self.params.secondary_structure.protein.helix,
        phil_sheets=self.params.secondary_structure.protein.sheet,
        pdb_hierarchy=self.whole_pdb_h)
    if self.ann.get_n_helices() + self.ann.get_n_sheets() == 0:
      self.ann = self.pdb_input.extract_secondary_structure()
    self.filtered_whole_ann = None
    if self.ann is not None:
      self.filtered_whole_ann = self.ann.deep_copy()

    filtered_ncs_restr_group_list = []
    if not self.params.ignore_ncs:
      ncs_obj = iotbx.ncs.input(
          hierarchy=self.whole_pdb_h,
          chain_max_rmsd=4.0,
          chain_similarity_threshold=0.99,
          residue_match_radius=999.0)
      print >> self.log, "Found NCS groups:"
      ncs_obj.show(format='phil', log=self.log)
      self.ncs_restr_group_list = ncs_obj.get_ncs_restraints_group_list(
          raise_sorry=False)
      total_ncs_selected_atoms = 0
      master_sel = flex.size_t([])
      filtered_ncs_restr_group_list = self.filter_ncs_restraints_group_list(
          self.whole_pdb_h, self.ncs_restr_group_list)

    if self.params.run_minimization_first:
      # running simple minimization and updating all
      # self.master, self.working, etc...
      self.whole_pdb_h.reset_atom_i_seqs()
      init_ref_map = self.prepare_init_reference_map(
        self.whole_pdb_h.extract_xray_structure(crystal_symmetry=self.cs),
        self.whole_pdb_h)
      self.minimize(
          hierarchy=self.whole_pdb_h,
          xrs=self.whole_pdb_h.extract_xray_structure(crystal_symmetry=self.cs),
          original_pdb_h=self.whole_pdb_h,
          grm=self.whole_grm,
          ncs_restraints_group_list=filtered_ncs_restr_group_list,
          excl_string_selection=None, # don't need if we have map
          ss_annotation=self.ann,
          reference_map=init_ref_map)
      self.init_gm_model_statistics = geometry_no_grm(
          pdb_hierarchy=self.whole_pdb_h,
          molprobity_scores=True)
      if self.params.debug:
        self.shift_and_write_result(
            hierarchy=self.whole_pdb_h,
            fname_suffix="init_gm",
            grm=self.whole_grm)

    if (self.init_gm_model_statistics is not None
        and self.init_gm_model_statistics.ramachandran_outliers == 0):
      print >> self.log, "Simple minimization was enough"
      # Early exit!!!
      self.shift_and_write_result(
          hierarchy=self.whole_pdb_h,
          fname_suffix="all_idealized",
          grm=self.whole_grm)
      self.final_model_statistics = geometry_no_grm(
          pdb_hierarchy=iotbx.pdb.input(
            source_info=None,
            lines=self.whole_pdb_h.as_pdb_string()).construct_hierarchy(),
          molprobity_scores=True)
      # self.original_boxed_hierarchy.write_pdb_file(file_name="original_boxed_end.pdb")
      self.time_for_run = time() - t_0
      return

    if not self.params.ignore_ncs:
      if len(filtered_ncs_restr_group_list) > 0:
        self.using_ncs = True
        master_sel = flex.bool(self.whole_pdb_h.atoms_size(), True)
        for ncs_gr in filtered_ncs_restr_group_list:
          for copy in ncs_gr.copies:
            master_sel.set_selected(copy.iselection, False)
        self.master_pdb_h = self.whole_pdb_h.select(master_sel).deep_copy()
        self.master_sel=master_sel
        self.master_pdb_h.reset_atom_i_seqs()

    if self.using_ncs:
      if self.params.debug:
        self.master_pdb_h.write_pdb_file("%s_master_h.pdb" % self.params.output_prefix)
      self.working_pdb_h = self.master_pdb_h
    else:
      self.working_pdb_h = self.whole_pdb_h
    self.working_pdb_h.reset_atom_i_seqs()


    self.working_xrs = self.working_pdb_h.extract_xray_structure(crystal_symmetry=self.cs)
    if self.using_ncs:
      self.whole_xrs = self.whole_pdb_h.extract_xray_structure(crystal_symmetry=self.cs)
    else:
      self.whole_xrs = self.working_xrs

    if self.reference_map is None and self.params.use_map_for_reference:
      # self.prepare_reference_map(xrs=self.whole_xrs, pdb_h=self.whole_pdb_h)
      # self.prepare_reference_map_2(xrs=self.whole_xrs, pdb_h=self.whole_pdb_h)
      self.prepare_reference_map_3(xrs=self.whole_xrs, pdb_h=self.whole_pdb_h)


    self.original_ann = None
    if self.ann is not None:
      self.original_ann = self.ann.deep_copy()
      print >> self.log, "Original SS annotation"
      print >> self.log, self.original_ann.as_pdb_str()
      self.ann.remove_short_annotations()
      self.filtered_whole_ann = self.ann.deep_copy()
      self.ann.remove_3_10_helices()
      self.filtered_whole_ann.remove_3_10_helices()
      self.ann.remove_empty_annotations(
          hierarchy=self.working_pdb_h)
      self.filtered_whole_ann.remove_empty_annotations(
          hierarchy=self.whole_pdb_h)
      self.ann.concatenate_consecutive_helices(
          hierarchy=self.whole_pdb_h)
      self.filtered_whole_ann.concatenate_consecutive_helices(
          hierarchy=self.whole_pdb_h)
      self.ann.split_helices_with_prolines(
          hierarchy=self.working_pdb_h,
          asc=None)
      self.filtered_whole_ann.split_helices_with_prolines(
          hierarchy=self.whole_pdb_h,
          asc=None)
      self.ann.filter_sheets_with_long_hbonds(
          hierarchy=self.working_pdb_h,
          asc=None)
      self.filtered_whole_ann.filter_sheets_with_long_hbonds(
          hierarchy=self.whole_pdb_h,
          asc=None)
      # print >> self.log, "Splitted SS annotation"
      # print >> self.log, ann.as_pdb_str()
      print >> self.log, "Filtered SS annotation"
      print >> self.log, self.ann.as_pdb_str()

    # getting grm with SS restraints
    self.get_grm()

    if (self.ann is None or
        self.ann.get_n_helices() + self.ann.get_n_sheets() == 0 or
        not self.params.ss_idealization.enabled):
      print >> self.log, "No secondary structure annotations found or SS idealization is disabled."
      print >> self.log, "Secondary structure substitution step will be skipped"
      self.log.flush()
      # here we want to do geometry minimization anyway!
      negate_selection = None
      if self.reference_map is None:
        outlier_selection_txt = mmtbx.building.loop_closure.utils. \
          rama_outliers_selection(self.working_pdb_h, self.rama_manager, 1)
        print >> self.log, "outlier_selection_txt", outlier_selection_txt
        negate_selection = "all"
        if outlier_selection_txt != "" and outlier_selection_txt is not None:
          negate_selection = "not (%s)" % outlier_selection_txt
      self.minimize(
          hierarchy=self.whole_pdb_h,
          xrs=self.whole_xrs,
          original_pdb_h=self.whole_pdb_h,
          grm=self.whole_grm,
          ncs_restraints_group_list=filtered_ncs_restr_group_list,
          excl_string_selection=negate_selection,
          ss_annotation=self.ann,
          reference_map=self.reference_map)
      # self.original_boxed_hierarchy.write_pdb_file(file_name="original_boxed_h_1.pdb")
    else:
      if self.params.debug:
        self.params.ss_idealization.file_name_before_regularization = \
            "%s_ss_before_reg.pdb" % self.params.output_prefix
      self.params.ss_idealization.skip_good_ss_elements = True
      ssb.substitute_ss(
          real_h=self.working_pdb_h,
          xray_structure=self.working_xrs,
          ss_annotation=self.ann,
          params=self.params.ss_idealization,
          grm=self.working_grm,
          fix_rotamer_outliers=True,
          cif_objects=self.cif_objects,
          verbose=True,
          reference_map=self.master_map,
          rotamer_manager=self.rotamer_manager,
          log=self.log)
      self.log.flush()

    for_stat_h = self.get_intermediate_result_hierarchy()
    self.after_ss_idealization = geometry_no_grm(
        pdb_hierarchy=for_stat_h,
        molprobity_scores=True)
    self.shift_and_write_result(
          hierarchy=for_stat_h,
          fname_suffix="ss_ideal_stat",
          grm=self.whole_grm)
    # self.after_ss_idealization = geometry_no_grm(
    #     pdb_hierarchy=iotbx.pdb.input(
    #       source_info=None,
    #       lines=self.working_pdb_h.as_pdb_string()).construct_hierarchy(),
    #     molprobity_scores=True)

    # Write resulting pdb file.
    if self.params.debug:
      self.shift_and_write_result(
          hierarchy=self.working_pdb_h,
          fname_suffix="ss_ideal",
          grm=self.working_grm)
    self.params.loop_idealization.minimize_whole = not self.using_ncs
    # self.params.loop_idealization.enabled = False
    # self.params.loop_idealization.variant_search_level = 0
    loop_ideal = loop_idealization(
        pdb_hierarchy=self.working_pdb_h,
        params=self.params.loop_idealization,
        secondary_structure_annotation=self.ann,
        reference_map=self.master_map,
        crystal_symmetry=self.working_xrs.crystal_symmetry(),
        grm=self.working_grm,
        rama_manager=self.rama_manager,
        rotamer_manager=self.rotamer_manager,
        log=self.log,
        verbose=True)
    self.log.flush()
    self.working_pdb_h = loop_ideal.resulting_pdb_h
    if self.params.debug:
      self.shift_and_write_result(
          hierarchy=self.working_pdb_h,
          fname_suffix="rama_ideal",
          grm=self.working_grm)
    for_stat_h = self.get_intermediate_result_hierarchy()
    # for_stat_h.write_pdb_file(file_name="compare_with_rama_ideal.pdb")
    self.after_loop_idealization = geometry_no_grm(
        pdb_hierarchy=for_stat_h,
        molprobity_scores=True)

    # self.after_loop_idealization = geometry_no_grm(
    #     pdb_hierarchy=iotbx.pdb.input(
    #       source_info=None,
    #       lines=loop_ideal.resulting_pdb_h.as_pdb_string()).construct_hierarchy(),
    #     molprobity_scores=True)

    # fixing remaining rotamer outliers
    # fixed_rot_pdb_h = loop_ideal.resulting_pdb_h.deep_copy()
    # fixed_rot_pdb_h.reset_atom_i_seqs()
    if (self.params.additionally_fix_rotamer_outliers and
        self.after_loop_idealization.rotamer_outliers > 0.004):
      print >> self.log, "Processing pdb file again for fixing rotamers..."
      self.log.flush()
      print >> self.log, "Fixing rotamers..."
      self.log.flush()
      if self.params.debug:
        self.shift_and_write_result(
          hierarchy=self.working_pdb_h,
          fname_suffix="just_before_rota")
      self.working_pdb_h = fix_rotamer_outliers(
          pdb_hierarchy=self.working_pdb_h,
          grm=self.working_grm.geometry,
          xrs=self.working_xrs,
          map_data=self.master_map,
          mon_lib_srv=self.mon_lib_srv,
          rotamer_manager=self.rotamer_manager,
          verbose=True)
    if self.params.debug:
      self.shift_and_write_result(
          hierarchy=self.working_pdb_h,
          fname_suffix="rota_ideal",
          grm=self.working_grm)
    cs_to_write = self.cs if self.shift_vector is None else None

    for_stat_h = self.get_intermediate_result_hierarchy()
    self.after_rotamer_fixing = geometry_no_grm(
        pdb_hierarchy=for_stat_h,
        molprobity_scores=True)

    # self.after_rotamer_fixing = geometry_no_grm(
    #     pdb_hierarchy=iotbx.pdb.input(
    #       source_info=None,
    #       lines=fixed_rot_pdb_h.as_pdb_string()).construct_hierarchy(),
    #     molprobity_scores=True)

    ref_hierarchy_for_final_gm = self.original_boxed_hierarchy
    if not self.params.use_starting_model_for_final_gm:
      ref_hierarchy_for_final_gm = self.whole_pdb_h
    ref_hierarchy_for_final_gm.reset_atom_i_seqs()
    # if self.params.additionally_fix_rotamer_outliers:
    #   ssb.set_xyz_smart(self.working_pdb_h, fixed_rot_pdb_h)

    # This massively repeats  self.get_intermediate_result_hierarchy
    self.whole_pdb_h = self.get_intermediate_result_hierarchy()
    # if self.using_ncs:
    #   print >> self.log, "Using ncs"
    #   # multiply back and do geometry_minimization for the whole molecule
    #   for ncs_gr in self.ncs_restr_group_list:
    #     # master_h = self.whole_pdb_h.select(ncs_gr.master_iselection)
    #     new_sites = self.working_pdb_h.atoms().extract_xyz()
    #     self.whole_pdb_h.select(ncs_gr.master_iselection).atoms().set_xyz(new_sites)
    #     for c in ncs_gr.copies:
    #       new_c_sites = c.r.elems * new_sites + c.t
    #       self.whole_pdb_h.select(c.iselection).atoms().set_xyz(new_c_sites)
    #   self.log.flush()
    # else:
    #   # still need to run gm if rotamers were fixed
    #   print >> self.log, "Not using ncs"
    if self.using_ncs:
      print >> self.log, "Using ncs"
    else:
      print >> self.log, "Not using ncs"

    # need to update SS manager for the whole model here.
    if self.params.use_ss_restraints:
      ss_params = sec_str_master_phil.fetch().extract()
      ss_params.secondary_structure.protein.remove_outliers = not self.params.ss_idealization.enabled
      ss_manager = manager(
          pdb_hierarchy=self.whole_pdb_h,
          geometry_restraints_manager=self.whole_grm.geometry,
          sec_str_from_pdb_file=self.filtered_whole_ann,
          params=ss_params.secondary_structure,
          mon_lib_srv=self.mon_lib_srv,
          verbose=-1,
          log=self.log)
      self.whole_grm.geometry.set_secondary_structure_restraints(
          ss_manager=ss_manager,
          hierarchy=self.whole_pdb_h,
          log=self.log)
    print >> self.log, "loop_ideal.ref_exclusion_selection", loop_ideal.ref_exclusion_selection
    print >> self.log, "Minimizing whole model"
    self.minimize(
        hierarchy=self.whole_pdb_h,
        xrs=self.whole_xrs,
        grm=self.whole_grm,
        ncs_restraints_group_list=filtered_ncs_restr_group_list,
        original_pdb_h=ref_hierarchy_for_final_gm,
        excl_string_selection=loop_ideal.ref_exclusion_selection,
        ss_annotation=self.ann,
        reference_map = self.reference_map)
    self.shift_and_write_result(
        hierarchy=self.whole_pdb_h,
        fname_suffix="all_idealized",
        grm=self.whole_grm)
    self.final_model_statistics = geometry_no_grm(
        pdb_hierarchy=iotbx.pdb.input(
          source_info=None,
          lines=self.whole_pdb_h.as_pdb_string()).construct_hierarchy(),
        molprobity_scores=True)
    # self.original_boxed_hierarchy.write_pdb_file(file_name="original_boxed_end.pdb")
    self.time_for_run = time() - t_0

  def minimize(self,
      hierarchy,
      xrs,
      original_pdb_h,
      grm,
      ncs_restraints_group_list,
      excl_string_selection,
      ss_annotation,
      reference_map):
    if reference_map is None:
      minimize_wrapper_for_ramachandran(
          hierarchy=hierarchy,
          xrs=xrs,
          original_pdb_h=original_pdb_h,
          grm=None, # anyway need to reprocess just because of reference model restraints
          excl_string_selection=excl_string_selection,
          number_of_cycles=self.params.number_of_refinement_cycles,
          log=self.log,
          ncs_restraints_group_list=ncs_restraints_group_list,
          ss_annotation=ss_annotation,
          mon_lib_srv=self.mon_lib_srv,
          ener_lib=self.ener_lib,
          rotamer_manager=self.rotamer_manager)
    else:
      print >> self.log, "Using map as reference"
      self.log.flush()
      mwwm = minimize_wrapper_with_map(
          pdb_h=hierarchy,
          xrs=xrs,
          target_map=reference_map,
          grm=grm,
          mon_lib_srv=self.mon_lib_srv,
          rotamer_manager=self.rotamer_manager,
          ncs_restraints_group_list=ncs_restraints_group_list,
          ss_annotation=ss_annotation,
          number_of_cycles=self.params.number_of_refinement_cycles,
          log=self.log)

  def shift_and_write_result(self, hierarchy, fname_suffix, grm=None):
    cs_to_write = self.cs if self.shift_vector is None else None
    pdb_h_shifted = hierarchy.deep_copy()
    pdb_h_shifted.reset_atom_i_seqs()
    if self.shift_vector is not None:
      atoms = pdb_h_shifted.atoms()
      sites_cart = atoms.extract_xyz()
      atoms.set_xyz(new_xyz=sites_cart-self.shift_vector)
    if self.params.debug:
      write_whole_pdb_file(
          file_name="%s_%s_nosh.pdb" % (self.params.output_prefix, fname_suffix),
          pdb_hierarchy=hierarchy,
          crystal_symmetry=self.cs,
          ss_annotation=self.ann)
    write_whole_pdb_file(
        file_name="%s_%s.pdb" % (self.params.output_prefix, fname_suffix),
        pdb_hierarchy=pdb_h_shifted,
        crystal_symmetry=cs_to_write,
        ss_annotation=self.filtered_whole_ann)
    # if grm is not None:
    #   grm.write_geo_file(
    #       sites_cart=hierarchy.atoms().extract_xyz(),
    #       site_labels= [atom.id_str() for atom in hierarchy.atoms()],
    #       file_name="%s_%s.geo" % (self.params.output_prefix, fname_suffix))

  def get_rmsd_from_start(self):
    if self.rmsd_from_start is not None:
      return self.rmsd_from_start
    # calculate rmsd
    self.rmsd_from_start = ssb.calculate_rmsd_smart(
        self.original_boxed_hierarchy,
        self.whole_pdb_h,
        backbone_only=True)
    return self.rmsd_from_start

  def get_rmsd_from_start2(self):
    return ssb.calculate_rmsd_smart(
        self.original_boxed_hierarchy,
        self.whole_pdb_h,
        backbone_only=False)

  def print_stat_comparison(self):
    if self.params.run_minimization_first:
      print >> self.log, "                        Starting    Init GM   SS ideal    Rama      Rota     Final"
      stat_obj_list = [self.init_model_statistics,
          self.init_gm_model_statistics,
          self.after_ss_idealization,
          self.after_loop_idealization,
          self.after_rotamer_fixing,
          self.final_model_statistics,]
    else:
      print >> self.log, "                        Starting    SS ideal    Rama      Rota     Final"
      stat_obj_list = [self.init_model_statistics,
          self.after_ss_idealization,
          self.after_loop_idealization,
          self.after_rotamer_fixing,
          self.final_model_statistics,]
    #                         Starting    SS ideal    Rama      Rota     Final
    # Molprobity Score     :      4.50      3.27      2.66      2.32      2.54
    for val_caption, val_name, val_format in [
        ("Molprobity Score", "mpscore", "{:10.2f}"),
        ("Clashscore", "clashscore", "{:10.2f}"),
        ("CBeta deviations", "c_beta_dev", "{:10d}"),
        ("Ramachandran outliers", "ramachandran_outliers", "{:10.2f}"),
        ("Ramachandran allowed", "ramachandran_allowed", "{:10.2f}"),
        ("Rotamer outliers", "rotamer_outliers", "{:10.2f}"),
        ("Cis-prolines", "n_cis_proline", "{:10d}"),
        ("Cis-general", "n_cis_general", "{:10d}"),
        ("Twisted prolines", "n_twisted_proline", "{:10d}"),
        ("Twisted general", "n_twisted_general", "{:10d}")]:
      l = "%-21s:" % val_caption
      for stat_obj in stat_obj_list:
        value = 99999
        if stat_obj is not None:
          l += val_format.format(getattr(stat_obj, val_name, 99999))
        else:
          l += val_format.format(0)
      print >> self.log, l

  def print_runtime(self):
    print >> self.log, "Time taken for idealization: %s" % str(
        datetime.timedelta(seconds=int(self.time_for_init + self.time_for_run)))

  @staticmethod
  def filter_ncs_restraints_group_list(whole_h, ncs_restr_group_list):
    def whole_chain_in_ncs(whole_h, master_iselection):
      m_c = whole_h.select(master_iselection)
      m_c_id = m_c.only_model().chains()[0].id
      for chain in whole_h.only_model().chains():
        if chain.id == m_c_id:
          if chain.atoms_size() == master_iselection.size():
            return True
          else:
            return False
    n_gr_to_remove = []
    for i, ncs_gr in enumerate(ncs_restr_group_list):
      if not whole_chain_in_ncs(whole_h, ncs_gr.master_iselection):
        n_gr_to_remove.append(i)
    result = deepcopy(ncs_restr_group_list)
    for i in reversed(n_gr_to_remove):
      del result[i]
    return result

def run(args):


  # processing command-line stuff, out of the object
  log = multi_out()
  log.register("stdout", sys.stdout)
  if len(args) == 0:
    format_usage_message(log)
    return
  inputs = mmtbx.utils.process_command_line_args(args=args,
      master_params=master_params())
  work_params = inputs.params.extract()
  inputs.params.show(prefix=" ", out=log)
  pdb_file_names = list(inputs.pdb_file_names)
  if work_params.file_name is not None:
    pdb_file_names += work_params.file_name
  if len(pdb_file_names) == 0:
    raise Sorry("No PDB file specified")
  if work_params.output_prefix is None:
    work_params.output_prefix = os.path.basename(pdb_file_names[0])
  log_file_name = "%s.log" % work_params.output_prefix
  logfile = open(log_file_name, "w")
  log.register("logfile", logfile)
  if work_params.loop_idealization.output_prefix is None:
    work_params.loop_idealization.output_prefix = "%s_rama_fixed" % work_params.output_prefix
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_file_names)
  pdb_input = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  map_data = None
  if inputs.ccp4_map is not None:
    print >> log, "Processing input CCP4 map file..."
    map_data = inputs.ccp4_map.data.as_double()
    print >> log, "Input map min,max,mean: %7.3f %7.3f %7.3f"%\
      map_data.as_1d().min_max_mean().as_tuple()
    map_data = map_data - flex.mean(map_data)
    sd = map_data.sample_standard_deviation()
    map_data = map_data/sd
    print >> log, "Rescaled map min,max,mean: %7.3f %7.3f %7.3f"%\
      map_data.as_1d().min_max_mean().as_tuple()
    inputs.ccp4_map.show_summary(prefix="  ")

  mi_object = model_idealization(
      pdb_input=pdb_input,
      cif_objects=inputs.cif_objects,
      map_data = map_data,
      crystal_symmetry = inputs.crystal_symmetry,
      params=work_params,
      log=log,
      verbose=True)
  mi_object.run()
  mi_object.print_stat_comparison()
  print >> log, "RMSD from starting model (backbone, all): %.4f, %.4f" % (
      mi_object.get_rmsd_from_start(), mi_object.get_rmsd_from_start2())
  mi_object.print_runtime()
  # add hydrogens if needed ?
  print >> log, "All done."
  log.close()

if __name__ == "__main__":
  run(sys.argv[1:])
