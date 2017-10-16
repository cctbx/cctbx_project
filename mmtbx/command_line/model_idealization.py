from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.model_idealization

import sys, os
import datetime
from time import time
from libtbx.utils import Sorry, multi_out, null_out
from libtbx import Auto, easy_pickle, group_args
from scitbx.array_family import flex
from cStringIO import StringIO

from cctbx import crystal
from iotbx.pdb import secondary_structure as ioss
from iotbx.pdb import write_whole_pdb_file
from iotbx.file_reader import any_file
from iotbx import reflection_file_utils
from iotbx.phil import process_command_line_with_files
import iotbx.ncs
import iotbx.phil
from cctbx import maptbx, miller

from mmtbx.secondary_structure import build as ssb
from mmtbx.secondary_structure import manager, sec_str_master_phil
import mmtbx.utils
from mmtbx.utils import fix_rotamer_outliers
from mmtbx.building.loop_idealization import loop_idealization
import mmtbx.building.loop_closure.utils
from mmtbx.refinement.geometry_minimization import minimize_wrapper_for_ramachandran
from mmtbx.ncs.ncs_utils import filter_ncs_restraints_group_list
from mmtbx.command_line.geometry_minimization import get_geometry_restraints_manager
from mmtbx.refinement.real_space.individual_sites import minimize_wrapper_with_map
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from mmtbx.rotamer.rotamer_eval import RotamerEval
from mmtbx_validation_ramachandran_ext import rama_eval
from mmtbx.validation.clashscore import check_and_add_hydrogen
import mmtbx.model

turned_on_ss = ssb.ss_idealization_master_phil_str
turned_on_ss = turned_on_ss.replace("enabled = False", "enabled = True")
master_params_str = """
model_file_name = None
  .type = path
  .multiple = True
  .short_caption = Model file
  .style = file_type:pdb bold input_file
  .expert_level = 0
map_file_name = None
  .type = path
  .help = User-provided map that will be used as reference
  .expert_level = 0
hkl_file_name = None
  .type = path
  .help = User-provided X-ray data to generate 2mFo-DFc map that will be used \
    as reference
  .expert_level = 0
data_labels = None
  .type = str
  .short_caption = Data labels
  .help = Labels for experimental data.
r_free_flags_labels = None
  .type = str
  .short_caption = Rfree labels
  .help = Labels for free reflections.
ligands_file_name = None
  .type = path
  .multiple = True
  .help = User-provided ligand restraints
  .expert_level = 0
mask_and_he_map = False
  .type = bool
  .help = Mask and histogram equalization of the input map
trim_alternative_conformations = False
  .type = bool
  .help = Leave only atoms with empty altloc
  .expert_level = 2
additionally_fix_rotamer_outliers = True
  .type = bool
  .help = At the late stage if rotamer is still outlier choose another one \
    with minimal clash with surrounding atoms
  .expert_level = 2
use_ss_restraints = True
  .type = bool
  .help = Use Secondary Structure restraints
  .expert_level = 2
add_hydrogens = False
  .type = bool
  .help = Add hydrogens to the model before rotamer fixing
  .expert_level = 2
use_starting_model_for_final_gm = False
  .type = bool
  .help = Use supplied model for final geometry minimization. Otherwise just \
    use self.
  .expert_level = 3
output_prefix = None
  .type = str
  .expert_level = 0
output_pkl = False
  .type = bool
  .expert_level = 3
use_map_for_reference = True
  .type = bool
  .expert_level = 1
run_minimization_first = True
  .type = bool
  .expert_level = 2
reference_map_resolution = 5
  .type = float
  .expert_level = 2
number_of_refinement_cycles = Auto
  .type = int
  .expert_level = 1
ignore_ncs = False
  .type = bool
  .help = Don't use NCS even if it is present in model.
  .expert_level = 2
filter_input_ss = True
  .type = bool
  .help = Filter input annotation
  .expert_level = 3
debug = False
  .type = bool
  .help = Output all intermediate files
  .expert_level = 3
verbose = False
  .type = bool
  .help = More output to log
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
               model, # shifted, with shift_manager
               # pdb_h, # shifted model
               # cif_objects=None,
               map_data = None, # shifted map_data
               # shift_manager = None, # manager used to shift and cut map/model
               # crystal_symmetry = None, # original consensus symmetry (map, model)
               # ss_annotation=None,
               params=None,
               log=sys.stdout,
               verbose=True):
    t_0 = time()
    self.model = model
    # self.cif_objects = cif_objects
    self.params = params
    self.log = log
    self.verbose = verbose

    # self.shift_manager = self.model.get_shift_manager()

    self.rmsd_from_start = None
    self.init_model_statistics = None
    self.init_gm_model_statistics = None
    self.after_ss_idealization = None
    self.after_loop_idealization = None
    self.after_rotamer_fixing = None
    self.final_model_statistics = None
    self.user_supplied_map = map_data
    self.reference_map = None # Whole map for all NCS copies
    self.master_map = None # Map for only one NCS copy, or == reference_map if no NCS
    self.init_ref_map = None # separate map for initial GM. Should be tighter than the 2 above

    self.whole_grm = None
    self.master_grm = None
    self.working_grm = None

    # self.mon_lib_srv = None
    # self.ener_lib = None
    # self.rotamer_manager = None
    # self.rama_manager = rama_eval()

    params_line = grand_master_phil_str
    params = iotbx.phil.parse(
        input_string=params_line, process_includes=True).extract()
    params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
    params.pdb_interpretation.peptide_link.ramachandran_restraints = True
    params.pdb_interpretation.peptide_link.oldfield.weight_scale=3
    params.pdb_interpretation.peptide_link.oldfield.plot_cutoff=0.03
    params.pdb_interpretation.peptide_link.apply_peptide_plane = True
    if self.params.loop_idealization.make_all_trans:
      params.pdb_interpretation.peptide_link.apply_all_trans = True
    params.pdb_interpretation.nonbonded_weight = 10000
    params.pdb_interpretation.c_beta_restraints=True
    params.pdb_interpretation.max_reasonable_bond_distance = None
    params.pdb_interpretation.ncs_search.enabled = True
    params.pdb_interpretation.ncs_search.chain_max_rmsd=4.0,
    params.pdb_interpretation.ncs_search.chain_similarity_threshold=0.99,
    params.pdb_interpretation.ncs_search.residue_match_radius=999.0
    params.pdb_interpretation.restraints_library.rdl = True
    self.model.set_pdb_interpretation_params(params)



    self.original_hierarchy = self.model.get_hierarchy() # original pdb_h, without any processing
    self.original_boxed_hierarchy = None # original and boxed (if needed)
    self.whole_pdb_h = None # boxed with processing (AC trimming, H trimming,...)
    self.master_pdb_h = None # master copy in case of NCS
    self.working_pdb_h = None # one to use for fixing (master_pdb_h or working_pdb_h)

    self.using_ncs = False
    self.ncs_restr_group_list = []
    self.filtered_ncs_restr_group_list = []

    self.init_ss_annotation = self.model.get_ss_annotation()

    # various checks, shifts, trims
    self.cs = self.original_cs = self.model.crystal_symmetry()
    if self.model.get_shift_manager() is not None:
      self.cs = self.model.get_shift_manager().box_crystal_symmetry

    # check self.cs (copy-paste from secondary_sturcure_restraints)
    corrupted_cs = False
    if self.cs is not None:
      if [self.cs.unit_cell(), self.cs.space_group()].count(None) > 0:
        corrupted_cs = True
        self.cs = None
      elif self.cs.unit_cell().volume() < 10:
        corrupted_cs = True
        self.cs = None

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

    self.original_boxed_hierarchy = self.model.get_hierarchy().deep_copy()
    # self.original_boxed_hierarchy.reset_atom_i_seqs()
    self.shift_vector = None
    if self.cs is None:
      assert self.model.get_shift_manager() is None
      # should it happen here?
      if corrupted_cs:
        print >> self.log, "Symmetry information is corrupted, "
      else:
        print >> self.log, "Symmetry information was not found, "
      print >> self.log, "putting molecule in P1 box."
      self.log.flush()
      from cctbx import uctbx
      # atoms = self.original_boxed_hierarchy.atoms()
      box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
        sites_cart=self.model.get_sites_cart(),
        buffer_layer=3)
      self.model.set_sites_cart(sites_cart=box.sites_cart)
      self.original_boxed_hierarchy = self.model.get_hierarchy().deep_copy()
      # atoms.set_xyz(new_xyz=box.sites_cart)
      self.cs = box.crystal_symmetry()
      self.model.set_crystal_symmetry(self.cs)
      self.shift_vector = box.shift_vector

    # self.original_boxed_hierarchy.write_pdb_file(file_name="original_boxed_h.pdb")
    if self.shift_vector is not None and self.params.debug:
      write_whole_pdb_file(
          file_name="%s_boxed.pdb" % self.params.output_prefix,
          pdb_hierarchy=self.model.get_hierarchy(),
          crystal_symmetry=self.model.crystal_symmetry(),
          ss_annotation=self.init_ss_annotation)

    # asc = self.original_boxed_hierarchy.atom_selection_cache()
    if self.params.trim_alternative_conformations:
      self.model.remove_alternative_conformations(always_keep_one_conformer=True)
      # sel = asc.selection("altloc ' '")
      # self.whole_pdb_h = self.original_boxed_hierarchy.select(sel).deep_copy()
      # print >> self.log, "Atoms in original/working model: %d/%d" % (
      #     self.original_boxed_hierarchy.atoms_size(), self.whole_pdb_h.atoms_size())
      # else:
      # self.whole_pdb_h = self.original_boxed_hierarchy.deep_copy()
      # self.whole_pdb_h.reset_atom_i_seqs()
    # self.whole_pdb_h = self.model.get_hierarchy()
    # Trimming hydrogens
    # Many intermediate variables are needed due to strange behavior of
    # selections described in
    # iotbx/pdb/tst_hierarchy.py:exercise_selection_and_deep_copy()
    # asc2 = self.whole_pdb_h.atom_selection_cache()
    # h_sel = asc2.selection("not (element H or element D)")
    # temp_h = self.whole_pdb_h.select(h_sel)
    # self.whole_pdb_h = temp_h.deep_copy()
    # self.whole_pdb_h.reset_atom_i_seqs()

    self.model = self.model.remove_hydrogens()
    self.whole_pdb_h = self.model.get_hierarchy()

    # for_stat_h = self.get_intermediate_result_hierarchy()

    # self.init_model_statistics = mmtbx.model.statistics(pdb_hierarchy=for_stat_h)
    # self.init_model_statistics = self.model.geometry_statistics()
    # print dir(self.init_model_statistics)
    self.time_for_init = time()-t_0

  def get_statistics(self, model):
    # should we shift here?
    # should we multiply NCS here?
    # print "model.restraints_manager", model.restraints_manager
    return model.geometry_statistics().result_for_model_idealization()

  def get_intermediate_result_hierarchy(self):
    result_h = self.model.get_hierarchy()
    result_h.reset_atom_i_seqs()
    if self.using_ncs:
      # multiply back and do geometry_minimization for the whole molecule
      for ncs_gr in self.filtered_ncs_restr_group_list:
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

  def prepare_user_map(self):
    print >> self.log, "Preparing user map..."
    self.map_shift_manager = mmtbx.utils.shift_origin(
      map_data         = self.user_supplied_map,
      xray_structure   = self.whole_pdb_h.extract_xray_structure(crystal_symmetry=self.cs),
      crystal_symmetry = self.cs)
    if(self.map_shift_manager.shift_cart is not None):
      # Need to figure out way to save the shift to shift back
      # and apply it to whole_pdb, master_pdb, etc. Don't forget about
      # boxing hierarchy when symmetry is not available or corrupted...
      raise Sorry("Map origin is not at (0,0,0). This is not implemented for model_idealization")
    map_data = self.map_shift_manager.map_data
    self.reference_map = map_data
    self.master_map = self.reference_map.deep_copy()
    if self.using_ncs and self.master_pdb_h is not None:
      # here we are negating non-master part of the model
      # self.master_sel=master_sel
      # self.master_map = self.reference_map.deep_copy()
      mask = maptbx.mask(
              xray_structure=self.whole_xrs.select(self.master_sel),
              n_real=self.master_map.focus(),
              mask_value_inside_molecule=1,
              mask_value_outside_molecule=-1,
              solvent_radius=0,
              atom_radius=1.)
      self.master_map = self.reference_map * mask
      if self.params.debug:
        iotbx.ccp4_map.write_ccp4_map(
            file_name="%s_3_master.map" % self.params.output_prefix,
            unit_cell=self.cs.unit_cell(),
            space_group=self.cs.space_group(),
            map_data=self.master_map,
            labels=flex.std_string([""]))
      self.master_map = map_data

  def prepare_init_reference_map(self, xrs, pdb_h):
    if self.user_supplied_map is not None:
      print >> self.log, "Using user-supplied map for initial GM."
      self.init_ref_map = self.reference_map
      return
    print >> self.log, "Preparing map for initial GM..."
    asc = pdb_h.atom_selection_cache()
    outlier_selection_txt = mmtbx.building.loop_closure.utils. \
          rama_score_selection(pdb_h, self.model.get_ramachandran_manager(), "outlier",1)
    # print >> self.log, "rama outlier selection:", outlier_selection_txt
    rama_out_sel = asc.selection(outlier_selection_txt)

    allowed_selection_txt = mmtbx.building.loop_closure.utils. \
          rama_score_selection(pdb_h, self.model.get_ramachandran_manager(), "allowed",0)
    # print >> self.log, "rama allowed selection:", allowed_selection_txt
    rama_allowed_sel = asc.selection(allowed_selection_txt)


    # side_chain_no_cb_selection = ~ xrs.main_chain_selection()
    side_chain_no_cb_selection = ~ xrs.backbone_selection()
    sc_rama_out = rama_out_sel & side_chain_no_cb_selection
    sc_rama_allowed =rama_allowed_sel & side_chain_no_cb_selection
    xrs=xrs.set_b_iso(value=10)
    xrs = xrs.set_b_iso(value=20, selection=side_chain_no_cb_selection)
    xrs = xrs.set_b_iso(value=25, selection=rama_allowed_sel)
    xrs = xrs.set_b_iso(value=50, selection=rama_out_sel)
    xrs = xrs.set_b_iso(value=40, selection=sc_rama_allowed)
    xrs = xrs.set_b_iso(value=70, selection=rama_out_sel)

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
    self.init_ref_map = init_reference_map

  def prepare_reference_map(self, xrs, pdb_h):
    print >> self.log, "Preparing reference map"
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
    xrs=xrs.set_b_iso(value=50)
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
          rama_score_selection(pdb_h, self.model.get_ramachandran_manager(), "outlier",1)
    asc = pdb_h.atom_selection_cache()
    # print >> self.log, "rama outlier selection:", outlier_selection_txt
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

  def update_ss_in_grm(self, ss_annotation):
    self.set_ss_restraints(ss_annotation)
    self.update_grms()

  def set_ss_restraints(self, ss_annotation):
    log = self.log
    if not self.verbose:
      log = null_out()
    if self.params.use_ss_restraints and ss_annotation is not None:
      ss_manager = manager(
          pdb_hierarchy=self.whole_pdb_h,
          geometry_restraints_manager=self.whole_grm.geometry,
          sec_str_from_pdb_file=ss_annotation,
          params=None,
          mon_lib_srv=self.model.get_mon_lib_srv(),
          verbose=-1,
          log=log)
      # self.whole_pdb_h.write_pdb_file(file_name="for_ss.pdb")
      self.whole_pdb_h.reset_atom_i_seqs()
      self.whole_grm.geometry.set_secondary_structure_restraints(
          ss_manager=ss_manager,
          hierarchy=self.whole_pdb_h,
          log=log)

  def update_grms(self):
    if self.using_ncs:
      self.master_grm = self.whole_grm.select(self.master_sel)
      self.working_grm = self.master_grm
    else:
      self.working_grm = self.whole_grm

  def get_grm(self):
    # first make whole grm using self.whole_pdb_h

    self.model.get_restraints_manager()
    # log = self.log
    # if not self.verbose:
    #   log = null_out()
    # if self.params.ignore_ncs:
    #   params.pdb_interpretation.ncs_search.enabled = False
    # processed_pdb_files_srv = mmtbx.utils.\
    #     process_pdb_file_srv(
    #         crystal_symmetry= self.whole_xrs.crystal_symmetry(),
    #         pdb_interpretation_params = params.pdb_interpretation,
    #         stop_for_unknowns         = False,
    #         log=log,
    #         cif_objects=self.cif_objects)
    # processed_pdb_file, junk = processed_pdb_files_srv.\
    #     process_pdb_files(
    #         raw_records=flex.split_lines(self.whole_pdb_h.as_pdb_string()))

    # self.model.get_mon_lib_srv() = processed_pdb_files_srv.mon_lib_srv
    # self.ener_lib = processed_pdb_files_srv.ener_lib
    # self.rotamer_manager = RotamerEval(mon_lib_srv=self.mon_lib_srv)

    # self.whole_grm = get_geometry_restraints_manager(
    #     processed_pdb_file, self.whole_xrs, params=params)

    self.whole_grm = self.model.get_restraints_manager()

    # set SS restratins
    self.set_ss_restraints(self.filtered_whole_ann)

    # now select part of it for working with master hierarchy
    self.update_grms()


  def get_filtered_ncs_group_list(self):
    if not self.params.ignore_ncs:
      # ncs_obj = self.model.get_ncs_obj()
      ncs_obj = iotbx.ncs.input(
          hierarchy=self.model.get_hierarchy(),
          chain_max_rmsd=4.0,
          chain_similarity_threshold=0.99,
          residue_match_radius=999.0)
      if ncs_obj is not None:
        print >> self.log, "Found NCS groups:"
        ncs_obj.show(format='phil', log=self.log)
        self.ncs_restr_group_list = ncs_obj.get_ncs_restraints_group_list(
            raise_sorry=False)
        total_ncs_selected_atoms = 0
        master_sel = flex.size_t([])
        self.filtered_ncs_restr_group_list = filter_ncs_restraints_group_list(
            self.whole_pdb_h, self.ncs_restr_group_list)
        if len(self.filtered_ncs_restr_group_list) > 0:
          self.using_ncs = True
          master_sel = flex.bool(self.whole_pdb_h.atoms_size(), True)
          for ncs_gr in self.filtered_ncs_restr_group_list:
            for copy in ncs_gr.copies:
              master_sel.set_selected(copy.iselection, False)
          self.master_pdb_h = self.whole_pdb_h.select(master_sel).deep_copy()
          self.master_sel=master_sel
          self.master_pdb_h.reset_atom_i_seqs()
          self.master_model = self.model.select(master_sel)


    if self.using_ncs:
      if self.params.debug:
        self.master_pdb_h.write_pdb_file("%s_master_h.pdb" % self.params.output_prefix)
      self.working_pdb_h = self.master_pdb_h
      self.working_model = self.master_model
    else:
      self.working_pdb_h = self.whole_pdb_h
      self.working_model = self.model.deep_copy()
    self.working_pdb_h.reset_atom_i_seqs()


    self.working_xrs = self.working_pdb_h.extract_xray_structure(crystal_symmetry=self.cs)
    if self.using_ncs:
      self.whole_xrs = self.whole_pdb_h.extract_xray_structure(crystal_symmetry=self.cs)
    else:
      self.whole_xrs = self.working_xrs


  def run(self):
    t_0 = time()
    self.ann = ioss.annotation.from_phil(
        phil_helices=self.params.secondary_structure.protein.helix,
        phil_sheets=self.params.secondary_structure.protein.sheet,
        pdb_hierarchy=self.whole_pdb_h)
    if self.ann.get_n_helices() + self.ann.get_n_sheets() == 0:
      self.ann = self.init_ss_annotation
      if self.ann is not None:
        self.ann.remove_empty_annotations(self.whole_pdb_h)
    self.filtered_whole_ann = None
    if self.ann is not None:
      self.filtered_whole_ann = self.ann.deep_copy()

    self.get_grm()
    self.get_filtered_ncs_group_list()

    self.init_model_statistics = self.get_statistics(self.model)

    # Here we are preparing maps if needed.
    if self.user_supplied_map is not None:
      self.prepare_user_map()

    if self.reference_map is None and self.params.use_map_for_reference:
      # self.prepare_reference_map(xrs=self.whole_xrs, pdb_h=self.whole_pdb_h)
      # self.prepare_reference_map_2(xrs=self.whole_xrs, pdb_h=self.whole_pdb_h)
      self.prepare_reference_map_3(xrs=self.whole_xrs, pdb_h=self.whole_pdb_h)

    if self.params.run_minimization_first:
      # running simple minimization and updating all
      # self.master, self.working, etc...
      self.whole_pdb_h.reset_atom_i_seqs()
      if self.init_ref_map is None:
        self.prepare_init_reference_map(
            self.whole_pdb_h.extract_xray_structure(crystal_symmetry=self.cs),
            self.whole_pdb_h)
      print >> self.log, "Minimization first"
      self.minimize(
          model=self.model,
          # xrs=self.whole_pdb_h.extract_xray_structure(crystal_symmetry=self.cs),
          original_pdb_h=self.whole_pdb_h,
          # grm=self.whole_grm,
          ncs_restraints_group_list=self.filtered_ncs_restr_group_list,
          excl_string_selection=None, # don't need if we have map
          # ss_annotation=self.ann,
          reference_map=self.init_ref_map,
          # reference_map=self.reference_map,
          )
      self.init_gm_model_statistics = self.get_statistics(self.model)
      if self.params.debug:
        self.shift_and_write_result(
            model = self.model,
            fname_suffix="init_gm")

    if (self.init_gm_model_statistics is not None
        and self.init_gm_model_statistics.ramachandran_outliers == 0
        and self.init_gm_model_statistics.twisted_general <= 0.01
        and self.init_gm_model_statistics.twisted_proline <= 0.01
        and self.init_gm_model_statistics.cis_general <= 0.01
        and self.init_gm_model_statistics.cis_proline <= 0.01):
      print >> self.log, "Simple minimization was enough"
      # Early exit!!!
      self.shift_and_write_result(
          model=self.model,
          fname_suffix="all_idealized")
      self.final_model_statistics = self.get_statistics(self.model)
      # self.original_boxed_hierarchy.write_pdb_file(file_name="original_boxed_end.pdb")
      self.time_for_run = time() - t_0
      if self.params.output_pkl:
        easy_pickle.dump(
            file_name="%s.pkl" % self.params.output_prefix,
            obj = self.get_stats_obj())
      return

    self.original_ann = None
    if self.ann is not None:
      self.original_ann = self.ann.deep_copy()
      print >> self.log, "Original SS annotation"
      print >> self.log, self.original_ann.as_pdb_str()
      if self.params.filter_input_ss:
        self.filtered_whole_ann = self.ann.filter_annotation(
            hierarchy=self.whole_pdb_h,
            asc=None)
        self.ann = self.ann.filter_annotation(
            hierarchy=self.working_pdb_h,
            asc=None)
      else:
        self.filtered_whole_ann = self.original_ann.deep_copy()
      # print >> self.log, "Splitted SS annotation"
      # print >> self.log, ann.as_pdb_str()
      print >> self.log, "Filtered SS annotation"
      print >> self.log, self.filtered_whole_ann.as_pdb_str()
      self.model.set_ss_annotation(self.filtered_whole_ann)

    # getting grm with SS restraints
    self.update_ss_in_grm(self.filtered_whole_ann)

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
          rama_score_selection(self.working_pdb_h, self.model.get_ramachandran_manager(), "outlier",1)
        print >> self.log, "outlier_selection_txt", outlier_selection_txt
        negate_selection = "all"
        if outlier_selection_txt != "" and outlier_selection_txt is not None:
          negate_selection = "not (%s)" % outlier_selection_txt
      # if self.params.run_minimization_first:
      self.minimize(
          model=self.model,
          # xrs=self.whole_xrs,
          original_pdb_h=self.whole_pdb_h,
          # grm=self.whole_grm,
          ncs_restraints_group_list=self.filtered_ncs_restr_group_list,
          excl_string_selection=negate_selection,
          # ss_annotation=self.ann,
          reference_map=self.reference_map)
      # self.original_boxed_hierarchy.write_pdb_file(file_name="original_boxed_h_1.pdb")
    else:
      if self.params.debug:
        self.params.ss_idealization.file_name_before_regularization = \
            "%s_ss_before_reg.pdb" % self.params.output_prefix
      self.params.ss_idealization.skip_good_ss_elements = True
      ssb.substitute_ss(
          # real_h=self.working_pdb_h,
          # xray_structure=self.working_xrs,
          # ss_annotation=self.ann,
          model = self.working_model,
          params=self.params.ss_idealization,
          # grm=self.workiworking_grm,
          fix_rotamer_outliers=True,
          # cif_objects=self.cif_objects,
          verbose=self.params.verbose,
          reference_map=self.master_map,
          # rotamer_manager=self.model.get_rotamer_manager(),
          log=self.log)
      self.log.flush()

    # for_stat_h = self.get_intermediate_result_hierarchy()
    self.after_ss_idealization = self.get_statistics(self.model)
    self.shift_and_write_result(
          model=self.working_model,
          fname_suffix="ss_ideal_stat")
    # self.after_ss_idealization = geometry_no_grm(
    #     pdb_hierarchy=iotbx.pdb.input(
    #       source_info=None,
    #       lines=self.working_pdb_h.as_pdb_string()).construct_hierarchy(),
    #     molprobity_scores=True)

    # Write resulting pdb file.
    if self.params.debug:
      self.shift_and_write_result(
          model=self.working_model,
          # hierarchy=self.working_pdb_h,
          fname_suffix="ss_ideal",
          # grm=self.working_grm,
          )
    self.params.loop_idealization.minimize_whole = not self.using_ncs
    self.params.loop_idealization.debug = self.params.debug or self.params.loop_idealization.debug
    # self.params.loop_idealization.enabled = False
    # self.params.loop_idealization.variant_search_level = 0
    print >> self.log, "Starting loop idealization"
    loop_ideal = loop_idealization(
        pdb_hierarchy=self.working_pdb_h,
        params=self.params.loop_idealization,
        secondary_structure_annotation=self.ann,
        reference_map=self.master_map,
        crystal_symmetry=self.working_xrs.crystal_symmetry(),
        grm=self.working_grm,
        rama_manager=self.model.get_ramachandran_manager(),
        rotamer_manager=self.model.get_rotamer_manager(),
        log=self.log,
        verbose=True)
    self.log.flush()
    self.working_pdb_h = loop_ideal.resulting_pdb_h
    if self.params.debug:
      self.shift_and_write_result(
          model = self.working_model,
          fname_suffix="rama_ideal")
    # for_stat_h = self.get_intermediate_result_hierarchy()
    # for_stat_h.write_pdb_file(file_name="compare_with_rama_ideal.pdb")
    self.after_loop_idealization = self.get_statistics(self.working_model)

    if self.params.add_hydrogens:
      print >> self.log, "Adding hydrogens"
      self.add_hydrogens()

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
          model = self.working_model,
          fname_suffix="just_before_rota")
      self.working_pdb_h = fix_rotamer_outliers(
          model = working_model,
          map_data=self.master_map,
          verbose=True)
    if self.params.debug:
      self.shift_and_write_result(
          model = self.working_model,
          fname_suffix="rota_ideal")
    cs_to_write = self.cs if self.shift_vector is None else None

    # for_stat_h = self.get_intermediate_result_hierarchy()
    self.after_rotamer_fixing = self.get_statistics(self.working_model)
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
    # assert 0
    if self.using_ncs:
      print >> self.log, "Using ncs"
      # assert 0
    else:
      print >> self.log, "Not using ncs"
      # assert 0

    # need to update SS manager for the whole model here.
    # print "self.filtered_ncs_restr_group_list", self.filtered_ncs_restr_group_list
    # STOP()
    if self.params.use_ss_restraints:
      ss_params = sec_str_master_phil.fetch().extract()
      ss_params.secondary_structure.protein.remove_outliers = not self.params.ss_idealization.enabled
      ss_manager = manager(
          pdb_hierarchy=self.whole_pdb_h,
          geometry_restraints_manager=self.whole_grm.geometry,
          sec_str_from_pdb_file=self.filtered_whole_ann,
          params=ss_params.secondary_structure,
          mon_lib_srv=self.model.get_mon_lib_srv(),
          verbose=-1,
          log=self.log)
      self.whole_grm.geometry.set_secondary_structure_restraints(
          ss_manager=ss_manager,
          hierarchy=self.whole_pdb_h,
          log=self.log)
    print >> self.log, "loop_ideal.ref_exclusion_selection", loop_ideal.ref_exclusion_selection
    print >> self.log, "Minimizing whole model"
    self.minimize(
        model = self.model,
        # hierarchy=self.whole_pdb_h,
        # xrs=self.whole_xrs,
        # grm=self.whole_grm,
        ncs_restraints_group_list=self.filtered_ncs_restr_group_list,
        original_pdb_h=ref_hierarchy_for_final_gm,
        excl_string_selection=loop_ideal.ref_exclusion_selection,
        # ss_annotation=self.ann,
        reference_map = self.reference_map)
    self.shift_and_write_result(
        model = self.model,
        fname_suffix="all_idealized")
    self.final_model_statistics = self.get_statistics(self.model)
    # mmtbx.model.statistics(
    #     pdb_hierarchy=iotbx.pdb.input(
    #       source_info=None,
    #       lines=self.whole_pdb_h.as_pdb_string()).construct_hierarchy())
    # self.original_boxed_hierarchy.write_pdb_file(file_name="original_boxed_end.pdb")
    self.time_for_run = time() - t_0
    if self.params.output_pkl or self.params.debug:
      easy_pickle.dump(
          file_name="%s.pkl" % self.params.output_prefix,
          obj = self.get_stats_obj())

  def add_hydrogens(self):
    # self.wo
    self.whole_pdb_h = self.get_intermediate_result_hierarchy()
    cs = self.whole_xrs.crystal_symmetry()
    pdb_str, changed = check_and_add_hydrogen(
        pdb_hierarchy=self.whole_pdb_h,
        file_name=None,
        nuclear=False,
        keep_hydrogens=False,
        verbose=False,
        model_number=0,
        n_hydrogen_cut_off=0,
        time_limit=120,
        allow_multiple_models=True,
        crystal_symmetry=self.cs,
        do_flips=True,
        log=self.log)
    self.whole_pdb_h = iotbx.pdb.input(lines=pdb_str, source_info=None).construct_hierarchy()
    # we need to renew everythig here: grm, whole model, xrs!
    self.whole_xrs = self.working_pdb_h.extract_xray_structure(crystal_symmetry=cs)
    self.get_grm()
    self.get_filtered_ncs_group_list()
    # temp workaround, do something here to handle NCS cases
    self.working_pdb_h = self.whole_pdb_h
    self.working_xrs = self.working_pdb_h.extract_xray_structure(crystal_symmetry=cs)
    self.working_grm = self.whole_grm



  def minimize(self,
      model,
      # hierarchy,
      # xrs,
      original_pdb_h,
      # grm,
      ncs_restraints_group_list,
      excl_string_selection,
      # ss_annotation,
      reference_map):
    # print "ncs_restraints_group_list", ncs_restraints_group_list
    # assert 0
    if reference_map is None:
      minimize_wrapper_for_ramachandran(
          hierarchy=model.get_hierarchy(),
          xrs=model.get_xray_structure(),
          original_pdb_h=original_pdb_h,
          grm=model.get_restraints_manager(), # anyway need to reprocess just because of reference model restraints
          processed_pdb_file = model.processed_pdb_file,
          excl_string_selection=excl_string_selection,
          number_of_cycles=self.params.number_of_refinement_cycles,
          log=self.log,
          ncs_restraints_group_list=ncs_restraints_group_list,
          ss_annotation=model.get_ss_annotation(),
          mon_lib_srv=model.get_mon_lib_srv(),
          ener_lib=model.get_ener_lib(),
          rotamer_manager=model.get_rotamer_manager())
    else:
      print >> self.log, "Using map as reference"
      self.log.flush()
      mwwm = minimize_wrapper_with_map(
          pdb_h=model.get_hierarchy(),
          xrs=model.get_xray_structure(),
          target_map=reference_map,
          grm=model.get_restraints_manager(),
          mon_lib_srv=model.get_mon_lib_srv(),
          rotamer_manager=model.get_rotamer_manager(),
          ncs_restraints_group_list=ncs_restraints_group_list,
          ss_annotation=model.get_ss_annotation(),
          number_of_cycles=self.params.number_of_refinement_cycles,
          log=self.log)

  def shift_and_write_result(self, model, fname_suffix=""):
    pdb_str = model.model_as_pdb()
    fname = "%s_%s.pdb" % (self.params.output_prefix, fname_suffix)
    with open(fname, 'w') as f:
      f.write(pdb_str)
    if self.params.debug:
      fname = "%s_%s_nosh.pdb" % (self.params.output_prefix, fname_suffix)
      pdb_str = model.model_as_pdb(do_not_shift_back=True)
      with open(fname, 'w') as f:
        f.write(pdb_str)
    # cs_to_write = self.original_cs
    # pdb_h_shifted = hierarchy.deep_copy()
    # if self.shift_manager is not None:
    #   self.shift_manager.shift_back(pdb_hierarchy=pdb_h_shifted)
    # pdb_h_shifted.reset_atom_i_seqs()
    # if self.params.debug:
    #   write_whole_pdb_file(
    #       file_name="%s_%s_nosh.pdb" % (self.params.output_prefix, fname_suffix),
    #       pdb_hierarchy=hierarchy,
    #       crystal_symmetry=self.cs,
    #       ss_annotation=self.ann)
    # write_whole_pdb_file(
    #     file_name="%s_%s.pdb" % (self.params.output_prefix, fname_suffix),
    #     pdb_hierarchy=pdb_h_shifted,
    #     crystal_symmetry=cs_to_write,
    #     ss_annotation=self.filtered_whole_ann)
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
        self.model.get_hierarchy(),
        backbone_only=True)
    return self.rmsd_from_start

  def get_rmsd_from_start2(self):
    return ssb.calculate_rmsd_smart(
        self.original_boxed_hierarchy,
        self.whole_pdb_h,
        backbone_only=False)

  def get_stats_obj(self):
    if self.params.run_minimization_first:
      stat_obj_list = [self.init_model_statistics,
          self.init_gm_model_statistics,
          self.after_ss_idealization,
          self.after_loop_idealization,
          self.after_rotamer_fixing,
          self.final_model_statistics,]
    else:
      stat_obj_list = [self.init_model_statistics,
          self.after_ss_idealization,
          self.after_loop_idealization,
          self.after_rotamer_fixing,
          self.final_model_statistics,]
    return group_args(
        geoms=stat_obj_list,
        rmsds=(self.get_rmsd_from_start(), self.get_rmsd_from_start2()),
        runtime=self.time_for_init + self.time_for_run)


  def print_stat_comparison(self):
    stat_obj_list = self.get_stats_obj()
    if self.params.run_minimization_first:
      print >> self.log, "                        Starting    Init GM   SS ideal    Rama      Rota     Final"
    else:
      print >> self.log, "                        Starting    SS ideal    Rama      Rota     Final"
    #                         Starting    SS ideal    Rama      Rota     Final
    # Molprobity Score     :      4.50      3.27      2.66      2.32      2.54
    for val_caption, val_name, val_format in [
        ("Molprobity Score", "mpscore", "{:10.2f}"),
        ("Clashscore", "clashscore", "{:10.2f}"),
        ("CBeta deviations", "c_beta_dev_percent", "{:10.2f}"),
        ("Ramachandran outliers", "ramachandran_outliers", "{:10.2f}"),
        ("Ramachandran allowed", "ramachandran_allowed", "{:10.2f}"),
        ("Rotamer outliers", "rotamer_outliers", "{:10.2f}"),
        ("Cis-prolines", "cis_proline", "{:10.2f}"),
        ("Cis-general", "cis_general", "{:10.2f}"),
        ("Twisted prolines", "twisted_proline", "{:10.2f}"),
        ("Twisted general", "twisted_general", "{:10.2f}"),
        # Until enabled in model.statistics
        # ("CaBLAM outliers", "cablam_outliers", "{:10.2f}"),
        # ("CaBLAM disfavored", "cablam_disfavored", "{:10.2f}"),
        # ("CaBLAM CA outliers", "cablam_ca_outliers", "{:10.2f}"),
        ]:
      l = "%-21s:" % val_caption
      for stat_obj in stat_obj_list.geoms:
        value = 99999
        if stat_obj is not None:
          l += val_format.format(getattr(stat_obj, val_name, 99999))
        else:
          l += val_format.format(0)
      print >> self.log, l

  def print_runtime(self):
    print >> self.log, "Time taken for idealization: %s" % str(
        datetime.timedelta(seconds=int(self.time_for_init + self.time_for_run)))

def get_map_from_hkl(hkl_file_object, params, xrs, log):
  print >> log, "Processing input hkl file..."
  crystal_symmetry = hkl_file_object.crystal_symmetry()
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry   = True,
    reflection_files = [hkl_file_object.file_content],
    err              = StringIO())


  parameters = mmtbx.utils.data_and_flags_master_params().extract()
  if (params.data_labels is not None):
    parameters.labels = params.data_labels
  if (params.r_free_flags_labels is not None):
    parameters.r_free_flags.label = params.r_free_flags_labels
  determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
    reflection_file_server = rfs,
    parameters             = parameters,
    keep_going             = True,
    working_point_group = crystal_symmetry.space_group().build_derived_point_group(),
    log                    = StringIO(),
    symmetry_safety_check  = True)
  f_obs = determined_data_and_flags.f_obs

  if (params.data_labels is None):
    params.data_labels = f_obs.info().label_string()
  r_free_flags = determined_data_and_flags.r_free_flags
  assert f_obs is not None
  print >> log,  "Input data:"
  print >> log, "  Iobs or Fobs:", f_obs.info().labels
  if (r_free_flags is not None):
    print >> log, "  Free-R flags:", r_free_flags.info().labels
    params.r_free_flags_labels = r_free_flags.info().label_string()
  else:
    print >> log, "  Free-R flags: Not present"

  fmodel = mmtbx.f_model.manager(
      f_obs        = f_obs,
      r_free_flags = r_free_flags,
      xray_structure = xrs)
  fmodel.update_all_scales()

  fft_map = fmodel.electron_density_map().fft_map(
    resolution_factor = 0.25,
    map_type          = "2mFo-DFc",
    use_all_data      = False) # Exclude free reflections
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded(in_place=False)
  if params.debug:
    fft_map.as_xplor_map(file_name="%s_21.map" % params.output_prefix)
    iotbx.ccp4_map.write_ccp4_map(
        file_name="%s_21.ccp4" % params.output_prefix,
        unit_cell=crystal_symmetry.unit_cell(),
        space_group=crystal_symmetry.space_group(),
        map_data=map_data,
        labels=flex.std_string([""]))
  return map_data, crystal_symmetry

def get_map_from_map(map_file_object, params, xrs, log):
  print >> log, "Processing input CCP4 map file..."
  map_data = map_file_object.file_content.data.as_double()
  try:
    # map_cs = map_content.file_object.crystal_symmetry()
    map_cs = map_file_object.crystal_symmetry()
  except NotImplementedError as e:
    pass
  print >> log, "Input map min,max,mean: %7.3f %7.3f %7.3f"%\
      map_data.as_1d().min_max_mean().as_tuple()
  if map_cs.space_group().type().number() not in [0,1]:
    print map_cs.space_group().type().number()
    raise Sorry("Only P1 group for maps is supported.")
  map_data = map_data - flex.mean(map_data)
  sd = map_data.sample_standard_deviation()
  map_data = map_data/sd
  print >> log, "Rescaled map min,max,mean: %7.3f %7.3f %7.3f"%\
    map_data.as_1d().min_max_mean().as_tuple()
  map_file_object.file_content.show_summary(prefix="  ")
  shift_manager = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure = xrs,
      map_data       = map_data.deep_copy(),
      box_cushion    = 5)
  sys.stdout.flush()
  xray_structure = shift_manager.xray_structure_box
  crystal_symmetry = xray_structure.crystal_symmetry()
  map_data = shift_manager.map_box

  if params.mask_and_he_map:
    print >> log, "Masking and histogram equalizing..."
    import boost.python
    cctbx_maptbx_ext = boost.python.import_ext("cctbx_maptbx_ext")
    xrs_p1 = xray_structure.expand_to_p1(sites_mod_positive=True)
    radii = flex.double(xrs_p1.scatterers().size(), 5.0)
    mask = cctbx_maptbx_ext.mask(
      sites_frac                  = xrs_p1.sites_frac(),
      unit_cell                   = xrs_p1.unit_cell(),
      n_real                      = map_data.all(),
      mask_value_inside_molecule  = 1,
      mask_value_outside_molecule = 0,
      radii                       = radii)
    map_data = mask*map_data
    from phenix.command_line.real_space_refine import write_ccp4_map
    write_ccp4_map(o=xray_structure.crystal_symmetry(), file_name="junk_mask.map",
     map_data=mask)
    del mask
    map_data = maptbx.volume_scale(map = map_data, n_bins = 10000).map_data()
    write_ccp4_map(o=xray_structure.crystal_symmetry(), file_name="junk_map.map",
     map_data=map_data)
  return map_data, map_cs, shift_manager

def run(args):
  # processing command-line stuff, out of the object
  log = multi_out()
  log.register("stdout", sys.stdout)
  if len(args) == 0:
    format_usage_message(log)
    return
  input_objects = process_command_line_with_files(
      args=args,
      master_phil=master_params(),
      pdb_file_def="model_file_name",
      map_file_def="map_file_name",
      reflection_file_def="hkl_file_name",
      cif_file_def="ligands_file_name")
  work_params = input_objects.work.extract()
  if [work_params.map_file_name, work_params.hkl_file_name].count(None) < 1:
    raise Sorry("Only one source of map could be supplied.")
  input_objects.work.show(prefix=" ", out=log)
  if len(work_params.model_file_name) == 0:
    raise Sorry("No PDB file specified")
  if work_params.output_prefix is None:
    work_params.output_prefix = os.path.basename(work_params.model_file_name[0])
  log_file_name = "%s.log" % work_params.output_prefix
  logfile = open(log_file_name, "w")
  log.register("logfile", logfile)
  err_log = multi_out()
  err_log.register(label="log", file_object=log)
  # err_log.register(label="stderr", file_object=sys.stderr)
  sys.stderr = err_log

  if work_params.loop_idealization.output_prefix is None:
    work_params.loop_idealization.output_prefix = "%s_rama_fixed" % work_params.output_prefix

  # Here we start opening files provided,
  # collect crystal symmetries
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=work_params.model_file_name)
  pdb_input = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))

  model = mmtbx.model.manager(
      model_input = pdb_input,
      restraint_objects = input_objects.cif_objects,
      process_input = False,
      log=log)

  pdb_cs = model.crystal_symmetry()
  pdb_h = model.get_hierarchy()
  map_cs = None
  crystal_symmetry = None
  map_data = None
  shift_manager = None
  map_content = input_objects.get_file(work_params.map_file_name)

  if map_content is not None:
    map_data, map_cs, shift_manager = get_map_from_map(
        map_content,
        work_params,
        xrs=model.get_xray_structure(),
        log=log)
    model.set_shift_manager(shift_manager)
    model.get_hierarchy().write_pdb_file("junk_shift.pdb")

  hkl_content = input_objects.get_file(work_params.hkl_file_name)
  if hkl_content is not None:
    map_data, map_cs = get_map_from_hkl(
        hkl_content,
        work_params,
        xrs=model.get_xray_structure(), # here we don't care about atom order
        log=log)

  # Crystal symmetry: validate and finalize consensus object
  if shift_manager is not None:
    crystal_symmetry = map_cs
  else:
    try:
      crystal_symmetry = crystal.select_crystal_symmetry(
          from_command_line     = None,
          from_parameter_file   = None,
          from_coordinate_files = [pdb_cs],
          from_reflection_files = [map_cs],
          enforce_similarity    = True)
    except AssertionError as e:
      if len(e.args)>0 and e.args[0].startswith("No unit cell and symmetry information supplied"):
        pass
      else:
        raise e
  # not sure this is right cs to set here...
  model.set_crystal_symmetry(crystal_symmetry)
  mi_object = model_idealization(
      model = model,
      # pdb_h=pdb_h,
      # cif_objects=input_objects.cif_objects,
      map_data = map_data,
      # shift_manager = shift_manager,
      # crystal_symmetry = crystal_symmetry,
      # ss_annotation = pdb_input.extract_secondary_structure(),
      params=work_params,
      log=log,
      verbose=False)
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
