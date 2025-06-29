"""Substitute secondary structure elements with ideal ones"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.model_idealization

import sys, os
import datetime
from time import time
from libtbx.utils import Sorry, multi_out, null_out
from libtbx import easy_pickle, group_args
import libtbx.load_env
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO

from cctbx import crystal
from cctbx import xray
from iotbx import reflection_file_utils
from iotbx.phil import process_command_line_with_files
import iotbx.ncs
import iotbx.phil
from cctbx import maptbx, miller

from mmtbx.secondary_structure.build import ss_idealization as ssb
from mmtbx.secondary_structure import manager, sec_str_master_phil
import mmtbx.utils
from mmtbx.building.loop_idealization import loop_idealization
from mmtbx.building.cablam_idealization import cablam_idealization
import mmtbx.building.loop_closure.utils
from mmtbx.refinement.geometry_minimization import minimize_wrapper_for_ramachandran
from mmtbx.refinement.real_space.individual_sites import minimize_wrapper_with_map
import mmtbx.model
import mmtbx.refinement.real_space.fit_residues
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager
from iotbx import extract_xtal_data

from elbow.command_line.ready_set import model_interface as ready_set_model_interface

from phenix.programs import phi_psi_2
import six

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
output_model_h = True
  .type = bool
  .expert_level = 2
use_map_for_reference = True
  .type = bool
  .expert_level = 1
run_minimization_first = True
  .type = bool
  .expert_level = 2
run_minimization_last = True
  .type = bool
  .expert_level = 2
use_hydrogens_in_minimization = False
  .type = bool
  .expert_level = 3
reference_map_resolution = 5
  .type = float
  .expert_level = 2
number_of_refinement_cycles = 5
  .type = int
  .expert_level = 1
cycles_to_converge = 2
  .type = int
  .help = Nuber of cycles of geometry optimization without considerable stat \
    change before stopping
  .expert_level=1
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
nonbonded_weight=10000
  .type = float
apply_all_trans = True
  .type = bool

%s
include scope mmtbx.geometry_restraints.ramachandran.old_master_phil
include scope mmtbx.secondary_structure.sec_str_master_phil_str
include scope mmtbx.building.loop_idealization.loop_idealization_master_phil_str
include scope mmtbx.building.cablam_idealization.master_phil_str
""" % turned_on_ss

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def format_usage_message(log):
  print("-"*79, file=log)
  msg = """\
phenix.model_idealization: Idealize model geometry.

Usage examples:
 phenix.model_idealization model.pdb
"""
  print(msg, file=log)
  print("-"*79, file=log)
  print(master_params().show(), file=log)

class model_idealization():
  def __init__(self,
               model, # shifted, with shift_manager
               map_data = None, # shifted map_data
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

    params = mmtbx.model.manager.get_default_pdb_interpretation_params()
    params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None

    params.pdb_interpretation.peptide_link.ramachandran_restraints = True
    params.pdb_interpretation.peptide_link.restrain_rama_outliers = self.params.restrain_rama_outliers
    params.pdb_interpretation.peptide_link.restrain_rama_allowed = self.params.restrain_rama_allowed
    params.pdb_interpretation.peptide_link.restrain_allowed_outliers_with_emsley = self.params.restrain_allowed_outliers_with_emsley
    params.pdb_interpretation.peptide_link.rama_weight = self.params.rama_weight
    params.pdb_interpretation.peptide_link.oldfield.weight_scale=self.params.oldfield.weight_scale
    params.pdb_interpretation.peptide_link.oldfield.plot_cutoff=self.params.oldfield.plot_cutoff

    params.pdb_interpretation.peptide_link.apply_peptide_plane = True
    if self.params.loop_idealization.make_all_trans:
      params.pdb_interpretation.peptide_link.apply_all_trans = self.params.apply_all_trans
    params.pdb_interpretation.nonbonded_weight = self.params.nonbonded_weight
    params.pdb_interpretation.c_beta_restraints=True
    params.pdb_interpretation.max_reasonable_bond_distance = None
    params.pdb_interpretation.ncs_search.enabled = True
    params.pdb_interpretation.ncs_search.chain_max_rmsd=4.0
    params.pdb_interpretation.ncs_search.chain_similarity_threshold=0.99
    params.pdb_interpretation.ncs_search.residue_match_radius=999.0
    params.pdb_interpretation.restraints_library.rdl = True
    params.pdb_interpretation.secondary_structure = self.params.secondary_structure
    self.params_for_model = params
    self.model.process(
      make_restraints=True, pdb_interpretation_params=params)


    self.original_hierarchy = self.model.get_hierarchy().deep_copy() # original pdb_h, without any processing
    self.original_boxed_hierarchy = None # original and boxed (if needed)

    self.filtered_ncs_restr_group_list = []

    self.init_ss_annotation = self.model.get_ss_annotation()

    # various checks, shifts, trims
    self.cs = self.original_cs = self.model.crystal_symmetry()

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
    # o_c.raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary()
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
    self.shift_vector = None
    if self.cs is None:
      assert self.model.get_shift_manager() is None
      # should it happen here?
      if corrupted_cs:
        print("Symmetry information is corrupted, ", file=self.log)
      else:
        print("Symmetry information was not found, ", file=self.log)
      print("putting molecule in P1 box.", file=self.log)
      self.log.flush()
      from cctbx import uctbx
      box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
        sites_cart=self.model.get_sites_cart(),
        buffer_layer=3)
      # Creating new xrs from box, inspired by extract_box_around_model_and_map
      sp = crystal.special_position_settings(box.crystal_symmetry())
      sites_frac = box.sites_frac()
      xrs_box = self.model.get_xray_structure().replace_sites_frac(box.sites_frac())
      xray_structure_box = xray.structure(sp, xrs_box.scatterers())
      self.model.set_xray_structure(xray_structure_box)
      self.cs = box.crystal_symmetry()
      self.shift_vector = box.shift_vector

    if self.shift_vector is not None and self.params.debug:
      self.model.pdb_or_mmcif_string_info(
        target_filename="%s_boxed.pdb" % self.params.output_prefix,
        write_file=True)

    if self.params.trim_alternative_conformations:
      self.model.remove_alternative_conformations(always_keep_one_conformer=True)

    self.model = self.model.remove_hydrogens()
    self.model_h = None

    self.time_for_init = time()-t_0

  def get_statistics(self, model):
    # should we shift here? No
    # should we multiply NCS here? No
    geometry = model.geometry_statistics().result()
    motifs = phi_psi_2.results_manager(model=model, log=null_out()).get_motif_count()
    mcounts = motifs.get_counts()
    res = {}
    # TODO is mcounts a dict ? If not consider changing back
    for key, value in six.iteritems(mcounts):
      res[key] = value.percent
    geometry.merge(group_args(**res))
    return geometry

  def prepare_user_map(self):
    print("Preparing user map...", file=self.log)
    self.map_shift_manager = mmtbx.utils.shift_origin(
      map_data         = self.user_supplied_map,
      xray_structure   = self.model.get_xray_structure(),
      crystal_symmetry = self.cs)
    if(self.map_shift_manager.shift_cart is not None):
      # Need to figure out way to save the shift to shift back
      # and apply it to whole_pdb, master_pdb, etc. Don't forget about
      # boxing hierarchy when symmetry is not available or corrupted...
      raise Sorry("Map origin is not at (0,0,0). This is not implemented for model_idealization")
    map_data = self.map_shift_manager.map_data
    self.reference_map = map_data
    self.master_map = self.reference_map.deep_copy()
    if self.model.ncs_constraints_present():
      mask = maptbx.mask(
              xray_structure=self.model.get_xray_structure().select(self.model.get_master_selection()),
              n_real=self.master_map.focus(),
              mask_value_inside_molecule=1,
              mask_value_outside_molecule=-1,
              solvent_radius=0,
              atom_radius=1.)
      self.master_map = self.reference_map * mask
      if self.params.debug:
        iotbx.mrcfile.write_ccp4_map(
            file_name="%s_3_master.map" % self.params.output_prefix,
            unit_cell=self.cs.unit_cell(),
            space_group=self.cs.space_group(),
            map_data=self.master_map,
            labels=flex.std_string([""]))
        iotbx.mrcfile.write_ccp4_map(
            file_name="%s_reference.map" % self.params.output_prefix,
            unit_cell=self.cs.unit_cell(),
            space_group=self.cs.space_group(),
            map_data=self.reference_map,
            labels=flex.std_string([""]))
      self.master_map = map_data

  def prepare_init_reference_map(self):
    xrs = self.model.get_xray_structure().deep_copy_scatterers()
    pdb_h = self.model.get_hierarchy().deep_copy()
    if self.user_supplied_map is not None:
      print("Using user-supplied map for initial GM.", file=self.log)
      self.init_ref_map = self.reference_map
      return
    print("Preparing map for initial GM...", file=self.log)
    asc = self.model.get_atom_selection_cache()
    outlier_selection_txt = mmtbx.building.loop_closure.utils. \
          rama_score_selection(pdb_h, self.model.get_ramachandran_manager(), "outlier",1)
    rama_out_sel = asc.selection(outlier_selection_txt)

    allowed_selection_txt = mmtbx.building.loop_closure.utils. \
          rama_score_selection(pdb_h, self.model.get_ramachandran_manager(), "allowed",0)
    rama_allowed_sel = asc.selection(allowed_selection_txt)


    side_chain_no_cb_selection = ~ self.model.sel_backbone()
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

  def prepare_reference_map_3(self):
    """ with ramachandran outliers """
    xrs = self.model.get_xray_structure().deep_copy_scatterers()
    pdb_h = self.model.get_hierarchy()
    print("Preparing reference map, method 3", file=self.log)
    outlier_selection_txt = mmtbx.building.loop_closure.utils. \
          rama_score_selection(pdb_h, self.model.get_ramachandran_manager(), "outlier",1)
    asc = self.model.get_atom_selection_cache()
    # print >> self.log, "rama outlier selection:", outlier_selection_txt
    rama_out_sel = asc.selection(outlier_selection_txt)
    xrs=xrs.set_b_iso(value=50)

    side_chain_no_cb_selection = ~ self.model.sel_backbone()
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
    if self.model.ncs_constraints_present():
      # here we are negating non-master part of the model
      # self.master_sel=master_sel
      # self.master_map = self.reference_map.deep_copy()
      mask = maptbx.mask(
              xray_structure=xrs.select(self.model.get_master_selection()),
              n_real=self.master_map.focus(),
              mask_value_inside_molecule=1,
              mask_value_outside_molecule=-1,
              solvent_radius=0,
              atom_radius=1.)
      self.master_map = self.reference_map * mask
      if self.params.debug:
        iotbx.mrcfile.write_ccp4_map(
            file_name="%s_3_master.map" % self.params.output_prefix,
            unit_cell=xrs.unit_cell(),
            space_group=xrs.space_group(),
            map_data=self.master_map,
            labels=flex.std_string([""]))

  def update_ss_in_grm(self, ss_annotation):
    self.set_ss_restraints(ss_annotation)

  def set_ss_restraints(self, ss_annotation, params=None):
    log = self.log
    if not self.verbose:
      log = null_out()
    if self.params.use_ss_restraints and ss_annotation is not None:
      ss_manager = manager(
          pdb_hierarchy=self.model.get_hierarchy(),
          geometry_restraints_manager=self.model.get_restraints_manager().geometry,
          sec_str_from_pdb_file=ss_annotation,
          params=None,
          mon_lib_srv=self.model.get_mon_lib_srv(),
          verbose=-1,
          log=log)
      self.model.get_restraints_manager().geometry.set_secondary_structure_restraints(
          ss_manager=ss_manager,
          hierarchy=self.model.get_hierarchy(),
          log=log)

  def _setup_model_h(self):
    if self.model_h is not None:
      return
    if not self.model.has_hd():
      # runs reduce internally
      assert (libtbx.env.has_module(name="reduce"))
      assert (libtbx.env.has_module(name="elbow"))
      self.model_h = ready_set_model_interface(
          model=self.model,
          params=["add_h_to_water=False",
                  "optimise_final_geometry_of_hydrogens=False"],
          )
    else:
      self.model_h = self.model.deep_copy()
    params_h = mmtbx.model.manager.get_default_pdb_interpretation_params()
    params_h.pdb_interpretation = self.model.get_current_pdb_interpretation_params().pdb_interpretation
    # customization for model with H
    params_h.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
    params_h.pdb_interpretation.max_reasonable_bond_distance = None
    params_h.pdb_interpretation.use_neutron_distances=True
    params_h.pdb_interpretation.ncs_search = self.params_for_model.pdb_interpretation.ncs_search
    params_h.pdb_interpretation.ncs_search.exclude_selection="water"
    #self.model_h.set_pdb_interpretation_params(params_h)
    self.model_h.get_restraints_manager()
    self.model_h.idealize_h_riding()
    self.model_h.setup_ncs_constraints_groups(filter_groups=True)
    self.model_h._update_master_sel()
    if self.params.debug:
      self.shift_and_write_result(
        model = self.model_h,
        fname_suffix="model_h")

  def _update_model_h(self):
    if self.model_h is None:
      self._setup_model_h()
    # transfer coords model -> model_h
    sc = self.model_h.get_sites_cart()
    sc.set_selected(~self.model_h.get_hd_selection(), self.model.get_sites_cart())
    self.model_h.set_sites_cart(sc)
    self.model_h.idealize_h_riding()

  def _update_model_from_model_h(self):
    self.model.set_sites_cart(
      sites_cart = self.model_h.get_hierarchy().select(~self.model_h.get_hd_selection()).atoms().extract_xyz())
    self.model.set_sites_cart_from_hierarchy(multiply_ncs=True)

  def idealize_rotamers(self):
    print("Fixing rotamers...", file=self.log)
    self.log.flush()
    if self.params.debug:
      self.shift_and_write_result(
        model = self.model,
        fname_suffix="just_before_rota")

    self._update_model_h()
    rotman = mmtbx.idealized_aa_residues.rotamer_manager.load(
          rotamers="favored")
    self.model_h.process(make_restraints=True)
    o = mmtbx.refinement.real_space.side_chain_fit_evaluator(
      pdb_hierarchy      = self.model_h.get_hierarchy(),
      crystal_symmetry   = self.model.crystal_symmetry(),
      rotamer_evaluator  = rotman.rotamer_evaluator,
      map_data           = self.master_map)
    result = mmtbx.refinement.real_space.fit_residues.run(
        vdw_radii         = self.model_h.get_vdw_radii(),
        bselection        = o.sel_all(),
        pdb_hierarchy     = self.model_h.get_hierarchy(),
        crystal_symmetry  = self.model.crystal_symmetry(),
        map_data          = self.master_map,
        rotamer_manager   = rotman,
        rotatable_hd      = self.model_h.rotatable_hd_selection(iselection=False),
        sin_cos_table     = scitbx.math.sin_cos_table(n=10000),
        backbone_sample   = False,
        mon_lib_srv       = self.model_h.get_mon_lib_srv(),
        log               = self.log)
    self.model_h.set_sites_cart_from_hierarchy()
    self._update_model_from_model_h()
    if self.params.debug:
      self.shift_and_write_result(
          model = self.model,
          fname_suffix="rota_ideal")

  def run(self):
    t_0 = time()
    self.ann = self.model.get_ss_annotation()
    self._setup_model_h()
    self.model.set_restraint_objects(self.model_h.get_restraint_objects())

    self.model.process(make_restraints=True)
    # set SS restratins
    self.set_ss_restraints(self.ann)

    self.model.setup_ncs_constraints_groups()

    self.init_model_statistics = self.get_statistics(self.model)

    #
    # Cablam idealization
    #
    if self.params.debug:
      self.shift_and_write_result(
          model = self.model,
          fname_suffix="start")
      self.shift_and_write_result(
          model = self.model_h,
          fname_suffix="start_h")
    self.params.cablam_idealization.find_ss_after_fixes = False
    ci_results = cablam_idealization(
        model=self.model,
        params=self.params.cablam_idealization,
        log=self.log).get_results()
    self.model = ci_results.model
    self.after_cablam_statistics = self.get_statistics(self.model)
    if self.params.debug:
      self.shift_and_write_result(
          model = self.model,
          fname_suffix="cablam_id")


    # Here we are preparing maps if needed.
    if self.user_supplied_map is not None:
      self.prepare_user_map()

    if self.reference_map is None and self.params.use_map_for_reference:
      self.prepare_reference_map_3()

    if self.params.run_minimization_first:
      # running simple minimization and updating all
      # self.master, self.working, etc...
      # self.whole_pdb_h.reset_atom_i_seqs()
      if self.init_ref_map is None:
        self.prepare_init_reference_map()
      print("Minimization first", file=self.log)
      self.minimize(
          model=self.model,
          original_pdb_h=self.original_hierarchy,
          excl_string_selection=None, # don't need if we have map
          reference_map=self.init_ref_map,
          )
      self.init_gm_model_statistics = self.get_statistics(self.model)
      if self.params.debug:
        self.shift_and_write_result(
            model = self.model,
            fname_suffix="init_gm")

    if (self.init_gm_model_statistics is not None
        and self.init_gm_model_statistics.ramachandran.outliers == 0
        and self.init_gm_model_statistics.omega.twisted_general <= 0.01
        and self.init_gm_model_statistics.omega.twisted_proline <= 0.01
        and self.init_gm_model_statistics.omega.cis_general <= 0.01
        and self.init_gm_model_statistics.omega.cis_proline <= 0.01
        and self.init_gm_model_statistics.rotamer.outliers <= 0.01):
      print("Simple minimization was enough", file=self.log)
      # Early exit!!!
      self.shift_and_write_result(
          model=self.model,
          fname_suffix="all_idealized")
      if self.params.output_model_h:
        self.shift_and_write_result(
            model=self.model_h,
            fname_suffix="all_idealized_h")

      self.final_model_statistics = self.get_statistics(self.model)
      self.time_for_run = time() - t_0
      if self.params.output_pkl:
        easy_pickle.dump(
            file_name="%s.pkl" % self.params.output_prefix,
            obj = self.get_stats_obj())
      return

    self.filtered_whole_ann = None
    if self.ann is not None:
      self.filtered_whole_ann = self.ann.deep_copy()
      print("Original SS annotation", file=self.log)
      print(self.ann.as_pdb_str(), file=self.log)
      if self.params.filter_input_ss:
        self.filtered_whole_ann = self.ann.filter_annotation(
            hierarchy=self.model.get_hierarchy(),
            asc=self.model.get_atom_selection_cache())
      print("Filtered SS annotation", file=self.log)
      print(self.filtered_whole_ann.as_pdb_str(), file=self.log)
      self.model.set_ss_annotation(self.filtered_whole_ann)

    # getting grm with SS restraints
    self.update_ss_in_grm(self.filtered_whole_ann)

    if (self.ann is None or
        self.ann.get_n_helices() + self.ann.get_n_sheets() == 0 or
        not self.params.ss_idealization.enabled):
      print("No secondary structure annotations found or SS idealization is disabled.", file=self.log)
      print("Secondary structure substitution step will be skipped", file=self.log)
      self.log.flush()
      # here we want to do geometry minimization anyway!
      negate_selection = None
      if self.reference_map is None:
        outlier_selection_txt = mmtbx.building.loop_closure.utils. \
          rama_score_selection(self.model.get_hierarchy(), self.model.get_ramachandran_manager(), "outlier",1)
        print("outlier_selection_txt", outlier_selection_txt, file=self.log)
        negate_selection = "all"
        if outlier_selection_txt != "" and outlier_selection_txt is not None:
          negate_selection = "not (%s)" % outlier_selection_txt
      # if self.params.run_minimization_first:
      # self.minimize(
      #     model=self.model,
      #     original_pdb_h=self.whole_pdb_h,
      #     ncs_restraints_group_list=self.filtered_ncs_restr_group_list,
      #     excl_string_selection=negate_selection,
      #     reference_map=self.reference_map)
    else:
      if self.params.debug:
        self.params.ss_idealization.file_name_before_regularization = \
            "%s_ss_before_reg.pdb" % self.params.output_prefix # PDB OK
      self.params.ss_idealization.skip_good_ss_elements = True
      ssb.substitute_ss(
          model = self.model,
          params=self.params.ss_idealization,
          reference_map=self.master_map,
          log=self.log)
      self.log.flush()

    self.after_ss_idealization = self.get_statistics(self.model)
    self.shift_and_write_result(
          model=self.model,
          fname_suffix="ss_ideal_stat")

    # Write resulting pdb file.
    if self.params.debug:
      self.shift_and_write_result(
          model=self.model,
          fname_suffix="ss_ideal",
          )
    # self.params.loop_idealization.minimize_whole = not self.model.ncs_constraints_present() and self.params.loop_idealization.minimize_whole
    self.params.loop_idealization.debug = self.params.debug or self.params.loop_idealization.debug
    # self.params.loop_idealization.enabled = False
    # self.params.loop_idealization.variant_search_level = 0
    print("Starting loop idealization", file=self.log)
    loop_ideal = loop_idealization(
        self.model,
        params=self.params.loop_idealization,
        reference_map=self.master_map,
        log=self.log,
        verbose=True)
    self.log.flush()
    if self.params.debug:
      self.shift_and_write_result(
          model = self.model,
          fname_suffix="rama_ideal")
    self.after_loop_idealization = self.get_statistics(self.model)

    # fixing remaining rotamer outliers
    if (self.params.additionally_fix_rotamer_outliers and
        self.after_loop_idealization.rotamer.outliers > 0.004):
      self.idealize_rotamers()



    self.after_rotamer_fixing = self.get_statistics(self.model)
    ref_hierarchy_for_final_gm = self.original_boxed_hierarchy
    if not self.params.use_starting_model_for_final_gm:
      ref_hierarchy_for_final_gm = self.model.get_hierarchy().deep_copy()
    ref_hierarchy_for_final_gm.reset_atom_i_seqs()

    if self.model.ncs_constraints_present():
      print("Using ncs", file=self.log)
      # assert 0
    else:
      print("Not using ncs", file=self.log)
      # assert 0

    # need to update SS manager for the whole model here.
    if self.params.use_ss_restraints:
      ss_params = sec_str_master_phil.fetch().extract()
      ss_params.secondary_structure.protein.remove_outliers = not self.params.ss_idealization.enabled
      self.set_ss_restraints(
          ss_annotation=self.filtered_whole_ann,
          params=ss_params.secondary_structure)
    if self.params.run_minimization_last:
      print("loop_ideal.ref_exclusion_selection", loop_ideal.ref_exclusion_selection, file=self.log)
      print("Minimizing whole model", file=self.log)
      self.minimize(
          model = self.model,
          original_pdb_h=ref_hierarchy_for_final_gm,
          excl_string_selection=loop_ideal.ref_exclusion_selection,
          reference_map = self.reference_map)
    self.shift_and_write_result(
        model = self.model,
        fname_suffix="all_idealized")
    if self.params.output_model_h:
      self.shift_and_write_result(
          model=self.model_h,
          fname_suffix="all_idealized_h")

    self.final_model_statistics = self.get_statistics(self.model)
    self.time_for_run = time() - t_0
    if self.params.output_pkl or self.params.debug:
      easy_pickle.dump(
          file_name="%s.pkl" % self.params.output_prefix, # PDB OK
          obj = self.get_stats_obj())

  def minimize(self,
      model,
      original_pdb_h,
      excl_string_selection,
      reference_map):
    # print "ncs_restraints_group_list", ncs_restraints_group_list
    # assert 0
    if reference_map is None:
      minimize_wrapper_for_ramachandran(
          model=model,
          original_pdb_h=original_pdb_h,
          excl_string_selection=excl_string_selection,
          number_of_cycles=self.params.number_of_refinement_cycles,
          log=self.log,
          )
      self._update_model_h()
    else:
      print("Using map as reference", file=self.log)
      self.log.flush()
      if self.params.use_hydrogens_in_minimization:
        self._update_model_h()
        mwwm = minimize_wrapper_with_map(
            model=self.model_h,
            target_map=reference_map,
            number_of_cycles=self.params.number_of_refinement_cycles,
            cycles_to_converge=self.params.cycles_to_converge,
            log=self.log)
        self._update_model_from_model_h()
      else:
        mwwm = minimize_wrapper_with_map(
            model=model,
            target_map=reference_map,
            number_of_cycles=self.params.number_of_refinement_cycles,
            cycles_to_converge=self.params.cycles_to_converge,
            log=self.log)
        self._update_model_h()

  def shift_and_write_result(self, model, fname_suffix=""):
    model.pdb_or_mmcif_string_info(
      target_filename="%s_%s.pdb" % (self.params.output_prefix, fname_suffix),
      write_file=True)
    if self.params.debug:
      model.pdb_or_mmcif_string_info(
        target_filename="%s_%s_nosh.pdb" % (self.params.output_prefix, fname_suffix),
        write_file=True,
        do_not_shift_back=True)

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
        self.model.get_hierarchy(),
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
    if self.after_cablam_statistics is not None:
      stat_obj_list.insert(1, self.after_cablam_statistics)
    return group_args(
        geoms=stat_obj_list,
        rmsds=(self.get_rmsd_from_start(), self.get_rmsd_from_start2()),
        runtime=self.time_for_init + self.time_for_run)

  def print_stat_comparison(self):
    stat_obj_list = self.get_stats_obj()
    if self.after_cablam_statistics is None:
      if self.params.run_minimization_first:
        print("                        Starting    Init GM   SS ideal    Rama      Rota     Final", file=self.log)
      else:
        print("                        Starting    SS ideal    Rama      Rota     Final", file=self.log)
    else:
      if self.params.run_minimization_first:
        print("                        Starting     Cablam   Init GM   SS ideal    Rama      Rota     Final", file=self.log)
      else:
        print("                        Starting     Cablam   SS ideal    Rama      Rota     Final", file=self.log)
    #                         Starting    SS ideal    Rama      Rota     Final
    # Molprobity Score     :      4.50      3.27      2.66      2.32      2.54
    for val_caption, val_name, val_subname, val_format in [
        ("Molprobity Score", "molprobity_score", "", "{:10.2f}"),
        ("Clashscore", "clash", "score", "{:10.2f}"),
        ("CBeta deviations", "c_beta", "outliers", "{:10.2f}"),
        ("Ramachandran outliers", "ramachandran", "outliers", "{:10.2f}"),
        ("Ramachandran allowed", "ramachandran", "allowed", "{:10.2f}"),
        ("Rotamer outliers", "rotamer", "outliers", "{:10.2f}"),
        ("Cis-prolines", "omega", "cis_proline", "{:10.2f}"),
        ("Cis-general", "omega", "cis_general", "{:10.2f}"),
        ("Twisted prolines", "omega", "twisted_proline", "{:10.2f}"),
        ("Twisted general", "omega", "twisted_general", "{:10.2f}"),
        ("CaBLAM outliers", "cablam", "outliers", "{:10.2f}"),
        ("CaBLAM disfavored", "cablam", "disfavored", "{:10.2f}"),
        ("CaBLAM CA outliers", "cablam", "ca_outliers", "{:10.2f}"),
        ("phi-psy2: Motif(10)", "MOTIF", "", "{:10.2f}"),
        ("phi-psy2: Motif(20)", "MOTIF20", "", "{:10.2f}"),
        ("phi-psy2: Motif(->)", "MOTIF...", "", "{:10.2f}"),
        ("phi-psy2: General", "GENERAL", "", "{:10.2f}"),
        ("phi-psy2: Outlier", "OUTLIER", "", "{:10.2f}"),
        ]:
      l = "%-21s:" % val_caption
      for stat_obj in stat_obj_list.geoms:
        value = 99999
        if stat_obj is not None:
          sub_class = getattr(stat_obj, val_name, None)
          if sub_class is not None:
            if val_subname != "":
              value = getattr(sub_class, val_subname, None)
            else:
              value = sub_class
          l += val_format.format(value)
        else:
          l += val_format.format(0)
      print(l, file=self.log)

  def print_runtime(self):
    print("Time taken for idealization: %s" % str(
        datetime.timedelta(seconds=int(self.time_for_init + self.time_for_run))), file=self.log)

def get_map_from_hkl(hkl_file_object, params, xrs, log):
  print("Processing input hkl file...", file=log)
  crystal_symmetry = hkl_file_object.crystal_symmetry()
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry   = True,
    reflection_files = [hkl_file_object.file_content],
    err              = StringIO())


  parameters = extract_xtal_data.data_and_flags_master_params().extract()
  if (params.data_labels is not None):
    parameters.labels = params.data_labels
  if (params.r_free_flags_labels is not None):
    parameters.r_free_flags.label = params.r_free_flags_labels
  determined_data_and_flags = extract_xtal_data.run(
    reflection_file_server = rfs,
    parameters             = parameters,
    keep_going             = True,
    working_point_group = crystal_symmetry.space_group().build_derived_point_group())
  f_obs = determined_data_and_flags.f_obs

  if (params.data_labels is None):
    params.data_labels = f_obs.info().label_string()
  r_free_flags = determined_data_and_flags.r_free_flags
  assert f_obs is not None
  print("Input data:", file=log)
  print("  Iobs or Fobs:", f_obs.info().labels, file=log)
  if (r_free_flags is not None):
    print("  Free-R flags:", r_free_flags.info().labels, file=log)
    params.r_free_flags_labels = r_free_flags.info().label_string()
  else:
    print("  Free-R flags: Not present", file=log)

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
    iotbx.mrcfile.write_ccp4_map(
        file_name="%s_21.ccp4" % params.output_prefix,
        unit_cell=crystal_symmetry.unit_cell(),
        space_group=crystal_symmetry.space_group(),
        map_data=map_data,
        labels=flex.std_string([""]))
  return map_data, crystal_symmetry

def get_map_from_map(map_file_object, params, xrs, log):
  print("Processing input CCP4 map file...", file=log)
  map_data = map_file_object.file_content.map_data()
  try:
    # map_cs = map_content.file_object.crystal_symmetry()
    map_cs = map_file_object.crystal_symmetry()
  except NotImplementedError as e:
    pass
  print("Input map min,max,mean: %7.3f %7.3f %7.3f"%\
      map_data.as_1d().min_max_mean().as_tuple(), file=log)
  if map_cs.space_group().type().number() not in [0,1]:
    print(map_cs.space_group().type().number())
    raise Sorry("Only P1 group for maps is supported.")
  map_data = map_data - flex.mean(map_data)
  sd = map_data.sample_standard_deviation()
  map_data = map_data/sd
  print("Rescaled map min,max,mean: %7.3f %7.3f %7.3f"%\
    map_data.as_1d().min_max_mean().as_tuple(), file=log)
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
    print("Masking and histogram equalizing...", file=log)
    import boost_adaptbx.boost.python as bp
    cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")
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
  pdb_cs = pdb_input.crystal_symmetry()
  crystal_symmetry = None
  map_cs = None
  map_content = input_objects.get_file(work_params.map_file_name)
  if map_content is not None:
    try:
      map_cs = map_content.crystal_symmetry()
    except NotImplementedError as e:
      pass

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


  model = mmtbx.model.manager(
      model_input = pdb_input,
      restraint_objects = input_objects.cif_objects,
      crystal_symmetry = crystal_symmetry,
      log=log)

  map_data = None
  shift_manager = None

  if map_content is not None:
    map_data, map_cs, shift_manager = get_map_from_map(
        map_content,
        work_params,
        xrs=model.get_xray_structure(),
        log=log)
    model.shift_model_and_set_crystal_symmetry(
      shift_cart=shift_manager.shift_cart)
    # model.get_hierarchy().write_pdb_file("junk_shift.pdb")

  hkl_content = input_objects.get_file(work_params.hkl_file_name) # PDB OK
  if hkl_content is not None:
    map_data, map_cs = get_map_from_hkl(
        hkl_content,
        work_params,
        xrs=model.get_xray_structure(), # here we don't care about atom order
        log=log)

  mi_object = model_idealization(
      model = model,
      map_data = map_data,
      params=work_params,
      log=log,
      verbose=False)
  mi_object.run()
  mi_object.print_stat_comparison()
  print("RMSD from starting model (backbone, all): %.4f, %.4f" % (
      mi_object.get_rmsd_from_start(), mi_object.get_rmsd_from_start2()), file=log)
  mi_object.print_runtime()
  # add hydrogens if needed ?
  print("All done.", file=log)
  log.flush()
  sys.stderr = sys.__stderr__
  log.close()

if __name__ == "__main__":
  run(sys.argv[1:])

