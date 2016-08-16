from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.model_idealization

import sys, os
from iotbx.pdb import secondary_structure as ioss
from mmtbx.secondary_structure import build as ssb
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


master_params_str = """
file_name = None
  .type = path
  .multiple = True
  .short_caption = Model file
  .style = file_type:pdb bold input_file
trim_alternative_conformations = False
  .type = bool
  .help = Leave only atoms with empty altloc
additionally_fix_rotamer_outliers = True
  .type = bool
  .help = At the late stage if rotamer is still outlier choose another one \
    with minimal clash with surrounding atoms
use_starting_model_for_final_gm = False
  .type = bool
  .help = Use supplied model for final geometry minimization. Otherwise just \
    use self.
output_prefix = None
  .type = str
include scope mmtbx.secondary_structure.build.ss_idealization_master_phil_str
include scope mmtbx.secondary_structure.sec_str_master_phil_str
include scope mmtbx.building.loop_idealization.loop_idealization_master_phil_str
"""

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
    self.after_ss_idealization = None
    self.after_loop_idealization = None
    self.after_rotamer_fixing = None
    self.final_model_statistics = None
    # various checks, shifts, trims
    self.cs = self.pdb_input.crystal_symmetry()
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
    self.original_hierarchy.reset_atom_i_seqs()
    pdb_h_raw = self.original_hierarchy.deep_copy()
    pdb_h_raw.reset_atom_i_seqs()
    self.shift_vector = None
    if self.cs is None:
      if corrupted_cs:
        print >> self.log, "Symmetry information is corrupted, "
      else:
        print >> self.log, "Symmetry information was not found, "
      print >> self.log, "putting molecule in P1 box."
      self.log.flush()
      from cctbx import uctbx
      atoms = pdb_h_raw.atoms()
      box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
        sites_cart=atoms.extract_xyz(),
        buffer_layer=3)
      atoms.set_xyz(new_xyz=box.sites_cart)
      self.cs = box.crystal_symmetry()
      self.shift_vector = box.shift_vector
    self.original_boxed_hierarchy = pdb_h_raw.deep_copy()
    self.original_boxed_hierarchy.reset_atom_i_seqs()

    if self.shift_vector is not None:
      write_whole_pdb_file(
          file_name="%s_boxed.pdb" % self.params.output_prefix,
          pdb_hierarchy=pdb_h_raw,
          crystal_symmetry=self.cs,
          ss_annotation=self.pdb_input.extract_secondary_structure())
    self.pdb_h = None
    if self.params.trim_alternative_conformations:
      asc = pdb_h_raw.atom_selection_cache()
      sel = asc.selection("altloc ' '")
      self.pdb_h = pdb_h_raw.select(sel)
      print >> self.log, "Atoms in original/working model: %d/%d" % (
          pdb_h_raw.atoms_size(), self.pdb_h.atoms_size())
    else:
      self.pdb_h = pdb_h_raw

    self.pdb_h.reset_atom_i_seqs()

    # couple checks if combined pdb_h is ok
    o_c = self.pdb_h.overall_counts()
    o_c.raise_duplicate_atom_labels_if_necessary()
    o_c.raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary()
    o_c.raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary()
    o_c.raise_improper_alt_conf_if_necessary()
    self.init_model_statistics = geometry_no_grm(
        pdb_hierarchy=iotbx.pdb.input(
          source_info=None,
          lines=pdb_h_raw.as_pdb_string()).construct_hierarchy(),
        molprobity_scores=True)
    self.time_for_init = time()-t_0


  def run(self):
    t_0 = time()

    master_pdb_h = None
    ncs_obj = iotbx.ncs.input(
        hierarchy=self.pdb_h,
        chain_max_rmsd=2.0,
        chain_similarity_threshold=0.99,
        residue_match_radius=999.0)
    print >> self.log, "Found NCS groups:"
    ncs_obj.show(format='phil')
    ncs_restr_group_list = ncs_obj.get_ncs_restraints_group_list(
        raise_sorry=False)
    using_ncs = False
    total_ncs_selected_atoms = 0
    master_sel = flex.size_t([])
    filtered_ncs_restr_group_list = self.filter_ncs_restraints_group_list(
        self.pdb_h, ncs_restr_group_list)
    if len(filtered_ncs_restr_group_list) > 0:
      using_ncs = True
      master_sel = flex.bool(self.pdb_h.atoms_size(), True)
      for ncs_gr in filtered_ncs_restr_group_list:
        for copy in ncs_gr.copies:
          master_sel.set_selected(copy.iselection, False)
      master_pdb_h = self.pdb_h.select(master_sel)
      master_pdb_h.reset_atom_i_seqs()

    if using_ncs:
      master_pdb_h.write_pdb_file("%s_master_h.pdb" % self.params.output_prefix)

    self.ann = ioss.annotation.from_phil(
        phil_helices=self.params.secondary_structure.protein.helix,
        phil_sheets=self.params.secondary_structure.protein.sheet,
        pdb_hierarchy=self.pdb_h)

    xrs = (master_pdb_h if master_pdb_h is not None else self.pdb_h).extract_xray_structure(crystal_symmetry=self.cs)
    if self.ann.get_n_helices() + self.ann.get_n_sheets() == 0:
      self.ann = self.pdb_input.extract_secondary_structure()
    self.original_ann = None
    if self.ann is not None:
      self.original_ann = self.ann.deep_copy()
      print >> self.log, "Original SS annotation"
      print >> self.log, self.original_ann.as_pdb_str()
    if self.ann is not None:
      self.ann.concatenate_consecutive_helices()
      self.ann.remove_empty_annotations(
          hierarchy=master_pdb_h if master_pdb_h is not None else self.pdb_h)
      self.ann.split_helices_with_prolines(
          hierarchy=master_pdb_h if master_pdb_h is not None else self.pdb_h,
          asc=None)
      # print >> self.log, "Splitted SS annotation"
      # print >> self.log, ann.as_pdb_str()
      self.ann.remove_short_annotations()
      print >> self.log, "Filtered SS annotation"
      print >> self.log, self.ann.as_pdb_str()
      # STOP()
    if (self.ann is None or
        self.ann.get_n_helices() + self.ann.get_n_sheets() == 0 or
        not self.params.ss_idealization.enabled):
      print >> self.log, "No secondary structure annotations found or SS idealization is disabled."
      print >> self.log, "Secondary structure substitution step will be skipped"
      self.log.flush()
      # here we want to do geometry minimization anyway!
      outlier_selection_txt = mmtbx.building.loop_closure.utils. \
        rama_outliers_selection(self.pdb_h, None, 1)
      print >> self.log, "outlier_selection_txt", outlier_selection_txt
      negate_selection = "all"
      if outlier_selection_txt != "" and outlier_selection_txt is not None:
        negate_selection = "not (%s)" % outlier_selection_txt
      minimize_wrapper_for_ramachandran(
          hierarchy=master_pdb_h if master_pdb_h is not None else self.pdb_h,
          xrs=xrs,
          original_pdb_h=master_pdb_h if master_pdb_h is not None else self.pdb_h,
          excl_string_selection=negate_selection,
          log=None,
          ss_annotation=self.ann)
    else:
      self.params.ss_idealization.file_name_before_regularization = \
          "%s_ss_before_reg.pdb" % self.params.output_prefix
      ssb.substitute_ss(
          real_h=master_pdb_h if master_pdb_h is not None else self.pdb_h,
          xray_structure=xrs,
          ss_annotation=self.ann,
          params=self.params.ss_idealization,
          fix_rotamer_outliers=True,
          cif_objects=self.cif_objects,
          verbose=True,
          log=self.log)
      self.log.flush()

    self.after_ss_idealization = geometry_no_grm(
        pdb_hierarchy=iotbx.pdb.input(
          source_info=None,
          lines=(master_pdb_h if master_pdb_h is not None else self.pdb_h).\
              as_pdb_string()).construct_hierarchy(),
        molprobity_scores=True)

    # Write resulting pdb file.
    self.shift_and_write_result(
        hierarchy=master_pdb_h if master_pdb_h is not None else self.pdb_h,
        fname_suffix="ss_ideal",)

    self.params.loop_idealization.minimize_whole = not using_ncs
    # self.params.loop_idealization.enabled = False
    # self.params.loop_idealization.variant_search_level = 0
    loop_ideal = loop_idealization(
        pdb_hierarchy=master_pdb_h if master_pdb_h is not None else self.pdb_h,
        params=self.params.loop_idealization,
        secondary_structure_annotation=self.ann,
        log=self.log,
        verbose=False)
    self.log.flush()
    # STOP()
    self.shift_and_write_result(
        hierarchy=loop_ideal.resulting_pdb_h,
        fname_suffix="rama_ideal")
    self.after_loop_idealization = geometry_no_grm(
        pdb_hierarchy=iotbx.pdb.input(
          source_info=None,
          lines=loop_ideal.resulting_pdb_h.as_pdb_string()).construct_hierarchy(),
        molprobity_scores=True)

    # fixing remaining rotamer outliers
    fixed_rot_pdb_h = loop_ideal.resulting_pdb_h.deep_copy()
    fixed_rot_pdb_h.reset_atom_i_seqs()
    if (self.params.additionally_fix_rotamer_outliers and
        self.after_loop_idealization.rotamer_outliers > 0.004):
      print >> self.log, "Processing pdb file again for fixing rotamers..."
      self.log.flush()
      # again get grm... - need to optimize somehow
      processed_pdb_files_srv = mmtbx.utils.\
          process_pdb_file_srv(
              crystal_symmetry= self.cs,
              pdb_interpretation_params = None,
              stop_for_unknowns         = False,
              log=self.log,
              cif_objects=None)
      processed_pdb_file, pdb_inp = processed_pdb_files_srv.\
          process_pdb_files(raw_records=flex.split_lines(fixed_rot_pdb_h.as_pdb_string()))

      grm = get_geometry_restraints_manager(
          processed_pdb_file, xrs, params=None)

      # mon_lib_srv and rotamer_manager already was created multiple times
      # to this moment in different procedures :(
      print >> self.log, "Fixing rotamers..."
      self.log.flush()
      fixed_rot_pdb_h = fix_rotamer_outliers(
          pdb_hierarchy=fixed_rot_pdb_h,
          grm=grm.geometry,
          xrs=xrs,
          mon_lib_srv=None,
          rotamer_manager=None)

    self.shift_and_write_result(
        hierarchy=fixed_rot_pdb_h,
        fname_suffix="rota_ideal")
    cs_to_write = self.cs if self.shift_vector is None else None
    self.after_rotamer_fixing = geometry_no_grm(
        pdb_hierarchy=iotbx.pdb.input(
          source_info=None,
          lines=fixed_rot_pdb_h.as_pdb_string()).construct_hierarchy(),
        molprobity_scores=True)

    ref_hierarchy_for_final_gm = self.original_boxed_hierarchy
    if not self.params.use_starting_model_for_final_gm:
      ref_hierarchy_for_final_gm = self.pdb_h
    ref_hierarchy_for_final_gm.reset_atom_i_seqs()
    if self.params.additionally_fix_rotamer_outliers:
      ssb.set_xyz_smart((master_pdb_h if master_pdb_h is not None else self.pdb_h), fixed_rot_pdb_h)
    if using_ncs:
      print >> self.log, "Using ncs"
      # multiply back and do geometry_minimization for the whole molecule
      for ncs_gr in ncs_restr_group_list:
        master_h = self.pdb_h.select(ncs_gr.master_iselection)
        for c in ncs_gr.copies:
          new_sites = master_h.atoms().extract_xyz()
          new_c_sites = c.r.elems * new_sites + c.t
          self.pdb_h.select(c.iselection).atoms().set_xyz(new_c_sites)
      self.log.flush()
    else:
      # still need to run gm if rotamers were fixed
      print >> self.log, "Not using ncs"

    print >> self.log, "loop_ideal.ref_exclusion_selection", loop_ideal.ref_exclusion_selection
    print >> self.log, "Minimizing whole model"
    minimize_wrapper_for_ramachandran(
        hierarchy=self.pdb_h,
        xrs=xrs,
        original_pdb_h=ref_hierarchy_for_final_gm,
        excl_string_selection=loop_ideal.ref_exclusion_selection,
        log=self.log,
        ss_annotation=self.ann)
    self.shift_and_write_result(
        hierarchy=self.pdb_h,
        fname_suffix="all_idealized")
    self.final_model_statistics = geometry_no_grm(
        pdb_hierarchy=iotbx.pdb.input(
          source_info=None,
          lines=self.pdb_h.as_pdb_string()).construct_hierarchy(),
        molprobity_scores=True)
    self.time_for_run = time() - t_0

  def shift_and_write_result(self, hierarchy, fname_suffix):
    cs_to_write = self.cs if self.shift_vector is None else None
    pdb_h_shifted = hierarchy.deep_copy()
    pdb_h_shifted.reset_atom_i_seqs()
    if self.shift_vector is not None:
      atoms = pdb_h_shifted.atoms()
      sites_cart = atoms.extract_xyz()
      atoms.set_xyz(new_xyz=sites_cart-self.shift_vector)
    # write_whole_pdb_file(
    #     file_name="%s_%s_nosh.pdb" % (self.params.output_prefix, fname_suffix),
    #     pdb_hierarchy=hierarchy,
    #     crystal_symmetry=self.cs,
    #     ss_annotation=self.ann)
    write_whole_pdb_file(
        file_name="%s_%s.pdb" % (self.params.output_prefix, fname_suffix),
        pdb_hierarchy=pdb_h_shifted,
        crystal_symmetry=cs_to_write,
        ss_annotation=self.original_ann)

  def get_rmsd_from_start(self):
    if self.rmsd_from_start is not None:
      return self.rmsd_from_start
    # calculate rmsd
    self.rmsd_from_start = ssb.calculate_rmsd_smart(
        self.original_boxed_hierarchy,
        self.pdb_h,
        backbone_only=True)
    return self.rmsd_from_start

  def get_rmsd_from_start2(self):
    return ssb.calculate_rmsd_smart(
        self.original_boxed_hierarchy,
        self.pdb_h,
        backbone_only=False)

  def print_stat_comparison(self):
    print >> self.log, "                        Starting    SS ideal    Rama      Rota     Final"
    #                         Starting    SS ideal    Rama      Rota     Final
    # Molprobity Score     :      4.50      3.27      2.66      2.32      2.54
    for val_caption, val_name, val_format in [
        ("Molprobity Score", "mpscore", "{:10.2f}"),
        ("Clashscore", "clashscore", "{:10.2f}"),
        ("CBeta deviations", "c_beta_dev", "{:10d}"),
        ("Ramachandran outliers", "ramachandran_outliers", "{:10.2f}"),
        ("Rotamer outliers", "rotamer_outliers", "{:10.2f}"),
        ("Cis-prolines", "n_cis_proline", "{:10d}"),
        ("Cis-general", "n_cis_general", "{:10d}"),
        ("Twisted prolines", "n_twisted_proline", "{:10d}"),
        ("Twisted general", "n_twisted_general", "{:10d}")]:
      l = "%-21s:" % val_caption
      for stat_obj in [
          self.init_model_statistics,
          self.after_ss_idealization,
          self.after_loop_idealization,
          self.after_rotamer_fixing,
          self.final_model_statistics
          ]:
        value = 99999
        if stat_obj is not None:
          l += val_format.format(getattr(stat_obj, val_name, 99999))
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
  log = sys.stdout
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
  # Massaging params, should be in the object probably
  #
  # !!!!!!!!!!!
  #
  work_params.ss_idealization.enabled=True
  # work_params.ss_idealization.file_name_before_regularization="before.pdb"
  if work_params.output_prefix is None:
    work_params.output_prefix = os.path.basename(pdb_file_names[0])
  if work_params.loop_idealization.output_prefix is None:
    work_params.loop_idealization.output_prefix = "%s_rama_fixed" % os.path.basename(pdb_file_names[0])
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_file_names)
  pdb_input = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  mi_object = model_idealization(
      pdb_input=pdb_input,
      params=work_params,
      log=sys.stdout,
      verbose=True)
  mi_object.run()
  mi_object.print_stat_comparison()
  print >> log, "RMSD from starting model (backbone, all): %.4f, %.4f" % (
      mi_object.get_rmsd_from_start(), mi_object.get_rmsd_from_start2())
  mi_object.print_runtime()
  # add hydrogens if needed ?
  print >> log, "All done."

if __name__ == "__main__":
  run(sys.argv[1:])
