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
use_starting_model_for_final_gm = True
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

    # various checks, shifts, trims
    self.cs = self.pdb_input.crystal_symmetry()
    # check self.cs (copy-paste from secondary_sturcure_restraints)
    corrupted_cs = False
    if self.cs is not None:
      if [self.cs.unit_cell(), self.cs.space_group()].count(None) > 0:
        corrupted_cs = True
        self.cs = None
      elif cs.unit_cell().volume() < 10:
        corrupted_cs = True
        self.cs = None

    self.shift_vector = None
    if self.cs is None:
      if corrupted_cs:
        print >> self.log, "Symmetry information is corrupted, "
      else:
        print >> self.log, "Symmetry information was not found, "
      print >> self.log, "putting molecule in P1 box."
      self.log.flush()
      from cctbx import uctbx
      atoms = self.pdb_input.atoms()
      box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
        sites_cart=atoms.extract_xyz(),
        buffer_layer=3)
      atoms.set_xyz(new_xyz=box.sites_cart)
      self.cs = box.crystal_symmetry()
      self.shift_vector = box.shift_vector
    pdb_h_raw = self.pdb_input.construct_hierarchy()

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
    self.original_hierarchy = self.pdb_h.deep_copy()
    self.original_hierarchy.reset_atom_i_seqs()
    # couple checks if combined pdb_h is ok
    o_c = self.pdb_h.overall_counts()
    o_c.raise_duplicate_atom_labels_if_necessary()
    o_c.raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary()
    o_c.raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary()
    o_c.raise_improper_alt_conf_if_necessary()

    self.init_model_statistics = geometry_no_grm(
        pdb_hierarchy=pdb_h_raw,
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
    print "Found NCS groups:"
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

    ann = ioss.annotation.from_phil(
        phil_helices=self.params.secondary_structure.protein.helix,
        phil_sheets=self.params.secondary_structure.protein.sheet,
        pdb_hierarchy=self.pdb_h)

    xrs = (master_pdb_h if master_pdb_h is not None else self.pdb_h).extract_xray_structure(crystal_symmetry=self.cs)
    if ann.get_n_helices() + ann.get_n_sheets() == 0:
      ann = self.pdb_input.extract_secondary_structure()
    original_ann = None
    if ann is not None:
      print "Annotation in idealization"
      print ann.as_pdb_str()
      original_ann = ann.deep_copy()
    if ann is not None:
      ann.remove_empty_annotations(
          hierarchy=master_pdb_h if master_pdb_h is not None else self.pdb_h)
    if (ann is None or
        ann.get_n_helices() + ann.get_n_sheets() == 0 or
        not self.params.ss_idealization.enabled):
      print >> self.log, "No secondary structure annotations found or SS idealization is disabled."
      print >> self.log, "Secondary structure substitution step will be skipped"
      self.log.flush()
      # here we want to do geometry minimization anyway!
      outlier_selection_txt = mmtbx.building.loop_closure.utils. \
        rama_outliers_selection(self.pdb_h, None, 1)
      print "outlier_selection_txt", outlier_selection_txt
      negate_selection = "all"
      if outlier_selection_txt != "" and outlier_selection_txt is not None:
        negate_selection = "not (%s)" % outlier_selection_txt
      minimize_wrapper_for_ramachandran(
          hierarchy=master_pdb_h if master_pdb_h is not None else self.pdb_h,
          xrs=xrs,
          original_pdb_h=master_pdb_h if master_pdb_h is not None else self.pdb_h,
          excl_string_selection=negate_selection,
          log=None,
          ss_annotation=ann)
    else:
      self.params.ss_idealization.file_name_before_regularization = \
          "%s_ss_before_reg.pdb" % self.params.output_prefix
      ssb.substitute_ss(
          real_h=master_pdb_h if master_pdb_h is not None else self.pdb_h,
          xray_structure=xrs,
          ss_annotation=ann,
          params=self.params.ss_idealization,
          fix_rotamer_outliers=True,
          cif_objects=self.cif_objects,
          verbose=True,
          log=self.log)
      self.log.flush()

    # Write resulting pdb file.
    pdb_h_shifted = (master_pdb_h if master_pdb_h is not None else self.pdb_h).deep_copy()
    pdb_h_shifted.reset_atom_i_seqs()
    write_whole_pdb_file(
        file_name="%s_ss_substituted_nosh.pdb" % self.params.output_prefix,
        pdb_hierarchy=pdb_h_shifted,
        crystal_symmetry=self.cs,
        ss_annotation=ann)
    if self.shift_vector is not None:
      # print >> self.log, "Shifting molecule back"
      atoms = pdb_h_shifted.atoms()
      sites_cart = atoms.extract_xyz()
      atoms.set_xyz(new_xyz=sites_cart-self.shift_vector)

    write_whole_pdb_file(
        file_name="%s_ss_substituted.pdb" % self.params.output_prefix,
        pdb_hierarchy=pdb_h_shifted,
        crystal_symmetry=self.cs,
        ss_annotation=ann)

    self.params.loop_idealization.minimize_whole = not using_ncs
    # self.params.loop_idealization.enabled = False
    # self.params.loop_idealization.variant_search_level = 0
    loop_ideal = loop_idealization(
        pdb_hierarchy=master_pdb_h if master_pdb_h is not None else self.pdb_h,
        params=self.params.loop_idealization,
        secondary_structure_annotation=ann,
        log=self.log,
        verbose=True)
    self.log.flush()

    # fixing remaining rotamer outliers
    fixed_rot_pdb_h = loop_ideal.resulting_pdb_h.deep_copy()
    fixed_rot_pdb_h.reset_atom_i_seqs()

    write_whole_pdb_file(
      file_name="%s_before_rot_fixing.pdb" % self.params.output_prefix,
      pdb_hierarchy=fixed_rot_pdb_h,
      crystal_symmetry=self.cs,
      ss_annotation=ann)

    if self.params.additionally_fix_rotamer_outliers:
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

    cs_to_write = self.cs if self.shift_vector is None else None
    write_whole_pdb_file(
      file_name="%s_ss_all_idealized.pdb" % self.params.output_prefix,
      pdb_hierarchy=fixed_rot_pdb_h,
      crystal_symmetry=cs_to_write,
      ss_annotation=ann)


    ref_hierarchy_for_final_gm = self.original_hierarchy
    if not self.params.use_starting_model_for_final_gm:
      ref_hierarchy_for_final_gm = self.pdb_h
    ref_hierarchy_for_final_gm.reset_atom_i_seqs()
    if using_ncs:
      print >> self.log, "Using ncs"
      ssb.set_xyz_smart(self.pdb_h, fixed_rot_pdb_h)
      # multiply back and do geometry_minimization for the whole molecule
      for ncs_gr in ncs_restr_group_list:
        master_h = self.pdb_h.select(ncs_gr.master_iselection)
        for c in ncs_gr.copies:
          new_sites = master_h.atoms().extract_xyz()
          new_c_sites = c.r.elems * new_sites + c.t
          self.pdb_h.select(c.iselection).atoms().set_xyz(new_c_sites)
      # and do geometry_minimization
      print >> self.log, "Minimizing whole model"
      self.log.flush()
      minimize_wrapper_for_ramachandran(
          hierarchy=self.pdb_h,
          xrs=xrs,
          # original_pdb_h=self.original_hierarchy,
          original_pdb_h=ref_hierarchy_for_final_gm,
          excl_string_selection=loop_ideal.ref_exclusion_selection,
          log=self.log,
          ss_annotation=original_ann)
    else:
      # still need to run gm if rotamers were fixed
      print >> self.log, "Not using ncs"
      print "loop_ideal.ref_exclusion_selection", loop_ideal.ref_exclusion_selection
      if self.params.additionally_fix_rotamer_outliers:
        ssb.set_xyz_smart(self.pdb_h, fixed_rot_pdb_h) # get out of if?
        minimize_wrapper_for_ramachandran(
            hierarchy=self.pdb_h,
            xrs=xrs,
            original_pdb_h=ref_hierarchy_for_final_gm,
            excl_string_selection=loop_ideal.ref_exclusion_selection,
            log=self.log,
            ss_annotation=original_ann)

    # shifting back if needed
    if self.shift_vector is not None:
      print >> self.log, "Shifting molecule back"
      self.log.flush()
      atoms = self.pdb_h.atoms()
      sites_cart = atoms.extract_xyz()
      atoms.set_xyz(new_xyz=sites_cart-self.shift_vector)

    write_whole_pdb_file(
        file_name="%s_ss_all_idealized_multip.pdb" % self.params.output_prefix,
        pdb_hierarchy=self.pdb_h,
        crystal_symmetry=cs_to_write,
        ss_annotation=original_ann)


    self.final_model_statistics = geometry_no_grm(
        pdb_hierarchy=self.pdb_h,
        molprobity_scores=True)
    self.time_for_run = time() - t_0

  def print_stat_comparison(self):
    print >> self.log, "                          Before     After"
                                      # '      2.33     22.33'
    print >> self.log, "Molprobity Score     :%10.2f%10.2f" % (
        self.init_model_statistics.mpscore,
        self.final_model_statistics.mpscore)
    print >> self.log, "Clashscore           :%10.2f%10.2f" % (
        self.init_model_statistics.clashscore,
        self.final_model_statistics.clashscore)
    print >> self.log, "CBeta deviations     :%10d%10d" % (
        self.init_model_statistics.c_beta_dev,
        self.final_model_statistics.c_beta_dev)
    print >> self.log, "Ramachandran outliers:%10.2f%10.2f" % (
        self.init_model_statistics.ramachandran_outliers,
        self.final_model_statistics.ramachandran_outliers)
    print >> self.log, "Rotamer outliers     :%10.2f%10.2f" % (
        self.init_model_statistics.rotamer_outliers,
        self.final_model_statistics.rotamer_outliers)
    print >> self.log, "Cis-prolines         :%10d%10d" % (
        self.init_model_statistics.n_cis_proline,
        self.final_model_statistics.n_cis_proline)
    print >> self.log, "Cis-general          :%10d%10d" % (
        self.init_model_statistics.n_cis_general,
        self.final_model_statistics.n_cis_general)
    print >> self.log, "Twisted prolines     :%10d%10d" % (
        self.init_model_statistics.n_twisted_proline,
        self.final_model_statistics.n_twisted_proline)
    print >> self.log, "Twisted general      :%10d%10d" % (
        self.init_model_statistics.n_twisted_general,
        self.final_model_statistics.n_twisted_general)

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
  mi_object.print_runtime()








  # add hydrogens if needed



  print >> log, "All done."
if __name__ == "__main__":
  run(sys.argv[1:])
