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

from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.rotalyze import rotalyze

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

def whole_chain_in_ncs(whole_h, master_iselection):
  m_c = whole_h.select(master_iselection)
  m_c_id = m_c.only_model().chains()[0].id
  print "m_c_id", m_c_id
  for chain in whole_h.only_model().chains():
    if chain.id == m_c_id:
      print "sizes:", chain.atoms_size(), master_iselection.size()
      if chain.atoms_size() == master_iselection.size():
        return True
      else:
        return False

def filter_ncs_restraints_group_list(whole_h, ncs_restr_group_list):
  n_gr_to_remove = []
  for i, ncs_gr in enumerate(ncs_restr_group_list):
    if not whole_chain_in_ncs(whole_h, ncs_gr.master_iselection):
      n_gr_to_remove.append(i)
  result = deepcopy(ncs_restr_group_list)
  for i in reversed(n_gr_to_remove):
    del result[i]
  return result

def get_ramachandran_scores(pdb_hierarchy):
  ramalyze_obj = ramalyze(pdb_hierarchy=pdb_hierarchy, outliers_only=False)
  return ramalyze_obj.percent_outliers, ramalyze_obj.percent_allowed, ramalyze_obj.percent_favored

def get_rotamer_scores(pdb_hierarchy):
  rotalyze_obj = rotalyze(pdb_hierarchy=pdb_hierarchy, outliers_only=False)
  return rotalyze_obj.percent_outliers

def run(args):
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
  cs = pdb_input.crystal_symmetry()
  # check cs (copy-paste from secondary_sturcure_restraints)

  corrupted_cs = False
  if cs is not None:
    if [cs.unit_cell(), cs.space_group()].count(None) > 0:
      corrupted_cs = True
      cs = None
    elif cs.unit_cell().volume() < 10:
      corrupted_cs = True
      cs = None

  shift_vector = None
  if cs is None:
    if corrupted_cs:
      print >> log, "Symmetry information is corrupted, "
    else:
      print >> log, "Symmetry information was not found, "
    print >> log, "putting molecule in P1 box."
    log.flush()
    from cctbx import uctbx
    atoms = pdb_input.atoms()
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart=atoms.extract_xyz(),
      buffer_layer=3)
    atoms.set_xyz(new_xyz=box.sites_cart)
    cs = box.crystal_symmetry()
    shift_vector = box.shift_vector
  pdb_h_raw = pdb_input.construct_hierarchy()
  init_rama_out = get_ramachandran_scores(pdb_h_raw)
  init_rota_out = get_rotamer_scores(pdb_h_raw)
  print >> log, "Initial rama/rota outliers:", init_rama_out[0], init_rota_out
  if shift_vector is not None:
    write_whole_pdb_file(
        file_name="%s_boxed.pdb" % work_params.output_prefix,
        pdb_hierarchy=pdb_h_raw,
        crystal_symmetry=cs,
        ss_annotation=pdb_input.extract_secondary_structure())
  if work_params.trim_alternative_conformations:
    asc = pdb_h_raw.atom_selection_cache()
    sel = asc.selection("altloc ' '")
    pdb_h = pdb_h_raw.select(sel)
    print >> log, "Atoms in original/working model: %d/%d" % (
        pdb_h_raw.atoms_size(), pdb_h.atoms_size())
  else:
    pdb_h = pdb_h_raw

  pdb_h.reset_atom_i_seqs()
  original_hierarchy = pdb_h.deep_copy()
  original_hierarchy.reset_atom_i_seqs()
  # couple checks if combined pdb_h is ok
  o_c = pdb_h.overall_counts()
  o_c.raise_duplicate_atom_labels_if_necessary()
  o_c.raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary()
  o_c.raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary()
  o_c.raise_improper_alt_conf_if_necessary()


  master_pdb_h = None
  ncs_obj = iotbx.ncs.input(
      hierarchy=pdb_h,
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
  filtered_ncs_restr_group_list = filter_ncs_restraints_group_list(
      pdb_h, ncs_restr_group_list)
  if len(filtered_ncs_restr_group_list) > 0:
    using_ncs = True
    master_sel = flex.bool(pdb_h.atoms_size(), True)
    for ncs_gr in filtered_ncs_restr_group_list:
      for copy in ncs_gr.copies:
        master_sel.set_selected(copy.iselection, False)
    master_pdb_h = pdb_h.select(master_sel)
    master_pdb_h.reset_atom_i_seqs()
  # for ncs_gr in ncs_restr_group_list:
  #   total_ncs_selected_atoms += ncs_gr.master_iselection.size()
  #   master_sel.extend(ncs_gr.master_iselection)
  #   for c in ncs_gr.copies:
  #     total_ncs_selected_atoms += c.iselection.size()
  # print "total_ncs_selected_atoms", total_ncs_selected_atoms
  # print "Total atoms in model", pdb_h.atoms_size()
  # if total_ncs_selected_atoms == pdb_h.atoms_size():
  #   using_ncs = True

  # print "master_sel.size()", master_sel.size()
  # print "list(master_sel)", list(master_sel)
  # master_sel = flex.size_t(sorted(list(master_sel)))
  # master_pdb_h = pdb_h.select(master_sel)
  if using_ncs:
    master_pdb_h.write_pdb_file("%s_master_h.pdb" % work_params.output_prefix)

  ann = ioss.annotation.from_phil(
      phil_helices=work_params.secondary_structure.protein.helix,
      phil_sheets=work_params.secondary_structure.protein.sheet,
      pdb_hierarchy=pdb_h)

  xrs = (master_pdb_h if master_pdb_h is not None else pdb_h).extract_xray_structure(crystal_symmetry=cs)
  if ann.get_n_helices() + ann.get_n_sheets() == 0:
    ann = pdb_input.extract_secondary_structure()
  original_ann = None
  if ann is not None:
    print "Annotation in idealization"
    print ann.as_pdb_str()
    original_ann = ann.deep_copy()
  if ann is not None:
    ann.remove_empty_annotations(
        hierarchy=master_pdb_h if master_pdb_h is not None else pdb_h)
  if (ann is None or
      ann.get_n_helices() + ann.get_n_sheets() == 0 or
      not work_params.ss_idealization.enabled):
    print >> log, "No secondary structure annotations found or SS idealization is disabled."
    print >> log, "Secondary structure substitution step will be skipped"
    log.flush()
    # here we want to do geometry minimization anyway!
    outlier_selection_txt = mmtbx.building.loop_closure.utils. \
      rama_outliers_selection(pdb_h, None, 1)
    print "outlier_selection_txt", outlier_selection_txt
    negate_selection = "all"
    if outlier_selection_txt != "" and outlier_selection_txt is not None:
      negate_selection = "not (%s)" % outlier_selection_txt
    minimize_wrapper_for_ramachandran(
        hierarchy=master_pdb_h if master_pdb_h is not None else pdb_h,
        xrs=xrs,
        original_pdb_h=master_pdb_h if master_pdb_h is not None else pdb_h,
        excl_string_selection=negate_selection,
        log=None,
        ss_annotation=ann)
  else:
    work_params.ss_idealization.file_name_before_regularization = \
        "%s_ss_before_reg.pdb" % work_params.output_prefix
    ssb.substitute_ss(
        real_h=master_pdb_h if master_pdb_h is not None else pdb_h,
        xray_structure=xrs,
        ss_annotation=ann,
        params=work_params.ss_idealization,
        fix_rotamer_outliers=True,
        cif_objects=inputs.cif_objects,
        verbose=True,
        log=log)
    log.flush()

  # Write resulting pdb file.
  pdb_h_shifted = (master_pdb_h if master_pdb_h is not None else pdb_h).deep_copy()
  pdb_h_shifted.reset_atom_i_seqs()
  write_whole_pdb_file(
      file_name="%s_ss_substituted_nosh.pdb" % work_params.output_prefix,
      pdb_hierarchy=pdb_h_shifted,
      crystal_symmetry=cs,
      ss_annotation=ann)
  if shift_vector is not None:
    # print >> log, "Shifting molecule back"
    atoms = pdb_h_shifted.atoms()
    sites_cart = atoms.extract_xyz()
    atoms.set_xyz(new_xyz=sites_cart-shift_vector)

  write_whole_pdb_file(
      file_name="%s_ss_substituted.pdb" % work_params.output_prefix,
      pdb_hierarchy=pdb_h_shifted,
      crystal_symmetry=cs,
      ss_annotation=ann)

  work_params.loop_idealization.minimize_whole = not using_ncs
  # work_params.loop_idealization.enabled = False
  # work_params.loop_idealization.variant_search_level = 0
  loop_ideal = loop_idealization(
      pdb_hierarchy=master_pdb_h if master_pdb_h is not None else pdb_h,
      params=work_params.loop_idealization,
      secondary_structure_annotation=ann,
      log=log,
      verbose=True)
  log.flush()

  # fixing remaining rotamer outliers
  fixed_rot_pdb_h = loop_ideal.resulting_pdb_h.deep_copy()
  fixed_rot_pdb_h.reset_atom_i_seqs()

  write_whole_pdb_file(
    file_name="%s_before_rot_fixing.pdb" % work_params.output_prefix,
    pdb_hierarchy=fixed_rot_pdb_h,
    crystal_symmetry=cs,
    ss_annotation=ann)

  if work_params.additionally_fix_rotamer_outliers:
    print >> log, "Processing pdb file again for fixing rotamers..."
    log.flush()
    # again get grm... - need to optimize somehow
    processed_pdb_files_srv = mmtbx.utils.\
        process_pdb_file_srv(
            crystal_symmetry= cs,
            pdb_interpretation_params = None,
            stop_for_unknowns         = False,
            log=log,
            cif_objects=None)
    processed_pdb_file, pdb_inp = processed_pdb_files_srv.\
        process_pdb_files(raw_records=flex.split_lines(fixed_rot_pdb_h.as_pdb_string()))

    grm = get_geometry_restraints_manager(
        processed_pdb_file, xrs, params=None)

    # mon_lib_srv and rotamer_manager already was created multiple times
    # to this moment in different procedures :(
    print >> log, "Fixing rotamers..."
    log.flush()
    fixed_rot_pdb_h = fix_rotamer_outliers(
        pdb_hierarchy=fixed_rot_pdb_h,
        grm=grm.geometry,
        xrs=xrs,
        mon_lib_srv=None,
        rotamer_manager=None)

  cs_to_write = cs if shift_vector is None else None
  write_whole_pdb_file(
    file_name="%s_ss_all_idealized.pdb" % work_params.output_prefix,
    pdb_hierarchy=fixed_rot_pdb_h,
    crystal_symmetry=cs_to_write,
    ss_annotation=ann)


  ref_hierarchy_for_final_gm = original_hierarchy
  if not work_params.use_starting_model_for_final_gm:
    ref_hierarchy_for_final_gm = pdb_h
  ref_hierarchy_for_final_gm.reset_atom_i_seqs()
  if using_ncs:
    print >> log, "Using ncs"
    ssb.set_xyz_smart(pdb_h, fixed_rot_pdb_h)
    # multiply back and do geometry_minimization for the whole molecule
    for ncs_gr in ncs_restr_group_list:
      master_h = pdb_h.select(ncs_gr.master_iselection)
      for c in ncs_gr.copies:
        new_sites = master_h.atoms().extract_xyz()
        new_c_sites = c.r.elems * new_sites + c.t
        pdb_h.select(c.iselection).atoms().set_xyz(new_c_sites)
    # and do geometry_minimization
    print >> log, "Minimizing whole model"
    log.flush()
    minimize_wrapper_for_ramachandran(
        hierarchy=pdb_h,
        xrs=xrs,
        # original_pdb_h=original_hierarchy,
        original_pdb_h=ref_hierarchy_for_final_gm,
        excl_string_selection=loop_ideal.ref_exclusion_selection,
        log=log,
        ss_annotation=original_ann)
  else:
    # still need to run gm if rotamers were fixed
    print >> log, "Not using ncs"
    print "loop_ideal.ref_exclusion_selection", loop_ideal.ref_exclusion_selection
    if work_params.additionally_fix_rotamer_outliers:
      ssb.set_xyz_smart(pdb_h, fixed_rot_pdb_h) # get out of if?
      minimize_wrapper_for_ramachandran(
          hierarchy=pdb_h,
          xrs=xrs,
          original_pdb_h=ref_hierarchy_for_final_gm,
          excl_string_selection=loop_ideal.ref_exclusion_selection,
          log=log,
          ss_annotation=original_ann)

  # shifting back if needed
  if shift_vector is not None:
    print >> log, "Shifting molecule back"
    log.flush()
    atoms = pdb_h.atoms()
    sites_cart = atoms.extract_xyz()
    atoms.set_xyz(new_xyz=sites_cart-shift_vector)

  write_whole_pdb_file(
      file_name="%s_ss_all_idealized_multip.pdb" % work_params.output_prefix,
      pdb_hierarchy=pdb_h,
      crystal_symmetry=cs_to_write,
      ss_annotation=original_ann)

  final_rama_out = get_ramachandran_scores(pdb_h)
  final_rota_out = get_rotamer_scores(pdb_h)
  print >> log, "Initial rama/rota outliers: %.2f %.2f" % (init_rama_out[0], init_rota_out)
  print >> log, "Final rama/rota outliers: %.2f %.2f" % (final_rama_out[0], final_rota_out)


  # add hydrogens if needed



  print >> log, "All done."
if __name__ == "__main__":
  run(sys.argv[1:])
