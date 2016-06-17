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


master_params_str = """
file_name = None
  .type = path
  .multiple = True
  .short_caption = Model file
  .style = file_type:pdb bold input_file
trim_alternative_conformations = False
  .type = bool
  .help = Leave only atoms with empty altloc
include scope mmtbx.secondary_structure.build.model_idealization_master_phil_str
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
  work_params.model_idealization.enabled=True
  # work_params.model_idealization.file_name_before_regularization="before.pdb"

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
  if work_params.trim_alternative_conformations:
    asc = pdb_h_raw.atom_selection_cache()
    sel = asc.selection("altloc ' '")
    pdb_h = pdb_h_raw.select(sel)
    print >> log, "Atoms in original/working model: %d/%d" % (
        pdb_h_raw.atoms().size(), pdb_h.atoms().size())
  else:
    pdb_h = pdb_h_raw
  pdb_h.reset_atom_i_seqs()
  # couple checks if combined pdb_h is ok
  o_c = pdb_h.overall_counts()
  o_c.raise_duplicate_atom_labels_if_necessary()
  o_c.raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary()
  o_c.raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary()
  o_c.raise_improper_alt_conf_if_necessary()

  ann = ioss.annotation.from_phil(
      phil_helices=work_params.secondary_structure.protein.helix,
      phil_sheets=work_params.secondary_structure.protein.sheet,
      pdb_hierarchy=pdb_h)
  xrs = pdb_h.extract_xray_structure(crystal_symmetry=cs)
  if ann.get_n_helices() + ann.get_n_sheets() == 0:
    ann = pdb_input.extract_secondary_structure()
  if ann is None or ann.get_n_helices() + ann.get_n_sheets() == 0:
    print >> log, "No secondary structure annotations found."
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
        hierarchy=pdb_h,
        xrs=xrs,
        original_pdb_h=pdb_h,
        excl_string_selection=negate_selection,
        log=None,
        ss_annotation=None)
  else:
    ssb.substitute_ss(
        real_h=pdb_h,
        xray_structure=xrs,
        ss_annotation=ann,
        params=work_params.model_idealization,
        cif_objects=inputs.cif_objects,
        verbose=True,
        log=log)
    log.flush()

  # Write resulting pdb file.
  write_whole_pdb_file(
      file_name="%s_ss_substituted.pdb" % os.path.basename(pdb_file_names[0]),
      pdb_hierarchy=pdb_h,
      crystal_symmetry=cs,
      ss_annotation=ann)

  loop_ideal = loop_idealization(
      pdb_hierarchy=pdb_h,
      params=work_params.loop_idealization,
      secondary_structure_annotation=ann,
      log=log,
      verbose=True)
  log.flush()

  # shifting back if needed
  if shift_vector is not None:
    print >> log, "Shifting molecule back"
    atoms = loop_ideal.resulting_pdb_h.atoms()
    sites_cart = atoms.extract_xyz()
    atoms.set_xyz(new_xyz=sites_cart-shift_vector)

  cs_to_write = cs if shift_vector is None else None
  write_whole_pdb_file(
    file_name="%s_ss_all_idealized.pdb" % os.path.basename(pdb_file_names[0]),
    pdb_hierarchy=loop_ideal.resulting_pdb_h,
    crystal_symmetry=cs_to_write,
    ss_annotation=ann)

  print >> log, "All done."
if __name__ == "__main__":
  run(sys.argv[1:])
