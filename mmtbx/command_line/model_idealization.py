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

master_params_str = """
file_name = None
  .type = path
  .multiple = True
  .short_caption = Model file
  .style = file_type:pdb bold input_file
include scope mmtbx.secondary_structure.build.model_idealization_master_phil_str
include scope mmtbx.secondary_structure.sec_str_master_phil_str
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
  pdb_h = pdb_input.construct_hierarchy()
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
  if ann.get_n_helices() + ann.get_n_sheets() == 0:
    ann = pdb_input.extract_secondary_structure()
  if ann is None or ann.get_n_helices() + ann.get_n_sheets() == 0:
    raise Sorry("No secondary structure annotations found.")
  ssb.substitute_ss(
      real_h=pdb_h,
      xray_structure=pdb_input.xray_structure_simple(),
      ss_annotation=ann,
      params=work_params.model_idealization,
      cif_objects=inputs.cif_objects,
      verbose=True,
      log=log,
      )
  # Write resulting pdb file.
  write_whole_pdb_file(
      file_name="%s_idealized.pdb" % os.path.basename(pdb_file_names[0]),
      pdb_hierarchy=pdb_h,
      crystal_symmetry=pdb_input.crystal_symmetry(),
      ss_annotation=ann)
  print >> log, "All done."
if __name__ == "__main__":
  run(sys.argv[1:])
