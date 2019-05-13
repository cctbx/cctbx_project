from __future__ import division
from __future__ import print_function
import iotbx.pdb
import sys
from scitbx.array_family import flex
from mmtbx.building.loop_idealization import loop_idealization, master_phil

import mmtbx.utils
from libtbx.utils import Sorry



def run(args, log=sys.stdout):
  # print "args", args

  inputs = mmtbx.utils.process_command_line_args(args=args,
      master_params=master_phil)
  work_params = inputs.params.extract()
  inputs.params.show(prefix=" ", out=log)
  pdb_file_names = list(inputs.pdb_file_names)
  if len(pdb_file_names) == 0:
    raise Sorry("No PDB file specified")
  work_params.loop_idealization.enabled=True
  work_params.loop_idealization.number_of_ccd_trials=1
  work_params.loop_idealization.minimize_whole=False
  work_params.loop_idealization.variant_search_level=1
  # print work_params.loop_idealization.output_prefix
  # STOP()
  # work_params.ss_idealization.file_name_before_regularization="before.pdb"
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_file_names)
  pdb_input = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  pdb_h = pdb_input.construct_hierarchy()


  loop_ideal = loop_idealization(
      pdb_hierarchy=pdb_h,
      params=work_params.loop_idealization,
      log=log)
  loop_ideal.resulting_pdb_h.write_pdb_file(
      file_name="%s_very_final.pdb" % work_params.loop_idealization.output_prefix)
  print("Outlier percentages: initial, after ccd, after minimization, berkeley after ccd, berkeley after minimization:", file=log)
  print(loop_ideal.p_initial_rama_outliers, end=' ', file=log)
  print(loop_ideal.p_before_minimization_rama_outliers, end=' ', file=log)
  print(loop_ideal.p_after_minimiaztion_rama_outliers, end=' ', file=log)
  print(loop_ideal.berkeley_p_before_minimization_rama_outliers, end=' ', file=log)
  print(loop_ideal.berkeley_p_after_minimiaztion_rama_outliers, file=log)



if (__name__ == "__main__"):
  run(sys.argv[1:])
