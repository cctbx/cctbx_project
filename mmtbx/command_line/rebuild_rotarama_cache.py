from __future__ import absolute_import, division, print_function

import os
import sys

# NB:  this can be run from the command line as "mmtbx.rebuild_rotarama_cache"
def run():
  from libtbx.utils import format_cpu_times
  from mmtbx.rotamer import rotamer_eval
  from mmtbx.rotamer import ramachandran_eval

  initial_current_working_directory = os.getcwd()
  rotamer_data_dir = rotamer_eval.find_rotarama_data_dir(optional=True)
  if rotamer_data_dir is None:
    print('  Rebuilding rotarama library skipped. Needs rotamer library.')
    return
  target_db = rotamer_eval.open_rotarama_dlite(
    rotarama_data_dir=rotamer_data_dir)
# rebuild_pickle_files(data_dir=rotamer_data_dir,
#   file_prefix="rota500-",
#   target_db=target_db,
#   amino_acids=rotamer_eval.aminoAcids)
  rebuild_pickle_files(data_dir=rotamer_data_dir,
    file_prefix="rota8000-",
    target_db=target_db,
    amino_acids=rotamer_eval.aminoAcids)
  #
  ramachandran_data_dir = rotamer_eval.find_rotarama_data_dir()
  target_db = rotamer_eval.open_rotarama_dlite(
    rotarama_data_dir=ramachandran_data_dir)
  rebuild_pickle_files(data_dir=rotamer_data_dir,
    file_prefix="rama8000-",
    target_db=target_db,
    amino_acids=ramachandran_eval.aminoAcids_8000)
# rebuild_pickle_files(data_dir=rotamer_data_dir,
#   file_prefix="rama500-",
#   target_db=target_db,
#   amino_acids=ramachandran_eval.aminoAcids)
  os.chdir(initial_current_working_directory)
  print(format_cpu_times())

def rebuild_pickle_files(data_dir, file_prefix, target_db, amino_acids):
  from libtbx import easy_pickle
  from libtbx.str_utils import show_string
  from mmtbx.rotamer.n_dim_table import NDimTable
  os.chdir(data_dir)
  print("Processing data files in %s:" % show_string(data_dir))
  for aa, aafile in amino_acids.items():
    data_file = file_prefix+aafile+".data"
    pickle_file = file_prefix+aafile+".pickle"
    pair_info = target_db.pair_info(
      source_path=data_file,
      target_path=pickle_file,
      path_prefix=data_dir)
    print("  %s -> %s:" % (data_file, pickle_file), end=' ')
    if not pair_info.needs_update:
      print("already up to date.")
    else:
      print("converting ...", end=' ')
      sys.stdout.flush()
      pair_info.start_building_target()
      ndt = NDimTable.createFromText(data_file)
      easy_pickle.dump(file_name=pickle_file, obj=ndt)
      pair_info.done_building_target()
      print("done.")
    sys.stdout.flush()
  target_db.write()

if __name__ == "__main__":
  run()
