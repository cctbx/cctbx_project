from __future__ import absolute_import, division, print_function
import libtbx.load_env
from libtbx import easy_pickle, dlite
from libtbx.utils import format_cpu_times
from libtbx.str_utils import show_string
from mmtbx.rotamer.n_dim_table import NDimTable
import os, sys
# NB:  this can be run from the command line as "mmtbx.rebuild_cablam_cache"

file_suffixes = [
  'expected.general',
  'expected.gly',
  'expected.transpro',
  'expected.cispro',
  'expected.general_CA',
  'expected.gly_CA',
  'expected.transpro_CA',
  'expected.cispro_CA',
  'motif.loose_alpha',
  'motif.loose_beta',
  'motif.loose_threeten',
  'motif.regular_alpha',
  'motif.regular_beta',
  'motif.regular_threeten',
  'proline.cis',
  'proline.trans']

def run():
  starting_dir = os.getcwd()
  #---Find cablam_data dir---
  cablam_dir = libtbx.env.find_in_repositories(
    os.path.join('chem_data','cablam_data'))
  if cablam_dir is None:
    cablam_dir = libtbx.env.find_in_repositories('cablam_data')
    if cablam_dir is None:
      cablam_dir = libtbx.env.find_in_repositories(
        os.path.join('ext_ref_files','cablam_data'))
      if cablam_dir is None:
        print('  Rebuilding CaBLAM contours skipped. Needs chem_data/cablam_data directory.')
        return
  #---end find cablam_data_dir---

  target_db = dlite.try_loading_db(os.path.join(cablam_dir,'cablam.dlite'))
  rebuild_pickle_files(
    data_dir=cablam_dir,
    file_prefix='cablam.8000.',
    target_db=target_db,
    suffixes=file_suffixes)
  os.chdir(starting_dir)
  print(format_cpu_times())

def rebuild_pickle_files(data_dir, file_prefix, target_db, suffixes):
  os.chdir(data_dir)
  print('Processing data files in %s:' % show_string(data_dir))
  for suffix in suffixes:
    data_file =   file_prefix + suffix + '.stat'
    pickle_file = file_prefix + suffix + '.pickle'
    pair_info = target_db.pair_info(
      source_path=data_file,
      target_path=pickle_file,
      path_prefix=data_dir)
    print("  %s -> %s:" % (data_file, pickle_file), end=' ')
    if not pair_info.needs_update: print("already up to date.")
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
