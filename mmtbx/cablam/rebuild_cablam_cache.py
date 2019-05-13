from __future__ import division
from __future__ import print_function
import libtbx.load_env
from libtbx import easy_pickle, dlite
#from libtbx import dlite
from libtbx.utils import Sorry
from libtbx.utils import format_cpu_times
from libtbx.str_utils import show_string
from mmtbx.rotamer.n_dim_table import NDimTable
import os, sys

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
        raise Sorry("""\
Cannot find chem_data/cablam_data directory.
  CaBLAM contours not rebuilt.
""")
  #---end find cablam_data_dir---
  rebuild_pickle_files(
    data_dir=cablam_dir,
    file_prefix='cablam.8000.',
    target_db=dlite.target_db(os.path.join(cablam_dir,'cablam.dlite')),
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
    if not pair_info.needs_update: print("alreadt up to date.")
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



###from __future__ import division
###from libtbx import easy_pickle
###from libtbx.utils import format_cpu_times
###from libtbx.str_utils import show_string
###from mmtbx.rotamer.n_dim_table import NDimTable
###from mmtbx.rotamer import rotamer_eval
###from mmtbx.rotamer import ramachandran_eval
###import sys, os
###
###def run():
###  initial_current_working_directory = os.getcwd()
###  #
###  ramachandran_data_dir = rotamer_eval.find_rotarama_data_dir()
###  target_db = rotamer_eval.open_rotarama_dlite(
###    rotarama_data_dir=ramachandran_data_dir)
###  rebuild_pickle_files(data_dir=rotamer_data_dir,
###    file_prefix="rama8000-",
###    target_db=target_db,
###    amino_acids=ramachandran_eval.aminoAcids_8000)
###  rebuild_pickle_files(data_dir=rotamer_data_dir,
###    file_prefix="rama500-",
###    target_db=target_db,
###    amino_acids=ramachandran_eval.aminoAcids)
###  os.chdir(initial_current_working_directory)
###  print format_cpu_times()
###
###def rebuild_pickle_files(data_dir, file_prefix, target_db, amino_acids):
###  os.chdir(data_dir)
###  print "Processing data files in %s:" % show_string(data_dir)
###  for aa, aafile in amino_acids.items():
###    data_file = file_prefix+aafile+".data"
###    pickle_file = file_prefix+aafile+".pickle"
###    pair_info = target_db.pair_info(
###      source_path=data_file,
###      target_path=pickle_file,
###      path_prefix=data_dir)
###    print "  %s -> %s:" % (data_file, pickle_file),
###    if not pair_info.needs_update:
###      print "already up to date."
###    else:
###      print "converting ...",
###      sys.stdout.flush()
###      pair_info.start_building_target()
###      ndt = NDimTable.createFromText(data_file)
###      easy_pickle.dump(file_name=pickle_file, obj=ndt)
###      pair_info.done_building_target()
###      print "done."
###    sys.stdout.flush()
###  target_db.write()

if __name__ == "__main__":
    run()
