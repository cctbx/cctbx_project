import libtbx.load_env

from libtbx import easy_pickle
from libtbx.utils import format_cpu_times
from libtbx.str_utils import show_string
from mmtbx.rotamer.n_dim_table import NDimTable
from mmtbx.rotamer import rotamer_eval
import sys

def run(args):
    rotamer_data_dir = rotamer_eval.find_rotamer_data_dir()
    print "Processing data files in %s:" % show_string(rotamer_data_dir)
    target_db = rotamer_eval.open_rotarama_dlite(
      rotamer_data_dir=rotamer_data_dir)
    for aa, aafile in rotamer_eval.aminoAcids.items():
        data_file = "rota500-"+aafile+".data"
        pickle_file = "rota500-"+aafile+".pickle"
        pair_info = target_db.pair_info(
          source_path=data_file,
          target_path=pickle_file,
          path_prefix=rotamer_data_dir)
        print "  %s -> %s:" % (data_file, pickle_file),
        if not pair_info.needs_update:
            print "already up to date."
        else:
            print "converting ...",
            sys.stdout.flush()
            pair_info.start_building_target()
            ndt = NDimTable.createFromText(data_file)
            easy_pickle.dump(file_name=pickle_file, obj=ndt)
            pair_info.done_building_target()
            print "done."
        sys.stdout.flush()
    target_db.write()
    print format_cpu_times()

if __name__ == "__main__":
    run(sys.argv[1:])
