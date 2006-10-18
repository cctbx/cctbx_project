import libtbx.load_env

from libtbx import easy_pickle
from libtbx import dlite
from libtbx.utils import Sorry, format_cpu_times
from libtbx.str_utils import show_string
from mmtbx.rotamer.n_dim_table import NDimTable
import os
import sys

def run(args):
    aminoAcids = {
        'arg' : 'arg',
        'asn' : 'asn',
        'asp' : 'asp',
        'cys' : 'cys',
        'gln' : 'gln',
        'glu' : 'glu',
        'his' : 'his',
        'ile' : 'ile',
        'leu' : 'leu',
        'lys' : 'lys',
        'met' : 'met',
        'phe' : 'phetyr',
        'pro' : 'pro',
        'ser' : 'ser',
        'thr' : 'thr',
        'trp' : 'trp',
        'tyr' : 'phetyr',
        'val' : 'val',
    }
    rotamer_data_dir = libtbx.env.find_in_repositories("rotarama_data")
    if rotamer_data_dir is None:
        rotamer_data_dir = libtbx.env.find_in_repositories(os.path.join("ext_ref_files", "rotarama_data"))
    if rotamer_data_dir is None:
        raise Sorry(
            "Can't find ext_ref_files/rotarama_data/\n"
            "  Please run\n"
            "    svn co svn://quiddity.biochem.duke.edu:21/phenix/rotarama_data\n"
            "  to resolve this problem.\n")
    print "Processing data files in %s:" % show_string(rotamer_data_dir)
    os.chdir(rotamer_data_dir)
    target_db = dlite.target_db("rotarama.dlite")
    for aa, aafile in aminoAcids.items():
        data_file = "rota500-"+aafile+".data"
        pickle_file = "rota500-"+aafile+".pickle"
        pair_info = target_db.pair_info(data_file, pickle_file)
        print "  rota500-%s.data -> rota500-%s.pickle:" % (aa, aa),
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

