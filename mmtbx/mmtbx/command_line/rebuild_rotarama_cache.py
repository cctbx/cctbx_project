import libtbx.load_env # required by PHENIX to set environment

from libtbx import easy_pickle
from libtbx import dlite
from libtbx.utils import Sorry
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
    target_db = dlite.target_db(os.path.join(rotamer_data_dir, "rotarama.dlite"))
    print "Converting data files in %s ..." % rotamer_data_dir
    for aa, aafile in aminoAcids.items():
        data_file = os.path.join(rotamer_data_dir, "rota500-"+aafile+".data")
        pickle_file = os.path.join(rotamer_data_dir, "rota500-"+aafile+".pickle")
        pair_info = target_db.pair_info(data_file, pickle_file)
        if not os.path.isfile(pickle_file) or pair_info.needs_update:
            print "  rota500-%s.data -> rota500-%s.pickle ..." % (aa, aa)
            pair_info.start_building_target()
            ndt = NDimTable.createFromText(data_file)
            easy_pickle.dump(file_name=pickle_file, obj=ndt)
            pair_info.done_building_target()
    target_db.write()

if __name__ == "__main__":
    run(sys.argv[1:])

