import libtbx.load_env # required by PHENIX to set environment

from n_dim_table import NDimTable
from libtbx import easy_pickle
import os


class RotamerEval:

    # This is shared among all instances of RotamerEval -- a class variable.
    # It holds a LOT of read-only data, so this helps save memory.
    aaTables = {} # maps "his" to a NDimTable object for histidine, etc.

    def __init__(self):
        aaTables = RotamerEval.aaTables # for convenience
        if len(aaTables) == 0:
            #print "need to load"
            # maps aa name to file name
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
            for aa, aafile in aminoAcids.items():
                #print "Loading %s ..." % aa
                data_file = os.path.join(rotamer_data_dir, "rota500-"+aafile+".data")
                pickle_file = os.path.join(rotamer_data_dir, "rota500-"+aafile+".pickle")
                if not os.path.isfile(pickle_file):
                    ndt = NDimTable.createFromText(data_file)
                    easy_pickle.dump(file_name=pickle_file, obj=ndt)
                else:
                    ndt = easy_pickle.load(file_name=pickle_file)
                aaTables[aa] = ndt
        else:
            pass #print "already loaded"

    def evaluate(self, aaName, chiAngles):
        '''Evaluates the specified rotamer from 0.0 (worst) to 1.0 (best).

        Values below 0.01 are generally considered outliers.
        If the 3-letter amino acid name is not recognized, returns None.'''
        aaName = aaName.lower()
        if aaName in RotamerEval.aaTables:
            return RotamerEval.aaTables[aaName].valueAt(chiAngles)
        else:
            return None




def run():
    r = RotamerEval()
    return

    #tbl = r.aaTables['val']
    #print tbl.whereIs([0.5])
    #print tbl.bin2index([0])
    #for y in tbl.lookupTable[0:20]: print y
    #return

    #print r.evaluate("ser", [60]  )
    #print r.evaluate("ser", (60,) )
    #print r.evaluate("ser",  60   ) # this one doesn't work -- no surprise

    # Based off new (Oct 2006) NDFTs built from top500-angles Makefile
    # Remaining inaccuracies are due to dihedrals being rounded off to one decimal place!
    for aminoAcid, chiAngles, molpValue in [
    ("MET", [80.4, -172.2, 177.5], 17.9),
    ("GLN", [166.0, 178.0, -107.4], 17.6),
    ("ILE", [60.3, 162.4], 26.3),
    ("PHE", [-60.7, 97.9], 96.3),
    ("VAL", [-179.8], 60.3),
    ("LYS", [-175.6, 176.2, -172.0, -174.2], 84.8),
    ("THR", [76.7], 11.8),
    ("LEU", [-68.2, -165.8], 17.3),
    ("THR", [70.7], 29.1),
    ("LYS", [-179.3, -179.4, -151.1, -49.3], 35.4),
    ("THR", [-63.4], 68.6),
    ("ILE", [125.7, -175.4], 0.1),
    ("THR", [66.5], 44.9),
    ("LEU", [-117.8, 30.2], 0.2),
    ("GLU", [-75.1, -167.9, 139.8], 53.4),
    ("VAL", [-62.5], 35.6),
    ("GLU", [-73.9, -54.5, -18.4], 54.2),
    ("PRO", [-29.0], 84.9),
    ("SER", [35.7], 2.0),
    ("ASP", [-80.6, -19.8], 69.2),
    ("THR", [60.6], 82.6),
    ("ILE", [-60.9, -54.6], 39.2),
    ("GLU", [-169.6, -175.1, 72.8], 36.2),
    ("ASN", [177.5, 53.8], 39.0),
    ("VAL", [168.2], 47.7),
    ("LYS", [-71.7, -173.9, 179.2, 179.4], 96.6),
    ("LYS", [-60.8, 169.3, 148.9, -89.1], 18.8),
    ("ILE", [-70.9, 166.5], 75.3),
    ("GLN", [176.9, 171.9, 35.2], 57.2),
    ("ASP", [-150.1, 65.5], 2.7),
    ("LYS", [78.3, 138.2, 62.4, -165.4], 6.6),
    ("GLU", [-60.1, -76.8, -36.2], 58.3),
    ("ILE", [-54.4, 161.0], 36.6),
    ("PRO", [-31.6], 78.4),
    ("PRO", [-28.4], 87.4),
    ("ASP", [134.6, -61.7], 0.0),
    ("GLN", [-61.9, -179.0, -165.3], 24.9),
    ("GLN", [-53.1, -179.9, 28.0], 67.4),
    ("ARG", [161.7, 173.6, 174.2, -106.7], 40.7),
    ("LEU", [-68.3, 166.9], 78.6),
    ("ILE", [-48.9, -58.1], 30.2),
    ("PHE", [178.0, 78.2], 93.4),
    ("LYS", [-61.5, 173.7, -111.9, -58.8], 8.3),
    ("GLN", [-172.6, 177.3, 118.5], 18.4),
    ("LEU", [-50.8, -172.7], 25.8),
    ("GLU", [173.0, 141.4, 172.4], 6.9),
    ("ASP", [-78.0, 177.6], 76.2),
    ("ARG", [-55.9, -71.7, 114.2, -128.0], 0.1),
    ("THR", [59.8], 72.1),
    ("LEU", [-60.3, -179.2], 84.5),
    ("SER", [59.4], 78.7),
    ("ASP", [-73.0, 157.0], 90.7),
    ("TYR", [-63.3, 103.6], 92.2),
    ("ASN", [-159.1, -145.0], 2.2),
    ("ILE", [-69.6, 176.4], 68.0),
    ("GLN", [-79.4, -161.7, -148.4], 9.6),
    ("LYS", [49.7, 165.7, 154.3, 72.9], 15.3),
    ("GLU", [-72.2, 126.6, 36.7], 3.4),
    ("SER", [-73.2], 31.8),
    ("THR", [-60.6], 97.0),
    ("LEU", [-43.5, -170.9], 11.6),
    ("HIS", [-69.1, -88.8], 90.4),
    ("LEU", [172.8, 65.9], 50.8),
    ("VAL", [177.0], 91.3),
    ("LEU", [-108.1, 39.2], 0.4),
    ("ARG", [133.9, -155.8, 27.2, -152.9], 0.0),
    ("LEU", [-92.5, 37.5], 2.0),
    ("ARG", [-146.6, 157.6, 92.9, -95.5], 3.2),
    ]:
        print "%.1f =? %.1f  %s" % (100*r.evaluate(aminoAcid, chiAngles), molpValue, aminoAcid)

if (__name__ == "__main__"):
    run()
