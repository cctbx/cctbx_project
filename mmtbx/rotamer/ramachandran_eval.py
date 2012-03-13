import mmtbx.rotamer
from mmtbx.rotamer.n_dim_table import NDimTable
from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
from mmtbx.rotamer.rotamer_eval import open_rotarama_dlite
from libtbx import easy_pickle
from libtbx.utils import Sorry
import weakref
import sys, os

# maps programatic name to file name
aminoAcids = {
    'general' : 'general',
    'glycine' : 'gly-sym',
    'proline' : 'pro',
    'prepro' : 'prepro',
}
aminoAcids_8000 = {
    'general' : 'general-noGPIVpreP',
    'glycine' : 'gly-sym',
    'cis-proline' : 'cispro',
    'trans-proline' : 'transpro',
    'pre-proline' : 'prepro-noGP',
    'isoleucine or valine' : 'ileval-nopreP',
}

class RamachandranEval:

    # This is shared among all instances of RamachandranEval -- a class variable.
    # It holds a LOT of read-only data, so this helps save memory.
    aaTables = {} # maps "his" to a NDimTable object for histidine, etc.

    def __init__(self):
        main_aaTables = RamachandranEval.aaTables
        self.aaTables = {}
        for aa,ndt_weakref in main_aaTables.items():
            # convert existing weak references to strong references
            self.aaTables[aa] = ndt_weakref()
        rama_data_dir = find_rotarama_data_dir()
        target_db = open_rotarama_dlite(rotarama_data_dir=rama_data_dir)
        for aa, aafile in aminoAcids_8000.items():
                if (self.aaTables.get(aa) is not None): continue
                data_file = "rama8000-"+aafile+".data"
                pickle_file = "rama8000-"+aafile+".pickle"
                pair_info = target_db.pair_info(
                  source_path=data_file,
                  target_path=pickle_file,
                  path_prefix=rama_data_dir)
                if pair_info.needs_update:
                    raise Sorry(
                        "chem_data/rotarama_data/*.pickle files are missing or out of date.\n"
                        "  Please run\n"
                        "    mmtbx.rebuild_rotarama_cache\n"
                        "  to resolve this problem.\n")
                ndt = easy_pickle.load(file_name=os.path.join(
                  rama_data_dir, pair_info.target.path))
                self.aaTables[aa] = ndt
                main_aaTables[aa] = weakref.ref(ndt)

    def evaluate(self, aaName, phiPsi):
        '''Evaluates the protein backbone conformation from 0.0 (worst) to 1.0 (best).

        Values below 0.0005 are considered outliers for the general case,
        or values below 0.002 for glycine, proline, and pre-proline.
        The field aaName should be one of "general", "glycine", "proline", or "prepro".
        If the aaName is not recognized, returns None.'''
        ndt = self.aaTables.get(aaName.lower())
        if (ndt is None): return None
        return ndt.valueAt(phiPsi)

    def evaluate_sites (self, aaName, phi_psi_i_seqs, sites_cart) :
      assert (aaName in ["general",
                         "glycine",
                         "cis-proline",
                         "trans-proline",
                         "pre-proline",
                         "isoleucine or valine"])
      (phi, psi) = mmtbx.rotamer.phi_psi_from_sites(
        i_seqs=phi_psi_i_seqs,
        sites_cart=sites_cart)
      return self.evaluate(aaName, (phi,psi))

def exercise(args):
  if (find_rotarama_data_dir(optional=True) is None):
    print "Skipping exercise(): rotarama_data directory not available"
  else:
    from mmtbx.command_line import rebuild_rotarama_cache
    rebuild_rotarama_cache.run()
    #
    from libtbx.test_utils import approx_equal
    #
    verbose = ("--verbose" in args)
    #
    r = RamachandranEval()
    tbl = r.aaTables['glycine']
    assert RamachandranEval().aaTables['glycine'] is tbl
    #
    # Based off new (Oct 2006) NDFTs built from top500-angles Makefile
    # Remaining inaccuracies are due to dihedrals being rounded off to
    # one decimal place!
    for aminoAcid, phiPsi, molpValue in [ #{{{
        ("general", [-91.17, 131.33], 40.94),
        ("general", [-139.95, 81.92], 1.47),
        ("general", [51.22, 23.16], 2.52),
        ("glycine", [-151.55, 88.38], 0.07),
        ("general", [-82.35, 160.71], 23.71),
        ("general", [-101.76, 135.67], 45.14),
        ("general", [-110.13, 120.96], 58.46),
        ("general", [-81.78, 128.43], 39.44),
        ("glycine", [-133.39, 160.58], 21.73),
        ("general", [-133.33, 127.68], 42.30),
        ("general", [-108.38, 126.39], 61.01),
        ("general", [-113.08, 126.63], 62.09),
        ("general", [-143.87, 166.78], 25.67),
        ("prepro", [-124.82, 120.58], 26.77),
        ("proline", [-64.55, 145.64], 94.90),
        ("general", [-45.79, -40.17], 10.23),
        ("general", [-63.09, 7.01], 0.47),
        ("general", [-100.14, -45.65], 7.66),
        ("general", [-86.64, -15.77], 36.67),
        ("general", [-154.47, 93.16], 1.42),
        ("general", [-112.91, 154.29], 27.39),
        ("prepro", [-164.30, 153.09], 10.20),
        ("proline", [-40.03, -62.10], 1.14),
        ("general", [-88.23, 54.64], 2.88),
        ("glycine", [-150.33, 8.77], 0.99),
        ("general", [-106.19, 65.65], 0.93),
        ("general", [-79.57, 41.90], 0.63),
        ("general", [-11.57, -82.78], 0.00),
        ("general", [54.67, 49.83], 15.13),
        ("general", [-31.21, 156.35], 0.00),
        ("general", [-5.62, -53.35], 0.00),
        ("general", [-60.03, -22.48], 63.20),
        ("general", [-71.41, -17.83], 64.40),
        ("general", [-80.18, -5.88], 44.61),
        ("general", [-90.33, -159.55], 0.66),
        ("glycine", [64.34, -132.09], 34.98),
        ("general", [-107.31, 27.71], 8.67),
        ("general", [-69.18, -45.50], 74.08),
        ("general", [-131.17, 129.92], 53.87),
        ("general", [-125.57, 136.49], 60.53),
        ("general", [-118.28, 151.90], 34.88),
        ("general", [-136.27, 122.13], 24.86),
        ("general", [-72.52, 148.04], 44.34),
        ("general", [-123.41, 135.11], 61.62),
        ("general", [-81.31, 137.39], 36.95),
        ("general", [-94.11, -22.71], 18.73),
        ("general", [-104.08, -30.86], 8.84),
        ("general", [-106.27, 124.65], 59.71),
        ("general", [-113.44, 157.11], 24.96),
        ("general", [-80.13, -171.26], 2.60),
        ("general", [-67.90, -10.85], 53.15),
        ("general", [-98.82, -6.37], 26.15),
        ("glycine", [92.77, -6.55], 73.68),
        ("general", [-88.65, 139.40], 33.51),
        ("glycine", [-53.72, 145.75], 24.34),
        ("general", [-69.26, 144.42], 48.55),
        ("general", [-107.68, 133.12], 57.79),
        ("glycine", [-77.98, -161.75], 22.90),
        ("general", [-88.43, 159.48], 18.09),
        ("general", [-69.18, 118.06], 13.32),
        ("general", [-118.20, 130.15], 63.55),
        ("general", [-101.63, 154.20], 20.70),
        ("general", [-60.43, 139.20], 53.08),
        ("general", [-98.65, -39.53], 9.62),
        ("glycine", [-170.85, 172.86], 45.91),
        ("general", [-83.58, 114.39], 22.07),
        ("general", [-86.93, 127.35], 36.67),
        ("general", [-133.54, 58.50], 1.43),
        ("general", [-104.10, -54.11], 3.80),
        ("prepro", [-95.84, 135.89], 33.59),
        ("proline", [-61.91, -44.54], 28.37),
        ("general", [-62.64, 83.24], 0.07),
        ("general", [-171.76, 111.15], 0.41),
        ("glycine", [-136.70, -139.20], 3.57),
        ("general", [42.42, 42.27], 3.32),
        ("glycine", [55.77, 28.48], 55.42),
        ("general", [-90.29, 157.93], 16.97),
        ("general", [-135.18, 139.86], 47.95),
        ("general", [-102.72, 137.51], 41.91),
        ("general", [-127.29, 134.66], 60.92),
        ("general", [-83.54, 142.97], 31.55),
        ("general", [-155.14, 146.33], 20.07),
        ("general", [-128.99, 131.88], 59.20),
        ("general", [-107.06, 146.26], 31.12),
        ("general", [-92.24, 97.79], 11.98),
        ("general", [-134.18, 101.62], 6.01),
        ("general", [-78.17, 152.18], 31.34),
        ("glycine", [-50.36, -27.45], 22.99),
        ("general", [-66.02, -20.32], 67.44),
        ("general", [-77.76, -2.10], 27.17),
        ("glycine", [81.90, -156.46], 41.39),
        ("prepro", [-107.81, 128.95], 28.69),
        ("proline", [-70.18, 133.33], 26.27),
        ("glycine", [-126.36, 3.87], 6.77),
        ("general", [-167.92, 174.97], 6.79),
        ("general", [-159.77, 150.09], 17.32),
        ("general", [-99.29, 152.73], 20.45),
        ("general", [-145.08, 134.33], 21.62),
        ("general", [-125.58, 126.16], 58.63),
        ("general", [-100.31, 110.18], 27.07),
        ("general", [-88.81, -7.47], 49.70),
        ("general", [-75.67, 160.16], 28.16),
    ]: #}}}
      r_eval = 100*r.evaluate(aminoAcid, phiPsi)
      if (verbose):
        print aminoAcid, "%4.1f %4.1f %4.1f" % (
          r_eval, molpValue, r_eval-molpValue)
      assert approx_equal(r_eval, molpValue, eps=0.9)
    #
    # check if tables are cleared from memory if all RamachandranEval instances
    # are gone
    for aa,ndt_weakref in RamachandranEval.aaTables.items():
      assert ndt_weakref() is not None
    del r
    del tbl
    for aa,ndt_weakref in RamachandranEval.aaTables.items():
      assert ndt_weakref() is None
    #
  print "OK"

if (__name__ == "__main__"):
    exercise(sys.argv[1:])
