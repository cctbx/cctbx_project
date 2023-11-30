from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("dirname", help="still process output folder", type=str)
parser.add_argument("--updatePhil", default=None, help="name of an exisiting stage 1 phil file  to update (just the init.Ncells portion)", type=str)
parser.add_argument("--expSuffix", help="extension of refined experiments", type=str, default="_refined.expt")
#parser.add_argument("--njobs", type=int, default=5, help="number of jobs (only runs on single node, no MPI)")
parser.add_argument("--plot", action="store_true", help="show a histogram at the end")
args = parser.parse_args()
# LIBTBX_SET_DISPATCHER_NAME diffBragg.estimate_Ncells_Eta
from mpi4py import MPI
COMM = MPI.COMM_WORLD
#from joblib import Parallel, delayed
import json
import numpy as np
from cctbx import uctbx
from scitbx.matrix import sqr
import os
import glob
from dxtbx.model import ExperimentList

glob_s = os.path.join(args.dirname, "*%s" % args.expSuffix)
fnames = glob.glob(glob_s)

#def main(jid):
all_Ns = []
all_mos_spreads = []
for i, f in enumerate(fnames):
    if i % COMM.size != COMM.rank:
        continue
    print(f)
    #Cs = ExperimentList.from_file(f, False).crystals()
    Cs = json.load(open(f, 'r'))['crystal']
    dom_sizes = np.array([C['ML_domain_size_ang'] for C in Cs])
    mos_spreads = [2*C['ML_half_mosaicity_deg'] for C in Cs]
    uc_vols = []
    for C in Cs:
        a = C['real_space_a']
        b = C['real_space_b']
        c = C['real_space_c']
        uc = uctbx.unit_cell(orthogonalization_matrix=sqr(a + b + c).transpose())
        uc_vols.append(uc.volume())

    Ns = dom_sizes / np.power(uc_vols, 1 / 3.)
    all_Ns += list(Ns)
    all_mos_spreads += list(mos_spreads)
#    return all_Ns, all_mos_spreads

all_Ns = COMM.reduce(all_Ns)
all_mos_spreads = COMM.reduce(all_mos_spreads)
#results = Parallel(n_jobs=args.njobs)(delayed(main)(j) for j in range(args.njobs))
#all_Ns = []
#all_mos_spreads = []
#for N,mos in results:
#    all_Ns += N
#    all_mos_spreads += mos

if COMM.rank==0:
    import pandas
    import pylab as plt
    df = pandas.DataFrame({"Ncells": all_Ns, "mos_spread_deg": all_mos_spreads})
    print(df.Ncells.describe())
    print(df.mos_spread_deg.describe())
    mean_N = df.Ncells.mean()
    mean_mos = df.mos_spread_deg.mean()
    print("mean Ncells=%f" % mean_N)
    print("mean mos_spread=%f (deg.)" % mean_mos)

    # mosaic spread should be at least 1e-5 if we want to try refining it
    mean_mos = max(mean_mos, 1e-5)

    phil = """\ninit {{
      Nabc = [{n},{n},{n}]
      eta_abc = [{m},{m},{m}]
    }}\n""".format(n=round(mean_N,4), m=mean_mos)

    if args.updatePhil is not None:
        with open(args.updatePhil, "r+") as o:
            s = o.read()
            s += phil
            o.write(s)
    if args.plot:
        df.hist(bins=100, log=True)
        plt.show()
