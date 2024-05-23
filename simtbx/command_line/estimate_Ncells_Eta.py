from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.estimate_Ncells_Eta

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("dirname", help="still process output folder", type=str, nargs="+")
parser.add_argument("--updatePhil", default=None, help="name of an exisiting stage 1 phil file  to update (just the init.Ncells portion)", type=str)
parser.add_argument("--expSuffix", help="extension of refined experiments (default: _refined.expt)", type=str, default="_refined.expt")
parser.add_argument("--thresh", type=float, default=7, help="MAD score for outliers (default=7 standard deviation above the median)")
parser.add_argument("--useMean", action="store_true", help="set Eta and Nabc using the mean (default is median)")
parser.add_argument("--NabcMax",  type=float, default=70, help="If estaimated Nabc is above this value, it will set to this value")
parser.add_argument("--NabcMin",  type=float, default=5, help="If estaimated Nabc is BELOW this value, it will be set to this value")
parser.add_argument("--EtaMax", type=float, default=0.5, help="If estimated Eta is above this range, it will be set to this value")
parser.add_argument("--EtaMin", type=float, default=1e-3, help="If estimated Eta is BELOW this range, it will be set to this value")

#parser.add_argument("--njobs", type=int, default=5, help="number of jobs (only runs on single node, no MPI)")
parser.add_argument("--plot", action="store_true", help="show a histogram at the end")
args = parser.parse_args()
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD
#from joblib import Parallel, delayed
import json
import numpy as np
from cctbx import uctbx
from scitbx.matrix import sqr
import os
import glob
from dxtbx.model import ExperimentList

fnames = []
for dirname in args.dirname:
    glob_s = os.path.join(dirname, "*%s" % args.expSuffix)
    fnames += glob.glob(glob_s)
if not fnames:
    if COMM.rank==0:
        print("no fnames")
    exit()

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
# template of the additional phil:
phil = """\ninit {{
  Nabc = [{n},{n},{n}]
  eta_abc = [{m},{m},{m}]
}}\n
"""

if COMM.rank==0:
    print("Obtained %d estimates ..." % len(all_Ns))
    import pandas
    import pylab as plt
    from simtbx.diffBragg import utils
    all_Ns = np.array(all_Ns)
    all_mos_spreads = np.array(all_mos_spreads)
    bad_Ns = utils.is_outlier(all_Ns, args.thresh)
    bad_mos_spreads = utils.is_outlier(all_mos_spreads, args.thresh)
    is_bad = np.logical_or(bad_Ns, bad_mos_spreads)
    print("Removing %d outlier estiamtes" % is_bad.sum())
    all_Ns = all_Ns[~is_bad]
    all_mos_spreads = all_mos_spreads[~is_bad]

    df = pandas.DataFrame({"Ncells": all_Ns, "mos_spread_deg": all_mos_spreads})
    print(df.Ncells.describe())
    print(df.mos_spread_deg.describe())
    if args.useMean:
        mean_N = df.Ncells.mean()
        mean_mos = df.mos_spread_deg.mean()
    else:
        mean_N = df.Ncells.median()
        mean_mos = df.mos_spread_deg.median()
    print("mean Ncells=%f" % mean_N)
    print("mean mos_spread=%f (deg.)" % mean_mos)

    if mean_mos > args.EtaMax or mean_mos < args.EtaMin:
        temp = mean_mos
        mean_mos = args.EtaMax if mean_mos > args.EtaMax else args.EtaMin
        print("Estimated Eta=%f, setting it to %f" % (temp, mean_mos))
    if mean_N > args.NabcMax or mean_N < args.NabcMin:
        temp = mean_N
        mean_N = args.NabcMax if mean_N > args.NabcMax else args.NabcMin
        print("Estimated N=%f, setting it to %f" %(temp, mean_N))

    phil = phil.format(n=round(mean_N,4), m=mean_mos)

    if args.updatePhil is not None:
        with open(args.updatePhil, "r") as o:
            s = o.read()
        s += phil
        with open(args.updatePhil, "w") as o:
            o.write(s)
    if args.plot:
        df.hist(bins=100, log=True)
        plt.show()
