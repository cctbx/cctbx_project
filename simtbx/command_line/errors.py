from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.errors

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("phil", type=str, help="the phil file used for per spot refinement")
parser.add_argument("indirs", type=str, nargs="+", help="hopper output folders with expers inside")
parser.add_argument("outdir", type=str, help="output folder where integration tables will be saved")
parser.add_argument("--ndev", type=int, default=1, help="number of gpu devices")
args = parser.parse_args()

from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

import os
from simtbx.diffBragg.parameter_errors import get_errors
import glob

if COMM.rank==0:
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

fnames = []
for dirname in args.indirs:
    fnames += glob.glob(dirname + "/expers/rank*/*expt")
if COMM.rank==0:
    print("Found %d expts total" %len(fnames))

failed_assert = []
devid = COMM.rank % args.ndev

for i,exp_f in enumerate(fnames):

    if i% COMM.size != COMM.rank:
            continue

    ref_f = exp_f.replace("expers/", "refls/").replace(".expt", ".refl")
    assert os.path.exists(ref_f)
    pkl_f = exp_f.replace("expers/", "pandas/").replace(".expt", ".pkl")
    assert os.path.exists(pkl_f)
    out_f = os.path.join(args.outdir + "/shot%d"% i)
    print("expt %s, shot %d on GPU device %d" % (out_f,i, devid  ))

    try:
        get_errors(args.phil, exp_f, ref_f, pkl_f, out_f, verbose=False, devid=devid)
    except AssertionError:
        failed_assert.append(exp_f)

failed_assert = COMM.reduce(failed_assert)
if COMM.rank==0:
        print("the following expts failed(!):", failed_assert)
