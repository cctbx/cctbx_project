from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.filt_refls

import argparse as ap

parser = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
parser.add_argument("input", type=str, help="hopper output dir (simtbx.diffBragg.hopper) or a pandas dataframe (combined pickle files from hopper output folder)")
parser.add_argument("out", type=str, help="output pickle filename, which will then be a suitable input for diffBragg.geometry_refiner")
parser.add_argument("--thresh", type=float, default=10, help="MAD score threshold for outliers. Lower to remove more outliers")
parser.add_argument("--tag", type=str, default="filt", help="New refl tables will be written with this tag added to end of the filename")
parser.add_argument("--ers", type=str, default=None, help="if provided, write an exp-ref-spec file suitable input for simtbx.diffBragg.hopper")
parser.add_argument("--mind", type=float, default=None, help="if provided, filter shots whose median prediction offset is above this number (units=pixels)")
parser.add_argument("--minN", type=int, default=None, help="filter expts if number of refls is below this number")
parser.add_argument("--reflMaxD", type=float, default=None)

args = parser.parse_args()


import pandas
import glob
import h5py
from dials.array_family import flex
import numpy as np
from simtbx.diffBragg import utils

try:
    fnames = glob.glob(args.input + "/pandas/rank*/*.pkl")
    assert fnames
    comb_df = pandas.concat([pandas.read_pickle(f) for f in fnames])
except Exception:
    comb_df = pandas.read_pickle(args.input)
comb_df.reset_index(inplace=True, drop=True)

R2names = []

def get_dist_from_R(R):
    """ returns prediction offset, R is reflection table"""
    x,y,_ = R['xyzobs.px.value'].parts()
    x2,y2,_ = R['xyzcal.px'].parts()
    dist = np.sqrt((x-x2)**2 + (y-y2)**2)
    return dist

from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

# one ranks dataframe
df = np.array_split(comb_df, COMM.size)[COMM.rank]

keep = []
n = 0
n2 = 0
all_pred_off = []
for i in range(len(df)):
    #if i % COMM.size != COMM.rank:
    #    continue
    row = df.iloc[i]
    h = h5py.File(row.stage1_output_img, 'r')
    Rname = row.opt_exp_name.replace("/expers/", "/refls/").replace(".expt", ".refl")
    R = flex.reflection_table.from_file(Rname)
    all_dists = get_dist_from_R(R)
    bad_dists = np.zeros(len(R), bool)
    if args.reflMaxD is not None:
        bad_dists = np.array(all_dists) > args.reflMaxD

    d = -1
    if len(all_dists) > 0:
        d = np.median(all_dists)
    K = True
    if args.mind is not None:
        K =  d < args.mind
    keep.append(K)

    vals = h['sigmaZ_vals'][()]
    vals[np.isnan(vals)] = np.inf

    bad = list(np.where(utils.is_outlier(vals, args.thresh))[0])
    sel = [R[i_r]['h5_roi_idx'] not in bad for i_r in range(len(R))]

    sel = np.logical_and(sel, ~bad_dists)
    R2 = R.select(flex.bool(sel))
    d2 = -1
    if len(R2)> 0:
        d2 = np.median(get_dist_from_R(R2))
    if args.minN is not None and  len(R2) < args.minN:
        keep[-1] = False

    if d < -1 or d2 < -1:
        keep[-1] = False

    if args.tag is not None and K:
        R2name = Rname.replace(".refl", "_%s.refl" % args.tag)
        R2.as_file(R2name)
        if COMM.rank==0:
            print(i, len(df), "New refl table written=%s" % R2name)
    else:
        R2name = Rname.replace(".refl", "_KEEP=False.refl")
    R2names.append(R2name)
    all_pred_off.append(d2)
    if keep[-1]:
        n += len(R)
        n2 += len(R2)


df['filtered_refls'] = R2names
df['pred_offsets'] = all_pred_off
df = df.loc[keep]
rank_dfs = COMM.gather(df)
n = COMM.reduce(n)
n2 = COMM.reduce(n2)
if COMM.rank==0:
    df = pandas.concat(rank_dfs)
    df.reset_index(inplace=True, drop=True)
    df.to_pickle(args.out)
    print("\nSummary\n<><><><><>")
    print("%d expts tot" % len(df))
    print("Filtered refls have median pred offset=%.3f pix" % df.pred_offsets.median())
    print("Wrote %s which can be passed into diffBragg.geometry_refiner  input_pickle=%s" % (args.out, args.out))
    print("Kept %d / %d refls. Removed %.2f %% "
        % (n2, n, (n-n2)/float(n)*100. ))
    if args.ers is not None:
        with open(args.ers, "w") as ersFile:
            for e, r, s in df[["exp_name", "filtered_refls", "spectrum_filename"]].values:
                if s is None:
                    s = ""
                ersFile.write("%s %s %s\n" % (e,r,s))
        print("Wrote %s which can be passed into simtbx.diffBragg.hopper exp_ref_spec_file=%s" % (args.ers, args.ers))
    with open(args.out +".exectution.txt", "w") as o:
        import sys,os
        o.write("filt_refls was run from folder: %s\n" % os.getcwd())
        o.write("The command line input was:\n")
        o.write(" ".join(sys.argv))
