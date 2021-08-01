# coding: utf-8
import glob
from joblib import Parallel,delayed

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--input", type=str, help="input pickle")
parser.add_argument("--glob", type=str, help="optional glob for input pickles", default=None)
parser.add_argument("--n", type=int, help="pixel cutoff", default=2)
parser.add_argument("--j", type=int, help="number of jobs", default=1)
parser.add_argument("--save", type=str, help="optional file name for saving figure output", default=None)
parser.add_argument("--signalcut", type=float, default=None, help="optional signal to backgroud cutoff")
parser.add_argument("--symbol", type=str, default="P1", help="space group symbol for mapping HKL list")
args = parser.parse_args()
import h5py
import numpy as np
import pandas
from scipy.ndimage.filters import gaussian_filter as GF
from scipy.spatial import distance, cKDTree
import os
import re
from simtbx.diffBragg import utils
from dials.array_family import flex
from scipy.ndimage import label, maximum_filter,center_of_mass, generate_binary_structure

def plot_overall_delta(df, savename=None):
    import pylab as plt
    all_ds = []
    all_ds_dials = []
    for name, X, Y, Xd, Yd in zip(df.exp_name, df.vec_x, df.vec_y, df.dials_vec_x, df.dials_vec_y):
        X = np.array(X)
        Y = np.array(Y)
        dist = np.sqrt(X ** 2 + Y ** 2)

        Xd = np.array(Xd)
        Yd = np.array(Yd)
        dist_d = np.sqrt(Xd ** 2 + Yd ** 2)
        dd = np.median(dist)
        all_ds.append(dd)
        dd2 = np.median(dist_d)
        all_ds_dials.append(dd2)
        print(name, "%d spots" % len(X),dd2-dd )
        if dd2 > dd:
            print("EXP %s shows improvement! (diff=%f)" % (name,dd2-dd) )
    ds = all_ds # [np.median(d) for d in all_ds]
    ds_dials = all_ds_dials #$[np.median(d) for d in all_ds_dials]
    diffs = np.sort(np.array(ds_dials) - np.array(ds))[::-1]
    pos = diffs >= 0
    x = np.arange(len(diffs))
    npos = sum(pos)
    frac_pos = 0 if npos==0 else npos / len(diffs) * 100.
    frac_neg = 100 - frac_pos
    plt.figure()
    plt.plot(x, diffs, 'kd', ms=2)
    plt.fill_between(x[pos], diffs[pos], color='C0', label="%.2f %%" % frac_pos)
    plt.fill_between(x[~pos], diffs[~pos], color='tomato', label="%.2f %%" % frac_neg)
    ax = plt.gca()
    plt.legend(prop={"size":11})
    ax.tick_params(labelsize=12)
    plt.xlabel("shot index (sorted by $\Delta_d$)", fontsize=14)
    plt.ylabel("$\Delta_d = d_{dials} - d_{diffBragg}$ (pixels)", fontsize=13)
    plt.grid(1, alpha=0.5)
    if savename is not None:
        plt.savefig(savename)
    else:
        plt.show()

    #res_bins = [30.0051,
    #            4.0169,
    #            3.1895,
    #            2.7867,
    #            2.5320,
    #            2.3506,
    #            2.2121,
    #            2.1013,
    #            2.0099,
    #            1.9325,
    #            1.8658,
    #            1.8075,
    #            1.7558,
    #            1.7096,
    #            1.6679,
    #            1.6300]

#def get_centroid(img):
#    bs = generate_binary_structure(2,3)
#    lab, nlab = label(img > np.percentile(img,90), bs)
#    # take the central label if 2 exist
#    if nlab > 1:
#        y, x = img.shape
#        coms = [center_of_mass(img, lab, i+1) for i in range(1, nlab)]
#        cent = y*.5, x*.5
#        com = coms[np.argmin([distance.euclidean(xy, cent) for xy in coms])]
#        nlab = 1
#    else:
#        com = center_of_mass(img, lab, 1)
#    return com, nlab
#
#def get_centroid2(img):
#    peakmask = maximum_filter(img, size=2)

print("pandas input")
if args.glob is not None:
    fnames = glob.glob(args.glob+"/pandas/rank*/*pkl")
    print("found %d files in glob" % len(fnames))
    df = pandas.concat([pandas.read_pickle(f) for f in fnames])
else:
    df = pandas.read_pickle(args.input)

#def get_strong_refl(name):
#    s = re.search("run[0-9]+_shot[0-9]+",name)
#    rs = name[s.start():s.end()]
#    strong_name = "/global/cfs/cdirs/m3562/der/indexed/%s_indexed.refl" % rs
#    return strong_name


if "stage1_output_img" in list(df):
    df["imgs"] = df.stage1_output_img

if "imgs" not in list(df):
    df['imgs'] = [f.replace("expers", "imgs").replace(".expt", ".h5") for f in df.opt_exp_name]

# NOTE new way
#df["refl_names"] = ["/global/cfs/cdirs/m3562/der/indexed/%s_indexed.refl" % utils.get_rs(f) for f in df.exp_name]
df["refl_names"] = df.stage1_refls #["/global/cfs/cdirs/m3562/der/indexed/%s_indexed.refl" % utils.get_rs(f) for f in df.exp_name]
assert os.path.exists(df.refl_names.values[0])

#if "stage1_refls" in list(df):
#    print("NEWWAY refls")
#    df['refl_names'] = df.stage1_refls
#
#elif "refl_names" not in list(df):
#    df['refl_names'] = [f.replace(".expt", ".refl") for f in df.exp_name]
#    #df['refl_names'] = [get_strong_refl(f) for f in df.exp_name]
#    #refls = [f.replace("_pathmod.expt", "_idx.refl") for f in df.exp_name]
#    #df['refl_names'] = [f.replace(".expt", "_indexed2.refl") for f in df.exp_name]
#    #df['refl_names'] = [f.replace(".expt", "_indexed3.refl") for f in df.exp_name]
#    #from IPython import embed
#    #embed()
#    #df['refl_names'] = [f.replace(".expt", "_expanded.refl") for f in df.exp_name]

def main(jid):
    dev_res =[]
    per_img_dists = []
    per_img_dials_dists =[]
    per_img_Z = []
    per_img_Z2 = []
    per_img_signal = []
    per_img_ref_index =[]
    per_img_vec_dists = []
    per_img_dials_vec_dists =[]
    img_names =[]
    per_img_shot_roi = []
    per_img_hkl = []
    Nno_sig = 0
    Nposs =0
    Nroi = 0
    for ii, (imgf,reff) in  enumerate(df[["imgs","refl_names"]].values):
        if ii % args.j != jid:
            continue
        if jid==0:
            print("Processing %d / %d" % (ii+1, len(df)))
        h = h5py.File(imgf, "r")
        dat = h['data']
        #if 'bragg' in list(h.keys()):
        mod = h['bragg']
        #else:
        mod_with_bg = h['model']
        nroi = len(dat.keys())
        R = flex.reflection_table.from_file(reff)
        R["refl_index"] = flex.int(range(len(R)))
        Rpp = utils.refls_by_panelname(R)
        pids = h['pids']
        rois = h['rois']
        all_dists = []
        all_dials_dists = []
        all_signal = []
        all_ref_index=[]
        all_vec_dists =[]
        all_dials_vec_dists = []
        all_Z =[]
        all_Z2 =[]
        all_hkl = []
        all_shot_roi = []
        sigma_rdout = h["sigma_rdout"][()]
        Nroi+= nroi
        Nposs += len(R)
        for img_i_roi in range(nroi):
            pid = pids[img_i_roi]
            if pid not in Rpp:
                continue

            ddd = dat["roi%d" % img_i_roi][()]
            mmm = mod["roi%d" % img_i_roi][()]
            mmm_with_bg = mod_with_bg["roi%d" % img_i_roi][()]
            if np.all(mmm==0):
                #print("Empty spot! continue")
                Nno_sig += 1
                continue
            #signal = mmm.max() / np.median(mmm)
            signal = mmm.max() #/ np.mean(mmm)
            if args.signalcut is not None and signal < args.signalcut:
                print("signal cut!")
                continue
            #noise = np.sqrt(ddd + sigma_rdout**2)
            noise = np.sqrt(mmm_with_bg + sigma_rdout**2)
            Z = (ddd-mmm_with_bg) / noise
            #Z2 = (ddd-mmm) / noise2
            roi_d = GF(ddd, 1) 
            roi_m = GF(mmm, 0 )
            #com_d, nlab_d = get_centroid(roi_d)
            #com_m, nlab_m = get_centroid(roi_m)


            ref_p = Rpp[pid]
            x1,x2,y1,y2 = rois[img_i_roi]

            Ycoor, Xcoor = np.indices(roi_m.shape)
            I = roi_m.ravel()
            Isum = I.sum()
            Ix = (I*Xcoor.ravel()).sum() / Isum
            Iy = (I*Ycoor.ravel()).sum() / Isum
            com_m = Iy, Ix
            nlab_m = 1




            com_d = .5*(y1+y2),.5*(x1+x2)
            if np.any(np.isnan(com_m)):#or np.any(np.isnan(com_d)):
                print("nan in com")
                continue
            if nlab_m != 1:
                print("multiple lab model!")
                #from IPython import embed
                #embed();exit()
                continue
            #if nlab_d != 1:
            #    print("multiple lab in data!")
            #    #from IPython import embed
            #    #embed();exit()
            #    continue
            xyz = np.array(ref_p["xyzobs.px.value"])
            xy = xyz[:,:2]
            tree = cKDTree(xy)
            y_com,x_com = com_d
            #y_com += y1+0.5
            #x_com += x1+0.5
            y_com +=0.5
            x_com+= 0.5
            y_nelder, x_nelder = com_m
            y_nelder += y1+0.5
            x_nelder += x1+0.5
            res = tree.query_ball_point((x_com, y_com), r=7)
            if len(res) == 0:
                #print("weird tree query multiple or no res")
                continue
            dists = [distance.euclidean(np.array(ref_p[i_r]['xyzobs.px.value'])[:2], (x_com,y_com)) for i_r in res]
            close = np.argmin(dists)
            i_roi = res[close]
            r = ref_p[i_roi]
            xcal,ycal,_ = r['xyzcal.px']
            xobs,yobs,_ = r['xyzobs.px.value']
            rlp = r['rlp']
            hkl = r['miller_index']
            reso = 1 / np.linalg.norm(rlp)
            refl_idx = r['refl_index']
            #d_dials = distance.euclidean((x_com, y_com), (xcal, ycal))
            d_dials = distance.euclidean((xobs, yobs), (xcal, ycal))
            d = distance.euclidean((xobs, yobs), (x_nelder, y_nelder))
            
            # output
            dev_res.append( (d, d_dials,imgf, i_roi, nroi))

            # distances
            all_dists.append(d)
            all_dials_dists.append(d_dials)
            # Z-score sigmas
            all_Z.append(np.std(Z))
            all_Z2.append(np.mean(Z))
            
            # vector differences
            vec_dials_d = np.array((xobs, yobs)) - np.array((xcal, ycal))
            vec_d = np.array((xobs, yobs)) - np.array((x_nelder, y_nelder))
            all_vec_dists.append(vec_d)
            all_dials_vec_dists.append(vec_dials_d)
            all_shot_roi.append(refl_idx)  # NOTE hijacking container for refls idx
            assert r['panel'] == pid
            all_ref_index.append(pid)  # NOTE hijacking refl index to be pid
            #all_signal.append(signal)
            all_signal.append(reso) # NOTE hijacking this container, rename to reso later
            all_hkl.append(hkl)

        img_names.append(imgf)
        per_img_hkl.append(tuple(utils.map_hkl_list(all_hkl, True, args.symbol)))
        per_img_ref_index.append(tuple(all_ref_index))
        per_img_signal.append(tuple(all_signal))
        per_img_dists.append(tuple(all_dists))
        per_img_dials_dists.append(tuple(all_dials_dists))
        per_img_vec_dists.append(tuple(all_vec_dists))
        per_img_dials_vec_dists.append(tuple(all_dials_vec_dists))
        per_img_Z.append(tuple(all_Z))
        per_img_Z2.append(tuple(all_Z2))
        per_img_shot_roi.append( all_shot_roi)

    return dev_res, per_img_dists, per_img_dials_dists, per_img_Z, per_img_Z2, per_img_vec_dists\
                ,per_img_dials_vec_dists, per_img_shot_roi, per_img_signal, img_names, per_img_ref_index,\
           (Nno_sig, Nposs, Nroi), per_img_hkl

results = Parallel(n_jobs=args.j)(delayed(main)(j) for j in range(args.j))

dev_res =[]
per_img_dists = []
per_img_dials_dists =[]
per_img_Z = []
per_img_Z2 = []
per_img_vec_dists =[]
per_img_signal = []
per_img_ref_index =[]
per_img_dials_vec_dists =[]
per_img_shot_roi = []
per_img_hkl = []
img_names =[]
n1=n2=n3 = 0
for r in results:
    dev_res +=r[0]
    per_img_dists += r[1]
    per_img_dials_dists +=r[2]
    per_img_Z += r[3]
    per_img_Z2 +=r[4]
    per_img_vec_dists += r[5]
    per_img_dials_vec_dists += r[6]
    per_img_shot_roi += r[7]
    per_img_signal += r[8]
    img_names += r[9]
    per_img_ref_index += r[10]
    nno_sig, nposs, nroi = r[11]
    per_img_hkl += r[12]
    n1 += nno_sig
    n2 += nposs
    n3 += nroi

df_process = pandas.DataFrame(
    {"sigmaZ_PoissonDat" : per_img_Z,
    "Z_PoissonMod" : per_img_Z2,
   # "signal_to_background" : per_img_signal,
    "resolution" : per_img_signal,  # NOTE hijacked the signal container to store reso, rename later
    "pred_offsets" : per_img_dists,
    #"refl_index": per_img_ref_index,
     "panel": per_img_ref_index,
    "pred_offsets_dials" : per_img_dials_dists,
     "hkl": per_img_hkl,
    "refls_idx" : per_img_shot_roi,"imgs": img_names})

df = pandas.merge(df, df_process, on="imgs", how='inner')
# NOTE uncomment for LS49
#df['tstamp'] = [ls49_utils.get_tstamp(f) for f in df.imgs]

vecx = [ tuple([i for i,j in v]) for v in per_img_vec_dists]
vecy = [ tuple([j for i,j in v]) for v in per_img_vec_dists]
df["vec_x"] = vecx
df["vec_y"] = vecy
dials_vecx = [ tuple([i for i,j in v]) for v in per_img_dials_vec_dists]
dials_vecy = [ tuple([j for i,j in v]) for v in per_img_dials_vec_dists]
df["dials_vec_x"] = dials_vecx
df["dials_vec_y"] = dials_vecy
#===================================================================

d,d2,imgf,_,_ = zip(*dev_res)
from pylab import *
bins = linspace(min(d+d2), max(d+d2), 100)
hist(d, bins=bins, histtype='step', lw=2, label="after stage1", color="C0")
hist(d2, bins=bins, histtype='step', lw=2, label="from dials", color="tomato")
d = np.array(d)
d2 = np.array(d2)
good_d = d[d < args.n]
before_good_d = d2[ d < args.n]

from itertools import groupby
dimg = list(zip(d, imgf))
gb = groupby(sorted(dimg, key=lambda x: x[1]), key=lambda x:x[1])

gb_results = {k:list(v) for k,v in gb}

#better_imgf = [v for sublist in [[ss[1] for ss in s] for s in stuff if len(s) > 1] for v in sublist]

better_d = []
better_imgf =[]
imgnames = list(gb_results.keys())
med_dists = []
for name in imgnames:
    img_dists = [v[0] for v in gb_results[name]]
    n_within_2 = sum([dist < args.n for dist in img_dists])
    if n_within_2 >= 2:
        better_d += img_dists
        better_imgf.append(name)
    med_dists.append(np.median(img_dists))

df["good_img"] = df.imgs.isin(better_imgf)
df_better = df.loc[df.imgs.isin(better_imgf)]
df_better.reset_index(inplace=True, drop=True)
df_better["predictions"] = df_better.refl_names
#df_better.to_pickle(args.input.replace(".pkl", "_stg2.pkl"))

nwithin = len(good_d)
npx = np.median(good_d)
npx_before = np.median(before_good_d)
print("good", np.median(good_d))
print("Before good", np.median(before_good_d))
hist(good_d, bins=bins, histtype='stepfilled', lw=2, label="%d spots within %d pix" % (nwithin,args.n), alpha=0.5, color="C0")
hist(before_good_d, bins=bins, histtype='stepfilled', lw=2, label="%d spots prior to refinement" % nwithin, alpha=0.5, color="tomato")
#hist(better_d, bins=bins, histtype='stepfilled', lw=2, label="%d spots on images containing 2 or more \npredictions within %.1f pix of observations (%d images)" % (len(better_d),args.n, len(set(better_imgf))), alpha=0.5, color="k")
legend()
xlabel("pixels")
xlabel("|xobs - xcal| pixels")
xlabel("|xobs - xcal| pixels", fontsize=12)
title("xcal and xobs for %d refls (med %.2f px from %.2f px)\n (%s)\n%d NoSig,%d GroupA,%d Stage1Input" % (len(d), npx, npx_before, args.input if args.input is not None else args.glob,
                                                                      n1,n2,n3), fontsize=12)
ax = gca()
ax.tick_params(labelsize=12)
print("Number of modeled refls within %d pix of obs: %d " % (args.n,nwithin))
if args.save is not None:
    plt.savefig(args.save)
    pkl_name = os.path.splitext(args.save)[0] + "_pandas.pkl"
    df.to_pickle(pkl_name)
    save_name = os.path.splitext(args.save)[0] + "_overallDelta.png"
    plot_overall_delta(df,save_name)
else:
    plot_overall_delta(df)
    #show()



def make_plot(d, d2, imgf, n=2):
    bins = linspace(min(d + d2), max(d + d2), 100)
    hist(d, bins=bins, histtype='step', lw=2, label="after basin-hopping", color="C0")
    hist(d2, bins=bins, histtype='step', lw=2, label="from dials", color="tomato")
    d = np.array(d)
    d2 = np.array(d2)
    good_d = d[d < n]
    before_good_d = d2[d < n]

    from itertools import groupby
    dimg = list(zip(d, imgf))
    gb = groupby(sorted(dimg, key=lambda x: x[1]), key=lambda x: x[1])

    gb_results = {k: list(v) for k, v in gb}

    # better_imgf = [v for sublist in [[ss[1] for ss in s] for s in stuff if len(s) > 1] for v in sublist]

    better_d = []
    better_imgf = []
    imgnames = list(gb_results.keys())
    med_dists = []
    for name in imgnames:
        img_dists = [v[0] for v in gb_results[name]]
        n_within_2 = sum([dist < n for dist in img_dists])
        if n_within_2 >= 2:
            better_d += img_dists
            better_imgf.append(name)
        med_dists.append(np.median(img_dists))

    nwithin = len(good_d)
    hist(good_d, bins=bins, histtype='stepfilled', lw=2, label="%d spots within %d pix" % (nwithin, n), alpha=0.5,
         color="C0")
    hist(before_good_d, bins=bins, histtype='stepfilled', lw=2, label="%d spots prior to refinement" % nwithin,
         alpha=0.5, color="tomato")
    #hist(better_d, bins=bins, histtype='stepfilled', lw=2,
    #     label="%d spots on images containing 2 or more \npredictions within %.1f pix of observations (%d images)" % (
    #     len(better_d), n, len(set(better_imgf))), alpha=0.5, color="k")
    legend()
    xlabel("pixels")
    xlabel("|xobs - xcal| pixels")
    xlabel("|xobs - xcal| pixels", fontsize=12)
    title("xcal and xobs for %d refls\n" % (len(d)), fontsize=12)
    ax = gca()
    ax.tick_params(labelsize=12)
    show()

