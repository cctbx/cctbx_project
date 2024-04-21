from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.pred_offsets



def main(fnames, skip_empty=True):
    from dials.array_family import flex
    import numpy as np
    """fnames is a list of stage1 reflection table outputs (if hopper was run using debug_mode=True)"""
    all_d = []
    all_shotd = []
    nref_per_shot = []
    for i_f,f in enumerate(fnames):
        R = flex.reflection_table.from_file(f)
        if len(R)==0:
            if not skip_empty:
                all_d.append(None)
                all_shotd.append(np.inf)
                nref_per_shot.append(0)
            continue
        x,y,_=R['xyzobs.px.value'].parts()
        x2,y2,_=R['xyzcal.px'].parts()
        assert x2 > 0
        assert y2 > 0
        d = np.sqrt( (x-x2)**2 + (y-y2)**2)
        med_d = np.median(d)
        print('(%d/%d) diffBragg: %.3f (%s) Nref=%d'
              % (i_f+1, len(fnames), med_d, f, len(d)))
        all_d.append(d)
        all_shotd.append(med_d)
        nref_per_shot .append( len(d))

    all_d = np.hstack([val for val in all_d if val is not None])
    print("median over %d shots=%f pixels (%d refls)" % (len(nref_per_shot), np.median(all_d), len(all_d)))
    print("Min refls per shot=%d, max refls per shot = %d, ave refls per shot=%.1f"
          % (np.min(nref_per_shot), np.max(nref_per_shot), np.mean(nref_per_shot)))
    return all_d, all_shotd, nref_per_shot


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("globs", nargs="+", help="input globs  (in quotes) for groups of reflection tables.",
                        type=str)
    parser.add_argument("--noPlot", action="store_true", help="dont display graphic")
    parser.add_argument("--linear", action="store_true", help="use linear histogram")
    parser.add_argument("--nbins", type=int, nargs=2, default=[200, 40])
    args = parser.parse_args()

    import glob
    import pylab as plt
    plt.style.use("ggplot")

    for glob_s in args.globs:
        fnames = glob.glob(glob_s)
        all_d, all_shotd, nref_per_shot = main(fnames)
        plt.subplot(121)
        plt.hist( all_d, bins=args.nbins[0], histtype='step', lw=2, label=glob_s)
        if not args.linear:
            plt.gca().set_yscale("log")
        plt.xlabel("|calc-obs| (pixels)", fontsize=17)
        plt.ylabel("num refls", fontsize=15)
        #plt.legend(prop={'size':12})
        plt.gca().tick_params(direction='in', which='both', labelsize=15)
        plt.subplot(122)
        plt.hist( all_shotd, bins=args.nbins[1], histtype='step', lw=2, label=glob_s)
        if not args.linear:
            plt.gca().set_yscale("log")
        plt.xlabel("median |calc-obs| (pixels)", fontsize=17)
        plt.ylabel("num shots", labelpad=0, fontsize=17)
        #legend(prop={'size':11})
        plt.gca().tick_params(direction='in', which='both', labelsize=15)

    plt.gcf().set_size_inches((7,4))
    #plt.suptitle(" %d shots ; %d refls" % (len(nref_per_shot), len(all_d)))
    plt.subplots_adjust(bottom=0.2, right=0.97, top=0.92, left=0.1)

    if args.noPlot:
        exit()

    plt.show()
