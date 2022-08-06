from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME diffBragg.pred_offsets

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("globs", nargs="+", help="input globs  (in quotes) for groups of reflection tables.",
                    type=str)
parser.add_argument("--noPlot", action="store_true", help = "dont display graphic")
parser.add_argument("--linear", action="store_true", help="use linear histogram")
parser.add_argument("--nbins", type=int, nargs=2, default=[200,40])
args = parser.parse_args()

from pylab import *
from dials.array_family import flex
import glob


style.use("ggplot")
for glob_s in args.globs:
    fnames = glob.glob(glob_s)
    all_d = []
    all_shotd = []
    nref_per_shot = []
    for f in fnames:
        R = flex.reflection_table.from_file(f)
        if len(R)==0:
            continue
        x,y,_=R['xyzobs.px.value'].parts()
        x2,y2,_=R['xyzcal.px'].parts()
        assert x2 > 0
        assert y2 > 0
        d = np.sqrt( (x-x2)**2 + (y-y2)**2)
        med_d = np.median(d)
        print('diffBragg: %.3f (%s) Nref=%d' % (np.median(d), f, len(d) ) )
        all_d.append(d)
        all_shotd.append( np.median(d))
        nref_per_shot .append( len(d))

    all_d = hstack(all_d)
    print("median over %d shots=%f pixels (%d refls)" % (len(nref_per_shot), median(all_d), len(all_d)))
    print("Min refls per shot=%d, max refls per shot = %d, ave refls per shot=%.1f" % (min(nref_per_shot), max(nref_per_shot), mean(nref_per_shot)))
    subplot(121)
    hist( all_d, bins=args.nbins[0], histtype='step', lw=2, label=glob_s)
    if not args.linear:
        gca().set_yscale("log")
    xlabel("|calc-obs| (pixels)", fontsize=17)
    ylabel("num refls", fontsize=15)
    #legend(prop={'size':12})
    gca().tick_params(direction='in', which='both', labelsize=15)
    subplot(122)
    hist( all_shotd, bins=args.nbins[1], histtype='step', lw=2, label=glob_s)
    if not args.linear:
        gca().set_yscale("log")
    xlabel("median |calc-obs| (pixels)", fontsize=17)
    ylabel("num shots", labelpad=0, fontsize=17)
    #legend(prop={'size':11})
    gca().tick_params(direction='in', which='both', labelsize=15)

gcf().set_size_inches((7,4))
suptitle(" %d shots ; %d refls" % (len(nref_per_shot), len(all_d)))
subplots_adjust(bottom=0.2, right=0.97, top=0.92, left=0.1)

if args.noPlot:
    exit()

show()
