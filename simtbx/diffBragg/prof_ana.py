

import sys
import glob

def get_color(r, g, b):
    return r/255., g/255., b/255.
funcs = ["self._update_Fcell()",
         "self.D.add_diffBragg_spots(pfs, nom_h)",
         "self._MPI_barrier()"]
funcs2 =[
    "expt_list = ExperimentListFactory.from_json_file(exper_name, check_format=True)",
    "shot_modeler.GatherFromExperiment(expt, refls, sg_symbol=self.symbol)",
    "COMM.Barrier()",
    #"self._gather_Hi_information()"
]

def main(input_glob):
    logfiles = glob.glob(input_glob)
    if not logfiles:
        raise IOError("no files found in glob!")

    log = []
    for logfile in logfiles:
        print("opening log %s" % logfile)
        log += open(logfile, 'r').readlines()
    log = [l for l in log if l.startswith('RANK')]
    key = lambda x: int( x.strip().split('|')[0].strip().split("RANK")[1])
    from itertools import groupby
    gb = groupby(sorted(log, key=key), key)
    logs_by_rank = {k:list(v) for k,v in gb}
    ranks = list(logs_by_rank.keys())

    total = ["def _compute_functional_and_gradients"]

    rank_times = {f:[] for f in funcs+funcs2+total}
    total_times = {r:0 for r in ranks}
    setup_times = {r:0 for r in ranks}

    for rank in ranks:
        log = logs_by_rank[rank]
        for i in range(1,len(log)):
            l = log[i]
            if "def _setup" in l:

                ts = float(log[i+1].strip().split()[5])/1000.
                setup_times[rank] += ts

            for tt in total:
                if tt in l:
                    ts = float(log[i+1].strip().split()[5])/1000.
                    total_times[rank] += ts

            for ff in funcs[:2] + funcs2[:2]: #+ funcs2[3:]:
                if ff in l:
                    ts = float(l.strip().split()[3])/1000.
                    print(rank, ff, ts)
                    rank_times[ff].append((rank, ts))
            if funcs[2] in l and "aggregate barrier" in log[i-1]:
                ff = funcs[2]
                ts = float(l.strip().split()[3])/1000.
                print(rank, ff, ts)
                rank_times[ff].append((rank, ts))
            if funcs2[2] in l and " ENTER BARRIER" in log[i-1]:
                ff = funcs2[2]
                ts = float(l.strip().split()[3])/1000.
                print(rank, ff, ts)
                rank_times[ff].append((rank, ts))
    return ranks, rank_times, setup_times, total_times

#=============================================================
#out10 = main("10node_3/*prof*.log")
out10 = main("perl_10node_1per/*prof*.log")
out15 = main("perl_15node_1per/*prof*.log")
COMPARE_ONE = False
if COMPARE_ONE:
    out20 = main("1per/*prof*.log")
else:
    #out20 = main("20node_3/*prof*.log")
    out20 = main("perl_20node_1per/*prof*.log")

#NUM_ITER = 101
#STARTUP_TIME_10 = 493.358889
#STARTUP_TIME_20 = 309.496389
STARTUP_TIME_10 = 37.28
STARTUP_TIME_15 = 31.564
STARTUP_TIME_20 = 28.145

out10_summit = main("10node_1per/*prof*.log")
#out15_summit = main("perl_15node_1per/*prof*.log")
out20_summit = main("20node_1per/*prof*.log")

STARTUP_TIME_10_SUMMIT = 169
STARTUP_TIME_15_SUMMIT = None
STARTUP_TIME_20_SUMMIT = 156.43

out10 += (STARTUP_TIME_10, 10,)
out15 += (STARTUP_TIME_15, 15,)
out20 += (STARTUP_TIME_20, 20,)

out10_summit += (STARTUP_TIME_10_SUMMIT, 10,)
out20_summit += (STARTUP_TIME_20_SUMMIT, 20,)

colors = [get_color(255,221,113),
          "tomato",
          get_color(23,190,207)]
colors2 = [
    "#777777",
    get_color(152,223,138),
    get_color(23,190,207),
    get_color(255,128,14)]

from pylab import *

fig, (ax1,ax2) = subplots(nrows=2, ncols=1, figsize=(6,4.1))
from textwrap import wrap
def wrap_t(txt, n=22):
    return "\n".join(wrap( txt, n))

def add_bars_to_ax_startup(ax,main_out):
    ranks, rank_times, setup_times, total_times, STARTUP_TIME, NODES = main_out
    #prep_time = 90
    prep_time = 9
    ax.bar(np.array(ranks)+0.5, [prep_time]*len(ranks), color=get_color(162,200,236), width=1,
            label=wrap_t("prep_dataframe (diffBragg): %.1f$\pm$0.0 sec." % prep_time))

    bottom = prep_time

    for i_f, f in enumerate(funcs2):
        r,t = map(np.array, zip(*rank_times[f]))
        ave_t = np.mean(t)
        sig_t = np.std(t)
        if i_f==0:
            label = wrap_t("from_json_file (dxtbx): %.1f$\pm$%.1f sec." % (ave_t, sig_t), 18)
        elif i_f==1:
            label = wrap_t("GatherFromExp (diffBragg): %.1f$\pm$%.1f sec." % (ave_t, sig_t))
        else:
            label= wrap_t("barrier (mpi4py): %.1f$\pm$%.1f sec." % (ave_t, sig_t), 20)
        ax.bar(r+0.5, t, bottom=bottom, width=1, color=colors2[i_f], label=label)
        if bottom is None:
            bottom = t
        else:
            bottom += t
    handles, leg_labels = ax.get_legend_handles_labels()
    leg1 = ax.legend(handles[::-1], leg_labels[::-1], prop={'size': 7}, loc="center left",
                      bbox_to_anchor=(1, 0.47), labelspacing=1)
    fr = leg1.get_frame()
    fr.set_alpha(1)
    fr.set_facecolor('w')
    fr.set_edgecolor('w')
    fr.set_linewidth(0.5)
    ax.text(1.05, 0.97, "%d node job" % NODES, fontsize=10, transform=ax.transAxes, fontweight='bold')
    setup_t = np.array([setup_times[r] for r in ranks])
    #bar(ranks, setup_t, bottom=bottom, width=1, color=get_color(247,182,210))
    other = np.ones(len(ranks))*STARTUP_TIME - bottom
    #bar(ranks, other, bottom=bottom, width=1, color='lightgray')
    ax.vlines(np.arange(0, NODES*18, 18), 0, STARTUP_TIME, ls='--', lw=0.5, color='k', alpha=0.5)
    #ax.hlines(STARTUP_TIME, min(ranks),max(ranks), color='k', ls='--', lw=0.5)
    ax.set_ylim(0, 493.358889+70)
    ax.text(3, STARTUP_TIME+20, s="total startup time: %.2f sec." % STARTUP_TIME, color='k', fontsize=8, ha='left')
    for i in range(NODES):
        ax.text(i*18+9-7,STARTUP_TIME-20, s="N%02d" % (i+1), fontsize=6, color='#777777', ha='left')
    ax.set_xlim(0,360)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.tick_params(direction='in')
    ax.tick_params(labelsize=8,pad=1)

#add_bars_to_ax_startup(ax1,out10)
#add_bars_to_ax_startup(ax2,out20)
##ax1.set_xticklabels([])
#ax2.set_xlabel("rank", fontsize=12)
#
##fig.set_size_inches((W+0.1*W, H))
#subplots_adjust(left=0.12, top=0.95, bottom=0.1, right=0.78, hspace=0.15)
##fig.set_size_inches((5.9, 4.1))
#ax1.text(-0.15, -0.2, 'seconds', fontsize=12, rotation=90,transform=ax1.transAxes )
#show()
#exit()
# ====================================================================================

#ax1.legend()

def add_bars_to_ax_iters(ax, main_out, NUM_ITER):
    ranks, rank_times, setup_times, total_times, STARTUP_TIME, NODES = main_out
    bottom = None
    for i_f, f in enumerate(funcs):
        r,t = map(np.array, zip(*rank_times[f]))
        t /= NUM_ITER
        ave_t = np.mean(t)
        sig_t = np.std(t)
        if i_f==0:
            label = wrap_t("update_Fcell (diffBragg): %.1f$\pm$%.1f sec." % (ave_t, sig_t), 20)
        elif i_f==1:
            label = wrap_t("add_diffBragg_spots (diffBragg): %.1f$\pm$%.1f sec." % (ave_t, sig_t))
        else:
            label= wrap_t("barrier (mpi4py): %.1f$\pm$%.1f sec." % (ave_t, sig_t), 20)
        ax.bar(r+0.5, t, bottom=bottom, width=1, color=colors[i_f], label=label)
        if bottom is None:
            bottom = t
        else:
            bottom += t

    handles, leg_labels = ax.get_legend_handles_labels()
    leg2 = ax.legend(handles[::-1], leg_labels[::-1], prop={'size': 7}, loc="center left",
                      bbox_to_anchor=(1, 0.35), labelspacing=1)
    fr2 = leg2.get_frame()
    fr2.set_alpha(1)
    fr2.set_facecolor('w')
    fr2.set_edgecolor('w')
    fr2.set_linewidth(0.5)
    ax.text(1.05, 0.73, "%d node job" % NODES, fontsize=10, transform=ax.transAxes, fontweight='bold')

    total_t = np.mean(list(total_times.values()))/ NUM_ITER
    #ax.hlines(total_t, min(ranks),max(ranks), color='k', ls='--', lw=0.5)
    #ax.text(5, total_t+0.7, s="average iteration time: %.2f sec." % total_t, color='k', fontsize=8, ha='left')
    ax.text(5, 8, s="average iteration time: %.2f sec." % total_t, color='k', fontsize=8, ha='left')
    k = get_color(65,68,81)
    other = np.ones(len(ranks))*total_t - bottom
    #bar(r, other, bottom=bottom, width=1, color='lightgray')
    ax.vlines(np.arange(0, NODES*18, 18), 0, total_t, ls='--', lw=0.5, color='k', alpha=0.5)
    ax.set_ylim(0, total_t+1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0,11)
    ax.set_xlim(0,360)
    for i in range(NODES):
        ax.text(i*18+9-7,total_t-0.5, s="N%02d" % (i+1), fontsize=6, color='#777777', ha='left')
    ax.tick_params(labelsize=8, pad=1)

NN = 502
#NN = 101
add_bars_to_ax_iters(ax1,out10, 502)
#add_bars_to_ax_iters(ax1, out15, NN)
add_bars_to_ax_iters(ax2, out10_summit, 101)

#add_bars_to_ax_iters(ax1, out10_summit, NN)
#add_bars_to_ax_iters(ax2, out20_summit, NN)
#ax1.set_xticklabels([])
ax2.set_xlabel("rank", fontsize=12)

ax1.text(-0.11, -0.2, 'seconds', fontsize=12, rotation=90,transform=ax1.transAxes)
#fig.set_size_inches((W+0.1*W, H))
subplots_adjust(left=0.08, top=0.95, bottom=0.1, right=0.76, hspace=0.15)
#fig.set_size_inches((5.9, 4.1))
show()



