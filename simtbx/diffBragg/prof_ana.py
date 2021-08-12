

import sys
import glob

logfiles = glob.glob(sys.argv[1])
STARTUP_TIME=309.496389
if not logfiles:
    raise IOError("no files found in glob!")
NUM_ITER = 101


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


def get_color(r, g, b):
    return r/255., g/255., b/255.


funcs = ["self._update_Fcell()",
         "self.D.add_diffBragg_spots(pfs, nom_h)",
         "self._MPI_barrier()"]
colors = [get_color(255,221,113),
          "tomato",
          get_color(23,190,207)]
total = ["def _compute_functional_and_gradients"]

funcs2 =[
    "expt_list = ExperimentListFactory.from_json_file(exper_name, check_format=True)",
    "shot_modeler.GatherFromExperiment(expt, refls, sg_symbol=self.symbol)",
    "COMM.Barrier()",
    #"self._gather_Hi_information()"
]

colors2 = [
    "#777777",
    get_color(152,223,138),
    get_color(23,190,207),
    get_color(255,128,14)]
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



from pylab import *
figure(1)
subplot(211)
ax1 = gca()
prep_time = 90
bar(ranks, [prep_time]*len(ranks), color=get_color(162,200,236), width=1)
bottom = prep_time
for i_f, f in enumerate(funcs2):
    r,t = map(np.array, zip(*rank_times[f]))
    bar(r, t, bottom=bottom, width=1, color=colors2[i_f], label="%d" % i_f)
    if bottom is None:
        bottom = t
    else:
        bottom += t
setup_t = np.array([setup_times[r] for r in ranks])
#bar(ranks, setup_t, bottom=bottom, width=1, color=get_color(247,182,210))
other = np.ones(len(ranks))*STARTUP_TIME - bottom
#bar(ranks, other, bottom=bottom, width=1, color='lightgray')
ax1.vlines(np.arange(0, 20*18, 18), 0, STARTUP_TIME, ls='--', lw=0.5, color='k', alpha=0.5)
#ax1.hlines(STARTUP_TIME, min(ranks),max(ranks), color='k', ls='--', lw=0.5)
ax1.set_ylim(0, STARTUP_TIME+38)
ax1.text(0, STARTUP_TIME+10, s="total startup time: %.2f sec." % STARTUP_TIME, color='k', fontsize=8, ha='left')
for i in range(20):
    ax1.text(i*18+9-7,STARTUP_TIME-20, s="N%02d" % (i+1), fontsize=6, color='#777777', ha='left')


#ax1.legend()


subplot(212)
ax2 = gca()
bottom = None
for i_f, f in enumerate(funcs):
    r,t = map(np.array, zip(*rank_times[f]))
    t /= NUM_ITER
    bar(r, t, bottom=bottom, width=1, color=colors[i_f])
    if bottom is None:
        bottom = t
    else:
        bottom += t
total_t = np.mean(list(total_times.values()))/ NUM_ITER
#ax2.hlines(total_t, min(ranks),max(ranks), color='k', ls='--', lw=0.5)
mid_rank = (max(ranks) - min(ranks)) / 2.
ax2.text(0, total_t+0.4, s="average iteration time: %.2f sec." % total_t, color='k', fontsize=8, ha='left')
k = get_color(65,68,81)
other = np.ones(len(ranks))*total_t - bottom
#bar(r, other, bottom=bottom, width=1, color='lightgray')
ax2.vlines(np.arange(0, 20*18, 18), 0, total_t, ls='--', lw=0.5, color='k', alpha=0.5)
ax2.set_ylim(0, total_t+1)
#plot([min(ranks), max(ranks)], [total_t, total_t], color=k, lw=1)
#vlines([min(ranks), max(ranks)], 0, total_t, color=k, lw=1)
for i in range(20):
    ax2.text(i*18+9-7,total_t-0.5, s="N%02d" % (i+1), fontsize=6, color='#777777', ha='left')

for ax in ax1,ax2:
    ax.tick_params(labelsize=8)
    #ymin, ymax = ax.get_ylim()

ax2.set_xlabel('rank', fontsize=13)
ax1.set_xticklabels([])
show()



