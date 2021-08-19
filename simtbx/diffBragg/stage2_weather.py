
import sys
import glob
import numpy as np
from pylab import *

key_phrases = {
    0: 'MPI aggregation of func and grad',
    1: 'run diffBragg',
    2: 'finished diffBragg',
    3: 'BEGIN FUNC GRAD',
    4: "start update Fcell",
    5: "done update Fcell",
    6: "EVENT: BEGIN loading experiment list",
    7: "EVENT: DONE loading experiment list",
    8: 'read input pickle',
    9: 'prep dataframe',
    10: "aggregate barrier",
    11: 'Functional',
    12: "EVENT: launch refiner",
    13: "EVENT: Gathering global HKL information",
    14: "Setup begins!",
    15: "DATA; ENTER BARRIER",
    16: "DATA; EXIT BARRIER",
    17: "Setup begins",
    18: "Setup ends",
    19: "FINISHED gather global HKL information",
    20: "EVENT: BEGIN prep dataframe",
    21: "EVENT: DONE prep dataframe",
    22: "DUMP param and Zscore data",
    23: "Time to dump param and Zscore data",
    24: "EVENT: LOADING ROI DATA",
    25: "EVENT: DONE LOADING ROI",
    26: "Time for MPIaggregation",
    27: "_launcher runno setup"
}

def get_color(r,g,b):
    return r/255.,g/255.,b/255.

durations = {(key_phrases[4], key_phrases[5]): {"color": get_color(255,221,113)},
             (key_phrases[1], key_phrases[2]): {"color": 'tomato'},
             (key_phrases[10], key_phrases[11]): {"color": get_color(23,190,207)},
             (key_phrases[15], key_phrases[16]): {"color": get_color(23,190,207)},
             (key_phrases[17], key_phrases[18]): {"color": get_color(247,182,210)},
             (key_phrases[13], key_phrases[19]): {"color": get_color(255,128,14)},
             (key_phrases[20], key_phrases[21]): {"color": get_color(162,200,236)},
             #(key_phrases[22], key_phrases[23]): {"color": get_color(255,221,113)},
             (key_phrases[24], key_phrases[25]): {"color": get_color(103,191,92)},
             #(key_phrases[11], key_phrases[26]): {"color": get_color(158,218,229)},
             (key_phrases[6], key_phrases[7]): {"color": '#777777'}}

major_events = {
    key_phrases[3]: {'marker': 'd', 'ms': 4, 'color': get_color(152,  223,   138)},
    #key_phrases[27]: {'marker': 's', 'ms': 4, 'color': 'w'},
}


def main(logfiles):
    log = []
    for logfile in logfiles:
        print("opening log %s" % logfile)
        log += open(logfile, 'r').readlines()

    print("getting rank")
    ranks = [int(l.split("|")[0].strip().split("RANK")[1].split(":")[0]) for l in log if l.startswith("RANK")]
    nranks = len(set(ranks))
    print("Found logs for %d ranks" % nranks)

    logs_per_rank = {rnk: [] for rnk in set(ranks)}
    for line, rank in zip(log, ranks):
        logs_per_rank[rank].append(line)
    print("Made log for each rank")
    from datetime import datetime
    import logging
    F = logging.Formatter()
    dateformat = F.default_time_format+",%f"
    # TODO figure out why loggers format of 0-padded millisecond isnt reproducible
    # convert timestamps
    # initial time in ranks log
    # find earliest time in all ranks
    ranks = sorted(logs_per_rank.keys())
    all_dt0 = []
    for rank in ranks:
        rank_log = logs_per_rank[rank]
        rank_t0 = rank_log[0].split("|")[1].strip()
        rank_dt0 = datetime.strptime(rank_t0 + "000", dateformat)
        all_dt0.append(rank_dt0)
    dt0 = min(all_dt0)

    print("Converting timestamps")
    import time
    t = time.time()
    all_delt = {rank:[] for rank in ranks}
    for rank in logs_per_rank:
        rank_log = logs_per_rank[rank]

        phrase_delt = {}
        for phrase in key_phrases.values():
            phrase_log = [l for l in rank_log if phrase in l]
            if not phrase_log:
                print("WARNING phrase %s not in log" % phrase)
                continue
            ts = [l.split("|")[1].strip() for l in phrase_log]
            dt = [datetime.strptime(t+"000", dateformat) for t in ts]
            seconds_from_start = [(t - dt0).total_seconds() for t in dt]
            phrase_delt[phrase] = seconds_from_start
        print(rank)
        all_delt[rank] = phrase_delt
    t = time.time()-t
    print("Took %f sec" % t)

    duration_timers = {}
    startup_times = {}
    collections = {}
    for rank in ranks:
        D = all_delt[rank]

        refine_start = D[key_phrases[27]][0]
        startup_times[rank] = refine_start
        #for event in major_events:
        #    tvals = D[event]
        #    offset = [rank] * len(D[event])
        #    ax.plot(np.array(tvals)-refine_start, offset, lw=0, mew=0.3, mec='k', **major_events[event])

        collections[rank] = []
        for i, dur in enumerate(list(durations)):
            start,stop = dur

            tvals_start = np.array(D[start]) - refine_start
            tvals_stop = np.array(D[stop]) - refine_start
            patches = []
            if dur not in duration_timers:
                duration_timers[dur] = []
            for t1,t2 in zip(tvals_start, tvals_stop):
                print(i, t1, t2)
                duration_timers[dur].append(t2-t1)
                xy = np.array([(t1,rank-0.5), (t1, rank+0.5), (t2, rank+0.5), (t2, rank-0.5)])
                color = durations[(start, stop)]["color"]
                patch = mpl.patches.Polygon(xy=xy, color=color, closed=True) #, ec='k', lw=0.5)
                patches.append(patch)
            C = mpl.collections.PatchCollection(patches, match_original=True)
            collections[rank].append(C)
    return ranks, all_delt, collections


#logfiles10 = glob.glob("10node_3/h37n16-main_stage2.log")
#logfiles20 = glob.glob("20node_3/h50n01-main_stage2.log")
#logfiles1 = glob.glob("1per/h34n16-main_stage2.log")
#out1 = main(logfiles1)
#out10 = main(logfiles10)
#out20 = main(logfiles20)
#np.savez("out1", a=out1[0], b=out1[1], c=out1[2])
#np.savez("out10", a=out10[0], b=out10[1], c=out10[2])
#np.savez("out20", a=out20[0], b=out20[1], c=out20[2])


O = np.load("out10.npz", allow_pickle=True)
out10 = O['a'][()], O['b'][()], O['c'][()]
#O2 = np.load("out20.npz", allow_pickle=True)
#out20 = O2['a'][()], O2['b'][()], O2['c'][()]

O1 = np.load("out1.npz", allow_pickle=True)
out20 = O1['a'][()], O1['b'][()], O1['c'][()]

BEFORE = False
if BEFORE:
    major_events = {}

ONE_COMP = True

def add_weather_to_ax(ax, main_out):
    ranks, all_delt, collections = main_out
    ranks = list(ranks)
    refine_start = all_delt[ranks[0]][key_phrases[27]][0]

    #duration_timers = {}
    #startup_times = {}
    for rank in ranks:
        D = all_delt[rank]

    #    refine_start = D[key_phrases[27]][0]
    #    startup_times[rank] = refine_start
        for event in major_events:
            tvals = D[event]
            offset = [rank] * len(D[event])
            ax.plot(np.array(tvals)-refine_start, offset, '.', color='k', mew=0.3, ms=3, mec='k')

    #    for i, dur in enumerate(list(durations)):
    #        start,stop = dur

    #        tvals_start = np.array(D[start]) - refine_start
    #        tvals_stop = np.array(D[stop]) - refine_start
    #        patches = []
    #        if dur not in duration_timers:
    #            duration_timers[dur] = []
    #        for t1,t2 in zip(tvals_start, tvals_stop):
    #            print(i, t1, t2)
    #            duration_timers[dur].append(t2-t1)
    #            xy = np.array([(t1,rank-0.5), (t1, rank+0.5), (t2, rank+0.5), (t2, rank-0.5)])
    #            color = durations[(start, stop)]["color"]
    #            patch = mpl.patches.Polygon(xy=xy, color=color, closed=True) #, ec='k', lw=0.5)
    #            patches.append(patch)
    #        C = mpl.collections.PatchCollection(patches, match_original=True)
    for rank in ranks:
        for C in collections[rank]:
            ax.add_collection(C)

    #for rank in ranks:
    #    D = all_delt[rank]
    #    times = D[key_phrases[5]]
    #    times = np.sort(times)
    #    tper = np.diff(times)
    #    print(rank, "min, max, mean, median", np.min(tper), np.max(tper), np.mean(tper), np.median(tper))
    #    print("first 5:", tper[:5])

    #for dur in duration_timers:
    #    print(dur)
    #    print("Average duration: %f sec" % np.mean(duration_timers[dur]))
    #    print("Total duration: %f sec" %  ( np.sum(duration_timers[dur]) / len(ranks)))
    #    print()


    #xlabel("post-startup runtime (sec.)", fontsize=12)
    xl = ax.get_xlim()
    for r in ranks + [min(ranks)-1, max(ranks)+1]:
        if ONE_COMP:
            ax.plot([-refine_start,240], [r-0.5, r-0.5], color='k', lw=0.25)
        else:
            ax.plot([-refine_start,60], [r-0.5, r-0.5], color='k', lw=0.25)

    ax.set_facecolor('w') #get_color(171,171,171))
    ax.set_xlim(-0.1, 42) #-refine_start, 0)
    #xlim(-refine_start, 0)
    ax.set_ylim(min(ranks)-0.5, max(ranks)+0.5)
    nranks = len(ranks)
    ax.set_yticklabels(list(map(str, np.arange(nranks)[::3])))
    ax.set_yticks(ranks[::3])
    ax.tick_params(labelsize=8, pad=1)

    ax.vlines(-refine_start, min(ranks)-0.5,max(ranks)+0.5, color='k',ls='--', lw=1.25)
    #print("Mean startup time=%f sec" % np.mean(list(startup_times.values())))
    #gca().set_xlim((-493.38, 0.0))
    print(ax.get_xlim())

fig, (ax1,ax2) = subplots(nrows=2, ncols=1, figsize=(6,3))

add_weather_to_ax(ax1, out10)
add_weather_to_ax(ax2, out20)
#ax2.set_xlabel("startup time (seconds before first iteration)", fontsize=12)
ax2.set_xlabel("time after first iteration (seconds)", fontsize=12)
ax1.text(-0.1, -0.2, 'rank', fontsize=12, rotation=90,transform=ax1.transAxes )


patches = [
    mpl.patches.Patch(color=get_color(255,221,113), ec='k', lw=0.5,label="update_Fcell"),
    mpl.patches.Patch(color='tomato', ec='k', lw=0.5, label='add_diffBragg_spots'),
    mpl.patches.Patch(color=get_color(23,190,207),ec='k', lw=0.5, label='MPI barrier'),
    mpl.patches.Patch(color='w', ec='k', lw=0.5, label='untracked'),
    plt.plot([], [], ls="", marker='.', color='k',ms=5,mec=None,  label="iteration begin")[0]
]

patches_before = [
    mpl.patches.Patch(color=get_color(162,200,236), ec='k', lw=0.5,label="prep_dataframe"),
    mpl.patches.Patch(color='#777777', ec='k', lw=0.5, label='from_json_file'),
    mpl.patches.Patch(color=get_color(103,191,92), ec='k', lw=0.5, label='GatherFromExp'),
    mpl.patches.Patch(color=get_color(23,190,207), ec='k', lw=0.5, label='MPI barrier'),
    mpl.patches.Patch(color=get_color(255,128,14), ec='k', lw=0.5, label='gather_Hi_info'),
    mpl.patches.Patch(color=get_color(247,182,210), ec='k', lw=0.5, label='setup'),
    mpl.patches.Patch(color='w', ec='k', lw=0.5, label='untracked'),
    plt.plot([],[], ls="--", color='k',label="program start")[0]
]


#leg = ax1.legend(handles=patches, markerscale=0.1,
#                 bbox_to_anchor=(.99,-0.13),
#                 prop={'size':7.5},
#                 loc="center left")
#fr.set_linewidth(0.5)

subplots_adjust(left=0.08, top=0.97, bottom=0.14, right=0.76, hspace=0.25)

if BEFORE:
    ax1.set_xlim(-500,0)
    ax2.set_xlim(-500,0)
    ax2.set_xlabel("time before first iteration (seconds)", fontsize=12)
    leg = ax1.legend(handles=patches_before,
                     bbox_to_anchor=(1,-0.13),
                     prop={'size':7.5},
                     loc="center left")
else:
    leg = ax1.legend(handles=patches, markerscale=1,
                     bbox_to_anchor=(.99,-0.13),
                     prop={'size':7.5},
                     loc="center left")
fr = leg.get_frame()
fr.set_alpha(1)
fr.set_facecolor('w')
fr.set_edgecolor('w')
if BEFORE:
    ax1.text(1.01, 0.90, "10 node job", fontsize=10, transform=ax1.transAxes, fontweight='bold')
    if ONE_COMP:
        ax1.text(1.01, 0.05, "10 node job", fontsize=10, transform=ax2.transAxes, fontweight='bold')
    else:
        ax1.text(1.01, 0.05, "20 node job", fontsize=10, transform=ax2.transAxes, fontweight='bold')

else:
    ax1.text(1.01, 0.95, "10 node job", fontsize=10, transform=ax1.transAxes, fontweight='bold')
    if ONE_COMP:
        ax1.text(1.01, -0.05, "10 node job", fontsize=10, transform=ax2.transAxes, fontweight='bold')
    else:
        ax1.text(1.01, -0.05, "20 node job", fontsize=10, transform=ax2.transAxes, fontweight='bold')
#ax2.spines['left'].set_position(('data', -310))
if not BEFORE and ONE_COMP:
    ax1.set_xlim(0,120)
    ax2.set_xlim(0,120)
    yl = ax2.get_ylim()
    ax2.set_ylim(yl[0], yl[1]+12)
    ax2.set_xlim(0,140)
    ax1.set_xlim(0,140)
show()
