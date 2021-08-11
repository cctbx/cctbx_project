
import sys
import glob
logfiles = glob.glob(sys.argv[1])

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
    26: "Time for MPIaggregation"
}

def get_color(r,g,b):
    return r/255.,g/255.,b/255.

durations = {(key_phrases[4], key_phrases[5]): {"color": 'C4'},
             (key_phrases[1], key_phrases[2]): {"color": 'tomato'},
             (key_phrases[10], key_phrases[11]): {"color": 'C0'},
             (key_phrases[15], key_phrases[16]): {"color": 'C0'},
             (key_phrases[17], key_phrases[18]): {"color": get_color(247,182,210)},
             (key_phrases[13], key_phrases[19]): {"color": get_color(255,195,86)},
             (key_phrases[20], key_phrases[21]): {"color": get_color(162,200,236)},
             (key_phrases[22], key_phrases[23]): {"color": get_color(255,221,113)},
             (key_phrases[24], key_phrases[25]): {"color": get_color(103,191,92)},
             #(key_phrases[11], key_phrases[26]): {"color": get_color(158,218,229)},
             (key_phrases[6], key_phrases[7]): {"color": '#777777'}}

major_events = {
    key_phrases[3]: {'marker': 'd', 'ms': 4, 'color': get_color(152,  223,   138)},
}

log = []
for logfile in logfiles:
    print("opening log %s" % logfile)
    log += open(logfile, 'r').readlines()

print("getting rank")
ranks = [int(l.split("|")[0].strip().split("RANK")[1].split(":")[0]) for l in log if l.startswith("RANK")]
nranks = len(set(ranks))
print("Found logs for %d ranks" % nranks)

logs_per_rank = {rnk:[] for rnk in set(ranks)}
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

from pylab import *
for rank in ranks:
    D = all_delt[rank]
    #from IPython import embed
    #embed()
    #for phrase in key_phrases:
    #    tvals = D[phrase]
    #    offset = [rank] * len(D[phrase])
    #    scatter(tvals, offset, **key_phrases[phrase])
    for event in major_events:
        tvals = D[event]
        offset = [rank] * len(D[event])
        plot(tvals, offset, lw=0, mew=0.3, mec='k', **major_events[event])

    colors = ['C5', 'tomato', 'C0']
    for i, (start, stop) in enumerate(list(durations)):
        tvals_start = D[start]
        tvals_stop = D[stop]
        patches = []

        for t1,t2 in zip(tvals_start, tvals_stop):
            print(i, t1, t2)
            xy = np.array([(t1,rank-0.5), (t1, rank+0.5), (t2, rank+0.5), (t2, rank-0.5)])
            color = durations[(start, stop)]["color"]
            patch = mpl.patches.Polygon(xy=xy, color=color, closed=True) #, ec='k', lw=0.5)
            patches.append(patch)
        C = mpl.collections.PatchCollection(patches, match_original=True)
        gca().add_collection(C)
        #    plot([t1,t2], [rank,rank], color='C%d' % i)
xlabel("runtime (sec)", fontsize=15)
ylabel("rank")
xl = gca().get_xlim()
for r in ranks:
    plot(xl, [r-0.5, r-0.5], color='k', lw=0.4)
gca().set_facecolor('lightgray')
show()
