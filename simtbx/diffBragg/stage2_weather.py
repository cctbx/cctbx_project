
import sys
import glob
logfiles = glob.glob(sys.argv[1])
key_phrases = [
    'MPI aggregation of func and grad',
    'run diffBragg',
    'BEGIN FUNC GRAD',
    "update Fcell",
    "EVENT: BEGIN loading experiment list",
    'read input pickle',
    'prep dataframe']

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

    phrase_delt = []
    for phrase in key_phrases:
        phrase_log = [l for l in rank_log if phrase in l]
        if not phrase_log:
            print("WARNING phrase %s not in log" % phrase)
            continue
        ts = [l.split("|")[1].strip() for l in phrase_log]
        dt = [datetime.strptime(t+"000", dateformat) for t in ts]
        seconds_from_start = [(t - dt0).total_seconds() for t in dt]
        phrase_delt.append(seconds_from_start)
    print(rank)
    all_delt[rank] = phrase_delt
t = time.time()-t
print("Took %f sec" % t)

from pylab import *
for rank in ranks:
    D = all_delt[rank]
    plot(np.array(D[0])/1., [rank]*len(D[0]), 'o', ms=3.5, color='C0')
    plot(np.array(D[1])/1., [rank+0]*len(D[1]), '.', ms=2.5, color='tomato')
    plot(np.array(D[2])/1., [rank+0]*len(D[2]), 'd', ms=3.5,color='C2')
    plot(np.array(D[3]) / 1., [rank + 0] * len(D[3]), '>', ms=3.5, color='C4')
    plot(np.array(D[4]) / 1., [rank + 0] * len(D[4]), 's', ms=2.5, color='#777777')
    plot(np.array(D[5]) / 1., [rank + 0] * len(D[5]), '^', ms=2.5, color='C1')
    plot(np.array(D[6]) / 1., [rank + 0] * len(D[6]), '*', ms=2.5, color='C5')
xlabel("runtime (sec)", fontsize=15)
ylabel("rank")
show()
