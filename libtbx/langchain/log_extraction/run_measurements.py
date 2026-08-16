#!/usr/bin/env python
"""P0 measurements: the frozen-screen denominator, performance, and the
pre-registered baselines every later phase is scored against.

Run with: PHENIX_LOG_CORPUS=... python3 run_measurements.py
"""
from __future__ import absolute_import, division, print_function
import os, sys, time, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from log_extraction.log_structure_extractor import scan, SCREEN_RULES, UNPARSED_STATUSES

CORPUS = os.environ["PHENIX_LOG_CORPUS"]
logs = []
for base, _, files in os.walk(CORPUS):
    for n in sorted(files):
        if n.endswith(".log") and not n.startswith("._"):
            logs.append(os.path.join(base, n))

rules = collections.Counter(); per_log = []
screen_total = 0; claimed_total = 0; outside = 0
statuses = collections.Counter(); refusals = collections.Counter()
total_bytes = 0; t0 = time.time(); slowest = (0, None)
for p in logs:
    raw = open(p, "rb").read(); total_bytes += len(raw)
    text = raw.decode("utf-8", "replace")
    t1 = time.time(); s = scan(text); dt = time.time() - t1
    if dt > slowest[0]: slowest = (dt, os.path.basename(p))
    from log_extraction.log_structure_extractor import screen_line, split_lines
    lines = split_lines(text)
    screen_total += sum(1 for l in lines if screen_line(l))
    for u in s.unparsed:
        rules[u.screen_rule] += 1
        statuses[u.status] += 1
        if u.excluded_by: refusals[u.excluded_by] += 1
    outside += len(s.claimed_outside_screen)
    claimed_total += len(s.sections)
    per_log.append((len(s.unparsed), s.n_lines, os.path.basename(p)))
elapsed = time.time() - t0

cand = screen_total
lines = sum(l for _, l, _ in per_log)
print("corpus2/work n=%d logs, %d lines, %.1f MB" % (len(logs), lines, total_bytes/1e6))
print()
print("FROZEN SCREEN -- the denominator every later phase is scored against")
print("  candidates            %d  (%.1f%% of all lines)" % (cand, 100.0*cand/lines))
for r in SCREEN_RULES:
    print("    %-14s %8d  (%.1f%% of what remains unparsed)"
          % (r, rules[r], 100.0*rules[r]/max(sum(rules.values()), 1)))
print("  logs with zero screened candidates  %d"
      % sum(1 for c, _, _ in per_log if c == 0))
q = sorted(c for c, _, _ in per_log)
print("  candidates per log: median %d, 95th %d, max %d (%s)"
      % (q[len(q)//2], q[int(.95*len(q))], q[-1],
         max(per_log)[2]))
print()
print("PARSER LEDGER")
print("  sections claimed        %d" % claimed_total)
print("  lines claimed outside the frozen screen  %d" % outside)
print("  unparsed: %s" % ", ".join("%s=%d" % (k, statuses[k])
                                   for k in UNPARSED_STATUSES))
print("  refusals: %s" % ", ".join("%s=%d" % kv
                                   for kv in sorted(refusals.items())))
print()
print("PERFORMANCE (budget: >=5 MB/s corpus-wide, <2 s any single log)")
print("  %.1f MB/s corpus-wide; slowest single log %.3f s (%s)"
      % (total_bytes/1e6/elapsed, slowest[0], slowest[1]))
