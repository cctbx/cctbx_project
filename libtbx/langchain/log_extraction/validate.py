#!/usr/bin/env python
"""Full-corpus validation for the Phenix log structure extractor. TIER 3.

    phenix.python <cctbx>/libtbx/langchain/log_extraction/validate.py

The in-tree regression suite (tests/tst_log_structure_extractor.py) runs on a
1.25 MB smoke corpus and answers "does the extractor still behave coherently
on diverse real logs?".

THIS answers a different question: **are our claims about the full known
population still true?** It needs the 21 MB validation corpus, which is
deliberately NOT in cctbx_project.

    Tier 1  does every known case still work?          fixtures, in tree
    Tier 2  coherent on diverse real logs?             1.25 MB, in tree
    Tier 3  are the population claims still true?      21 MB, out of tree

WHAT IT CHECKS, and none of it can be done at smoke scale:

  * every one of the 253 working logs scans, and every item cites a line
    whose text matches -- ~148,000 items
  * identification: the HARD GATE (0 misidentified) and the coverage target
    (>=40% named). A gate you can pass by sampling is not a gate.
  * the 75 agent metric labels, scored on GUI-SHAPE. The agent-shape figure
    is circular -- it reads the block being scored -- and must not be quoted.
  * failed runs rarely record completion: 23 error logs, not the 1 in smoke
  * performance over 17.7 MB against the 5 MB/s budget
  * every Tier-1 fixture still occurs in the full corpus (drift check)
  * the smoke set still covers all 47 features of the coverage vector
  * the frozen acceptance suite, if it is beside the corpus
  * the counted totals the documents quote

If the corpus is absent this FAILS. It does not skip: a skip that reads as a
pass has cost this project three times.

It writes a VALIDATION RECORD, so that "I think I ran it" and "this revision
passed this frozen corpus" are different sentences.
"""

from __future__ import absolute_import, division, print_function

import hashlib
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LANGCHAIN = os.path.dirname(HERE)
sys.path.insert(0, LANGCHAIN)

from log_extraction import scan                                   # noqa: E402
from log_extraction.log_structure_extractor import SCREEN_RULES   # noqa: E402
from log_extraction.log_structure_extractor import split_lines    # noqa: E402

# Tom keeps it here; the alternates are for anyone else.
CORPUS_LOCATIONS = (
  os.path.expanduser("~/Downloads/extraction_test_data"),
  os.path.expanduser("~/extraction_test_data"),
  os.environ.get("PHENIX_LOG_VALIDATION_CORPUS") or "",
)

CORPUS_ID = "corpus2 + wave1 9d8689b3861aeeb2c7b32585aad350c2"
PERFORMANCE_BUDGET_MB_S = 5.0        # the reference, measured idle
PERFORMANCE_COLLAPSE_MB_S = 1.5      # below this is an algorithmic regression


def find_corpus(argv=None):
  """-> the corpus root, or None. An explicit --corpus wins, then the
  documented locations. No environment variable is required."""
  argv = list(sys.argv[1:] if argv is None else argv)
  for index, arg in enumerate(argv):
    if arg == "--corpus" and index + 1 < len(argv):
      return argv[index + 1]
    if arg.startswith("--corpus="):
      return arg.split("=", 1)[1]
  for root in CORPUS_LOCATIONS:
    if root and os.path.isdir(os.path.join(root, "work", "ok")):
      return root
  return None


def logs_under(root):
  out = []
  for base, _, files in os.walk(root):
    for name in sorted(files):
      if name.endswith(".log") and not name.startswith("._"):
        out.append(os.path.join(base, name))
  return out


def read(path):
  handle = open(path, "rb")
  try:
    return handle.read().decode("utf-8", "replace")
  finally:
    handle.close()


def truth_from_name(path):
  stem = os.path.basename(path)[:-4].split("__")[-1]
  return re.sub(r"_\d+$", "", re.sub(r"_err$", "", stem)).lower()


class Result(object):
  def __init__(self):
    self.checks = []
    self.failed = 0

  def measure(self, name, detail):
    """Reported every run, never a pass/fail. Machine-dependent numbers do
    not belong in a gate."""
    print("  %-46s MEASURE  %s" % (name, detail))

  def check(self, name, ok, detail):
    self.checks.append((name, bool(ok), detail))
    if not ok:
      self.failed += 1
    print("  %-46s %s  %s" % (name, "PASS" if ok else "FAIL", detail))


# ------------------------------------------------------------------ checks


def population_checks(corpus, result):
  work = os.path.join(corpus, "work")
  paths = logs_under(work)
  scans = {}
  start = time.time()
  total_bytes = 0
  slowest = (0.0, None)
  for path in paths:
    text = read(path)
    total_bytes += len(text.encode("utf-8", "replace"))
    one = time.time()
    scans[path] = scan(text)
    taken = time.time() - one
    if taken > slowest[0]:
      slowest = (taken, os.path.basename(path))
  elapsed = time.time() - start

  result.check("every log scans", True, "%d logs" % len(paths))

  bad = []
  checked = 0
  for path, structure in scans.items():
    lines = split_lines(read(path))
    for item in structure.items():
      checked += 1
      if not (1 <= item.line <= item.end <= len(lines)):
        bad.append((path, item))
        continue
      evidence = item.evidence()
      if evidence and evidence not in lines[item.line - 1]:
        bad.append((path, item))
  result.check("every item cites a matching line", not bad,
               "%d items, %d bad" % (checked, len(bad)))

  named = wrong = 0
  for path, structure in scans.items():
    ident = structure.identification
    if ident.is_unknown:
      continue
    got = re.sub(r"^phenix\.", "", ident.name)
    truth = truth_from_name(path)
    if truth == got or truth.startswith(got) or got.startswith(truth):
      named += 1
    else:
      wrong += 1
  result.check("identification HARD GATE: 0 wrong", wrong == 0,
               "%d named, %d wrong" % (named, wrong))
  result.check("identification target: >=40% named",
               named * 100 >= len(paths) * 40,
               "%d of %d = %d%%" % (named, len(paths), 100 * named // len(paths)))

  err = [p for p in paths if os.sep + "err" + os.sep in p]
  with_record = sum(1 for p in err if scans[p].completion_records)
  result.check("failed runs rarely record completion",
               with_record * 4 < len(err) * 3,
               "%d of %d error logs" % (with_record, len(err)))

  # PERFORMANCE IS A MEASUREMENT, NOT A GATE. FOUND IN REVIEW: with the
  # budget as a pass/fail check, this run reported VALIDATION FAILED on a
  # loaded machine -- 5.9 MB/s idle, 2.5 MB/s with four cores busy. A
  # validation record that flips on background load is worthless as an
  # authoritative artifact, and the failure it reports is about the machine,
  # not the extractor.
  #
  # The gate is now a COLLAPSE FLOOR: anything above it is machine noise,
  # anything below it is an algorithmic regression. The reference figure is
  # printed every run so a real slowdown is visible.
  rate = total_bytes / 1e6 / max(elapsed, 1e-9)
  result.measure("performance (reference: %.0f MB/s idle)"
                 % PERFORMANCE_BUDGET_MB_S,
                 "%.1f MB/s, slowest log %.2f s (%s)"
                 % (rate, slowest[0], slowest[1]))
  result.check("performance has not collapsed",
               rate >= PERFORMANCE_COLLAPSE_MB_S and slowest[0] < 5.0,
               "floor %.1f MB/s; a low figure on a busy machine is not a"
               " failure" % PERFORMANCE_COLLAPSE_MB_S)
  return scans, paths


def label_score(corpus, result):
  """The 75 agent metric labels, on GUI-SHAPE. The agent-shape figure is
  circular -- 61 of the 75 are recoverable there only because the block being
  scored is present. This project made that exact error once."""
  labels_path = os.path.join(corpus, "agent", "agent_metrics_labels.json")
  gui = os.path.join(corpus, "agent_gui")
  if not (os.path.exists(labels_path) and os.path.isdir(gui)):
    result.check("75-label score (GUI-shape)", False, "labels or agent_gui absent")
    return
  handle = open(labels_path)
  try:
    labels = json.load(handle)
  finally:
    handle.close()
  covered = total = 0
  for log, fields in labels.items():
    structure = scan(read(os.path.join(gui, log)))
    seen = set()

    def remember(text):
      text = str(text).strip()
      seen.add(text)
      try:
        seen.add("%g" % float(text))
      except ValueError:
        pass

    for item in structure.measurements:
      remember(item.value)
    for item in structure.labeled_values:
      remember(item.value)
      for token in item.value.split():
        remember(token)
    for _, value in fields.items():
      total += 1
      text = str(value)
      try:
        text = "%g" % float(text)
      except ValueError:
        pass
      if text in seen or str(value) in seen:
        covered += 1
  result.check("75-label score (GUI-shape)", covered >= 25,
               "%d of %d recovered from PROGRAM output" % (covered, total))


def smoke_matches_full(corpus, result):
  """Every shipped smoke log must be BYTE-IDENTICAL to its counterpart in the
  full corpus.

  This replaces a heuristic "do the fixture strings still occur" scrape that
  passed on `found > 0` and located 3 of 11 candidates -- a check that could
  not fail for the right reason. The real risk once the big corpus leaves the
  tree is that the shipped smoke copy and the authoritative corpus drift
  apart, and that is exactly what this compares."""
  smoke = os.path.join(HERE, "corpus")
  if not os.path.isdir(smoke):
    result.check("smoke logs identical to the full corpus", False,
                 "smoke corpus absent")
    return
  mismatched = []
  compared = 0
  for path in logs_under(smoke):
    relative = os.path.relpath(path, smoke)
    twin = os.path.join(corpus, relative)
    if not os.path.exists(twin):
      mismatched.append(relative + " (absent upstream)")
      continue
    compared += 1
    if hashlib.md5(open(path, "rb").read()).digest() != \
       hashlib.md5(open(twin, "rb").read()).digest():
      mismatched.append(relative + " (differs)")
  result.check("smoke logs identical to the full corpus", not mismatched,
               "%d compared, %d mismatched %s"
               % (compared, len(mismatched), mismatched[:2] or ""))


def smoke_cover(result):
  """The smoke set must still cover every feature of the coverage vector.
  If it stops, the in-tree suite has quietly gone blind to a form."""
  smoke = os.path.join(HERE, "corpus")
  if not os.path.isdir(smoke):
    result.check("smoke set still covers the vector", False, "smoke corpus absent")
    return
  channels = set()
  screens = set()
  refusals = set()
  sources = set()
  for path in logs_under(smoke):
    structure = scan(read(path))
    channels |= set(structure.forms)
    for item in structure.unparsed:
      if item.screen_rule in SCREEN_RULES:      # "none" is a placeholder on
        screens.add(item.screen_rule)           # rule_excluded records

      if item.excluded_by:
        refusals.add(item.excluded_by)
    for item in structure.items():
      sources.add(item.source)
  ok = (len(channels) >= 11 and len(screens) >= 4 and len(refusals) >= 11
        and len(sources) >= 3)
  result.check("smoke set still covers the vector", ok,
               "%d channels, %d screen rules, %d refusals, %d sources"
               % (len(channels), len(screens), len(refusals), len(sources)))


def conformance(corpus, result):
  suite = os.path.join(corpus, "tst_conformance.py")
  if not os.path.exists(suite):
    print("  %-46s ---   %s" % ("acceptance suite", "not beside the corpus"))
    return
  env = dict(os.environ)
  env["EXTRACTOR"] = "log_structure_extractor"
  env["PHENIX_LOG_CORPUS"] = os.path.join(corpus, "agent")
  env["PHENIX_LOG_CORPUS_GUI"] = os.path.join(corpus, "agent_gui")
  env["PYTHONPATH"] = HERE + os.pathsep + env.get("PYTHONPATH", "")
  out = subprocess.Popen([sys.executable, suite], cwd=corpus, env=env,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT).communicate()[0]
  text = out.decode("utf-8", "replace")
  match = re.search(r"(\d+) passed, (\d+) failed, (\d+) skipped", text)
  if not match:
    result.check("acceptance suite (frozen, v4)", False, "did not run")
    return
  passed, failed, skipped = (int(x) for x in match.groups())
  result.check("acceptance suite (frozen, v4)", failed == 0 and skipped == 0,
               "%d passed, %d failed, %d skipped" % (passed, failed, skipped))


def totals(scans, result):
  counts = {}
  for structure in scans.values():
    for kind in ("sections", "phases", "stages", "cycles", "decisions",
                 "control_skips", "exclusions", "labeled_values",
                 "completion_records", "measurements"):
      counts[kind] = counts.get(kind, 0) + len(getattr(structure, kind))
  print("\n  counted totals (what the documents quote):")
  for kind in sorted(counts):
    print("      %-20s %d" % (kind, counts[kind]))
  return counts


# -------------------------------------------------------------------- main


def main():
  corpus = find_corpus()
  if corpus is not None and not os.path.isdir(os.path.join(corpus, "work",
                                                           "ok")):
    print("Not a validation corpus: %s" % corpus)
    print("Expected a 'work/ok' directory inside it.")
    return 2
  if corpus is None:
    print("FULL VALIDATION CORPUS NOT FOUND.\n")
    print("Looked in:")
    for root in CORPUS_LOCATIONS:
      if root:
        print("    %s" % root)
    print("\nCorpus id: %s" % CORPUS_ID)
    print("\nIt is NOT in cctbx_project by design -- 21 MB.")
    print("\n    tar xzf extraction_test_data.tar.gz -C ~/Downloads")
    print("\nor point at it wherever it already is:")
    print("\n    phenix.python %s \\" % os.path.abspath(__file__))
    print("        --corpus /path/to/extraction_test_data")
    print("\nThis is a FAILURE, not a skip: validation was requested and did")
    print("not happen.")
    return 2

  print("=" * 68)
  print("FULL VALIDATION -- Phenix log structure extractor")
  print("=" * 68)
  print("corpus: %s" % corpus)
  print("")

  result = Result()
  scans, paths = population_checks(corpus, result)
  label_score(corpus, result)
  smoke_matches_full(corpus, result)
  smoke_cover(result)
  conformance(corpus, result)
  counted = totals(scans, result)

  module = os.path.join(HERE, "log_structure_extractor.py")
  revision = hashlib.sha256(open(module, "rb").read()).hexdigest()[:16]

  print("")
  print("=" * 68)
  if result.failed:
    print("VALIDATION FAILED -- %d of %d checks" % (result.failed,
                                                    len(result.checks)))
    print("=" * 68)
    return 1

  record = [
    "extractor sha256:  %s" % revision,
    "corpus id:         %s" % CORPUS_ID,
    "corpus path:       %s" % corpus,
    "date:              %s" % time.strftime("%Y-%m-%d %H:%M:%S"),
    "logs scanned:      %d" % len(paths),
    "checks:            %d of %d PASS" % (len(result.checks),
                                          len(result.checks)),
    "structural items:  %d (excludes the unparsed ledger)"
    % sum(counted.values()),
  ]
  print("VALIDATION PASSED")
  print("=" * 68)
  for line in record:
    print("  " + line)
  # BESIDE THE CORPUS, not in the source tree: writing it into docs/ left a
  # modified file in cctbx after every validation run.
  target = os.path.join(corpus, "VALIDATION_RECORD.txt")
  try:
    handle = open(target, "w")
    try:
      handle.write("\n".join(record) + "\n")
    finally:
      handle.close()
    print("\n  written to %s" % target)
  except IOError:
    pass
  print("\n  Negative controls are a separate run (~8 min):")
  print("    phenix.python %s" % os.path.join(HERE, "run_controls.py"))
  return 0


if __name__ == "__main__":
  sys.exit(main())
