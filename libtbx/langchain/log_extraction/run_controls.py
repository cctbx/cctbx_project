#!/usr/bin/env python
"""Negative controls for the P0 harness.

Run with:
  PHENIX_LOG_CORPUS=... PHENIX_LOG_CORPUS_AGENT=... python3 run_controls.py

Each control disables ONE feature, re-runs the whole suite in a fresh
interpreter, and reports pass/fail. Without the second number, "30 passed"
only shows that the tests agree with themselves.

A control that changes NOTHING is itself a finding: it means no test actually
exercises that feature.
"""

from __future__ import absolute_import, division, print_function

import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))          # .../log_extraction
LANGCHAIN = os.path.dirname(HERE)                          # .../langchain
TESTS = os.path.join(LANGCHAIN, "tests")

# (name, what it breaks, (pattern, replacement) applied to log_structure_extractor.py)
CONTROLS = [
  ("C1 screen disabled",
   "screen_line always returns None -- nothing is ever a candidate",
   ('def screen_line(text):',
    'def screen_line(text):\n  return None\n\ndef _screen_line_disabled(text):')),

  ("C2 regions disabled",
   "find_regions returns [] -- wrapper text is indistinguishable",
   ('def find_regions(lines):\n  """-> list of Region, 1-based inclusive. Pure; takes the split lines."""',
    'def find_regions(lines):\n  """disabled"""\n  return []\n\n\ndef _find_regions_disabled(lines):')),

  ("C3 off-by-one line numbers",
   "items cite index instead of index+1 -- the span verifier must catch it",
   ('  for index, raw in enumerate(lines):\n    line_no = index + 1\n'
    '    rule = screen_line(raw)',
    '  for index, raw in enumerate(lines):\n    line_no = max(1, index)\n'
    '    rule = screen_line(raw)')),

  ("C4 source always program",
   "wrapper regions no longer mark their items -- uncertainty becomes a claim",
   ('def _source_at(regions, line):\n  for r in regions:',
    'def _source_at(regions, line):\n  return SOURCE_PROGRAM\n\n\n'
    'def _source_at_disabled(regions, line):\n  for r in regions:')),

  ("C13 sections disabled",
   "no section is ever found -- P1's whole output",
   ('def find_sections(lines):', 'def find_sections(lines):\n'
    '  return [], {}, []\n\n\ndef _find_sections_disabled(lines):')),

  ("C14 form B disabled",
   "only titles embedded in a rule count -- the prototype's blind spot",
   ('    if index + 1 < n and _SECTION_RULE_RE.match(lines[index + 1]):',
    '    if False and index + 1 < n:')),

  ("C15 banner rules underline titles",
   "* delimiters treated as underlines -- 824 junk titles",
   ('SECTION_UNDERLINE_CHARS = set("-=")',
    'SECTION_UNDERLINE_CHARS = set("-=*#_")')),

  ("C16 labeled values ignore other parsers' claims",
   "a section title is captured a second time as a labeled value",
   ('    if line_no in claimed:\n      continue\n'
    '    if screen_line(raw) != SCREEN_KEY_VALUE:',
    '    if False:\n      continue\n'
    '    if screen_line(raw) != SCREEN_KEY_VALUE:')),

  ("C23 Starting is treated as a phase",
   "the reinstated defect: `Starting CC of ligand` becomes a phase",
   ('_PHASE_RE = re.compile(r"^\\s*Running\\s+(\\S.*?)\\s*$")',
    '_PHASE_RE = re.compile(r"^\\s*(?:Running|Starting)\\s+(\\S.*?)\\s*$")')),

  ("C24 control skips merged into exclusions",
   "the complaint case turns into a phase that never ran",
   ('    control = _CONTROL_SKIP_RE.match(rest)', '    control = None')),

  ("C25 exclusions require a name",
   "3 of find_reference's 18 vanish",
   ('_ITEM_SKIP_RE = re.compile(r"^\\s*(\\S.*?)?\\s*-\\s+(\\S.*?)\\s*$")',
    '_ITEM_SKIP_RE = re.compile(r"^\\s*(\\S.*?)\\s*-\\s+(\\S.*?)\\s*$")')),

  ("C26 unrecognised skips are dropped",
   "234 lines vanish instead of being reported",
   ('    refusals.append((line_no, raw, REFUSAL_SKIP_UNRECOGNISED))',
    '    pass')),

  ("C27 the end summary row counts as a stage",
   "29 stages and a final r_free of 0.4942 -- the number that hides the finding",
   ('      if name == STAGE_SUMMARY_ROW:', '      if False:')),

  ("C28 stage rows are matched globally",
   "rows from outside the table join it -- 40 instead of 29",
   ('    position = index + 1\n    while position < n:\n'
    '      row = _STAGE_ROW_RE.match(lines[position])\n      if not row:\n'
    '        break',
    '    position = index + 1\n    while position < n:\n'
    '      row = _STAGE_ROW_RE.match(lines[position])\n      if not row:\n'
    '        position += 1\n        continue')),

  ("C29 the 999.90 sentinel is quoted as an R-factor",
   "a nonsense R-factor a consumer could quote",
   ('      if work == CYCLE_SENTINEL or free == CYCLE_SENTINEL:\n'
    '        sentinel = True', '      if False:\n        sentinel = True')),

  ("C30 R/Rfree is not split",
   "one key two values collapsed to none",
   ('    pair = _CYCLE_PAIR_RE.search(rest)', '    pair = None')),

  ("C31 a bare cycle counter becomes a record",
   "`Cycle 3 of morphing chain trace` reported as a cycle with no metrics",
   ('    if not metrics and not sentinel:', '    if False:')),

  ("C32 decisions disabled",
   "the program's own branch announcements vanish",
   ('def find_decisions(lines):', 'def find_decisions(lines):\n'
    '  return [], {}, []\n\n\ndef _find_decisions_disabled(lines):')),

  ("C33 reasonless settings become decisions",
   "1,075 values with no branch enter the highest-value channel",
   ('    decision = _DECISION_RE.match(rest)\n    if not decision:',
    '    decision = _DECISION_RE.match(rest)\n    if False:')),

  ("C34 reasonless settings are dropped",
   "1,075 lines vanish instead of being reported",
   ('        reason = REFUSAL_SETTING_NARRATION\n'
    '      refusals.append((line_no, raw, reason))\n      continue',
    '        reason = REFUSAL_SETTING_NARRATION\n      continue')),

  ("C35 all refused settings share one reason",
   "the misleading diagnosis found in the P5 review, reinstated",
   ('      if _SETTING_UP_RE.match(rest):\n        reason = REFUSAL_SETTING_UP',
    '      if False:\n        reason = REFUSAL_SETTING_UP')),

  ("C36 completion records disabled",
   "the only content-derived outcome signal vanishes",
   ('def find_completion_records(lines):', 'def find_completion_records(lines):\n'
    '  return [], {}\n\n\ndef _find_completion_records_disabled(lines):')),

  ("C37 no clustering",
   "a three-line ending counts as three records",
   ('    joins = (events\n             and not names_child',
    '    joins = (False\n             and not names_child')),

  ("C38 top_level without evidence",
   "every last record claims to be the run's own ending",
   ('    elif position == len(events) - 1 and tail_is_blank:',
    '    elif position == len(events) - 1:')),

  ("C39 Finished with X absorbs what follows",
   "the P6 defect: a child ending swallows the run's own",
   ('             and not _FINISHED_WITH_RE.match(events[-1][0][1])',
    '             and True')),

  ("C40 identification disabled",
   "abstains on every log -- passes the hard gate, useless",
   ('def identify_program(lines):', 'def identify_program(lines):\n'
    '  return Identification()\n\n\ndef _identify_program_disabled(lines):')),

  ("C41 degenerate banners are trusted",
   "`Starting phenix` and `Starting libtbx.start_process` become names",
   ('BANNER_BLACKLIST = frozenset(["phenix", "libtbx.start_process", "job"])',
    'BANNER_BLACKLIST = frozenset()')),

  ("C42 the banner takes the rest of the line",
   "the P7 latent defect: `Starting CC of ligand ...` becomes a program",
   ('_BANNER_RE = re.compile(r"^[\\s*]*Starting\\s+(\\S+?)\\s*\\**\\s*$")',
    '_BANNER_RE = re.compile(r"^[\\s*]*Starting\\s+(\\S.*?)\\s*$")')),

  ("C43 the agent prefix hides the banner",
   "9 of 20 wave-1 logs lose their identification",
   ('    match = _BANNER_RE.match(strip_agent_prefix(raw))',
    '    match = _BANNER_RE.match(raw)')),

  ("C44 agent measurements disabled",
   "the labeled test set we already own is not extracted",
   ('def find_agent_measurements(lines):',
    'def find_agent_measurements(lines):\n  return [], {}\n\n\n'
    'def _find_agent_measurements_disabled(lines):')),

  ("C45 the opening rule closes the block",
   "the P8 defect: the channel silently emits nothing at all",
   ('      if seen_any:\n        inside = False', '      inside = False')),

  ("C46 agent measurements hardcode source=agent",
   "a mid-file quoted block is asserted to be the agent's",
   ('      source=_source_at(structure.regions, line_no)))\n\n'
    '  # attached numbers',
    '      source=SOURCE_AGENT))\n\n  # attached numbers')),

  ("C47 attached measurements claim literal evidence",
   "the span verifier is made to assert a derived name is on the line",
   ('    return self.name if self.source == SOURCE_AGENT else ""',
    '    return self.name')),

  ("C48 terminal failures disabled",
   "an aborted run reports only a child completion, as it used to",
   ('def find_terminal_failures(lines):', 'def find_terminal_failures(lines):\n'
    '  return [], {}\n\n\ndef _find_terminal_failures_disabled(lines):')),

  ("C49 the Sorry remedy is dropped",
   "the most useful text in a failed log is thrown away",
   ('      remedy = "\\n".join(l for l in lines[index + 1:end + 1]).strip() or None',
    '      remedy = None')),

  ("C50 the preamble is not a region",
   "`Texas Engineering Experiment Station` is a section again",
   ('def find_preamble(lines):', 'def find_preamble(lines):\n  return None\n\n\n'
    'def _find_preamble_disabled(lines):')),

  ("C51 cap truncation is silent",
   "a Sorry block cut by the cap no longer says so",
   ('      truncated = (end - index) >= SORRY_MAX_LINES',
    '      truncated = False')),

  ("C52 compound pairs are not split",
   "`R: 0.42  Rfree: 0.48` buries the second number again",
   ('    pairs = split_compound_pairs(raw.strip())', '    pairs = []')),

  ("C53 the whole-line test for pairs is dropped",
   "`Space group: C 1 2 1 (No. 5)` becomes `C`",
   ('  if re.sub(r"[\\s:]", "", text) != packed:\n    return []',
    '  if False:\n    return []')),

  ("C54 grep searches the ledger, not the raw lines",
   "the very lines --grep exists to find are the ones it would miss",
   ('      hits = [(number, text)\n'
    '              for number, text in enumerate(lines_of(structure_text), 1)\n'
    '              if lowered in text.lower()]',
    '      hits = [(u.line, u.text) for u in structure.unparsed\n'
    '              if lowered in u.text.lower()]')),

  ("C22 the ledger stops crediting parsers",
   "claimed lines are reported as unparsed again -- unparsed never falls",
   ('    if line_no in claimed:\n      continue\n    structure.unparsed.append',
    '    if False:\n      continue\n    structure.unparsed.append')),

  ("C17 labeled values disabled",
   "the long tail loses its only channel",
   ('def find_labeled_values(lines, claimed, limit=LABEL_DISTINCT_LIMIT):',
    'def find_labeled_values(lines, claimed, limit=LABEL_DISTINCT_LIMIT):\n'
    '  return [], {}, []\n\n\ndef _find_labeled_values_disabled('
    'lines, claimed, limit=LABEL_DISTINCT_LIMIT):')),

  ("C18 identical pairs do not collapse",
   "the same fact stated twice is reported twice",
   ('      if value in seen:', '      if False:')),

  ("C19 runaway labels are not capped",
   "`Target left` x1999 enters the channel unchecked",
   ('    if len(distinct) > limit:', '    if False:')),

  ("C20 generic_only folded back into unclaimed",
   "the anti-gaming split collapses to one number again",
   ('        text=raw, line=line_no, screen_rule=rule, status=GENERIC_ONLY,',
    '        text=raw, line=line_no, screen_rule=rule, status=UNCLAIMED,')),

  ("C21 section lookup uses a stateful cursor",
   "the P2 defect: lookups out of ascending order all return None",
   ('    position = bisect.bisect_right(section_starts, line_no) - 1',
    '    position = len(section_starts) - 1')),

  ("C5 is_flat counts unparsed",
   "a log that only yielded 'I could not read this' reads as structured",
   ('    return not [k for k in ITEM_KINDS if k != "unparsed" and getattr(self, k)]',
    '    return not [k for k in ITEM_KINDS if getattr(self, k)]')),

  ("C6 identification gates extraction",
   "a program_name hint suppresses items -- the layering invariant broken",
   ('  lines = split_lines(text)\n  structure = LogStructure(n_lines=len(lines))',
    '  lines = split_lines(text)\n  if program_name:\n    lines = lines[:1]\n'
    '  structure = LogStructure(n_lines=len(lines))')),

  ("C7 items() unsorted",
   "items no longer walk in line order -- spatial grouping breaks",
   ('    return sorted(out, key=lambda it: (it.line, ITEM_KINDS.index(it.kind)))',
    '    return out')),

  ("C8 line validation removed",
   "Item accepts line=0, negatives and bools",
   ('    if isinstance(line, bool) or not isinstance(line, int) or line < 1:\n'
    '      raise ValueError("line must be a 1-based int, got %r" % (line,))',
    '    line = line if isinstance(line, int) and line >= 1 else 1')),

  ("C10 agent prefix not stripped",
   "strip_agent_prefix returns the line unchanged, prefix and all",
   ('  if text.startswith(AGENT_LINE_PREFIX):\n'
    '    return text[len(AGENT_LINE_PREFIX):].lstrip(" ")\n  return text',
    '  return text')),

  ("C11 agent prefix line swallowed by the header",
   "the program's own banner is marked agent again -- the P0 review defect",
   ('    regions.append(Region("agent_prefix", i + 1, i + 1, SOURCE_PROGRAM))',
    '    regions.append(Region("agent_prefix", i + 1, i + 1, SOURCE_AGENT))')),

  ("C12 splitlines instead of newline splitting",
   "line numbers drift on the one corpus log carrying bare CR characters",
   ('  lines = text.split("\\n")\n  if lines and lines[-1] == "":\n    lines.pop()',
    '  lines = text.splitlines()\n  if False:\n    lines.pop()')),

  ("C9 source validation removed",
   "Item accepts any string as a source",
   ('    if source not in SOURCES:\n'
    '      raise ValueError("source must be one of %r, got %r" % (SOURCES, source))',
    '    pass')),
]

RESULT = re.compile(r"(\d+) passed, (\d+) failed, (\d+) total")

# Controls run against every STRIDEth log. A control has to DISCRIMINATE, not
# to measure -- every reported number comes from the full corpus. Set to 1 to
# run controls at full size.
STRIDE = int(os.environ.get("PHENIX_CONTROL_STRIDE", "1"))  # the smoke corpus is 19 logs; sampling it is meaningless


def run(directory, stride=None):
  env = dict(os.environ)
  if stride:
    env["PHENIX_LOG_CORPUS_STRIDE"] = str(stride)
  # --count, NOT the house fail-fast entry point. A control result is the
  # DIFFERENCE a disabled feature makes, so the magnitude matters: "24 passed,
  # 10 failed" says more than "failed". Fail-fast stops at the first one.
  process = subprocess.Popen(
    [sys.executable,
     os.path.join(directory, "tests", "tst_log_structure_extractor.py"),
     "--count"],
    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
  out = process.communicate()[0].decode("utf-8", "replace")
  match = RESULT.search(out)
  if not match:
    return None, None, out
  return int(match.group(1)), int(match.group(2)), out


def _run_one(job):
  """One control in its own copy of the tree. Controls are independent, so
  they run in parallel -- serial, 21 controls x a 16 s suite is over five
  minutes and the runner stops being something anyone runs."""
  name, why, (pattern, replacement) = job
  tmp = tempfile.mkdtemp()
  try:
    # mirror the langchain/ layout so the test module's PROJECT_ROOT
    # resolves and `from log_extraction...` imports the PATCHED copy
    target = os.path.join(tmp, "langchain")
    os.makedirs(target)
    shutil.copytree(HERE, os.path.join(target, "log_extraction"))
    shutil.copytree(TESTS, os.path.join(target, "tests"))
    path = os.path.join(target, "log_extraction",
                        "log_structure_extractor.py")
    handle = open(path)
    try:
      source = handle.read()
    finally:
      handle.close()
    if pattern not in source:
      return name, None, None, "patch target not found"
    handle = open(path, "w")
    try:
      handle.write(source.replace(pattern, replacement, 1))
    finally:
      handle.close()
    passed, failed, out = run(target, stride=STRIDE)
    if passed is None:
      return name, None, None, "suite did not run"
    return name, passed, failed, None
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def main():
  # FOUND IN REVIEW: this required an env var that nothing else needs any
  # more. The smoke corpus ships with the package and corpus_paths finds it,
  # so the check is now "is there a corpus at all", not "did you export".
  sys.path.insert(0, LANGCHAIN)
  from log_extraction import corpus_paths
  if not corpus_paths.working_corpus():
    print("No corpus found -- controls would run against fixtures only,"
          " which is exactly the weak measurement they exist to avoid.\n")
    print(corpus_paths.missing_corpus_message())
    return 1

  baseline_pass, baseline_fail, baseline_out = run(LANGCHAIN, stride=STRIDE)
  print("baseline (corpus stride %d)   %2d passed, %d failed"
        % (STRIDE, baseline_pass, baseline_fail))
  if baseline_fail:
    # FOUND IN P8: at stride 12 the baseline was 125/2 and the runner
    # reported every control anyway -- 47 numbers measured against a broken
    # reference. A control result is a DIFFERENCE from the baseline, so a red
    # baseline makes all of them meaningless. Refuse rather than mislead.
    print("\nBASELINE IS NOT CLEAN -- %d test(s) already failing.\n"
          "Control results are differences from the baseline, so every"
          " number below would be measured against a broken reference.\n"
          "Refusing to run. Fix the baseline first:" % baseline_fail)
    for line in baseline_out.splitlines():
      if "FAILED" in line or "ERROR" in line:
        print("   " + line.strip())
    return 1
  print("-" * 72)

  weak = []
  try:
    from multiprocessing import Pool
    pool = Pool(4)
    results = pool.map(_run_one, CONTROLS)
    pool.close()
    pool.join()
  except Exception:                       # fall back to serial
    results = [_run_one(job) for job in CONTROLS]

  order = dict((name, position) for position, (name, _, _)
               in enumerate(CONTROLS))
  for name, passed, failed, problem in sorted(results,
                                              key=lambda r: order[r[0]]):
    if problem:
      print("%-44s SKIPPED -- %s" % (name, problem))
      weak.append(name)
      continue
    flag = "" if failed > baseline_fail else "   <-- NO EFFECT"
    if flag:
      weak.append(name)
    print("%-44s %2d passed, %d failed%s" % (name, passed, failed, flag))

  print("-" * 72)
  if weak:
    print("CONTROLS WITH NO EFFECT (untested features): %s" % ", ".join(weak))
    return 1
  print("every control moved the numbers: no P0 feature is untested")
  return 0


if __name__ == "__main__":
  sys.exit(main())
