"""Deterministic reader for Phenix log files.

Lives at libtbx/langchain/log_structure_extractor.py. NOT `log_extractor`:
that name collides with the v113 suite already registered in
tests/run_all_tests.py, and two modules of one name shadow each other on
sys.path with no warning.

  from log_structure_extractor import scan
  structure = scan(open(path, "rb").read().decode("utf-8", "replace"))

WHAT THIS ANSWERS: what did this run actually do? It returns structure and
measurements -- sections, phases, stages, cycles, decisions, skips,
exclusions, labeled values, completion records -- each carrying the line it
came from. It reads what the program says about itself. It does not
interpret, rank, diagnose, or call a model; that is the consumer's job.

THREE LAYERS, KEPT APART BY CONSTRUCTION:

  A  generic observation -- needs no header, no program identity, no section
  B  program identification -- reads whatever evidence exists, may abstain
  C  program-specific interpretation -- out of scope

scan() NEVER REQUIRES PROGRAM IDENTIFICATION. Failure to identify a program
cannot suppress any otherwise extractable structure. Headerlessness is a
normal state: a 40-line headerless log yielding three cited facts and no
program name is a successful extraction. This is enforced by a corpus-level
test and a negative control, not by intention.

STATUS: phase P0. The types, the line discipline, the wrapper-region map and
the frozen candidate screen are implemented. NO CONTENT PARSER EXISTS YET, so
every screened candidate is reported as unparsed/unclaimed. That is the
intended P0 output: it establishes the denominator before any parser can
shrink it.

Standard library only. Runs under python and libtbx.python. scan() is pure:
no file I/O, no global state, no logging. File handling lives in main().
"""

from __future__ import absolute_import, division, print_function

import bisect
import re
import sys

__version__ = "0.1.0-P0"

# --------------------------------------------------------------- item types
#
# Every item carries a 1-based `line`. The consumer cites evidence by span so
# a human can check it; an item without provenance is unusable. `source` says
# who wrote the text and `unknown` is a real value -- see REGIONS below.
# `section_id` is optional EVERYWHERE and None is normal, so that headerless
# programs are not second-class inputs.

SOURCE_PROGRAM = "program"
SOURCE_AGENT = "agent"
SOURCE_UNKNOWN = "unknown"
SOURCES = (SOURCE_PROGRAM, SOURCE_AGENT, SOURCE_UNKNOWN)


class Item(object):
  """Base for everything scan() returns. Comparable and repr-able so tests can
  assert on whole structures rather than on field-by-field spot checks."""

  kind = "item"
  fields = ()

  def __init__(self, line, source=SOURCE_PROGRAM, section_id=None, end=None):
    # bool is a subclass of int, so True would silently become line 1
    if isinstance(line, bool) or not isinstance(line, int) or line < 1:
      raise ValueError("line must be a 1-based int, got %r" % (line,))
    if source not in SOURCES:
      raise ValueError("source must be one of %r, got %r" % (SOURCES, source))
    self.line = line
    self.end = end if end is not None else line
    self.source = source
    self.section_id = section_id

  def _key(self):
    return (self.kind, self.line, self.end, self.source, self.section_id,
            tuple(getattr(self, f) for f in self.fields))

  def __eq__(self, other):
    return isinstance(other, Item) and self._key() == other._key()

  def __ne__(self, other):
    return not self.__eq__(other)

  def __hash__(self):
    return hash(self._key())

  def __repr__(self):
    parts = ["line=%d" % self.line]
    if self.end != self.line:
      parts.append("end=%d" % self.end)
    for f in self.fields:
      parts.append("%s=%r" % (f, getattr(self, f)))
    if self.source != SOURCE_PROGRAM:
      parts.append("source=%r" % self.source)
    return "%s(%s)" % (type(self).__name__, ", ".join(parts))

  def evidence(self):
    """The token that must appear on the cited line. The span verifier checks
    this against the log text; it is what catches an item citing the wrong
    line, or the wrong log, and validating cleanly anyway."""
    return ""


class Section(Item):
  kind = "sections"
  fields = ("title",)

  def __init__(self, title, line, end=None, **kw):
    Item.__init__(self, line, end=end, **kw)
    self.title = title

  @property
  def start(self):
    """Alias for `line`. Requirements 3 names a Section's first line `start`;
    every other item calls it `line`. Both are provided so the frozen
    conformance suite and the uniform item API agree without either moving."""
    return self.line

  def evidence(self):
    return self.title


class Phase(Item):
  kind = "phases"
  fields = ("name",)

  def __init__(self, name, line, **kw):
    Item.__init__(self, line, **kw)
    self.name = name

  def evidence(self):
    return self.name


class Stage(Item):
  kind = "stages"
  fields = ("name", "metrics", "description")

  def __init__(self, name, line, metrics=None, description=None, **kw):
    Item.__init__(self, line, **kw)
    self.name = name
    self.metrics = dict(metrics or {})
    self.description = description

  def _key(self):
    return (self.kind, self.line, self.end, self.source, self.section_id,
            self.name, tuple(sorted(self.metrics.items())), self.description)

  def evidence(self):
    return self.name


class Cycle(Item):
  kind = "cycles"
  fields = ("number", "metrics", "sentinel")

  def __init__(self, number, line, metrics=None, sentinel=False, **kw):
    Item.__init__(self, line, **kw)
    self.number = number
    self.metrics = dict(metrics or {})
    self.sentinel = bool(sentinel)

  def _key(self):
    return (self.kind, self.line, self.end, self.source, self.section_id,
            self.number, tuple(sorted(self.metrics.items())), self.sentinel)


class Decision(Item):
  kind = "decisions"
  fields = ("setting", "value", "reason")

  def __init__(self, setting, value, reason, line, **kw):
    Item.__init__(self, line, **kw)
    self.setting = setting
    self.value = value
    self.reason = reason

  def evidence(self):
    return self.setting


class ControlSkip(Item):
  """A phase that did not run. NOT the same thing as an Exclusion -- merging
  the two either invents a phase that never existed or loses the item."""

  kind = "control_skips"
  fields = ("what", "reason")

  def __init__(self, what, reason, line, **kw):
    Item.__init__(self, line, **kw)
    self.what = what
    self.reason = reason


class Exclusion(Item):
  """A candidate the program rejected. The item name is OPTIONAL: requiring
  one silently dropped 3 of find_reference's 18."""

  kind = "exclusions"
  fields = ("item", "reason")

  def __init__(self, item, reason, line, **kw):
    Item.__init__(self, line, **kw)
    self.item = item
    self.reason = reason

  def evidence(self):
    return self.reason


class Measurement(Item):
  """A number arriving ATTACHED to a structure already parsed -- a stage row,
  a cycle line, the agent metrics block. Generic `Key: value` goes to
  LabeledValue instead, and a line is never emitted as both."""

  kind = "measurements"
  fields = ("name", "value", "unit", "context")

  # `unit` IS PERMANENTLY None IN v1, and is retained rather than removed
  # because the consuming project built against the stub. Populating it
  # means parsing units out of values, which is interpretation (non-goal
  # 2.10). Documented in requirements v2 3.3 rather than quietly dropped.

  def __init__(self, name, value, line, unit=None, context=None, **kw):
    Item.__init__(self, line, **kw)
    self.name = name
    self.value = value
    self.unit = unit
    self.context = context

  def evidence(self):
    """Only the agent block writes a literal `Name: value` line. An ATTACHED
    measurement is a re-emission of a stage row or cycle line whose own span
    is already verified, and its metric name (`r_free` from `R/Rfree=a/b`) is
    DERIVED and appears nowhere on the line. Claiming it as evidence would
    make the span verifier assert something untrue.

    DEFECT FOUND IN P8: the first version returned the name unconditionally
    and the corpus span verifier failed on `r_free` from a cycle line."""
    return self.name if self.source == SOURCE_AGENT else ""


class LabeledValue(Item):
  """A `Key: value` record, uninterpreted. No unit normalisation, no type
  inference, no mapping of `R-free`/`R Free`/`r_free` onto one name -- that
  mapping is the maintenance treadmill this tool exists to escape.

  repeat_count > 1 means the same label appeared repeatedly in this log and
  was collapsed; `line` is the first occurrence and `end` the last."""

  kind = "labeled_values"
  fields = ("label", "value", "repeat_count")

  def __init__(self, label, value, line, repeat_count=1, end=None, **kw):
    Item.__init__(self, line, end=end, **kw)
    self.label = label
    self.value = value
    self.repeat_count = repeat_count

  def evidence(self):
    return self.label


class CompletionRecord(Item):
  """A record that something finished. NEVER a claim that the run succeeded:
  of 30 failed runs in corpus2, 27 contain no error keyword at all in what
  the program itself wrote, so absence of error is not evidence of success.

  applies_to is `top_level` only with positive evidence; otherwise `unknown`.
  One corpus log carries 59 of these (phaser writes EXIT STATUS per module),
  so a singular `termination` field would pick one of 59 and call it the
  outcome of the run."""

  kind = "completion_records"
  fields = ("text", "applies_to")

  TOP_LEVEL = "top_level"
  UNKNOWN = "unknown"

  def __init__(self, text, line, applies_to=UNKNOWN, **kw):
    Item.__init__(self, line, **kw)
    self.text = text
    self.applies_to = applies_to

  def evidence(self):
    return self.text.strip()


class TerminalFailure(Item):
  """The program's own account of why it stopped: a traceback, a `Sorry:`
  diagnosis with the remedy that follows it, or an explicit abort marker.

  **NOT a verdict that the run failed** (non-goal 2.8). It reports text the
  program wrote, with line numbers, and a log may carry both this and a
  completion record -- in which case both are shown and the reader sees the
  conflict.

  ADDED after an AutoSol run ended `STOPWIZARD IN AutoSol_run_1_/STOPWIZARD`
  on line 3877 and the report said only `Finished with search...` from line
  1994, 1,883 lines earlier. A reader could have concluded it finished.

  MEASURED on 291 previously unseen raw captures: 25 carry a `Traceback`,
  25 a `Sorry:`, 19 end `Please try again.`, 2 `STOPWIZARD`, 3
  `EndOfResolve` -- 27 logs with some terminal marker (9%). In the GUI-shape
  corpus the same markers appear in 5 of 253 (2%), because there the failure
  went to stderr and was stripped into the agent wrapper. The long-standing
  claim that "27 of 30 failed runs say nothing" is true of GUI-shape logs
  and FALSE of raw captures."""

  kind = "terminal_failures"
  fields = ("failure_kind", "text", "remedy", "truncated")

  TRACEBACK = "traceback"
  SORRY = "sorry"
  ABORT = "abort_marker"

  def __init__(self, failure_kind, text, line, remedy=None, end=None,
               truncated=False, **kw):
    Item.__init__(self, line, end=end, **kw)
    self.failure_kind = failure_kind
    self.text = text
    self.remedy = remedy
    # True when the block hit SORRY_MAX_LINES rather than reaching its own
    # end. Reported, never silent.
    self.truncated = truncated

  def evidence(self):
    return self.text.strip()[:40]


# unparsed statuses
UNCLAIMED = "unclaimed"
GENERIC_ONLY = "generic_only"
RULE_EXCLUDED = "rule_excluded"
UNPARSED_STATUSES = (UNCLAIMED, GENERIC_ONLY, RULE_EXCLUDED)


class Unparsed(Item):
  """Structure seen but not understood. A REQUIREMENT, not a nicety: two
  parser bugs were caught in the prototype solely because unattached items
  were reported rather than dropped.

  Three statuses, always reported together, because a single number is
  gameable from both ends:
    unclaimed      -- no parser claimed it
    generic_only   -- claimed only by labeled_values: recorded, not understood
    rule_excluded  -- refused by a NAMED admission rule, reported with the rule
  """

  kind = "unparsed"
  fields = ("screen_rule", "status", "excluded_by", "text")

  # RENAMED FROM v1's `Unparsed(kind, line, text)`. `kind` means the CHANNEL
  # name on every other item ("sections", "cycles") and is what items()
  # sorts by, so reusing it for "which screen rule matched" is a collision.
  # Declared in requirements v2 3.2 rather than left as a silent deviation.

  def __init__(self, text, line, screen_rule, status=UNCLAIMED,
               excluded_by=None, **kw):
    Item.__init__(self, line, **kw)
    if status not in UNPARSED_STATUSES:
      raise ValueError("status must be one of %r" % (UNPARSED_STATUSES,))
    self.text = text
    self.screen_rule = screen_rule
    self.status = status
    self.excluded_by = excluded_by

  def evidence(self):
    return self.text.strip()[:40]


ITEM_KINDS = ("sections", "phases", "stages", "cycles", "decisions",
              "control_skips", "exclusions", "measurements", "labeled_values",
              "completion_records", "terminal_failures", "unparsed")

# kinds that say WHAT HAPPENED rather than that something began. Reported
# separately from bare reach, because adding section recognition buys 21
# points of headline coverage and little diagnostic value.
DIAGNOSTIC_KINDS = ("stages", "cycles", "decisions", "control_skips",
                    "exclusions", "terminal_failures")
BASIC_KINDS = ("labeled_values", "measurements", "completion_records")


# ------------------------------------------------------------ identification


class Candidate(object):
  """A ranked identification candidate.

  Two evidence kinds, kept apart because they license different confidence:

    self_identification -- a POSITIONAL banner the program wrote about itself
                           (`Starting phenix.xtriage` on line 1). High.
    embedded_program    -- a block or command CONSISTENT WITH some program.
                           Low, and never promoted. A wizard that launches a
                           child immediately can make the two disagree:
                           predict_and_build says `phenix.refine` 38 times and
                           never names itself.
  """

  SELF_IDENTIFICATION = "self_identification"
  EMBEDDED_PROGRAM = "embedded_program"

  def __init__(self, name, confidence, evidence_kind, line):
    self.name = name
    self.confidence = confidence
    self.evidence_kind = evidence_kind
    self.line = line

  def __eq__(self, other):
    return (isinstance(other, Candidate)
            and (self.name, self.confidence, self.evidence_kind, self.line)
            == (other.name, other.confidence, other.evidence_kind, other.line))

  def __ne__(self, other):
    return not self.__eq__(other)

  def __repr__(self):
    return "Candidate(%r, %.2f, %s, line=%d)" % (
      self.name, self.confidence, self.evidence_kind, self.line)


class Identification(object):
  """Ranked candidates plus the signals that fired, or an explicit abstention.

  Abstention is `candidates == []`, never a low-confidence guess. A wrong
  program name is worse than no program name: ai_analysis derived the program
  from the FILENAME, and when derivation failed it carried on and printed a
  confident report about a program that never ran. NEVER use the filename.

  signals_fired is reported because a fingerprint that rots silently is worse
  than one that fails loudly -- a run identifying from one signal where it
  used to use three is a warning that only shows if the signals are named."""

  def __init__(self, candidates=None, signals_fired=None):
    self.candidates = list(candidates or [])
    self.signals_fired = list(signals_fired or [])

  @property
  def name(self):
    return self.candidates[0].name if self.candidates else None

  @property
  def is_unknown(self):
    return not self.candidates

  def __eq__(self, other):
    return (isinstance(other, Identification)
            and self.candidates == other.candidates
            and self.signals_fired == other.signals_fired)

  def __ne__(self, other):
    return not self.__eq__(other)

  def __repr__(self):
    if self.is_unknown:
      return "Identification(unknown, signals=%r)" % (self.signals_fired,)
    return "Identification(%r, signals=%r)" % (self.candidates,
                                               self.signals_fired)


# ------------------------------------------------------------------- regions
#
# Two blocks in some logs are written by the PHENIX AI Agent, not by the
# program: the WORKING DIRECTORY / COMMAND THAT WAS RUN header and the
# ERROR TEXT / MTZ LABEL INFO / FINAL QUALITY METRICS footer. A log from the
# GUI -- the common case in production -- has neither.
#
# Regions are computed BEFORE any parser reads a line, so that no parser's
# behaviour can become accidentally dependent on wrapper text.

HEADER_KEYS = ("WORKING DIRECTORY:", "COMMAND THAT WAS RUN:")

# The agent does not add a line here -- it PREFIXES the program's own first
# line. `LOG TEXT: Starting AutoBuild` is nine characters of agent in front of
# the program's own banner.
#
# DEFECT FOUND IN P0 REVIEW: the first version swallowed this whole line into
# the agent header region. MEASURED: 15 of the 20 agent-shape logs carry the
# prefix and in 9 of them the payload is a `Starting <program>` banner -- so
# the rule buried the single best identification signal inside a region
# marked "written by the agent". P7 would then have had to either ignore the
# marking or re-derive the banner, and one of those is a silent misread.
#
# The line is now its own region whose SOURCE IS THE PROGRAM, because the
# content is the program's; only the prefix is not. Parsers call
# strip_agent_prefix() before reading such a line.
AGENT_LINE_PREFIX = "LOG TEXT:"
FOOTER_KEYS = ("ERROR TEXT:", "MTZ LABEL INFO:",
               "FINAL QUALITY METRICS REPORT:")
_RULE_ONLY = re.compile(r"^[-=*#_]{10,}\s*$")

# A wrapper block sitting in the MIDDLE of a file is ambiguous: it may be the
# agent wrapping this log, or this log QUOTING an agent (ai_agent_*.log files
# quote their children's wrapper blocks verbatim). We cannot tell, so we say
# unknown rather than assert `agent`. Turning absence of contrary evidence
# into a claim is the specific error this project has already made once.
#
# REVISION RECORDED (P0 review, two drafts). Draft 1 ended a footer block at
# the next wrapper key and judged mid-file-ness by how much text FOLLOWED the
# block -- which could never fire, because with no following key the block
# always ran to EOF and `remaining` was always 0. Draft 2 keyed on the start
# position alone and then called a short fixture's genuine footer ambiguous.
#
# The rule that holds on everything measured uses two signals:
#   1. POSITION -- a genuine wrapper footer sits in the back half of the file
#      (corpus2 error originals put ERROR TEXT at line 8041 of 8225);
#   2. REPETITION -- a genuine wrapper writes each key AT MOST ONCE, while a
#      log that QUOTES wrapper blocks repeats them (p9-sad__ai_agent_1 carries
#      8 FINAL QUALITY METRICS blocks, sec17-sad__ai_agent_2 carries 4).
# Either signal failing means we cannot tell, and we say unknown rather than
# assert `agent`. Erring toward "we do not know" is the point.
_FOOTER_TAIL_FRACTION = 0.5


# THE PHENIX PREAMBLE. Every Phenix log opens with a developer list,
# copyright, third-party components and the citation -- identical in every
# log, and not about the run. It is boilerplate, and it was producing
# sections like `Texas Engineering Experiment Station` (a copyright
# continuation line that happens to sit above a rule) and
# `http://www.phenix-online.org/`.
#
# Bounded, marked and excluded from sections and labeled values -- the same
# machinery as the agent wrapper, so the extent stays auditable and nothing
# is silently dropped. `source` stays `program`: the program did write it.
_PREAMBLE_START_RE = re.compile(r"^\s*Phenix developers include:")
_PREAMBLE_END_RE = re.compile(r"Acta Cryst")
PREAMBLE_MAX_LINES = 90


def find_preamble(lines):
  """-> a Region for the Phenix citation preamble, or None."""
  for index, raw in enumerate(lines[:20]):
    if not _PREAMBLE_START_RE.match(raw):
      continue
    end = index
    limit = min(len(lines), index + PREAMBLE_MAX_LINES)
    for position in range(index, limit):
      if _PREAMBLE_END_RE.search(lines[position]):
        end = position
        # the citation runs a few lines past `Acta Cryst`; stop at the blank
        while end + 1 < limit and lines[end + 1].strip():
          end += 1
        break
    else:
      return None
    return Region("phenix_preamble", index + 1, end + 1, SOURCE_PROGRAM)
  return None


class Region(object):
  """Diagnostic, NOT part of the frozen output contract. Exposed because it is
  useful to tools that want to strip or audit wrapper text, but a consumer
  should read `source` on items rather than depend on this."""

  def __init__(self, kind, start, end, source):
    self.kind = kind
    self.start = start
    self.end = end
    self.source = source

  def contains(self, line):
    return self.start <= line <= self.end

  def __repr__(self):
    return "Region(%s, %d-%d, source=%s)" % (self.kind, self.start, self.end,
                                             self.source)


def find_regions(lines):
  """-> list of Region, 1-based inclusive. Pure; takes the split lines."""
  regions = []
  n = len(lines)
  if n == 0:
    return regions

  i = 0
  while i < n and lines[i].startswith(HEADER_KEYS):
    i += 1
  if i > 0:
    regions.append(Region("agent_header", 1, i, SOURCE_AGENT))
  # The prefix does NOT require a header above it. DEFECT FOUND IN P0 REVIEW:
  # 15 of the 20 wave-1 GUI-shape logs begin with a bare `LOG TEXT:` and no
  # header at all -- an artifact of how that corpus was stripped -- and the
  # first version emitted no prefix region for any of them, so the commonest
  # case in the shipped GUI corpus went unmarked.
  if i < n and lines[i].startswith(AGENT_LINE_PREFIX):
    regions.append(Region("agent_prefix", i + 1, i + 1, SOURCE_PROGRAM))

  for j in range(i, n):
    if not lines[j].startswith(FOOTER_KEYS):
      continue
    if any(r.contains(j + 1) for r in regions):
      continue
    start = j + 1
    # a bare rule immediately above the block belongs to the block
    if j > 0 and _RULE_ONLY.match(lines[j - 1]) and not any(
        r.contains(j) for r in regions):
      start = j
    quoted = any(
      sum(1 for line in lines if line.startswith(key)) > 1
      for key in FOOTER_KEYS)
    # position is judged from the KEY line, not from the decorative rule
    # above it, or a short file's genuine footer reads as ambiguous
    in_back_half = (j + 1) > n * _FOOTER_TAIL_FRACTION
    source = SOURCE_AGENT if (in_back_half and not quoted) \
        else SOURCE_UNKNOWN
    regions.append(Region("agent_footer", start, n, source))
    break
  return regions


def strip_agent_prefix(text):
  """Remove the agent's `LOG TEXT:` prefix, leaving the program's own text.

  Returns the text unchanged when there is no prefix, so a parser can call it
  unconditionally. Note that the GUI-shape corpus shipped with wave 1 still
  carries this prefix on line 1 of 15 of its 20 logs -- an artifact of how it
  was stripped -- so tolerating it is not optional."""
  if text.startswith(AGENT_LINE_PREFIX):
    return text[len(AGENT_LINE_PREFIX):].lstrip(" ")
  return text


def split_lines(text):
  """Split on newlines ONLY, the way grep -n, sed -n and an editor do.

  DEFECT FOUND IN P0 REVIEW. str.splitlines() also breaks on bare CR, VT, FF,
  FS, GS, RS, NEL, U+2028 and U+2029. corpus2/work/ok/fem_407.log carries 100
  bare CR characters: splitlines() gives it 539 lines where `grep -n` gives
  440, so every cited line after the first CR pointed a human at the wrong
  place -- in a tool whose entire value is that a claim can be checked.

  A trailing CR from a CRLF file is dropped so the line text matches too."""
  if not text:
    return []
  lines = text.split("\n")
  if lines and lines[-1] == "":
    lines.pop()
  return [line[:-1] if line.endswith("\r") else line for line in lines]


def _source_at(regions, line):
  for r in regions:
    if r.contains(line):
      return r.source
  return SOURCE_PROGRAM


# ------------------------------------------------- the frozen candidate screen
#
# FROZEN IN P0, BEFORE ANY PARSER EXISTS, and deliberately over-selective.
#
# `unparsed` is worthless if the implementer also decides what counts as
# "structure seen but not understood" -- the easy way to make it zero is to
# see less. So the screen is written first and does not move. It will
# over-select; much prose matches `Key: value`. That is the design, and the
# cost is paid in a large `generic_only` count rather than in a flattering
# `unclaimed` one.
#
# Four shapes, and nothing is filtered out here. Refusals happen at ADMISSION,
# where each one must name the rule that refused it (Unparsed.excluded_by), so
# a reviewer can disagree with a rule instead of with a number that quietly
# moved.

SCREEN_RULE_LINE = "rule_line"
SCREEN_KEY_VALUE = "key_value"
SCREEN_VERB = "verb"
SCREEN_NUMERIC_ROW = "numeric_row"
SCREEN_RULES = (SCREEN_RULE_LINE, SCREEN_KEY_VALUE, SCREEN_VERB,
                SCREEN_NUMERIC_ROW)

_SCREEN_RULE_RE = re.compile(r"^\s*[-=*#_]{10,}\s*$")
# SCREEN VERSION 2 (signed off by Tom, P4). `{1,48}` after the leading letter
# meant a label of 2-49 characters, so a SINGLE-CHARACTER LABEL could not
# match. MEASURED on corpus2/work n=253: 124 such lines, and 115 of them --
# across 33 logs, mostly phaser and autosol -- were invisible to the screen
# entirely, not merely unclaimed:
#
#     R:   0.45  Rfree:   0.46
#
# R and R-free, unreadable. Now `{0,48}` = 1-49 characters, which ADDS short
# labels without tightening the upper bound. Measured cost: +124 candidates
# on 107,390 (+0.1%); every unparsed baseline recorded before P4 shifts by
# that much and has been re-recorded.
#
# Two readings were measured and rejected. `{1,48}` (frozen) misses short
# labels; `{0,47}` = 1-48, which is what the old comment CLAIMED, would have
# LOST 213 candidates by tightening the maximum. Correcting code to match its
# own comment would have made things worse.
#
# KNOWN LIMITATION, recorded rather than fixed: `R: 0.26  Rfree: 0.29` is two
# facts on one line and is captured as one -- label "R", value
# "0.26  Rfree: 0.29". Splitting compound lines is interpretation.
SCREEN_VERSION = 2
_SCREEN_KV_RE = re.compile(
  r"^ {0,3}([A-Za-z][A-Za-z0-9 _()/.'-]{0,48}):[ \t]+(\S.*?)\s*$")
# indentation is allowed deliberately: `            Running build_chains on
# segment 1` is a real phase announcement, and refine's `    Cycle  Chain
# plDDT ...` is a real table header. 35 such lines in the first 80 logs.
_SCREEN_VERB_RE = re.compile(r"^\s*(Skipping|Setting|Running|Cycle)\b")
_NUMBER_RE = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


_NUMERIC_START = set("0123456789+-.")


def _is_numeric_row(text, minimum=4):
  """At least `minimum` numeric tokens. NOT a majority of tokens.

  MEASURED IN P0 REVIEW, because the plan's wording ("a row of >=4
  whitespace-separated numeric columns") is ambiguous and a majority reading
  is the obvious one. On corpus2/work n=253 the three readings give:

    >=4 numeric tokens (this)               46,643 candidates
    >=4 and a numeric majority              35,280
    >=4, majority, ignoring separator cells 39,181

  Both stricter readings LOSE GENUINE TABLE ROWS -- `| 4_0 (c) | 0 | 0.00 |`
  has label cells, and `CELL from /path ...` is a real record. Including
  prose that happens to carry four numbers is the cheaper error for a screen
  whose whole job is to over-select, so the loose reading stands and the
  plan's wording is what gets corrected."""
  tokens = text.split()
  if len(tokens) < minimum:
    return False
  # cheap first pass: a numeric token must START with a digit or sign. This
  # is a pure optimisation -- the regex still decides -- and it exists
  # because this function is the hot path (42,426 calls on the largest log,
  # 38% of scan time before it was added).
  numeric = 0
  remaining = len(tokens)
  for token in tokens:
    if token[0] in _NUMERIC_START and _NUMBER_RE.match(token.rstrip(":,")):
      numeric += 1
      if numeric >= minimum:
        return True
    remaining -= 1
    if numeric + remaining < minimum:
      return False
  return numeric >= minimum


def screen_line(text):
  """-> the screen rule this line matches, or None. Pure, order-independent,
  and the first rule that matches wins. FROZEN: changing this changes every
  unparsed number ever measured, so it moves only with a recorded measurement
  and a suite version bump."""
  if _SCREEN_RULE_RE.match(text):
    return SCREEN_RULE_LINE
  if _SCREEN_VERB_RE.match(text):
    return SCREEN_VERB
  if _SCREEN_KV_RE.match(text) and not text.rstrip().endswith(":"):
    return SCREEN_KEY_VALUE
  if _is_numeric_row(text):
    return SCREEN_NUMERIC_ROW
  return None


# ------------------------------------------------------------------ sections
#
# TWO FORMS, both measured on corpus2/work n=253 before this was written:
#
#   A  a title embedded in a rule       158 logs, 219 distinct titles
#      ============= Collecting inputs =============
#      ----------Processing PDB file(s)----------
#
#   B  a title on its own line over a rule   225 logs, 2109 distinct titles
#      Processing files:
#      -------------------------------------------
#
# Form B is the Phenix Program Template's own shape and the prototype missed
# it entirely, which is most of why the prototype found sections in only 15
# of the 20 wave-1 logs.
#
# A BARE RULE WITH NO TITLE IS A SEPARATOR, NOT A SECTION -- the corpus has
# thousands.
#
# Three refusals, each NAMED so that a reviewer can disagree with a rule
# rather than with a number that moved (see Unparsed.rule_excluded):
#   too_long      a title over 60 characters. Requirements 6 names this
#                 explicitly as something to report rather than accept.
#   table_header  a title containing `|`. Measured: `| d_spacing | Expected
#                 rel. I | Data Z-score |` sits above a rule 38 times and is
#                 a table header, not a section.
#   numeric_row   a title that is mostly numbers, for the same reason.
# On corpus2/work these refuse 2,200 / 270 / 97 lines respectively.

# UNDERLINE CHARACTERS. Form B accepts only `-` and `=`.
#
# MEASURED on corpus2/work n=253: of 8,644 form-B candidates, the rule
# beneath is `-` 5,967 times, `=` 1,818, `*` 824, `#` 35. The `*` rules are
# BANNER DELIMITERS, not underlines -- they wrap a block above and below --
# and the titles they would create are the junk: `Unset PHENIX_OVERWRITE_ALL
# if overwriting is not desired`, `hydrogens.refine=riding`. Restricting to
# `-` and `=` refuses 859 of 8,644 (10%).
#
# A blank-line-above rule was measured first and REJECTED: it would refuse
# 4,603 of 8,644 including `OUTPUT FILES` (240), `TRANSLATIONAL NCS` (192),
# `WILSON DISTRIBUTION` (126) and `REFINEMENT` (69), which are unmistakably
# headings. It costs more than it saves.
SECTION_UNDERLINE_CHARS = set("-=")
SECTION_MAX_TITLE = 60
_SECTION_INLINE_RE = re.compile(
  r"^\s*[-=*#_]{3,}\s*(\S.*?)\s*[-=*#_]{3,}\s*$")
_SECTION_RULE_RE = re.compile(r"^\s*[-=*#_]{10,}\s*$")
_ALL_RULE_CHARS_RE = re.compile(r"^[-=*#_\s]+$")

REFUSAL_TOO_LONG = "section_title_too_long"
REFUSAL_TABLE_HEADER = "section_title_is_a_table_header"
REFUSAL_NUMERIC = "section_title_is_a_numeric_row"
REFUSAL_BANNER = "rule_is_a_banner_not_an_underline"


def _refuse_title(title):
  """-> a named refusal, or None if the title is acceptable."""
  text = title.strip()
  if not text or _ALL_RULE_CHARS_RE.match(text):
    return REFUSAL_BANNER
  # ORDER MATTERS, and the more specific diagnosis wins. A 61-character
  # table header is both over-long and a table header; reporting `too_long`
  # would send a reader looking for a truncation problem that is not there.
  if "|" in text:
    return REFUSAL_TABLE_HEADER
  if _is_numeric_row(text):
    return REFUSAL_NUMERIC
  if len(text) > SECTION_MAX_TITLE:
    return REFUSAL_TOO_LONG
  return None


def find_sections(lines):
  """-> (sections, claimed_lines, refusals).

  `claimed_lines` maps a 1-based line number to the kind that claimed it, so
  the unparsed ledger can tell a line we understood from one we did not.
  `refusals` is a list of (line, text, reason)."""
  sections, claimed, refusals = [], {}, []
  n = len(lines)
  for index, raw in enumerate(lines):
    line_no = index + 1
    if _SECTION_RULE_RE.match(raw):
      continue     # a bare rule is a separator, and MUST be tested first:
                   # the inline pattern otherwise matches a pure rule, taking
                   # its middle dashes for a title (defect found in P1)
    match = _SECTION_INLINE_RE.match(raw)
    if match:
      reason = _refuse_title(match.group(1))
      if reason is None:
        sections.append([match.group(1).strip(), line_no, None])
        claimed[line_no] = "sections"
      else:
        refusals.append((line_no, raw, reason))
      continue
    if not raw.strip():
      continue     # a blank line above a rule is not a refused title, it is
                   # not a title at all. DEFECT FOUND IN P1 REVIEW: without
                   # this, every blank line above a rule was emitted as a
                   # rule_excluded record, inflating the refusal count and
                   # breaking the paired-shape invariant on one log.
    if index + 1 < n and _SECTION_RULE_RE.match(lines[index + 1]):
      if not set(lines[index + 1].strip()) <= SECTION_UNDERLINE_CHARS:
        refusals.append((line_no, raw, REFUSAL_BANNER))
        continue
      reason = _refuse_title(raw)
      if reason is None:
        sections.append([raw.strip(), line_no, None])
        claimed[line_no] = "sections"
        claimed[line_no + 1] = "sections"
      else:
        refusals.append((line_no, raw, reason))

  # a section runs until the next one begins, or to the end of the log
  for position, section in enumerate(sections):
    following = sections[position + 1][1] - 1 if position + 1 < len(sections) \
        else n
    section[2] = max(section[1], following)
  return sections, claimed, refusals


# --------------------------------------------------- phases, skips, exclusions
#
# `Running X` IS a phase announcement. **`Starting X` is NOT.** autobuild
# writes `Starting job`, `Starting mtz file`, `Starting model list`;
# ligandfit writes `Starting CC of ligand as input to map: 0.714`, which is a
# metric wearing a verb. An early draft of the requirements had this
# backwards and a survey of the actual lines corrected it. Do not reinstate
# `Starting` without measuring it.
#
# TWO KINDS OF SKIP, AND MERGING THEM COSTS THE PROJECT ITS FLAGSHIP CASE:
#
#   CONTROL SKIP   Skipping model-building as 'maps_only' is set
#                  -- a phase that did not run
#   ITEM EXCLUSION Skipping 9GSD.A - not x-ray and not computational
#                  Skipping - not protein or too few residues.
#                  -- a candidate the program rejected
#
# Eighteen item exclusions in two reason-groups ARE the entire content of the
# user complaint that started the consuming project. Merging the two either
# invents a phase that never existed or loses the complaint. The `as`/
# `because` form is control flow; the dash form is an item. THE ITEM NAME IS
# OPTIONAL -- requiring one silently dropped 3 of find_reference's 18.
#
# MEASURED on corpus2/work n=253 before writing this: 381 `Skipping` lines in
# 39 logs -- 99 as/because, 49 dash, and 233 MATCHING NEITHER, in shapes like
# `Skipping remainder region 1 (already written out)` and `Skipping CORR_RMS
# in CC-scoring...`. Those 233 are REPORTED as unparsed, exactly as
# requirements 6 asks, and NOT fitted with a third pattern: a control-versus-
# item distinction guessed from a corpus I have already read is the kind of
# wrong that is invisible until it costs a claim.

_PHASE_RE = re.compile(r"^\s*Running\s+(\S.*?)\s*$")
_SKIP_RE = re.compile(r"^\s*Skipping\b(.*)$")
_CONTROL_SKIP_RE = re.compile(r"^\s*(\S.*?)\s+\b(?:as|because)\b\s+(\S.*?)\s*$")
_ITEM_SKIP_RE = re.compile(r"^\s*(\S.*?)?\s*-\s+(\S.*?)\s*$")

REFUSAL_SKIP_UNRECOGNISED = "skip_matches_neither_form"


def find_phases_and_skips(lines):
  """-> (phases, control_skips, exclusions, claims, refusals).

  Each entry is a tuple ending in its 1-based line number.

  NOTE ON CLAIMS. These are SPECIFIC parsers and they ignore other parsers'
  claims: `Running refine_ca to refine and make full model for ...` is both a
  phase announcement and a form-B section title, and requirements 3 says
  nothing is deduplicated across kinds -- both are true and both are emitted.
  The GENERIC channel (labeled_values) does the opposite and skips claimed
  lines, because a heading also captured as a nameless fact is one fact
  counted twice."""
  phases, controls, items = [], [], []
  claims, refusals = {}, []
  for index, raw in enumerate(lines):
    line_no = index + 1
    match = _PHASE_RE.match(raw)
    if match:
      phases.append((match.group(1), line_no))
      claims[line_no] = "phases"
      continue
    match = _SKIP_RE.match(raw)
    if not match:
      continue
    rest = match.group(1)
    control = _CONTROL_SKIP_RE.match(rest)
    if control:
      controls.append((control.group(1), control.group(2), line_no))
      claims[line_no] = "control_skips"
      continue
    item = _ITEM_SKIP_RE.match(rest)
    if item:
      name = item.group(1)
      items.append((name.strip() if name else None, item.group(2), line_no))
      claims[line_no] = "exclusions"
      continue
    refusals.append((line_no, raw, REFUSAL_SKIP_UNRECOGNISED))
  return phases, controls, items, claims, refusals


# ------------------------------------------------------- stage tables, cycles
#
# phenix.refine tabulates its own progress, with a legend above it:
#
#   1_bss: bulk solvent correction and/or (anisotropic) scaling
#   1_xyz: refinement of coordinates
#   ------------------------------------------------------------------
#    stage r-work r-free bonds angles b_min b_max b_ave n_water shift
#         0    : 0.4747 0.5371 0.010  1.282  28.2 133.4  65.0   0   0.000
#         1_bss: 0.4064 0.4880 0.010  1.282  26.9 132.1  63.7   0   0.000
#
# ROWS ARE ANCHORED TO THEIR HEADER, not matched globally: the row shape
# matches 40 lines in refine_5_beta_blip where only 29 belong to the table.
#
# THE `end` ROW IS NOT A STAGE. The table's last row is named `end` and
# carries r_free 0.4942, while the last real stage (3_adp) carries 0.4935.
# 0.4942 is exactly the number the consuming project identified as HIDING the
# finding -- "the summary 0.4880 -> 0.4942 hides it". A parser that keeps
# `end` reports the summary and loses the trajectory.
#
# The legend names the stage. A stage with no legend entry keeps its code and
# reports description=None rather than inventing one.
#
# This table appears inside WIZARD logs too -- autobuild, autosol and
# predict_and_build embed refine runs -- so one parser serves several
# programs. 72 headers in 47 logs of corpus2/work.
#
# CYCLES:
#   Cycle 2  R/Rfree=0.20/0.23  Built=125  Placed=125 Resolution=2.5 A
#   Cycle 2  R/Rfree=999.90/999.90  Built=0  Placed=0 Resolution=2.1 A
#
# `R/Rfree=a/b` is ONE KEY AND TWO VALUES and must be split into r_work and
# r_free. **999.90 IS A SENTINEL meaning "no usable result"**: the cycle is
# flagged and the nonsense R-factors are NOT emitted as numbers a consumer
# could quote. A bare `Cycle 3 of morphing chain trace` carries no metrics --
# a counter, not a record -- and is reported as unparsed.

STAGE_SUMMARY_ROW = "end"
REFUSAL_STAGE_SUMMARY = "stage_row_is_the_table_summary"
CYCLE_SENTINEL = 999.90
REFUSAL_CYCLE_NO_METRICS = "cycle_line_carries_no_metrics"

_STAGE_HEADER_RE = re.compile(r"^\s*stage\s+r-work\s+r-free\b(.*)$", re.I)
_STAGE_ROW_RE = re.compile(r"^\s*(\S+?)\s*:\s+([-\d.]+(?:\s+[-\d.]+)+)\s*$")
_STAGE_LEGEND_RE = re.compile(r"^\s*(\S+?)\s*:\s+([A-Za-z]\S*.*)$")
_CYCLE_RE = re.compile(r"^\s*Cycle\s+(\d+)\b(.*)$")
_CYCLE_PAIR_RE = re.compile(r"R/Rfree\s*=\s*([-\d.]+)\s*/\s*([-\d.]+)")
_CYCLE_KV_RE = re.compile(r"(\b[A-Za-z][A-Za-z_]*)\s*=\s*([-\d.]+)")


def _to_number(text):
  try:
    return float(text)
  except ValueError:
    return None


def find_stage_tables(lines):
  """-> (stages, claims, refusals). Rows are taken only as a contiguous run
  after a header, and the trailing `end` summary row is not a stage."""
  stages, claims, refusals = [], {}, []
  n = len(lines)
  for index, raw in enumerate(lines):
    header = _STAGE_HEADER_RE.match(raw)
    if not header:
      continue
    columns = ["stage", "r_work", "r_free"] + header.group(1).split()
    claims[index + 1] = "stages"

    # the legend sits just above, before any rule line
    legend = {}
    back = index - 1
    while back >= 0:
      text = lines[back]
      if _SECTION_RULE_RE.match(text) or not text.strip():
        back -= 1
        if index - back > 12:
          break
        continue
      entry = _STAGE_LEGEND_RE.match(text)
      if not entry or _STAGE_ROW_RE.match(text):
        break
      legend[entry.group(1)] = entry.group(2).strip()
      back -= 1

    position = index + 1
    while position < n:
      row = _STAGE_ROW_RE.match(lines[position])
      if not row:
        break
      name = row.group(1)
      values = row.group(2).split()
      claims[position + 1] = "stages"
      if name == STAGE_SUMMARY_ROW:
        # REPORTED, not silently swallowed. DEFECT FOUND IN P4 REVIEW: the
        # row was claimed and emitted nowhere. It happened to surface anyway,
        # because it sits above a rule and the SECTION parser refused it as
        # `section_title_is_a_numeric_row` -- a coincidence, and a misleading
        # diagnosis for anyone reading the ledger. 77 such rows across the
        # two corpora, every one of them relying on that accident.
        refusals.append((position + 1, lines[position], REFUSAL_STAGE_SUMMARY))
        position += 1
        continue
      metrics = {}
      for column, value in zip(columns[1:], values):
        number = _to_number(value)
        if number is not None:
          metrics[column] = number
      stages.append((name, position + 1, metrics, legend.get(name)))
      position += 1
  return stages, claims, refusals


def find_cycles(lines):
  """-> (cycles, claims, refusals)."""
  cycles, claims, refusals = [], {}, []
  for index, raw in enumerate(lines):
    match = _CYCLE_RE.match(raw)
    if not match:
      continue
    line_no = index + 1
    rest = match.group(2)
    metrics = {}
    sentinel = False
    pair = _CYCLE_PAIR_RE.search(rest)
    if pair:
      work, free = _to_number(pair.group(1)), _to_number(pair.group(2))
      if work == CYCLE_SENTINEL or free == CYCLE_SENTINEL:
        sentinel = True               # no usable result -- do not quote it
      else:
        if work is not None:
          metrics["r_work"] = work
        if free is not None:
          metrics["r_free"] = free
    for key, value in _CYCLE_KV_RE.findall(rest):
      if key in ("R", "Rfree"):
        continue
      number = _to_number(value)
      if number is not None:
        metrics[key] = number
    if not metrics and not sentinel:
      refusals.append((line_no, raw, REFUSAL_CYCLE_NO_METRICS))
      continue
    cycles.append((int(match.group(1)), line_no, metrics, sentinel))
    claims[line_no] = "cycles"
  return cycles, claims, refusals


# ------------------------------------------------------------- measurements
#
# ATTACHED ONLY. A number reaches this channel because it arrived inside a
# structure already parsed -- a stage row, a cycle line, or the agent's
# metrics block. Generic `Key: value` stays in `labeled_values`,
# uninterpreted. Per-program metric patterns are the treadmill this tool
# exists to escape (non-goal 2.9), so v1 does not grow a vocabulary.
#
# WHAT THIS COSTS, MEASURED AND STATED RATHER THAN HIDDEN. Scored against
# agent_metrics_labels.json (75 field-instances, 14 logs):
#
#     agent-shape   attached 12 + labeled_values 61 = 73 of 75
#     GUI-shape     attached 12 + labeled_values 12 = 24 of 75
#
# THE HONEST NUMBER IS 24. The agent-shape 73 is circular -- 61 of those come
# from the agent's own metrics block, which is the answer sheet. This project
# has already made that exact error once (the "67 of 75 appear verbatim"
# claim was measured against text that contained the block), so the score is
# taken on GUI-shape.
#
# THE 51 UNCOVERED VALUES ARE IN THE LOG, in shapes the frozen screen refuses:
#
#     All-atom Clashscore : 3.56     indent 4, space before the colon
#   Clashscore            =   3.56   an `=` assignment
#
# Extending the screen to take both would add 64,451 candidates to 41,023
# (+157%) and move the label score from 24 to 30. That is a bad trade and it
# is NOT taken here; the arithmetic is in the handoff as a sign-off item.

AGENT_METRICS_MARKER = "FINAL QUALITY METRICS REPORT:"
_AGENT_METRIC_RE = re.compile(r"^\s*([A-Za-z][^:]{0,48}):\s+(\S.*?)\s*$")


def find_agent_measurements(lines):
  """-> (records, claims). The agent's own metrics block, marked as its own.

  Extracting it is useful; presenting it as the program's measurement is not
  (requirements 4.8.3), so every record carries source='agent' and the block
  is claimed so it never also appears as a labeled value."""
  records, claims = [], {}
  inside = False
  seen_any = False
  for index, raw in enumerate(lines):
    if raw.strip().startswith(AGENT_METRICS_MARKER):
      inside, seen_any = True, False
      continue
    if not inside:
      continue
    if _RULE_ONLY.match(raw):
      # DEFECT FOUND IN P8: the block opens `FINAL QUALITY METRICS REPORT:`
      # then a `-----` rule then the metrics, and closes with `*****`. The
      # first version treated ANY rule as the close, so the opening rule
      # ended the block immediately and the channel emitted nothing at all.
      if seen_any:
        inside = False
      continue
    match = _AGENT_METRIC_RE.match(raw)
    if not match:
      if raw.strip():
        inside = False                 # unrelated content ends the block
      continue
    seen_any = True
    records.append((match.group(1).strip(), match.group(2).strip(), index + 1))
    claims[index + 1] = "measurements"
  return records, claims


# -------------------------------------------------------- terminal failures

ABORT_MARKERS = ("STOPWIZARD", "EndOfResolve")
_TRACEBACK_RE = re.compile(r"^Traceback \(most recent call last\)")
_SORRY_RE = re.compile(r"^\s*Sorry[:,]\s*(.*)$")
# What ends a `Sorry:` block: the run resuming, or a box-drawing border.
# FOUND IN REVIEW: one block ran to the 40-line cap and stopped there by
# luck -- the ai_agent frames its output in box characters, and the real end
# was the frame's bottom border one line later. A cap that silently decides
# a boundary is the silent-truncation class this project keeps meeting, so
# the boundary rule is widened AND cap-truncation is now flagged rather than
# hidden.
_RESUMES_RE = re.compile(
  r"^\s*(Starting|Running|PHENIX|Traceback|Setting|Skipping|\[NOTICE\]"
  r"|[\u2514\u2518\u255a\u255d])")
SORRY_MAX_LINES = 40


def find_terminal_failures(lines):
  """-> (failures, claims). Reports what the program said about stopping."""
  failures, claims = [], {}
  index = 0
  n = len(lines)
  while index < n:
    raw = lines[index]
    if _TRACEBACK_RE.match(raw):
      end = index
      while end + 1 < n and not _SORRY_RE.match(lines[end + 1]) \
          and (lines[end + 1].startswith((" ", "\t"))
               or re.match(r"^\w+(Error|Exception):", lines[end + 1])):
        end += 1
      failures.append((TerminalFailure.TRACEBACK, raw.strip(), index + 1,
                       None, end + 1, False))
      for line_no in range(index + 1, end + 2):
        claims[line_no] = "terminal_failures"
      index = end + 1
      continue
    match = _SORRY_RE.match(raw)
    if match:
      # the remedy is what follows, to the end of the block. Blank lines are
      # inside it (`xray_data.r_free_flags.generate=True` sits between two),
      # so the block ends at TWO blanks or at a line that resumes the run.
      end = index
      while end + 1 < n and end - index < SORRY_MAX_LINES:
        following = lines[end + 1]
        if _RESUMES_RE.match(following):
          break
        if not following.strip() and (end + 2 >= n or not lines[end + 2].strip()):
          break
        end += 1
      remedy = "\n".join(l for l in lines[index + 1:end + 1]).strip() or None
      truncated = (end - index) >= SORRY_MAX_LINES
      failures.append((TerminalFailure.SORRY, match.group(1).strip() or
                       raw.strip(), index + 1, remedy, end + 1, truncated))
      for line_no in range(index + 1, end + 2):
        claims[line_no] = "terminal_failures"
      index = end + 1
      continue
    for marker in ABORT_MARKERS:
      if marker in raw:
        failures.append((TerminalFailure.ABORT, raw.strip(), index + 1,
                         None, index + 1, False))
        claims[index + 1] = "terminal_failures"
        break
    index += 1
  return failures, claims


# ---------------------------------------------------------- identification
#
# LAYER B. It reads whatever evidence exists and may return unknown, and it
# can NEVER gate layer-A extraction (see the module docstring): a 40-line
# headerless log that yields three cited facts and no program name is a
# successful extraction.
#
# ONE SIGNAL, deliberately. `Starting <program>` on the first lines, after
# the agent prefix is stripped. It is POSITIONAL self-identification -- the
# program writing about itself -- not a mention. A frequency rule over
# `phenix.X` mentions identifies every wizard as its child:
# predict_and_build_2_bromodomain says `phenix.refine` 38 times and never
# names itself.
#
# TWO DEGENERATE FORMS ARE BLACKLISTED because they identify nothing:
#   `Starting phenix`               15 logs -- the name is simply missing
#   `Starting libtbx.start_process` 11 logs -- the launcher, not the program
# `Starting phenix` with no name is worth reporting to Phenix: it is the
# self-identification fix §4.9 asks for, failing.
#
# TWO OTHER SIGNALS WERE MEASURED AND ARE NOT USED:
#
#   The STRUCTURAL HYPOTHESIS from §4.9 -- "the invoking program writes its
#   own parameter block first" -- FAILS. On the 107 banner-abstainers in
#   corpus2/work: 6 correct, 50 WRONG, 51 with no root at all. Roots are not
#   program names; `scaling` is xtriage's (27 wrong by itself) and
#   `data_manager` is the shared scope. Salvaging it needs a root->program
#   map from the PHIL masters, which is non-goal 2.3.
#
#   The `Found phil, .../<stem>.eff` record is REJECTED ON PRINCIPLE. It
#   scores 13 correct / 1 wrong / 93 absent on corpus2 abstainers, but the
#   GUI names the .eff after the job, which is also what names the log file.
#   Using it is THE FILENAME BY PROXY, and §4.9 forbids the filename because
#   ai_analysis derived the program that way and then printed a confident
#   report about a program that never ran. On wave-1 it is 0 correct and 3
#   wrong -- the same weakness, showing.

IDENTIFICATION_BANNER_LINES = 6

# SECOND SELF-IDENTIFICATION SIGNAL: the Phenix program header.
#
#   ****************************************************************
#             PHENIX predict_and_build  Wed Jun 17 10:06:16 2026
#   ****************************************************************
#
# or, in the wizard preamble, without the rules:
#
#                  PHENIX autosol  Sun Jul 19 10:20:17 2026
#
# The SHAPE comes from the Phenix source, not from guessing at logs: ~60
# programs print `"\n" + 60*"*" + "\n" + 10*" " + "PHENIX <name>"`, and
# ai_agent/ai_analysis print it with the name substituted at runtime. The
# name is whatever precedes the weekday, so multi-word names survive
# (`PHENIX Plan SAD experiment`), as do CamelCase (`AutoSol`) and dotted
# (`phenix.ai_agent`) forms -- normalisation lowercases, strips a `phenix.`
# prefix and joins words with underscores.
#
# MEASURED on all 657 distinct logs held: present in 36, and the name agrees
# with the filename in 35. The one "disagreement" is the signal being RIGHT
# and the filename being wrong -- `ligand_identification_68.log` is a GUI job
# that runs LigandFit, and says so three times (`Running LigandFit process
# 4`, `PHENIX ligandfit`, `Starting LigandFit with the command:`). It is the
# first demonstrated defect in the filename-derived truth labels, and it is
# exactly why requirements 4.9 forbids the filename as evidence.
#
# The window is larger than the banner's because the header sits after the
# citation preamble -- line 59 in predict_and_build.
IDENTIFICATION_HEADER_LINES = 80
SIGNAL_PHENIX_HEADER = "phenix_program_header"
_PHENIX_HEADER_RE = re.compile(
  r"^\s*PHENIX\s+(\S.*?)\s+(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun)\b")
BANNER_BLACKLIST = frozenset(["phenix", "libtbx.start_process", "job"])
SIGNAL_BANNER = "starting_banner"

# THE NAME MUST BE THE WHOLE REST OF THE LINE, one token, nothing after it.
#
# LATENT DEFECT FOUND IN P7 REVIEW (before it bit): the first version
# captured `(\S.*?)` -- everything to end of line. No current identification
# is affected, because within the six-line window every banner in both
# corpora is a bare name. But the corpus proves the unsafe shape exists:
# `Starting AutoBuild with the command:` appears at line ~55 in four logs,
# and ligandfit writes `Starting CC of ligand as input to map: 0.714` -- a
# metric wearing a verb. Under the old pattern the first would have produced
# `phenix.autobuild with the command:` and the second `phenix.cc of ligand
# as input to map: 0.714`, and a first-token rule would have produced
# `phenix.cc`, which is a WRONG NAME and fails the hard gate.
#
# Requiring a bare single token abstains on both. It costs nothing measurable
# today and closes the hole.
# Decoration around the banner is still a banner. FOUND WHEN THE HOLDOUT WAS
# OPENED: 18 of the 87 held-out logs write `**Starting phenix.molprobity **`,
# a form absent from the working corpus, and the extractor abstained on every
# one. Abstention, not a wrong name -- the hard gate held -- but it is the
# single largest coverage miss the holdout exposed.
#
# THE 32% HOLDOUT COVERAGE MEASURED BEFORE THIS FIX IS THE NUMBER OF RECORD.
# Any number after it is in-sample: the holdout is opened once, and a score
# taken after tuning to it measures nothing.
_BANNER_RE = re.compile(r"^[\s*]*Starting\s+(\S+?)\s*\**\s*$")


def _normalise_program(name):
  lowered = name.strip().lower()
  lowered = re.sub(r"^phenix\.", "", lowered)
  lowered = re.sub(r"^voyager\.", "", lowered)
  return re.sub(r"\s+", "_", lowered)


def identify_program(lines):
  """-> Identification. Abstains by returning empty candidates."""
  fired = []
  for index, raw in enumerate(lines[:IDENTIFICATION_BANNER_LINES]):
    match = _BANNER_RE.match(strip_agent_prefix(raw))
    if not match:
      continue
    fired.append(SIGNAL_BANNER)
    name = _normalise_program(match.group(1))
    if name in BANNER_BLACKLIST:
      return Identification(signals_fired=fired)     # saw it, cannot use it
    return Identification(
      [Candidate("phenix." + name, 0.95, Candidate.SELF_IDENTIFICATION,
                 index + 1)], signals_fired=fired)

  # the program header, second because the banner is earlier and already
  # measured; both are positional self-identification, neither is a word list
  for index, raw in enumerate(lines[:IDENTIFICATION_HEADER_LINES]):
    match = _PHENIX_HEADER_RE.match(strip_agent_prefix(raw))
    if not match:
      continue
    fired.append(SIGNAL_PHENIX_HEADER)
    name = _normalise_program(match.group(1))
    if name in BANNER_BLACKLIST:
      return Identification(signals_fired=fired)
    return Identification(
      [Candidate("phenix." + name, 0.95, Candidate.SELF_IDENTIFICATION,
                 index + 1)], signals_fired=fired)
  return Identification(signals_fired=fired)


# ------------------------------------------------------- completion records
#
# A record that SOMETHING FINISHED. **Never a claim that the run
# succeeded**, and an empty list is never a claim that it failed: of 30 failed
# runs in corpus2, 27 contain no error keyword at all in what the program
# itself wrote, so absence of error is not evidence of success. See non-goal
# 2.8. This channel reports observations; the verdict is the consumer's.
#
# WHY A LIST AND NOT A SCALAR. MEASURED on corpus2/work/ok n=230, after
# clustering adjacent markers: 70 logs have no completion event, 119 have
# exactly one, and 41 have MORE THAN ONE. One phaser log carries 59.
#
# CORRECTION RECORDED (P6 review): that "59" was attributed to
# `EXIT STATUS: SUCCESS` in the planning notes and repeated since. WRONG.
# Re-measured across both corpora: `EXIT STATUS` appears **exactly once** in
# each of the 24 logs that carry it, and the per-module marker is
# `Finished: <timestamp>` -- 59 in a2u-globulin-mr__phaser_3, 34 in
# AF_POMGNT2__phaser_6, 29 in three more. The conclusion (a list, not a
# scalar) was right; the evidence cited for it was not. That claim was
# propagated into the plan and three handoffs before anyone checked it.
#
# WHY THE THREE-LINE ENDING IS ONE RECORD. The Program Template closes with
# `Job complete` / `usr+sys time:` / `wall clock time:` on consecutive lines.
# Counting three would treble every completion count.
#
# APPLIES_TO NEEDS POSITIVE EVIDENCE, and the corpus offers exactly one form
# of it: `Finished with search...` / `Finished with predict_chain` NAMES ITS
# CHILD -- 113 such lines. Otherwise a record is `top_level` only when it is
# the LAST event in the file with nothing but blank lines after it, and
# `unknown` in every other case. A child completing is not the run
# completing, and in a wizard log we usually cannot tell which we are
# looking at.

_COMPLETION_PATTERNS = (
  (re.compile(r"^\s*Job complete\b"), "job_complete"),
  (re.compile(r"^\s*Finished\b(.*)$"), "finished"),
  (re.compile(r"^\s*usr\+sys time:"), "usr_sys_time"),
  (re.compile(r"wall clock time:"), "wall_clock"),
  (re.compile(r"^\s*EXIT STATUS:\s*(\S+)"), "exit_status"),
)
_FINISHED_WITH_RE = re.compile(r"^\s*Finished\s+with\s+(\S.*?)\s*$")
COMPLETION_CLUSTER_GAP = 6


def find_completion_records(lines):
  """-> (records, claims). Each record is (text, line, applies_to)."""
  hits = []
  for index, raw in enumerate(lines):
    for pattern, kind in _COMPLETION_PATTERNS:
      if pattern.search(raw):
        hits.append((index, raw, kind))
        break

  # Cluster by proximity AND by kind. DEFECT FOUND IN P6 (its own fixture):
  # proximity alone merged two genuinely separate `Job complete` lines two
  # apart into one event. A REPEATED KIND STARTS A NEW EVENT -- the Program
  # Template's ending is three DIFFERENT markers in a row, so requiring
  # distinct kinds keeps that as one record while keeping two completions
  # apart however close they sit.
  # Cluster by proximity AND by kind. DEFECT FOUND IN P6 (its own fixture):
  # proximity alone merged two genuinely separate `Job complete` lines two
  # apart into one event. A REPEATED KIND STARTS A NEW EVENT.
  #
  # SECOND DEFECT, found by a missed prediction: `Finished with X` is
  # SELF-CONTAINED and must not absorb what follows it. In 26 logs a child's
  # `Finished with process_predicted_model` sat within six lines of the
  # program's own `Job complete` / `usr+sys` / `wall clock` ending, and
  # because every kind differed they merged -- the run's own ending
  # disappearing into a record attributed to the child.
  events = []
  for index, raw, kind in hits:
    names_child = _FINISHED_WITH_RE.match(raw) is not None
    joins = (events
             and not names_child
             and not _FINISHED_WITH_RE.match(events[-1][0][1])
             and index - events[-1][-1][0] <= COMPLETION_CLUSTER_GAP
             and kind not in [k for _, _, k in events[-1]])
    if joins:
      events[-1].append((index, raw, kind))
    else:
      events.append([(index, raw, kind)])

  # nothing but blank lines after the final event is the only evidence the
  # corpus offers that a record belongs to the outermost run
  tail_is_blank = False
  if events:
    last = events[-1][-1][0]
    tail_is_blank = all(not lines[j].strip()
                        for j in range(last + 1, len(lines)))

  records, claims = [], {}
  for position, event in enumerate(events):
    index, raw, kind = event[0]
    child = _FINISHED_WITH_RE.match(raw)
    if child:
      applies_to = child.group(1)
    elif position == len(events) - 1 and tail_is_blank:
      applies_to = CompletionRecord.TOP_LEVEL
    else:
      applies_to = CompletionRecord.UNKNOWN
    records.append((raw.strip(), index + 1, applies_to))
    for member_index, _, _ in event:
      claims[member_index + 1] = "completion_records"
  return records, claims


# --------------------------------------------------------------- decisions
#
#   Setting rebuild_in_place=False as maps_only=True
#   Setting n_cycle_build=1 as nbatch >1 (nbatch =3)
#
# The program changing its own configuration AND STATING THE REASON. This is
# the highest-value structure in the corpus: it is the program announcing a
# branch. Setting, value and reason are captured separately.
#
# ONLY THE REASONED FORM COUNTS. MEASURED on corpus2/work n=253: 1,283
# `Setting` lines in 133 logs, of which only 208 carry a reason. The other
# 1,075 state a value and no branch --
#
#   Setting output.overwrite to True                      129
#   Setting macro_cycles to  3                             69
#   Setting default value of  True  for  thorough_denmod   41
#   Setting up prediction keywords...                      25  (a gerund)
#
# -- and are REPORTED as unparsed rather than fitted with a second pattern.
# Same precedent as P3: a control skip with no subject was refused rather
# than invented, and 234 unmatched skips were reported rather than guessed
# at. The 1,075 are a candidate v2 channel, recorded with their measurement
# so the decision is made deliberately rather than on first sighting.

# THREE NAMED REFUSALS, not one. DEFECT FOUND IN P5 REVIEW: every refused
# `Setting` line carried `setting_states_no_reason`, and for 490 of the 1,075
# that is the wrong diagnosis. A reader auditing the ledger would go looking
# for a missing reason on lines that are not assignments at all.
#
# MEASURED on corpus2/work n=253, of the 1,075 refused:
#   572  `Setting X to Y`            -- an assignment, genuinely no reason
#   310  `Setting default value of True for X`, `Setting parameters`
#                                    -- a narration, no `name = value` at all
#   180  `Setting up prediction keywords...`  -- a gerund
#    13  `Setting X=Y`               -- an assignment, genuinely no reason
#
# The same failure as D20 and D31: a refusal that does not say accurately
# WHY is a number that moved for unknown reasons.
REFUSAL_SETTING_NO_REASON = "setting_states_no_reason"
REFUSAL_SETTING_UP = "setting_up_is_not_an_assignment"
REFUSAL_SETTING_NARRATION = "setting_line_assigns_nothing"

_SETTING_UP_RE = re.compile(r"^\s*up\b")
_SETTING_ASSIGNS_RE = re.compile(r"^\s*\S+\s*(?:=|\bto\b)")

_SETTING_RE = re.compile(r"^\s*Setting\b(.*)$")
_DECISION_RE = re.compile(
  r"^\s*(\S+?)\s*=\s*(\S+?)\s+\b(?:as|because)\b\s+(\S.*?)\s*$")


def find_decisions(lines):
  """-> (decisions, claims, refusals)."""
  decisions, claims, refusals = [], {}, []
  for index, raw in enumerate(lines):
    match = _SETTING_RE.match(raw)
    if not match:
      continue
    line_no = index + 1
    rest = match.group(1)
    decision = _DECISION_RE.match(rest)
    if not decision:
      if _SETTING_UP_RE.match(rest):
        reason = REFUSAL_SETTING_UP
      elif _SETTING_ASSIGNS_RE.match(rest) or " to " in rest:
        reason = REFUSAL_SETTING_NO_REASON
      else:
        reason = REFUSAL_SETTING_NARRATION
      refusals.append((line_no, raw, reason))
      continue
    decisions.append((decision.group(1), decision.group(2),
                      decision.group(3), line_no))
    claims[line_no] = "decisions"
  return decisions, claims, refusals


# ------------------------------------------------------------ labeled values
#
# A `Key: value` record, UNINTERPRETED. No unit normalisation, no type
# inference, no mapping of `R-free` / `R Free` / `r_free` onto one name --
# that mapping is the maintenance treadmill this tool exists to escape, and
# it belongs to a consumer who can revise it without touching this file.
#
# This channel is what the long tail depends on. 30% of successful logs and
# 9% of failed ones carry a stage table, cycle or decision; most Phenix
# programs are 20-300 line tools that announce nothing and still report their
# inputs, counts, space group and resolution.
#
# TWO COLLAPSES, both measured on corpus2/work n=253 before being written.
# Raw admission gives 39,957 items.
#
#   1. IDENTICAL (label, value) pairs collapse ALWAYS: 39,957 -> 26,119
#      (-35%). The same fact stated twice is one fact. `repeat_count` holds
#      the number of occurrences, `line` the first and `end` the last, so
#      nothing about where it happened is lost.
#
#   2. A label with MORE THAN `LABEL_DISTINCT_LIMIT` DISTINCT VALUES in one
#      log collapses to a single item, and the lines it would have produced
#      are recorded as `rule_excluded` -- reported, never silently dropped.
#
# WHY 50, with the alternatives that were measured and rejected:
#
#   N=12  21,389 items, 141 label groups collapsed
#   N=20  22,683 items,  54 groups
#   N=50  23,711 items,  20 groups          <- chosen
#
# The low values look attractive on the item count and cost real series. In
# the 21-80 zone sit `CC` with 26 distinct values in a predict_and_build log
# and `Residual` with 22 in an ensembler log -- per-cycle series, which are
# exactly the finding this project exists to surface. The runaway cases are
# far above: `Target left` 1,999, `Chi-restraints excluded` 467, `REVERT`
# 194. 50 separates them cleanly, sits above the largest genuine series
# measured (36), and fires on 20 groups in the whole corpus.
#
# The per-log burden after collapse 1 alone is already median 47 / max 2,271;
# collapse 2 caps the worst log at 1,224. So the second collapse is a guard
# against runaway, not a filter -- which is why its limit is set high.

LABEL_DISTINCT_LIMIT = 50
REFUSAL_LABEL_RUNAWAY = "label_exceeds_distinct_value_limit"


# A line may carry SEVERAL pairs: `R:   0.42  Rfree:   0.48` is two numbers,
# and capturing it as label "R" value "0.42  Rfree:   0.48" buries the second
# one. 9% of the key:value lines across 554 logs carry two or more pairs.
#
# SPLIT ONLY WHEN THE WHOLE LINE IS PAIRS. `Space group: C 1 2 1 (No. 5)` has
# one colon and a value of several tokens; a naive pair split turns it into
# `C` and loses the rest. So the line is split only if a sequence of
# `label: single-token` pairs consumes it entirely -- which is a SHAPE, not a
# judgement about which numbers matter.
_PAIR_RE = re.compile(r"([A-Za-z][A-Za-z0-9_()/.'-]{0,30}):\s+([^\s:]+)")


def split_compound_pairs(text):
  """-> [(label, value), ...]. A single pair unless the line is all pairs."""
  pairs = _PAIR_RE.findall(text)
  if len(pairs) < 2:
    return []
  # the test is whether the pairs account for the line once whitespace is
  # removed -- counting characters, not guessing at a slack allowance
  packed = "".join(a + b for a, b in pairs)
  if re.sub(r"[\s:]", "", text) != packed:
    return []                      # there is other content; leave it alone
  return pairs


def find_labeled_values(lines, claimed, limit=LABEL_DISTINCT_LIMIT):
  """-> (records, claims, refusals).

  `records` are (label, value, first_line, last_line, repeat_count) tuples in
  first-appearance order. Lines already claimed by another parser are skipped
  so that a section title is never also a labeled value."""
  groups = []
  index = {}
  for position, raw in enumerate(lines):
    line_no = position + 1
    if line_no in claimed:
      continue
    if screen_line(raw) != SCREEN_KEY_VALUE:
      continue
    match = _SCREEN_KV_RE.match(raw)
    if not match:
      continue
    pairs = split_compound_pairs(raw.strip())
    if not pairs:
      pairs = [(match.group(1).strip(), match.group(2).strip())]
    for label, value in pairs:
      label = label.strip()
      value = value.strip()
      if label not in index:
        index[label] = len(groups)
        groups.append((label, []))
      groups[index[label]][1].append((value, line_no))

  records, claims, refusals = [], {}, []
  for label, occurrences in groups:
    distinct = []
    seen = {}
    for value, line_no in occurrences:
      if value in seen:
        entry = distinct[seen[value]]
        entry[2] = line_no
        entry[3] += 1
      else:
        seen[value] = len(distinct)
        distinct.append([value, line_no, line_no, 1])
    if len(distinct) > limit:
      first_value, first_line = occurrences[0]
      records.append((label, first_value, first_line, occurrences[-1][1],
                      len(occurrences)))
      claims[first_line] = "labeled_values"
      for _, line_no in occurrences[1:]:
        refusals.append((line_no, lines[line_no - 1], REFUSAL_LABEL_RUNAWAY))
      continue
    for value, first_line, last_line, count in distinct:
      records.append((label, value, first_line, last_line, count))
    # every occurrence of this label is understood, not just the first of
    # each distinct value. DEFECT FOUND IN REVIEW: this claim loop was nested
    # INSIDE the loop above, so it re-walked all occurrences once per distinct
    # value -- O(distinct x occurrences) for the same result. A label with 50
    # distinct values over 1,000 lines did 50,000 dictionary writes to claim
    # 1,000 lines.
    for _, line_no in occurrences:
      claims[line_no] = "labeled_values"
  return records, claims, refusals


# -------------------------------------------------------------- LogStructure


class LogStructure(object):
  """What scan() returns.

  Order is preserved within each kind -- stage and cycle sequence IS the
  finding. Nothing is deduplicated across kinds: a line may be both a phase
  announcement and a decision, and both are emitted."""

  def __init__(self, n_lines=0):
    self.n_lines = n_lines
    for kind in ITEM_KINDS:
      setattr(self, kind, [])
    self.identification = Identification()
    self.regions = []
    # lines a parser claimed that the frozen screen never proposed. A rising
    # number here says the screen is too narrow -- reported, not hidden.
    self.claimed_outside_screen = []

  # ---- views

  @property
  def forms(self):
    """The kinds that are non-empty, in a fixed order."""
    return [k for k in ITEM_KINDS if getattr(self, k)]

  @property
  def is_flat(self):
    """No structure found at all. `unparsed` does not count as structure: a
    log whose only output is 'I saw things I could not read' is flat."""
    return not [k for k in ITEM_KINDS if k != "unparsed" and getattr(self, k)]

  def items(self):
    """Every item, of every kind, in line order.

    When the program is unknown and section_id is None -- the normal state for
    a headerless log -- a consumer's only remaining context is SPATIAL
    PROXIMITY, and it can only use that if items can be walked in document
    order across kinds."""
    out = []
    for kind in ITEM_KINDS:
      out.extend(getattr(self, kind))
    return sorted(out, key=lambda it: (it.line, ITEM_KINDS.index(it.kind)))

  def unparsed_counts(self):
    """The three numbers, always together (see Unparsed)."""
    counts = dict((s, 0) for s in UNPARSED_STATUSES)
    for u in self.unparsed:
      counts[u.status] += 1
    return counts

  def reach(self):
    """Three metrics, always reported together.

    structural -- any item at all
    basic      -- at least one labeled value, measurement or completion record
    diagnostic -- at least one stage, cycle, decision, skip or exclusion

    A program with no stages BY DESIGN should not be scored as a failure when
    it correctly yields its inputs, counts and space group."""
    return dict(
      structural=not self.is_flat,
      basic=any(getattr(self, k) for k in BASIC_KINDS),
      diagnostic=any(getattr(self, k) for k in DIAGNOSTIC_KINDS))

  # ---- methods from the contract

  def metric_moves(self, metric, threshold=0.0):
    """Stages where `metric` changed by at least `threshold`. This is what
    turns a table into a finding."""
    moves = []
    previous = None
    for stage in self.stages:
      if metric not in stage.metrics:
        continue
      value = stage.metrics[metric]
      if previous is not None and abs(value - previous) >= threshold:
        moves.append(dict(stage=stage.name, before=previous, after=value,
                          delta=value - previous, line=stage.line,
                          description=stage.description))
      previous = value
    return moves

  def exclusion_groups(self):
    """Exclusions grouped by reason. Eighteen identical exclusions are one
    finding, not eighteen."""
    order, groups = [], {}
    for exclusion in self.exclusions:
      if exclusion.reason not in groups:
        groups[exclusion.reason] = []
        order.append(exclusion.reason)
      groups[exclusion.reason].append(exclusion)
    return [dict(reason=r, count=len(groups[r]),
                 lines=[e.line for e in groups[r]],
                 items=[e.item for e in groups[r]]) for r in order]

  def __repr__(self):
    return "LogStructure(%d lines, forms=%r, unparsed=%r)" % (
      self.n_lines, self.forms, self.unparsed_counts())


# --------------------------------------------------------------------- scan


def scan(text, program_name=None):
  """Read a Phenix log. Pure: text in, structure out.

  `program_name` is an OPTIONAL hint. The extractor must work without it --
  the consumer often does not know it, and one production defect came from a
  wrong guess at it. It may never gate layer-A extraction; a corpus test and a
  negative control enforce that.

  Encoding-safe by contract: callers decode with errors="replace" (see main).
  """
  if text is None:
    text = ""
  lines = split_lines(text)
  structure = LogStructure(n_lines=len(lines))
  structure.regions = find_regions(lines)
  preamble = find_preamble(lines)
  if preamble is not None:
    structure.regions.append(preamble)

  claimed = {}
  refusals = []

  def in_preamble(line_no):
    return preamble is not None and preamble.contains(line_no)

  raw_sections, section_claims, section_refusals = find_sections(lines)
  raw_sections = [x for x in raw_sections if not in_preamble(x[1])]
  section_claims = dict((k, v) for k, v in section_claims.items()
                        if not in_preamble(k))
  section_refusals = [r for r in section_refusals if not in_preamble(r[0])]
  claimed.update(section_claims)
  refusals.extend(section_refusals)
  for position, (title, start, end) in enumerate(raw_sections):
    source = _source_at(structure.regions, start)
    # A SECTION MUST NOT CROSS A PROVENANCE BOUNDARY. DEFECT FOUND IN P1
    # REVIEW: a section runs to the next heading or to EOF, so the last
    # program-written section of an agent-wrapped log swallowed the agent's
    # footer -- 16 cases in the corpora, including molprobity's `Summary`
    # "containing" the agent's FINAL QUALITY METRICS REPORT. A consumer
    # grouping by section would then attribute the agent's numbers to the
    # program, which is the exact failure requirements 4.8 exists to prevent.
    for region in structure.regions:
      if region.start > start and region.source != source:
        end = min(end, region.start - 1)
    structure.sections.append(Section(
      title, start, end=max(start, end), section_id=position, source=source))

  # A LINE -> SECTION LOOKUP. The first version scanned the whole section
  # list per line: 389 sections x 11,450 items on the largest log, and one
  # log carries 419 sections.
  #
  # DEFECT FOUND IN P2 BY THE EQUIVALENCE TEST: the second version used a
  # forward CURSOR, correct only if every call arrives in ascending line
  # order. P2 added a second pass over the file, so the cursor was already at
  # the end when the ledger loop restarted at line 1, and lookups returned
  # None for the whole file. Stateless bisect instead -- correct regardless
  # of call order, still O(log n). The optimisation-equivalence test caught
  # it; nothing else would have, because the wrong answer was a plausible
  # None rather than a crash.
  section_starts = [x.line for x in structure.sections]
  section_bounds = [(x.line, x.end, x.section_id) for x in structure.sections]

  def section_of(line_no):
    position = bisect.bisect_right(section_starts, line_no) - 1
    if position < 0:
      return None
    start, end, ident = section_bounds[position]
    return ident if start <= line_no <= end else None

  phases, controls, items, ps_claims, ps_refusals = find_phases_and_skips(
    lines)
  claimed.update(ps_claims)
  refusals.extend(ps_refusals)
  for name, line_no in phases:
    structure.phases.append(Phase(
      name, line_no, section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))
  for what, reason, line_no in controls:
    structure.control_skips.append(ControlSkip(
      what, reason, line_no, section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))
  for name, reason, line_no in items:
    structure.exclusions.append(Exclusion(
      name, reason, line_no, section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))

  stage_rows, stage_claims, stage_refusals = find_stage_tables(lines)
  claimed.update(stage_claims)
  refusals.extend(stage_refusals)
  for name, line_no, metrics, description in stage_rows:
    structure.stages.append(Stage(
      name, line_no, metrics=metrics, description=description,
      section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))

  cycle_rows, cycle_claims, cycle_refusals = find_cycles(lines)
  claimed.update(cycle_claims)
  refusals.extend(cycle_refusals)
  for number, line_no, metrics, sentinel in cycle_rows:
    structure.cycles.append(Cycle(
      number, line_no, metrics=metrics, sentinel=sentinel,
      section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))

  failure_rows, failure_claims = find_terminal_failures(lines)
  claimed.update(failure_claims)
  for failure_kind, text, line_no, remedy, end, truncated in failure_rows:
    structure.terminal_failures.append(TerminalFailure(
      failure_kind, text, line_no, remedy=remedy, end=end,
      truncated=truncated,
      section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))

  completion_rows, completion_claims = find_completion_records(lines)
  claimed.update(completion_claims)
  for text, line_no, applies_to in completion_rows:
    structure.completion_records.append(CompletionRecord(
      text, line_no, applies_to=applies_to, section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))

  decision_rows, decision_claims, decision_refusals = find_decisions(lines)
  claimed.update(decision_claims)
  refusals.extend(decision_refusals)
  for setting, value, reason, line_no in decision_rows:
    structure.decisions.append(Decision(
      setting, value, reason, line_no, section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))

  agent_metrics, agent_claims = find_agent_measurements(lines)
  claimed.update(agent_claims)
  for name, value, line_no in agent_metrics:
    structure.measurements.append(Measurement(
      name, value, line_no, context=AGENT_METRICS_MARKER,
      section_id=section_of(line_no),
      # DEFECT FOUND IN P8 REVIEW: this hardcoded source='agent'. A metrics
      # block quoted MID-FILE by the agent's own log sits in a region whose
      # source is `unknown` -- we cannot tell the agent wrapping this log
      # from this log quoting an agent -- and hardcoding turned that
      # uncertainty back into a claim, which is the one thing source=unknown
      # exists to prevent. 76 measurements in corpus2/work were affected.
      source=_source_at(structure.regions, line_no)))

  # attached numbers, re-emitted flat so a consumer can query one channel
  for stage in structure.stages:
    for metric_name in sorted(stage.metrics):
      structure.measurements.append(Measurement(
        metric_name, stage.metrics[metric_name], stage.line,
        context="stage:" + stage.name, section_id=stage.section_id,
        source=stage.source))
  for cycle in structure.cycles:
    for metric_name in sorted(cycle.metrics):
      structure.measurements.append(Measurement(
        metric_name, cycle.metrics[metric_name], cycle.line,
        context="cycle:%d" % cycle.number, section_id=cycle.section_id,
        source=cycle.source))

  label_records, label_claims, label_refusals = find_labeled_values(
    lines, claimed)
  label_records = [r for r in label_records if not in_preamble(r[2])]
  label_claims = dict((k, v) for k, v in label_claims.items()
                      if not in_preamble(k))
  claimed.update(label_claims)
  refusals.extend(label_refusals)
  generic_only = set(label_claims)
  for label, value, first, last, count in label_records:
    structure.labeled_values.append(LabeledValue(
      label, value, first, repeat_count=count, end=last,
      section_id=section_of(first),
      source=_source_at(structure.regions, first)))

  # THE UNPARSED LEDGER. A screened candidate no parser claimed is
  # `unclaimed`; a line a parser REFUSED by a named rule is `rule_excluded`
  # and carries the rule. Claiming a line the frozen screen never proposed is
  # counted separately -- a rising count there says the screen is too narrow,
  # which is the honest counterweight to freezing it.
  refused_lines = dict((line, reason) for line, _, reason in refusals)
  for index, raw in enumerate(lines):
    line_no = index + 1
    rule = screen_line(raw)
    if line_no in refused_lines:
      structure.unparsed.append(Unparsed(
        text=raw, line=line_no, screen_rule=rule or "none",
        status=RULE_EXCLUDED, excluded_by=refused_lines[line_no],
        section_id=section_of(line_no),
        source=_source_at(structure.regions, line_no)))
      continue
    if rule is None:
      continue
    if line_no in generic_only:
      # recorded, but not understood as anything in particular -- the honest
      # successor to a single `unparsed` number, and expected to be large
      structure.unparsed.append(Unparsed(
        text=raw, line=line_no, screen_rule=rule, status=GENERIC_ONLY,
        section_id=section_of(line_no),
        source=_source_at(structure.regions, line_no)))
      continue
    if line_no in claimed:
      continue
    structure.unparsed.append(Unparsed(
      text=raw, line=line_no, screen_rule=rule, status=UNCLAIMED,
      section_id=section_of(line_no),
      source=_source_at(structure.regions, line_no)))

  structure.claimed_outside_screen = sorted(
    line for line in claimed if screen_line(lines[line - 1]) is None)

  # LAYER B, LAST. Computed after every layer-A parser and used by none of
  # them. The order is the enforcement: a parser cannot consult a program
  # name that does not exist yet. REVISION (P7 review): the first version
  # sat mid-way through the parsers with a comment claiming it was last.
  structure.identification = identify_program(lines)

  return structure


# ---------------------------------------------------------------------- main


_LONG_PATH_RE = re.compile(r"(/[^\s/]+){3,}/?")


def _elide_paths(text):
  """`/Users/terwill/Documents/gene-5-mad/AutoSol_run_1_/TEMP0/build_2/...`
  is not a finding; the line it sits on is. Long paths collapse to
  `.../basename` so the finding stays on one line."""
  def shorten(match):
    whole = match.group(0)
    return ".../" + whole.rstrip("/").rsplit("/", 1)[-1] if len(whole) > 30 \
        else whole
  return _LONG_PATH_RE.sub(shorten, text)


def report(structure, path, out=None, limit=12, verbose=False,
           value_limit=30, only_label=None, pattern=None,
           structure_text=None, grep_limit=40):
  """A readable account of one log, for a person at a terminal.

  Every line carries the log line it came from, because the point of this
  tool is that a claim can be checked. Long channels are truncated with a
  count, so a 1,000-stage refine run does not scroll away the finding."""
  # `out=sys.stdout` as a DEFAULT ARGUMENT binds the stream at definition
  # time, so redirecting sys.stdout afterwards had no effect and the report
  # went to the real terminal while the caller captured nothing. Found by the
  # test written for this function, on its first run.
  stream = sys.stdout if out is None else out

  def lines_of(text):
    return split_lines(text) if text else []

  def emit(text=""):
    stream.write(text + "\n")

  emit("=" * 72)
  emit(path)
  emit("=" * 72)

  identification = structure.identification
  if identification.is_unknown:
    signals = ", ".join(identification.signals_fired) or "no signal found"
    emit("program:  UNKNOWN  (%s)" % signals)
  else:
    top = identification.candidates[0]
    emit("program:  %s   (%s, line %d)"
         % (identification.name, top.evidence_kind, top.line))
  reach = structure.reach()
  emit("lines:    %d" % structure.n_lines)
  emit("reach:    %s" % (", ".join(k for k in ("structural", "basic",
                                               "diagnostic") if reach[k])
                         or "none"))

  # HOW IT ENDED, first, because it is what a reader most needs and what the
  # report used to bury. An AutoSol run that ended in STOPWIZARD on line 3877
  # showed only a child completion from line 1994.
  emit("")
  if structure.terminal_failures:
    emit("HOW IT ENDED:  the program reported a failure -- see below")
  elif any(c.applies_to == CompletionRecord.TOP_LEVEL
           for c in structure.completion_records):
    emit("HOW IT ENDED:  a top-level completion record was found")
  elif structure.completion_records:
    emit("HOW IT ENDED:  completion records found, none of them top-level")
  else:
    emit("HOW IT ENDED:  NO completion record  (not the same as failed)")

  if structure.terminal_failures:
    emit("")
    emit("THE PROGRAM'S OWN ACCOUNT OF STOPPING (%d)"
         % len(structure.terminal_failures))
    for failure in structure.terminal_failures:
      emit("  %6d  [%s] %s" % (failure.line, failure.failure_kind,
                               failure.text[:70]))
      if failure.remedy:
        for remedy_line in failure.remedy.splitlines()[:8]:
          if remedy_line.strip():
            emit("          | %s" % remedy_line.strip()[:70])
      if failure.truncated:
        emit("          | ... block hit the %d-line cap; see the log from"
             " line %d" % (SORRY_MAX_LINES, failure.line))

  def block(title, items, render):
    if not items:
      return
    emit("")
    emit("%s (%d)" % (title, len(items)))
    for item in items[:limit]:
      text = _elide_paths(render(item))
      if len(text) > 96:
        # a chain sequence is one labeled value 1,000 characters long; the
        # line number is what matters, so point at it rather than print it
        text = text[:96] + "... (%d chars, see line %d)" % (len(text),
                                                            item.line)
      emit("  %6d  %s" % (item.line, text))
    if len(items) > limit:
      emit("  %6s  ... %d more" % ("", len(items) - limit))

  block("DECISIONS", structure.decisions,
        lambda x: "%s = %s   because %s" % (x.setting, x.value, x.reason))
  block("SKIPPED PHASES", structure.control_skips,
        lambda x: "%s   because %s" % (x.what, x.reason))

  if structure.exclusions:
    emit("")
    emit("EXCLUDED CANDIDATES (%d, in %d groups)"
         % (len(structure.exclusions), len(structure.exclusion_groups())))
    for group in structure.exclusion_groups():
      named = [i for i in group["items"] if i]
      emit("  %4d x  %s" % (group["count"], group["reason"]))
      emit("          lines %s%s"
           % (", ".join(str(n) for n in group["lines"][:8]),
              " ..." if len(group["lines"]) > 8 else ""))
      if named:
        emit("          e.g. %s" % ", ".join(named[:6]))

  def stage_line(x):
    # lead with the two numbers a crystallographer reads first; the rest are
    # on the cited line. Alphabetical order buried r_free behind b_max.
    lead = " ".join("%s=%s" % (k, x.metrics[k])
                    for k in ("r_work", "r_free") if k in x.metrics)
    return "%-14s %-30s  %s" % (x.name, lead, x.description or "")
  block("STAGES", structure.stages, stage_line)

  moves = structure.metric_moves("r_free", 0.001)
  worse = [m for m in moves if m["delta"] > 0]
  if worse:
    emit("")
    emit("R-FREE WENT UP AT (%d of %d stages that moved it)"
         % (len(worse), len(moves)))
    for move in worse[:limit]:
      emit("  %6d  %-14s %.4f -> %.4f  (+%.4f)%s"
           % (move["line"], move["stage"], move["before"], move["after"],
              move["delta"],
              "   " + move["description"] if move["description"] else ""))
  block("CYCLES", structure.cycles,
        lambda x: "cycle %-3d %s%s"
        % (x.number, " ".join("%s=%s" % kv for kv in sorted(x.metrics.items())),
           "   NO USABLE RESULT (sentinel)" if x.sentinel else ""))
  block("COMPLETION", structure.completion_records,
        lambda x: "%-40s applies to: %s" % (x.text[:40], x.applies_to))
  block("MEASUREMENTS", [m for m in structure.measurements
                         if m.source == SOURCE_AGENT],
        lambda x: "%s = %s   (added by the agent, not the program)"
        % (x.name, x.value))
  block("PHASES", structure.phases, lambda x: x.name)

  # background, last and collapsed unless asked for: on a 3,877-line AutoSol
  # run these are 57 sections and 342 labeled values, and they buried the
  # three findings above them.
  if verbose:
    block("SECTIONS", structure.sections, lambda x: x.title)
    block("LABELED VALUES", structure.labeled_values,
          lambda x: "%s: %s%s" % (x.label, x.value,
                                  "   (x%d, to line %d)"
                                  % (x.repeat_count, x.end)
                                  if x.repeat_count > 1 else ""))
  else:
    # LABELED VALUES ARE WHERE THE NUMBERS ARE. Collapsing them to a count
    # was a regression: `Map-model CC: 0.47` and `R: 0.42 Rfree: 0.48` were
    # captured all along and the report hid them behind one BACKGROUND line.
    #
    # 342 values on one AutoSol run is too many to list, so they are GROUPED
    # BY LABEL -- 259 distinct here -- showing the last value and the line to
    # look at. Grouping is factual, the same as exclusion_groups; it is not a
    # judgement about which numbers matter, which is non-goal 2.1. Use
    # --label to pull one label's whole series.
    if structure.labeled_values:
      groups = {}
      for item in structure.labeled_values:
        groups.setdefault(item.label, []).append(item)
      emit("")
      emit("VALUES THE RUN REPORTED (%d values, %d distinct labels)"
           % (len(structure.labeled_values), len(groups)))
      shown = 0
      # MOST RECENTLY REPORTED FIRST. Chronology, not a judgement about
      # importance -- but a run's last numbers are the ones a reader came
      # for, and ascending order buried `Map-model CC` and `R/Rfree` under
      # 30 lines of setup.
      for label in sorted(groups, key=lambda k: -groups[k][-1].line):
        items = groups[label]
        if shown >= value_limit:
          emit("  %6s  ... %d more labels (--all, or --label <name>)"
               % ("", len(groups) - shown))
          break
        last = items[-1]
        suffix = "   (%d values, from line %d)" % (len(items), items[0].line) \
            if len(items) > 1 else ""
        emit("  %6d  %s: %s%s" % (last.line, label,
                                  _elide_paths(last.value)[:60], suffix))
        shown += 1
    if structure.sections:
      emit("")
      emit("SECTIONS  %d   (--all to list them)" % len(structure.sections))

  if pattern is not None:
    # SEARCH THE RAW LINES, not the unparsed ledger.
    #
    # I claimed once that unread prose "sits in the unparsed ledger". IT DOES
    # NOT: `unparsed` holds only lines the frozen SCREEN proposed -- 737 of
    # this AutoSol log's 3,877 -- and `Model is in /path` matches no screen
    # rule at all, so it is in neither the channels nor the ledger. Anything
    # that only searched `unparsed` would have missed exactly the lines this
    # option exists to find.
    #
    # Each hit says whether anything claimed the line, so a reader can tell
    # "the tool read this and put it somewhere" from "the tool never saw it".
    lowered = pattern.lower()
    claimed_by = {}
    for kind in ITEM_KINDS:
      for item in getattr(structure, kind):
        if kind != "unparsed":
          claimed_by.setdefault(item.line, kind)
    for item in structure.unparsed:
      claimed_by.setdefault(item.line, "unparsed/" + item.status)
    if structure_text is None:
      emit("")
      emit("--grep needs the log text; call report(..., structure_text=text)")
      hits = []
    else:
      hits = [(number, text)
              for number, text in enumerate(lines_of(structure_text), 1)
              if lowered in text.lower()]
    emit("")
    emit("LINES MATCHING %r (%d)" % (pattern, len(hits)))
    for number, text in hits[:grep_limit]:
      emit("  %6d  [%s] %s"
           % (number, claimed_by.get(number, "not read"),
              _elide_paths(text.strip())[:64]))
    if len(hits) > grep_limit:
      emit("  %6s  ... %d more" % ("", len(hits) - grep_limit))

  if only_label:
    # EXACT label first, substring only as a fallback. `--label R` matched
    # every label containing the letter r -- 196 of them, which is not a
    # series, it is the whole log again.
    exact = [x for x in structure.labeled_values
             if x.label.lower() == only_label.lower()]
    wanted = exact or [x for x in structure.labeled_values
                       if only_label.lower() in x.label.lower()]
    emit("")
    emit("EVERY VALUE FOR %r (%d)" % (only_label, len(wanted)))
    for item in wanted:
      emit("  %6d  %s: %s" % (item.line, item.label,
                              _elide_paths(item.value)[:70]))

  counts = structure.unparsed_counts()
  emit("")
  emit("NOT UNDERSTOOD  unclaimed=%d  generic_only=%d  rule_excluded=%d"
       % (counts[UNCLAIMED], counts[GENERIC_ONLY], counts[RULE_EXCLUDED]))
  # SHAPES, not just totals. "unclaimed=175" is a number nobody can act on;
  # the same 175 lines shown as six masked shapes turn out to be xtriage
  # tables, PDB CRYST/SCALE records and bare rules -- reassuring rather than
  # alarming, and it took a measurement to know which.
  shapes = {}
  for item in structure.unparsed:
    if item.status != UNCLAIMED:
      continue
    masked = re.sub(r"[-+]?[0-9][0-9.]*", "#", item.text.strip())[:56]
    if masked not in shapes:
      shapes[masked] = [0, item.line]
    shapes[masked][0] += 1
  if shapes:
    emit("    the commonest unclaimed shapes (digits masked):")
    for masked, (count, first) in sorted(shapes.items(),
                                         key=lambda kv: -kv[1][0])[:8]:
      emit("      %4d x  line %-6d %s" % (count, first, masked))
  refusals = {}
  for item in structure.unparsed:
    if item.excluded_by:
      refusals[item.excluded_by] = refusals.get(item.excluded_by, 0) + 1
  for reason in sorted(refusals):
    emit("    refused: %-40s %d" % (reason, refusals[reason]))
  emit("")
  emit("Nothing above is interpreted. Every line number is real -- check it.")


def main(argv=None):
  """Thin I/O wrapper. scan() does no file I/O by contract, so this is where
  reading and decoding live. Never gates on the `.log` extension: a production
  defect traced to `if ext != ".log": return None`, after which the program
  name went unset and the pipeline carried on and printed a confident report
  about a program that never ran."""
  argv = list(sys.argv[1:] if argv is None else argv)
  if not argv or argv[0] in ("-h", "--help"):
    print(__doc__.strip())
    print("")
    print("  log_structure_extractor.py <file>...     a readable report")
    print("  log_structure_extractor.py --all <file>          everything")
    print("  log_structure_extractor.py --label R <file>      one label's series")
    print("  log_structure_extractor.py --grep 'Model is in' <file>"
          "   any line, read or not")
    print("  log_structure_extractor.py --summary <file>...   one line each")
    return 0
  summary = "--summary" in argv
  verbose = "--all" in argv
  # Options are consumed BY POSITION, not by value. FOUND IN REVIEW: the
  # first version dropped any path equal to the pattern, so
  # `--grep Sorry Sorry.log` searched nothing and said so about no file.
  only_label = None
  pattern = None
  paths = []
  position = 0
  while position < len(argv):
    argument = argv[position]
    if argument in ("--grep", "--label") and position + 1 < len(argv):
      if argument == "--grep":
        pattern = argv[position + 1]
      else:
        only_label = argv[position + 1]
      position += 2
      continue
    if argument.startswith("--grep="):
      pattern = argument.split("=", 1)[1]
    elif argument.startswith("--label="):
      only_label = argument.split("=", 1)[1]
    elif not argument.startswith("-"):
      paths.append(argument)
    position += 1
  if summary and (pattern or only_label):
    # silently ignoring a flag is how a user comes to distrust the output
    sys.stderr.write("--summary ignores --grep and --label; drop --summary"
                     " to use them\n")
    return 2
  status = 0
  for path in paths:
    try:
      handle = open(path, "rb")
      try:
        text = handle.read().decode("utf-8", "replace")
      finally:
        handle.close()
    except IOError as error:
      sys.stderr.write("%s: %s\n" % (path, error))
      status = 1
      continue
    structure = scan(text)
    if summary:
      counts = structure.unparsed_counts()
      print("%s: %d lines, forms=%s, reach=%s, unparsed=%s" % (
        path, structure.n_lines, ",".join(structure.forms) or "-",
        ",".join(k for k, v in sorted(structure.reach().items()) if v) or "-",
        ",".join("%s=%d" % (k, counts[k]) for k in UNPARSED_STATUSES)))
    else:
      report(structure, path, verbose=verbose, only_label=only_label,
             pattern=pattern, structure_text=text)
  return status


if __name__ == "__main__":
  try:
    sys.exit(main())
  except IOError:
    # a shipped CLI must survive `... | head`. FOUND IN THE FINAL REVIEW:
    # piping --help into head raised BrokenPipeError and printed a traceback.
    sys.exit(0)
