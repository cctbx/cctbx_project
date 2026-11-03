"""Check that the figures in a report are supported by the log.

Step 3 of IMPLEMENTATION.md. Runs **after** generation, appends a
visible flag, and never rewrites.

WHY IT MUST TOLERATE ROUNDING AND BOUNDS

Phase 0a's sweep found **zero fabricated figures in twenty reports**, so
this guards a problem we do not have -- which makes rarity the whole
value. A verifier that fires often teaches users which warnings to
ignore.

In the first real report examined, two figures were derived rather than
quoted:

    "Resolution (1.74 A)"              from 1.74434
    "Twin fraction estimates (<0.03)"  a bound over 0.00, 0.008,
                                       0.029, 0.022

**Both correct; a naive check flags both.** Two false positives on the
first report examined is the failure mode above, arriving immediately.

WHAT IT DELIBERATELY DOES NOT DO

* It does not check that a figure carries the right LABEL. If `0.26`
  appears as R-work when the log says R-free, the number is supported
  and the claim is wrong. **"The numbers are checked" overstates this**;
  the badge is *Data-source checks passed*.
* It does not reach interpretation. The wrong-phasing recommendation of
  defect 1, and "the enantiomorphic space groups I 4 and I 41" -- a
  wrong domain claim in a real shipped report -- pass every check here.
  No deterministic checker reaches them.
"""
from __future__ import absolute_import, division, print_function

import re

#: Figures below this are usually indices, counts or cycle numbers
#: rather than measurements, and appear in prose constantly.
TRIVIAL = 3

#: A report figure is supported if some log figure rounds to it at the
#: report's own precision. 1.74 accepts 1.74434; it does not accept 1.8.
#: The trailing guard rejects a following digit or word character but
#: MUST allow a sentence-ending period. The first version used
#: `(?![\w.])`, which made every figure at the end of a sentence
#: invisible -- "The R-free was 0.187." matched nothing. That is a false
#: NEGATIVE, and a checker that silently cannot see a whole class of
#: figures is worse than one that flags too much. Caught by the control
#: test on its first run.
#: Scientific notation must be matched FIRST, or "7.68e+00" yields the
#: exponent "00" as a figure and the mantissa is lost -- so a log value
#: printed as 7.68e+00 was invisible while a spurious 0.0 was created.
#: That was the third of the three forms the plan named (rounding,
#: thousands separators, scientific notation) and the third one missed.
NUMBER = re.compile(
    r"(?<![\w.,])("
    r"\d+(?:\.\d+)?[eE][+-]?\d+"          # 7.68e+00
    r"|\d{1,3}(?:,\d{3})+(?:\.\d+)?"      # 8,830
    r"|\d+(?:\.\d+)?"                     # 1.74434
    r")(?!\.?\d)(?!\w)")

#: Words that turn a figure into a bound rather than a quantity.
BOUND_BEFORE = re.compile(
    r"(?:<|>|under|below|above|less than|greater than|no more than|"
    r"at most|at least|approximately|about|~)\s*$", re.I)


class Flag(object):
    __slots__ = ("figure", "kind", "context")

    def __init__(self, figure, kind, context):
        self.figure = figure
        self.kind = kind
        self.context = context

    def __repr__(self):
        return "Flag(%s, %s)" % (self.figure, self.kind)

    def __str__(self):
        return "%s (%s): ...%s..." % (self.figure, self.kind, self.context)


def _numbers(text):
    """(value, string, position) for every figure in `text`.

    Thousands separators are normalised. Without this, "8,830" was read
    as two figures -- 8 and 830 -- and 830 is not in a log that says
    8830. **Every one of the 43 false positives in the first
    measurement was a comma**, which is why the plan specified
    "allowing for rounding, thousands separators and scientific
    notation" and why implementing two of the three was not enough.
    """
    out = []
    for m in NUMBER.finditer(text):
        raw = m.group(1)
        try:
            out.append((float(raw.replace(",", "")), raw, m.start()))
        except ValueError:
            pass
    return out


def _rounds_to(candidate, shown):
    """Does `candidate` round to `shown` at `shown`'s precision?

    The precision comes from the REPORT, not from a fixed tolerance:
    a report saying 1.74 is claiming two decimals, and 1.74434 satisfies
    that. A report saying 1.744 is claiming three, and 1.74 does not.
    """
    shown = shown.replace(",", "")
    if "." in shown and "e" not in shown.lower():
        places = len(shown.split(".")[1])
        if round(candidate, places) == round(float(shown), places):
            return True
        # Reports TRUNCATE as well as round: a log saying 30.9956
        # appears as "30.9-12.4 A" for a shell range. Refusing that
        # produced four false positives on two real reports.
        # Guard the truncation compare: a log carrying 1e+308 or an
        # inf overflows int(). Found by the fabrication control, which
        # injected a figure that pushed a scientific-notation value
        # down this path.
        import math
        if not math.isfinite(candidate) or abs(candidate) > 1e15:
            return False
        factor = 10.0 ** places
        return int(candidate * factor) == int(float(shown) * factor)
    return abs(candidate - float(shown)) < 0.5


def check_numbers(report, log, evidence=""):
    """Return a list of Flag for figures the log does not support.

    `evidence` is any additional supporting text (a deterministic
    evidence panel, when one is in use).
    """
    source = log + "\n" + (evidence or "")
    source_values = [v for v, _s, _p in _numbers(source)]
    source_strings = set(s for _v, s, _p in _numbers(source))

    flags = []
    for value, shown, pos in _numbers(report):
        if value <= TRIVIAL and "." not in shown:
            continue                      # counts, cycles, list indices
        if shown in source_strings or shown.replace(",", "") in \
                source_strings:
            continue                      # quoted verbatim

        if any(_rounds_to(v, shown) for v in source_values):
            continue                      # rounded, e.g. 1.74 <- 1.74434

        # A percentage in the report against a fraction in the log:
        # "73.5%" where the log says 0.735. Common, unambiguous, and
        # cheap to accept.
        following = report[pos + len(shown):pos + len(shown) + 2]
        if following.startswith("%"):
            if any(_rounds_to(v * 100.0, shown) for v in source_values):
                continue

        preceding = report[max(0, pos - 24):pos]
        if BOUND_BEFORE.search(preceding):
            # A bound is supported when log values actually sit under
            # (or over) it and close to it -- "<0.03" over a maximum of
            # 0.029. A bound nothing approaches is not support.
            near = [v for v in source_values
                    if 0 < value - v <= max(0.5 * value, 0.01)]
            if near:
                continue

        context = " ".join(
            report[max(0, pos - 30):pos + len(shown) + 20].split())
        flags.append(Flag(shown, "not found in the log", context))
    return flags


def check_program_names(report, registry):
    """Any `phenix.<name>` in the report must be a real program.

    One line, and it would have caught `phenix.royal_` on its first
    appearance. This check BLOCKS rather than flags.
    """
    bad = []
    for name in set(re.findall(r"phenix\.[a-z_0-9]+", report)):
        if name not in registry:
            bad.append(name)
    return sorted(bad)


def summarise(flags, blocked):
    """One line for the operator, in shadow mode."""
    if blocked:
        return "BLOCK: report names programs that do not exist: %s" % (
            ", ".join(blocked))
    if not flags:
        return "Data-source checks passed."
    return "Data-source checks: %d figure(s) not found in the log." % len(
        flags)
