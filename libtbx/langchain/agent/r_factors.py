"""Extract R-work and R-free from a Phenix log, reading the right column.

Every pattern that previously did this located the LABEL and took the
first number after it. In the paired forms that number is the
R-**work**, so `R/Rfree=0.21/0.27` yielded 0.21 as the R-free. The value
then drove stop decisions: at 2.5 A the success threshold is 0.23, so a
misread 0.21 stopped a run that a correct 0.27 would have continued.

Six sites extracted R-free from text; three misread the same way:

    knowledge/patterns.yaml  metrics.r_free   -> graph_nodes 611 and 1508
    knowledge/programs.yaml  autobuild r_free -> STEP 1 registry extraction
    agent/metrics_analyzer.py:270

The registry pattern requires a hyphen, so it misses `R/Rfree=` and
lands on the correct `Final R-work = ..., R-free = ...` line **by
accident**. Given `R/R-free:` it fails identically. It is not a safe
fallback, which is why this helper exists rather than a rule that
prefers one pattern over another.

FORMS HANDLED
-------------
Collected from 19 corpus logs plus two production runs:

    Cycle 3  R/Rfree=0.21/0.26  Built=121
    Model: .../refine_1.pdb  R/Rfree=0.21/0.27
    New values of R/Rfree:    0.21/   0.27
    with R/Rfree=   0.21/  0.27 and map-model CC =   0.75
    Best solution on cycle: 14    R/Rfree=0.38/0.42
    R/R-free: 0.5/0.54
    R and R-free: 0.20 0.26
    R:   0.21  Rfree:   0.27
    Start R-work = 0.1808, R-free = 0.2308
    Final R-work = 0.3802, R-free = 0.4498
    R Free: 0.4498
    R-free                =   0.2267

REJECTED
--------
    Cycle 2  R/Rfree=999.90/999.90     sentinel, not a measurement
    R-free likelihood based estimates  prose
    stage r-work r-free bonds angles   table header

An R factor is a fraction: anything at or above 1.0 is not one. 999.90
appears in no source file, so nothing depends on it being passed
through; it is returned as absent rather than as a number.

DELIBERATELY CONSERVATIVE
-------------------------
An unrecognised line yields None rather than a guess. The hand-labelled
expectations in tst_r_factor_extraction.py cover every shape found in 21
logs, but a form none of them contains would otherwise be read by
whichever pattern happened to match — which is the failure this replaces.
"""

from __future__ import absolute_import, division, print_function

import re


# An R factor is a fraction.  1.0 or above is a sentinel or a parse
# error, never a measurement.
_MAX_PLAUSIBLE = 1.0

_NUM = r"([0-9]*\.?[0-9]+)"

# Paired forms: one label covering two columns, R-work FIRST.
# Tried before the single-value forms, because the single-value pattern
# also matches inside these lines and would take the R-work.
_PAIRS = (
    # R/Rfree=0.21/0.27   R/R-free: 0.5/0.54   R/Rfree:  0.21/  0.27
    re.compile(r"\bR\s*/\s*R-?\s*free\s*[=:]?\s*" + _NUM + r"\s*/\s*" + _NUM,
               re.IGNORECASE),
    # R and R-free: 0.20 0.26
    re.compile(r"\bR\s+and\s+R-?\s*free\s*[=:]\s*" + _NUM + r"\s+" + _NUM,
               re.IGNORECASE),
    # R:   0.21  Rfree:   0.27      R :  0.21  Rfree :  0.27
    re.compile(r"\bR\s*[=:]\s*" + _NUM + r"\s+R-?\s*free\s*[=:]\s*" + _NUM,
               re.IGNORECASE),
)

# Both labelled explicitly on one line, in either order.
_LABELLED_WORK = re.compile(r"\bR-?\s*work\s*[=:]\s*" + _NUM, re.IGNORECASE)
_LABELLED_FREE = re.compile(r"\bR[-_\s]?free\s*[=:]\s*" + _NUM, re.IGNORECASE)


def _plausible(text):
    """Return the float if it could be an R factor, else None."""
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    if value < 0.0 or value >= _MAX_PLAUSIBLE:
        return None
    return value


def extract_r_factors_from_line(line):
    """Return (r_free, r_work) for one line; either may be None.

    Paired forms are tried first: the single-value pattern matches
    inside them and would return the R-work column.
    """
    if not line:
        return None, None
    for pattern in _PAIRS:
        match = pattern.search(line)
        if match:
            work = _plausible(match.group(1))
            free = _plausible(match.group(2))
            return free, work
    free_match = _LABELLED_FREE.search(line)
    work_match = _LABELLED_WORK.search(line)
    free = _plausible(free_match.group(1)) if free_match else None
    work = _plausible(work_match.group(1)) if work_match else None
    return free, work


def extract_last_r_factors(text):
    """Return (r_free, r_work) from the LAST line that supplies each.

    Refinement logs carry one R-free per macro-cycle and the final value
    is the one that matters.  R-free and R-work are tracked separately:
    a line may give one without the other, and taking them from
    different lines is how an R-work once ended up stored as an R-free.
    """
    if not text:
        return None, None
    last_free = None
    last_work = None
    for line in text.split("\n"):
        if "free" not in line.lower():
            continue
        free, work = extract_r_factors_from_line(line)
        if free is not None:
            last_free = free
        if work is not None:
            last_work = work
    return last_free, last_work


def r_factors_are_consistent(r_free, r_work, tolerance=0.005):
    """False when r_free is below r_work by more than `tolerance`.

    R-free is computed on reflections held out of refinement, so it
    cannot be lower than R-work.  An inversion means at least one value
    was read from the wrong column.

    The tolerance matters: a run recorded r_free 0.21 and r_work 0.2107,
    two values a rounding apart, and a strict comparison would call that
    an error.
    """
    if r_free is None or r_work is None:
        return True
    return r_free >= (r_work - tolerance)
