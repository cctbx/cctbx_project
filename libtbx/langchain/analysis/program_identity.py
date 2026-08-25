"""Resolve which Phenix program produced a log.

Step 1 of IMPLEMENTATION.md. Replaces two measured defects:

  defect 3  On the default server path the name is never derived from
            the file name -- the derivation at ai_analysis.py:1063 sits
            inside run_job_locally, which the server path does not take
            -- so content substring patterns decide it. Wrong on 6 of 20
            corpus logs; seen live identifying an xtriage log as
            phenix.phaser via 'molecular replacement'.

  defect 4  ai_analysis.py:1067 writes "Log file for phenix.<derived>"
            into the prompt. On 20 of 20 corpus logs that named a
            program which does not exist.

  caller -> banner -> file name -> refuse

Measured 20 of 20, zero wrong, zero refused, on both corpus shapes.

DESIGN NOTES, EACH PAID FOR

* **No LLM-derived value is consulted.** The shipped code takes the
  reported name from a scrape of the model's own summary. That scrape
  returns '' on the format the model actually emits, which is why
  retrieval has been running against a constant string. Nothing here
  reads a summary, a completion, or a processed_log_dict, and
  tst_program_identity asserts that against this file's source.

* **The registry is loaded, not hardcoded.** An earlier draft of this
  resolver held 15 names against programs.yaml's 23 and would have
  refused a `polder` log as 'not a Phenix program'. A curated list
  drifts -- the same argument this project makes about the PHIL dumps.

* **A file name is INFERRED, never authoritative.** Registry membership
  proves a name exists, not that it is correct, and a user can rename a
  file. Callers that need certainty read `.authoritative`.

* **Self-contained.** Part A does not ship the evidence layer or the log
  extractor, so this imports neither. Banner detection is ten lines
  here rather than a dependency on 1,400.
"""
from __future__ import absolute_import, division, print_function

import os
import re

#: Fallback used only when programs.yaml cannot be read. Deliberately
#: short: it exists so the resolver degrades rather than crashes, not
#: as a second source of truth.
_FALLBACK = (
    "phenix.xtriage", "phenix.phaser", "phenix.refine", "phenix.autobuild",
    "phenix.autosol", "phenix.ligandfit", "phenix.molprobity",
    "phenix.predict_and_build", "phenix.find_reference", "phenix.mtriage",
)

#: Leading punctuation is allowed before the keyword.  GUI-shape logs
#: emit "**Starting phenix.molprobity **" with markdown asterisks, and
#: an anchored "^\s*Starting" missed every one of them -- 12 of 20 logs
#: that STATE their program were being reported as `filename`, i.e. as
#: inferred.  The answer was right and the provenance was wrong, which
#: is worse than it sounds: `authoritative` is the field a caller is
#: meant to branch on.
_BANNER = (
    re.compile(r"COMMAND THAT WAS RUN:.*?(phenix\.[A-Za-z_0-9]+)"),
    re.compile(r"Starting\s+(phenix\.[A-Za-z_0-9]+)"),
    re.compile(r"PHENIX\s+(phenix\.[A-Za-z_0-9]+)"),
    #: Bare names: "LOG TEXT: Starting AutoBuild", "Starting
    #: find_reference".  Resolved STRICTLY -- exact stem only, no prefix
    #: matching -- because "Starting refinement" would otherwise
    #: prefix-match phenix.refine, and "Starting libtbx.start_process"
    #: heads three refine logs and names something that is not the
    #: program that ran.
    re.compile(r"Starting\s+([A-Za-z][A-Za-z_0-9]{2,})\b"),
)

_HEAD = 5000          # banner scan window; banners are at the top or absent


class Identity(object):
    """The answer, with where it came from.

    `authoritative` is the field a caller should branch on. A banner or
    a caller-supplied name states what ran. A file name suggests it.
    """

    __slots__ = ("name", "source", "authoritative", "note")

    def __init__(self, name, source, authoritative, note=""):
        self.name = name
        self.source = source
        self.authoritative = authoritative
        self.note = note

    def __repr__(self):
        return "Identity(%r, source=%r, authoritative=%r)" % (
            self.name, self.source, self.authoritative)


def load_registry(path=None):
    """Program names from knowledge/programs.yaml."""
    if path and os.path.exists(path):
        try:
            import yaml
            with open(path) as handle:
                names = tuple(sorted(yaml.safe_load(handle)))
            extra = tuple(n for n in _FALLBACK if n not in names)
            return names + extra
        except Exception:                                  # noqa: BLE001
            pass
    return _FALLBACK


def _canonical(raw, registry):
    """Map a raw string to a registry program, or None.

    Exact, then the stem with any run number stripped, then a UNIQUE
    prefix match in either direction. The last clause exists because
    'PHASER_beta' carries no digits to strip: without it the resolver
    scored 19 of 20 while the measurement it was built from said 20 --
    a reimplementation of tested logic not being the tested logic.
    Uniqueness is required so a prefix can never choose between two
    programs.
    """
    if not raw:
        return None
    name = str(raw).strip().lower()
    if name in registry:
        return name
    if not name.startswith("phenix."):
        name = "phenix." + name
    if name in registry:
        return name

    stem = re.sub(r"_\d+.*$", "", name[len("phenix."):])
    if not stem:
        return None
    exact = [p for p in registry if p[len("phenix."):] == stem]
    if len(exact) == 1:
        return exact[0]
    near = [p for p in registry
            if stem.startswith(p[len("phenix."):])
            or p[len("phenix."):].startswith(stem)]
    return near[0] if len(near) == 1 else None


def _strict(raw, registry):
    """Exact registry membership after adding the prefix and stripping a
    run number. No prefix matching: a bare word taken off a banner is
    much weaker evidence than a dotted name, and prefix matching would
    turn "Starting refinement" into phenix.refine."""
    if not raw:
        return None
    name = str(raw).strip().lower()
    if name in registry:
        return name
    if not name.startswith("phenix."):
        name = "phenix." + name
    if name in registry:
        return name
    stem = "phenix." + re.sub(r"_\d+.*$", "", name[len("phenix."):])
    return stem if stem in registry else None


def _from_banner(log_text, registry):
    if not log_text:
        return None
    head = log_text[:_HEAD]
    for index, pattern in enumerate(_BANNER):
        strict = (index == len(_BANNER) - 1)      # the bare-name pattern
        for line in head.split("\n"):
            match = pattern.search(line)
            if match:
                resolved = (_strict(match.group(1), registry) if strict
                            else _canonical(match.group(1), registry))
                if resolved:
                    return resolved
    return None


def _from_filename(log_path, registry):
    if not log_path:
        return None
    stem = os.path.splitext(os.path.basename(log_path))[0]
    return _canonical(stem, registry)


def resolve_program(caller=None, log_text="", log_path=None,
                    registry_path=None, registry=None):
    """Resolve the program that produced a log.

    Returns an :class:`Identity`. `name` is None when nothing
    trustworthy was available -- **refuse rather than guess**, because
    the wrong program's documentation is worse than none and an invented
    program name is worse than both.
    """
    registry = registry or load_registry(registry_path)

    resolved = _canonical(caller, registry)
    if resolved:
        return Identity(resolved, "caller", True)
    if caller:
        return Identity(None, "none", False,
                        "caller supplied %r, which is not a Phenix "
                        "program" % caller)

    resolved = _from_banner(log_text, registry)
    if resolved:
        return Identity(resolved, "banner", True)

    resolved = _from_filename(log_path, registry)
    if resolved:
        return Identity(resolved, "filename", False,
                        "inferred from the file name; the log does not "
                        "state the program")

    return Identity(None, "none", False,
                    "the log does not state the program and the file "
                    "name does not resolve to one")


def header_line(identity):
    """The line to prepend to the log text, or None.

    Defect 4 in one function. The shipped code writes
    "Log file for phenix.%s" % <derived>, unconditionally, producing
    phenix.xtriage_1 on 20 of 20 runs. Here an unresolved identity
    yields None and **the caller writes nothing at all** -- silence is
    correct where the old code fabricated.
    """
    if identity is None or identity.name is None:
        return None
    return "Log file for %s" % identity.name
