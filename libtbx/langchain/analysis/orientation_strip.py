#!/usr/bin/env python
"""
Orientation strip -- the < 1 s, zero-LLM first screen.

    phenix.python orientation_strip.py

THE MECHANISM, and why it is this and not a ranking.

The obvious approach is to score every extracted item and take the top
N. That needs a scoring function, which is judgment, which is the thing
that goes wrong -- and it would need tuning per program, which is 23
tuning problems rather than one.

Instead: FIXED SLOTS. The strip always has the same six slots in the
same order. Each slot is filled by a deterministic query against the
extractor output and has a hard line cap. Four consequences:

  * "fits one screen" is guaranteed by construction, not hoped for
  * the shape is identical on every log, so a user learns once where to
    look -- which is most of what orientation means
  * an unfillable slot is VISIBLE as an omission rather than silently
    replaced by whatever ranked next
  * there is nothing to tune, so nothing to get wrong per program

THE DIVISION OF LABOUR. A better answer arrives ~30 s later, so this
strip carries only what is QUOTED OR TRIVIALLY DERIVED from the log.
No inference, no interpretation, no advice. That guarantees the strip
can never contradict the report on judgment -- only on fact, where the
strip is the more reliable of the two.

ONE ORDERING RULE. Slots are in fixed order EXCEPT that PROBLEMS
promotes to the top when non-empty. If the run died, that is the only
thing the user needs in the first second.

ABSENCE IS REPORTED. "no errors reported" is information; a missing
Problems slot is ambiguous. Slots that are merely empty are omitted.
"""
from __future__ import absolute_import, division, print_function

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE = os.path.dirname(HERE)          # .../libtbx/langchain

# The extractor is already in the tree; do NOT vendor a second copy.
# Its own comment warns that "two modules of one name silently shadow
# each other on sys.path".
#
# Its location is not settled from what I have been given: the module's
# docstring says libtbx/langchain/log_structure_extractor.py, the
# development package ships it as a log_extraction/ package with an
# __init__.py (byte-identical file), and Tom states there are no python
# files directly in langchain/. So try the package forms in order and
# fall back to a path search, rather than hard-coding a guess.
X = None
for _mod in ("libtbx.langchain.log_extraction.log_structure_extractor",
             "libtbx.langchain.analysis.log_structure_extractor",
             "libtbx.langchain.log_structure_extractor"):
    try:
        X = __import__(_mod, fromlist=["scan"])
        break
    except ImportError:
        continue
if X is None:                                             # pragma: no cover
    for _cand in (os.path.join(PACKAGE, "log_extraction"),
                  HERE, PACKAGE, os.path.join(HERE, "vendor"),
                  os.path.join(os.path.dirname(PACKAGE), "vendor")):
        if os.path.exists(os.path.join(_cand,
                                       "log_structure_extractor.py")):
            if _cand not in sys.path:
                sys.path.insert(0, _cand)
            break
    import log_structure_extractor as X                   # noqa: E402

# The registry lives at libtbx/langchain/knowledge/, one level up from
# this module. An earlier draft looked for HERE/knowledge and would have
# shipped silently falling back to a ten-name list.
KNOWLEDGE = os.path.join(PACKAGE, "knowledge", "programs.yaml")
if not os.path.exists(KNOWLEDGE):
    KNOWLEDGE = os.path.join(HERE, "knowledge", "programs.yaml")

# Program identity comes from the tested resolver. This module used to
# carry its own load_registry() and resolve(); they agreed with
# program_identity on ten of ten probes, by luck rather than
# construction, and the two have already diverged once in this project
# at 19/20 against a measured 20/20.
try:
    from libtbx.langchain.analysis.program_identity import resolve_program
except ImportError:                                       # pragma: no cover
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    from program_identity import resolve_program

LINE_BUDGET = 30

# Curated allowlists. These are the one piece of judgment in the
# mechanism, they are explicit, and they are program-independent.
SETUP_LABELS = re.compile(
    r"space group|unit cell|resolution|high.?resolution|low.?resolution|"
    r"number of (scatterers|models|atoms|chains|residues|reflections)|"
    r"total number of atoms|completeness|wavelength|solvent content|"
    r"number of copies|sequence|matthews|anisotropy|twin", re.I)

HEADLINE_METRICS = re.compile(
    r"^(r.?work|r.?free|r.?value|cc|cc.?mask|cc.?work|cc.?free|fom|"
    r"tfz|rfz|llg|clashscore|molprobity.?score|ramachandran|rotamer|"
    r"map.?cc|resolution|completeness|d.?min|rmsd|bond|angle)", re.I)

FILEISH = re.compile(r"\.(pdb|cif|mtz|map|ccp4|seq|fasta|eff|dat|sca|hkl)\b",
                     re.I)


def _clean(text, width=88):
    """One line, no runaway whitespace, truncated on a word boundary."""
    t = " ".join(str(text or "").split())
    if len(t) <= width:
        return t
    cut = t[:width].rsplit(" ", 1)[0]
    return cut + " ..."


def _last_by_name(items):
    """Final value wins. A log states R-free thirty times; the run's
    R-free is the last one."""
    out = {}
    for it in items:
        name = getattr(it, "name", None)
        if name:
            out[" ".join(str(name).split())] = it
    return out


# The registry. A name not in it is not shown -- the first version of
# this file printed "program: refine_3" and "program: PHASER_beta",
# which are not programs: the same class-4 defect found in
# ai_analysis.py's injected header.
#
# It is LOADED, not hardcoded. The hardcoded version held 15 names
# against programs.yaml's 23, and missed pdbtools, polder, map_to_model
# and seven more. A curated list drifts; that is the same argument this
# project makes about the PHIL dumps, and it applied here within a day.
def slot_run(st, log_name, identity=None):
    """What ran, and where the name came from.

    Delegates entirely to program_identity: caller -> log banner -> file
    name -> refuse, every result validated against the registry.
    Measured 20 of 20, zero wrong, zero refused, on both corpus shapes.
    """
    if identity is None:
        identity = resolve_program(
            log_text="\n".join(l for l in _head_lines(st)) or "",
            log_path=log_name,
            registry_path=KNOWLEDGE if os.path.exists(KNOWLEDGE) else None)
    if identity.name is None:
        return ["program: not stated in the log, and the file name does "
                "not resolve to a Phenix program"]
    lines = ["program: %s  (%s)" % (identity.name,
                                    "stated in the log"
                                    if identity.authoritative
                                    else "inferred from the file name")]
    return lines


def _head_lines(st):
    """The identification evidence the extractor already found, so the
    resolver sees a banner without re-reading the whole log."""
    out = []
    for cand in getattr(st.identification, "candidates", [])[:3]:
        out.append("Starting %s" % cand.name)
    return out



    """What ran. A caller-supplied name would take precedence here; with
    no caller, the log's banner beats the file name, and a real
    disagreement is flagged rather than resolved."""
    cands = st.identification.candidates
    banner = resolve(cands[0].name) if cands else None
    base = os.path.splitext(os.path.basename(log_name))[0]
    fname = resolve("_".join(base.split("_")[:-1]))

    if banner:
        lines = ["program: %s  (stated in the log)" % banner]
        if fname and fname != banner:
            lines.append("  note: the file name suggests %s -- flagged, not "
                         "resolved" % fname)
        return lines
    if fname:
        return ["program: %s  (from the file name; the log does not say)"
                % fname]
    return ["program: not stated in the log, and the file name does not "
            "resolve to a Phenix program"]


def slot_problems(st):
    lines = []
    for tf in getattr(st, "terminal_failures", [])[:2]:
        lines.append("FAILED: %s" % _clean(getattr(tf, "text", tf)))
    for cs in getattr(st, "control_skips", [])[:2]:
        # ControlSkip carries `what` and `reason`; it has NO `text`
        # attribute, so the old fallback chain ended at the object
        # itself and printed a Python repr into a user-facing file --
        # "ControlSkip(line=13829, what='seq alignment...'". Visible on
        # 4 of 20 corpus logs and on the 649 KB one.
        what = getattr(cs, "what", None)
        reason = getattr(cs, "reason", None)
        if not what:
            continue
        text = _clean(what) if not reason else _clean(
            "%s (%s)" % (what, reason))
        lines.append("skipped: %s" % text)
    return lines


def slot_inputs(st):
    lines, seen = [], set()
    for lv in st.labeled_values:
        label = " ".join(str(getattr(lv, "label", "")).split())
        value = " ".join(str(getattr(lv, "value", "")).split())
        if not value or value in ('""', "''"):
            continue
        if not SETUP_LABELS.search(label) and not FILEISH.search(value):
            continue
        key = label.lower()
        if key in seen:
            continue
        seen.add(key)
        lines.append("%s: %s" % (_clean(label, 34), _clean(value, 46)))
    return lines


def slot_outline(st):
    names = [_clean(getattr(p, "title", None) or getattr(p, "name", ""), 40)
             for p in (st.phases or st.sections)]
    names = [n for n in names if n]
    if not names:
        return []
    if len(names) <= 6:
        return [" -> ".join(names[:3])] + ([" -> ".join(names[3:])]
                                           if names[3:] else [])
    return ["%s -> %s -> ... -> %s   (%d steps)"
            % (names[0], names[1], names[-1], len(names))]


def slot_numbers(st):
    lines = []
    for name, it in _last_by_name(st.measurements).items():
        if not HEADLINE_METRICS.match(name):
            continue
        val = " ".join(str(getattr(it, "value", "")).split())
        unit = " ".join(str(getattr(it, "unit", "") or "").split())
        lines.append("%s = %s%s" % (_clean(name, 30), _clean(val, 24),
                                    (" " + unit) if unit else ""))
    return lines


def slot_outputs(st):
    """How the run ended.

    This slot used to be labelled "wrote:" and printed
    `completion_records`, which are lines like "Job complete" and
    "Finished: Fri Jun 19 15:22:26 2026" -- **not files**. On 13 of 20
    corpus logs it told the reader a timestamp had been written. The
    extractor exposes no output-file structure, so the honest fix is to
    label the slot for what it actually holds and stop claiming files
    were written.
    """
    seen = []
    for cr in getattr(st, "completion_records", []):
        text = _clean(getattr(cr, "text", ""), 78)
        if text and text not in seen:
            seen.append(text)
        if len(seen) >= 3:
            break
    return seen


def build(text, log_name, pending=True):
    st = X.scan(text)
    problems = slot_problems(st)

    slots = [
        ("WHAT RAN", slot_run(st, log_name), 3),
        ("ON WHAT", slot_inputs(st), 7),
        ("WHAT IT DID", slot_outline(st), 3),
        ("NUMBERS (final values)", slot_numbers(st), 7),
        ("FINISHED", slot_outputs(st), 3),
    ]
    if problems:
        slots.insert(0, ("PROBLEMS", problems, 4))
    elif st.is_flat or st.n_lines < 40:
        # "none reported" on an empty or truncated log is false comfort:
        # absence of evidence presented as evidence of absence, which is
        # this project's most common error class. Say what we actually
        # know.
        slots.insert(0, ("PROBLEMS", [
            "cannot say -- the log has only %d lines and little structure; "
            "it may be truncated or still being written" % st.n_lines], 1))
    else:
        slots.append(("PROBLEMS", ["no recognised error messages found in the log"], 1))

    out = []
    for title, lines, cap in slots:
        lines = [l for l in lines if l]
        if not lines:
            continue
        out.append("[%s]" % title)
        for l in lines[:cap]:
            out.append("  " + l)
        if len(lines) > cap:
            out.append("  ... and %d more, in the full report"
                       % (len(lines) - cap))
    # No time promise. The measured wait is dominated by a queue we do
    # not control, so a first screen that opens by missing a deadline is
    # worse than one that never sets it.
    #
    # `pending` is False when this is written to summary_file, which is
    # produced AFTER the analysis exists -- promising a pending analysis
    # in a saved file is simply false.
    if pending:
        out.append("[the full analysis is being prepared]")
    return out, st


def summary_for_log(log_text, log_name, identity=None):
    """The deterministic contents of summary_file.

    Replaces `log_info.summary = analysis`, which made summary_file a
    byte-identical copy of the report on 40 of 40 outputs across two
    passes. This is a summary of the LOG -- what ran, on what, what came
    out, what failed -- built without an LLM, so it cannot inherit the
    formatting drift that broke the old scrape or the runaway that
    turned a 312 KB log into 2 MB.

    Never raises: a summary is a convenience, and losing a report that
    already exists to a summariser fault would be a worse bug than the
    one this fixes.
    """
    try:
        lines, _structure = build(log_text, log_name, pending=False)
        return "\n".join(lines) + "\n"
    except Exception as error:                             # noqa: BLE001
        return ("Deterministic summary unavailable: %s\n"
                "The full analysis is in the analysis file.\n"
                % type(error).__name__)


def main():
    G = os.path.join(HERE, "test_corpus", "corpus_gui")
    import time
    worst = 0
    print("%-32s %6s %6s" % ("log", "lines", "secs"))
    for f in sorted(os.listdir(G)):
        text = open(os.path.join(G, f), errors="replace").read()
        t0 = time.time()
        lines, _ = build(text, f)
        dt = time.time() - t0
        worst = max(worst, len(lines))
        flag = "" if len(lines) <= LINE_BUDGET else "  OVER BUDGET"
        print("%-32s %6d %6.3f%s" % (f[:-4], len(lines), dt, flag))
    print("\nworst case %d lines against a budget of %d" % (worst, LINE_BUDGET))

    print("\n" + "=" * 66)
    for name in ("molprobity_4_esterase", "refine_3_esterase",
                 "PHASER_beta_blip"):
        text = open(os.path.join(G, name + ".log"), errors="replace").read()
        lines, _ = build(text, name + ".log")
        print("\n--- %s (%d lines) ---" % (name, len(lines)))
        print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
