"""Raw advice file-detection scanner — v119.H4.1 (Step 1F).

Telemetry-gathering prerequisite for the deferred Phase 2B
preprocessor deprecation (see PHASE2B_PREPROCESSOR_DEPRECATION_PLAN.md
§5.1).  Used by the metric-logging block in
phenix_ai/run_ai_analysis.py::run_advice_preprocessing to record
(llm_files, regex_files) tuples for Phase 2B's regex_recall
threshold computation.

This module exposes one public function:

    scan_files_in_advice(raw_advice, file_list_hint=None) -> list[str]

Returns basenames of crystallography files found in `raw_advice`,
unioned with basenames extracted from `file_list_hint`.  Never
raises; pathological input returns [].

History:
  v119.H4    -- initial scanner with 21-extension allowlist and
                consuming boundary character classes.  Measured
                recall 0.806 on golden master corpus.
  v119.H4.1  -- expanded to 37 extensions based on telemetry
                analysis; switched to non-consuming lookarounds
                with period excluded from trailing lookahead.
                Measured recall 0.981 on golden master corpus.

When Phase 2B activates, this module will gain
`infer_experiment_type()` and `scan_raw_advice()` helpers and
become the canonical (LLM-free) preprocessor.  v119.H4.1 ships
only `scan_files_in_advice` since that's all the metric block
needs.

Regex shape is non-backtracking by construction (see plan rev 3
§3.2):
  - No nested quantifiers (single `+` over non-overlapping class)
  - No overlapping alternation (extensions are literal strings)
  - No `.*` flanking optionals
  - Lookarounds are fixed-width / zero-width assertions, O(1)
    per position in Python's `re` engine
K_H4 includes a required ReDoS perf test that asserts <50ms on
a 100KB pathological input.
"""
from __future__ import absolute_import, division, print_function

import os
import re


# Closed allowlist of crystallography file extensions.
#
# Initial H4 set covers the canonical 20 X-ray/cryo-EM extensions.
# H4.1 expansion (16 new) is data-driven: ranks of `in_llm_only`
# extension counts across a 567-entry telemetry corpus from the
# May 2026 tutorial sweep.  Each addition has empirical evidence
# of LLM extraction in real production tutorials.
#
# A noisier allowlist (e.g., adding `.txt` or `.log`) would
# inflate `in_regex_only` (false positives) but cannot improve
# recall against the LLM, since the LLM does not list those
# extensions either.  Stay narrow to crystallography-relevant
# files.
_FILE_EXTENSIONS = (
    # ----- H4 initial allowlist (21) -----
    # Crystallographic data
    'mtz', 'sca', 'hkl', 'cif', 'mmcif',
    # Map data (cryo-EM and X-ray)
    'map', 'mrc', 'ccp4',
    # Model files
    'pdb', 'ent',
    # Sequence files
    'fa', 'fasta', 'seq',
    # Phaser/refinement output
    'phs', 'eff', 'def',
    # Generic crystal data
    'xds', 'sad', 'mad', 'iso', 'peak',

    # ----- H4.1 additions (16) -----
    # Data-driven, ordered by miss-frequency on golden master.
    'dat',          # sequence/data files (+79 misses fixed)
    'inp',          # PHENIX legacy input parameter files (+72)
    'cv',           # cross-validation reflection data (+27)
    'com',          # script files (+27)
    'json',         # AlphaFold PAE data files (+22)
    'png',          # AlphaFold visualizations (+22)
    'ncs_spec',     # NCS specifications (+16)
    'sh', 'csh',    # shell scripts (+11, +7)
    'xplor',        # XPLOR map format (+9)
    'param',        # restraint parameter files (+9)
    'cxs',          # ChimeraX scene files (+7)
    'py',           # demo/runner Python scripts (+6)
    'gz',           # archives, matches `.tar.gz` via greedy prefix (+3)
    'pdf',          # documentation referenced in tutorials (+3)
    'list',         # list files like Facts.list (+2)
)


# Non-backtracking regex (plan rev 3 §3.2).
#
# Structure:
#   (?<![\w.\-])         -- leading boundary: NOT preceded by
#                           a filename character.  Negative
#                           lookbehind; fixed-width; O(1).
#   ([\w.\-]+\.(?:...))  -- captured filename: greedy filename
#                           body followed by literal allowlisted
#                           extension.
#   (?![\w\-])           -- trailing boundary: NOT followed by
#                           a word/hyphen character.  Note that
#                           `.` is INTENTIONALLY EXCLUDED here:
#                           this lets `file.eff.` (sentence-end
#                           period) match.  The trade-off is
#                           that `file.pdb.bak` will match
#                           `file.pdb`, but `.bak` files aren't
#                           real inputs and the LLM doesn't list
#                           them either, so this only adds
#                           tolerable noise to `in_regex_only`
#                           (precision), not `in_llm_only`
#                           (recall).
_FILE_PATTERN = re.compile(
    r'(?<![\w.\-])'
    r'([\w.\-]+\.(?:%s))'
    r'(?![\w\-])'
    % '|'.join(_FILE_EXTENSIONS),
    re.IGNORECASE,
)


def _safe_basename(p):
    """Return basename(normpath(p)), or '' for anything weird."""
    try:
        if not p:
            return ''
        if not isinstance(p, str):
            return ''
        bn = os.path.basename(os.path.normpath(p))
        return bn or ''
    except Exception:
        return ''


def scan_files_in_advice(raw_advice, file_list_hint=None):
    """Regex-based file detection over raw user advice.

    v119.H4.1 (Step 1F).  Phase 2B prerequisite.

    Args:
        raw_advice: User's raw advice text (pre-preprocessing).
            May be any type; non-string input returns [].
        file_list_hint: Optional iterable of file paths the agent
            already knows about (from the request's `files` field).
            Hint files are unioned into the result so the scanner
            matches Phase 2B's planned "regex + context" detector.
            None or empty iterable is fine.

    Returns:
        list[str]: Basenames of files found.  Sorted ASCII
            ascending (case-insensitive), deduped case-insensitively
            (first-seen casing wins).  Empty list if nothing
            matches.

    Never raises.  All exceptional inputs return [].
    """
    # Use a case-insensitive map (lowercased basename -> original
    # basename) so dedup is case-insensitive but emitted output
    # preserves the original casing as observed in the input.
    seen = {}

    # 1. Scan raw_advice for filename patterns.
    if isinstance(raw_advice, str) and raw_advice:
        try:
            for match in _FILE_PATTERN.finditer(raw_advice):
                token = match.group(1)
                bn = _safe_basename(token)
                if bn:
                    seen.setdefault(bn.lower(), bn)
        except Exception:
            pass

    # 2. Union with file_list_hint.
    if file_list_hint:
        try:
            for path in file_list_hint:
                bn = _safe_basename(path)
                if bn:
                    seen.setdefault(bn.lower(), bn)
        except Exception:
            pass

    # 3. Sort by lowercased basename for stable cross-platform order.
    return sorted(seen.values(), key=str.lower)
