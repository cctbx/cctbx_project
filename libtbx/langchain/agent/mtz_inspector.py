"""MTZ inspection and per-program label selection (v119.H16).

This module replaces the agent's LLM-guessed crystallographic column
labels with deterministic cctbx-derived selection.  Two components:

1. ``inspect_mtz(mtz_path)`` — reads an MTZ via
   ``iotbx.reflection_file_reader`` and returns a structured dict
   describing its column structure (which arrays are intensities vs.
   amplitudes, which are anomalous, what the R-free column is).

2. ``select_obs_labels_for(program, mtz_info)`` — applies a
   per-program preference policy (refine prefers intensities, phaser
   prefers amplitudes, autosol prefers anomalous pairs) and returns
   the label string to inject into the command.

These functions are designed to be safe defaults: any failure path
returns ``None`` and the caller (the ``_apply_invariants()`` branch
for ``auto_fill_obs_labels`` in command_builder.py) silently skips
injection.  Pre-existing PHENIX behavior is preserved when this
module's outputs are absent.

Trigger for shipping this: 88 TIER-1 failures in two batch scans
matching "Sorry: Multiple equally suitable arrays of observed xray
data found", concentrated in AF_exoV_MRSAD and lysozyme-MRSAD
tutorials.  Their MTZs carry both anomalous pairs (I(+)/I(-)) and
merged intensities (IMEAN) and merged amplitudes (FOBS); PHENIX
refuses to guess and aborts.

The deterministic policy here picks correctly without LLM mediation.
"""

from __future__ import absolute_import, division, print_function

import logging
import os

logger = logging.getLogger(__name__)


# =====================================================================
# MTZ inspection
# =====================================================================

def inspect_mtz(mtz_path):
    """Inspect an MTZ file via cctbx and return its column structure.

    Returns a dict shaped like::

        {
          "anomalous_intensities": ["I(+),SIGI(+),I(-),SIGI(-)"],
          "anomalous_amplitudes":  [],
          "merged_intensities":    ["IMEAN,SIGIMEAN"],
          "merged_amplitudes":     ["FOBS,SIGFOBS"],
          "rfree_label":           "FreeR_flag",
          "is_multi_array":        True,
        }

    Each of the four data-category fields holds a list of comma-joined
    label strings.  Each string is in the form PHENIX expects after
    ``xray_data.labels=`` — so ``[0]`` of any category can be injected
    verbatim.

    ``is_multi_array`` is True when the file contains 2+ suitable data
    array sets across all categories.  This is the condition under
    which PHENIX refuses to guess and the H16 fix needs to fire.

    Returns ``None`` when the file is missing/unreadable, when cctbx
    isn't importable, or when the file isn't an MTZ.  Never raises.

    Args:
      mtz_path: str or pathlib.Path.  The MTZ file to inspect.

    Returns:
      dict with the structure above, or None.
    """
    if not mtz_path:
        return None
    mtz_path = str(mtz_path)
    if not os.path.exists(mtz_path):
        return None

    try:
        from iotbx.reflection_file_reader import any_reflection_file
    except ImportError:
        # cctbx isn't available — graceful degradation.  Caller
        # treats None as "no info, skip injection".
        return None

    try:
        rfile = any_reflection_file(mtz_path)
    except Exception:
        return None

    try:
        file_type = rfile.file_type() if hasattr(rfile, "file_type") else None
    except Exception:
        file_type = None
    # v119.H16.1 fix: the previous check `if file_type and
    # file_type != "ccp4_mtz"` had a logic hole — when cctbx returns
    # a falsy file_type (None for unrecognized formats), the
    # short-circuit `and` skipped the bail-out and we fell through
    # to as_miller_arrays() which returned an empty list, producing
    # a structured-but-empty result dict instead of None.  Tom's
    # production test on a non-MTZ text file caught this:
    # AssertionError "Non-MTZ file should return None; got
    # {'anomalous_intensities': [], ...}".  The fix is to require
    # exact "ccp4_mtz" match.
    if file_type != "ccp4_mtz":
        return None

    try:
        arrays = rfile.as_miller_arrays(merge_equivalents=False)
    except Exception:
        return None

    result = {
        "anomalous_intensities": [],
        "anomalous_amplitudes": [],
        "merged_intensities": [],
        "merged_amplitudes": [],
        "rfree_label": None,
        "is_multi_array": False,
    }

    for arr in arrays:
        try:
            entry = _classify_array(arr)
        except Exception:
            continue
        if entry is None:
            continue
        category, labels_str = entry
        if category == "rfree":
            # Take the FIRST r-free we encounter; multiple r-free
            # columns are rare and one is enough for the typical
            # use case.
            if result["rfree_label"] is None:
                result["rfree_label"] = labels_str
            continue
        result[category].append(labels_str)

    # is_multi_array: 2+ data array sets across the four data
    # categories.  When 0 or 1 set is present, PHENIX can pick
    # unambiguously and no injection is needed.
    data_array_count = (
        len(result["anomalous_intensities"])
        + len(result["anomalous_amplitudes"])
        + len(result["merged_intensities"])
        + len(result["merged_amplitudes"])
    )
    result["is_multi_array"] = data_array_count > 1

    return result


def _classify_array(arr):
    """Classify a Miller array into one of our 5 categories.

    Returns (category_str, labels_str) or None when the array isn't
    one we care about.  Categories are:
      - "anomalous_intensities"
      - "anomalous_amplitudes"
      - "merged_intensities"
      - "merged_amplitudes"
      - "rfree"

    ``labels_str`` is the comma-joined label list — the form PHENIX
    expects in ``xray_data.labels="..."``.
    """
    # Pull the labels first; we need them for rfree heuristic
    # regardless of the observation type.
    try:
        info = arr.info()
        labels = list(info.labels) if info and info.labels else []
    except Exception:
        labels = []
    if not labels:
        return None
    labels_str = ",".join(labels)

    # Check rfree by label name (more reliable than type-based
    # detection, since R-free is an integer flag array not an
    # "observation" array).
    if _is_rfree_labels(labels):
        return ("rfree", labels[0])

    # Determine anomalous flag
    try:
        anom = bool(arr.anomalous_flag())
    except Exception:
        anom = False

    # Determine intensity vs. amplitude.  cctbx provides
    # is_xray_intensity_array() and is_xray_amplitude_array()
    # which check the observation type.  If neither returns True,
    # this isn't a data array we should classify.
    is_intensity = False
    is_amplitude = False
    try:
        is_intensity = bool(arr.is_xray_intensity_array())
    except Exception:
        pass
    try:
        is_amplitude = bool(arr.is_xray_amplitude_array())
    except Exception:
        pass

    if is_intensity:
        return (
            "anomalous_intensities" if anom else "merged_intensities",
            labels_str,
        )
    if is_amplitude:
        return (
            "anomalous_amplitudes" if anom else "merged_amplitudes",
            labels_str,
        )
    return None


# Token heuristics for R-free detection.  cctbx doesn't always
# classify R-free arrays uniformly (it's an integer flag, not an
# observation), so we lean on the label naming convention.  False
# positives here are mild — an unrelated array with "free" in its
# label gets misclassified — but that's a rare case and doesn't
# affect the obs_labels injection path (R-free injection is not
# part of H16; rfree_label is for future use).
_RFREE_TOKENS = ("rfree", "r-free", "r_free", "freer", "free_r")


def _is_rfree_labels(labels):
    for label in labels:
        low = str(label).lower()
        for token in _RFREE_TOKENS:
            if token in low:
                return True
    return False


# =====================================================================
# Per-program label selection policy (H16 Item 3)
# =====================================================================

# Preference order per program.  Each list is the search order in
# inspect_mtz()'s result categories; the first non-empty category is
# the one we pick from, and within that category we use the first
# entry.
#
# Domain rationale (encoded so future readers don't have to re-derive):
#  - refine: intensities are statistically optimal for refinement
#    weighting; amplitudes are a fallback when intensities aren't
#    archived.  Anomalous data merges to give merged values too —
#    but the merged set is preferred over the raw anomalous pair
#    because refinement doesn't use the anomalous signal.
#  - phaser (MR): amplitudes (FOBS-style) are the canonical input
#    for molecular replacement search.  Intensities fall back.
#  - autosol (SAD/MRSAD): anomalous pairs are REQUIRED — the
#    program's job is to use the anomalous signal.  Amplitudes
#    fallback only if intensities aren't available; merged-only
#    data won't actually let autosol succeed but we still inject
#    something sensible so PHENIX produces its own error rather
#    than the multi-array Sorry.
#  - autobuild (scaling step): same as refine — intensities are
#    preferred for the integrated refinement steps.
_PROGRAM_PREFERENCES = {
    "refine": [
        "merged_intensities",
        "merged_amplitudes",
        "anomalous_intensities",
    ],
    "phaser": [
        "merged_amplitudes",
        "merged_intensities",
    ],
    # NOTE: autosol and autobuild are intentionally NOT in this
    # policy table.  Both programs use label-specification formats
    # that the agent's comma-pair convention doesn't cover (e.g.
    # space-separated, or multi-key, or program-specific PHIL
    # scopes that need their own value handlers).  Adding them
    # without verifying the format risks injecting malformed
    # commands.  Their multi-array Sorries persist post-H16; address
    # in a follow-up micro-fix that handles each program's actual
    # label syntax.
}


def select_obs_labels_for(program, mtz_info):
    """Pick the obs_labels string for a given program.

    Applies the program-specific preference list against the inspected
    MTZ structure.  Returns the first label string from the first
    non-empty preferred category.

    Args:
      program: str, the program name (with or without ``phenix.``
        prefix).  Variants like ``phenix.autobuild_denmod`` are
        normalized to their root program for policy lookup.
      mtz_info: dict from ``inspect_mtz()``, or None.

    Returns:
      str: the label string to inject (e.g. ``"IMEAN,SIGIMEAN"``).
      None: if mtz_info is missing, the program isn't in the policy
        table, or no preferred category has any entries.
    """
    if not mtz_info or not isinstance(mtz_info, dict):
        return None
    if not program:
        return None

    # Normalize the program name.  Strip phenix. prefix and any
    # variant suffix (autobuild_denmod → autobuild) so the policy
    # table only needs the root names.
    prog_key = str(program).replace("phenix.", "")
    # Take the part before the first underscore as the root
    # program name.  E.g. "autobuild_denmod" → "autobuild".
    root = prog_key.split("_", 1)[0]

    prefs = _PROGRAM_PREFERENCES.get(prog_key)
    if prefs is None:
        prefs = _PROGRAM_PREFERENCES.get(root)
    if prefs is None:
        return None

    for category in prefs:
        candidates = mtz_info.get(category, [])
        if candidates:
            return candidates[0]

    return None


# =====================================================================
# Helper for callers
# =====================================================================

def has_ambiguous_arrays(mtz_info):
    """Return True if the inspected MTZ has multiple suitable data
    array sets and would trigger PHENIX's "Multiple equally suitable
    arrays" Sorry.

    Convenience wrapper around mtz_info["is_multi_array"] with a None
    guard.  Used by the builder to gate injection: only fire when
    ambiguity exists.
    """
    if not mtz_info or not isinstance(mtz_info, dict):
        return False
    return bool(mtz_info.get("is_multi_array"))
