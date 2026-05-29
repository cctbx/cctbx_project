"""K_H15_ITEM_3: Reasoning/command divergence detection (v119.H15).

Tom's bromodomain cycle 7 had a classic LLM hallucination:
the reasoning said "I will use the MTZ file from the last
refinement" but the command used `7qz0.mtz` (raw input data,
no R-free flags).  Refinement failed.

Per H15 plan with Gemini critique integrated, Item 3 ships
this as a DETECT-ONLY measurement instrument.  The proper
architectural fix (constrained decoding / two-pass
generation) is in the post-v119 backlog.  This function
gives us production telemetry to scope that work.

Per H15 plan question 3: detect-only EXCEPT when the
divergent file literally doesn't exist on disk and would
cause OS-level FileNotFoundError — that case blocks.

6 tests:
  §A: Tom's exact case → suspect_count=1
  §B: Reasoning and command agree → suspect_count=0
  §C: No categorical phrase → no detection
  §D: Negation cue suppresses detection
  §E: FileNotFound case → blocking_error populated
  §F: No active_files context → graceful no-op
"""
from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_HERE)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from agent.command_postprocessor import (
    cross_check_reasoning_vs_command,
)


# =====================================================================
# §A: Tom's exact case
# =====================================================================

def test_toms_case_detects_last_refinement_mismatch():
    """Reproduce Tom's cycle-7 divergence exactly.

    Reasoning mentions "the MTZ file from the last refinement"
    but the command uses 7qz0.mtz (raw input).  active_files
    includes both the raw input and the refine output.
    Expected: suspect_count == 1 (one divergence detected)."""
    reasoning = (
        "The ligand has been successfully fitted and combined "
        "with the protein model. The next essential step is to "
        "refine this new protein-ligand complex against the "
        "experimental data to optimize its geometry and improve "
        "the R-factors. The current R-free of 0.272 is reasonable "
        "but can be lowered further. I will use the combined "
        "model from the previous pdbtools step and the MTZ file "
        "from the last refinement to ensure consistency of the "
        "R-free flags."
    )
    command = (
        "phenix.refine "
        "/Users/x/sub_06_pdbtools/refine_001_001_modified.pdb "
        "/Users/x/7qz0.mtz "
        "ordered_solvent=True output.prefix=refine_002"
    )
    active_files = [
        "/Users/x/7qz0.mtz",
        "/Users/x/sub_03_refine/refine_001_001.mtz",
    ]
    logs = []
    suspect_count, blocking_error = cross_check_reasoning_vs_command(
        reasoning=reasoning,
        command=command,
        program="refine",
        active_files=active_files,
        best_files={},
        # Existing-file check returns True for both — block-on-FNF
        # won't fire.  We're only testing detection here.
        file_exists_check=lambda p: True,
        log=logs.append,
    )
    assert suspect_count == 1, (
        "Expected suspect_count=1 (one divergence), got %d.  Logs: %r"
        % (suspect_count, logs))
    assert blocking_error is None, (
        "File exists per the check stub — no block.  Got: %r"
        % blocking_error)
    # The log should have one [DIVERGENCE] line
    div_logs = [L for L in logs if "[DIVERGENCE]" in L]
    assert len(div_logs) >= 1, (
        "Expected at least 1 [DIVERGENCE] log line, got %d: %r"
        % (len(div_logs), logs))
    print("  PASS: test_toms_case_detects_last_refinement_mismatch")


# =====================================================================
# §B: Reasoning and command agree
# =====================================================================

def test_reasoning_and_command_agree():
    """When the command USES the file the reasoning describes,
    no divergence.  Same reasoning text, but command uses
    refine_001_001.mtz (a refine output)."""
    reasoning = (
        "I will use the MTZ file from the last refinement."
    )
    command = (
        "phenix.refine "
        "/Users/x/sub_03_refine/refine_001_001.pdb "
        "/Users/x/sub_03_refine/refine_001_001.mtz"
    )
    active_files = [
        "/Users/x/7qz0.mtz",
        "/Users/x/sub_03_refine/refine_001_001.mtz",
    ]
    logs = []
    suspect_count, blocking_error = cross_check_reasoning_vs_command(
        reasoning=reasoning,
        command=command,
        active_files=active_files,
        file_exists_check=lambda p: True,
        log=logs.append,
    )
    assert suspect_count == 0, (
        "Refine output IS in command — no divergence; got %d.  "
        "Logs: %r" % (suspect_count, logs))
    assert blocking_error is None
    print("  PASS: test_reasoning_and_command_agree")


# =====================================================================
# §C: No categorical phrase → no detection
# =====================================================================

def test_no_categorical_phrase_no_detection():
    """When reasoning doesn't use one of the recognized
    categorical phrases, the detector doesn't try to match.
    Generic phrases like "the MTZ file" without qualifier are
    too vague — don't flag."""
    reasoning = (
        "I will run phenix.refine with the standard inputs. "
        "The MTZ file is provided and the model is in place."
    )
    command = (
        "phenix.refine /path/to/model.pdb /path/to/data.mtz"
    )
    suspect_count, blocking_error = cross_check_reasoning_vs_command(
        reasoning=reasoning,
        command=command,
        active_files=["/path/to/data.mtz"],
        file_exists_check=lambda p: True,
        log=None,
    )
    assert suspect_count == 0, (
        "Generic reasoning has no categorical phrase — should "
        "not flag.  Got %d." % suspect_count)
    print("  PASS: test_no_categorical_phrase_no_detection")


# =====================================================================
# §D: Negation cue suppresses
# =====================================================================

def test_negation_cue_suppresses_detection():
    """When the reasoning explicitly says "instead of" or
    "rather than", the LLM is deliberately overriding the
    obvious choice.  Don't flag — false positives erode
    trust in the [DIVERGENCE] signal."""
    reasoning = (
        "I will use the raw experimental data instead of the "
        "last refinement MTZ, because I need to regenerate "
        "the R-free flags from scratch."
    )
    command = (
        "phenix.refine /path/to/model.pdb /path/to/7qz0.mtz "
        "xray_data.r_free_flags.generate=True"
    )
    active_files = [
        "/path/to/7qz0.mtz",
        "/path/to/sub_03_refine/refine_001_001.mtz",
    ]
    suspect_count, blocking_error = cross_check_reasoning_vs_command(
        reasoning=reasoning,
        command=command,
        active_files=active_files,
        file_exists_check=lambda p: True,
        log=None,
    )
    assert suspect_count == 0, (
        "Negation cue ('instead of') should suppress detection; "
        "got %d." % suspect_count)
    print("  PASS: test_negation_cue_suppresses_detection")


# =====================================================================
# §E: FileNotFound case → blocking_error populated
# =====================================================================

def test_filenotfound_blocks_command():
    """When the divergent file doesn't exist on disk,
    blocking_error is populated.  This is the one case Phase 1
    blocks on (per H15 plan question 3): running would cause
    OS-level FileNotFoundError, so refuse pre-emptively."""
    reasoning = (
        "I will use the MTZ file from the last refinement."
    )
    command = (
        "phenix.refine /Users/x/model.pdb "
        "/Users/x/does_not_exist.mtz"
    )
    logs = []
    suspect_count, blocking_error = cross_check_reasoning_vs_command(
        reasoning=reasoning,
        command=command,
        active_files=["/Users/x/sub_03_refine/refine_001_001.mtz"],
        # File doesn't exist (lambda always returns False)
        file_exists_check=lambda p: False,
        log=logs.append,
    )
    assert suspect_count >= 1, (
        "Should detect the divergence, got %d" % suspect_count)
    assert blocking_error is not None, (
        "Missing file should set blocking_error; got None.  Logs: %r"
        % logs)
    assert "does not exist" in blocking_error, (
        "blocking_error should describe the missing file; got %r"
        % blocking_error)
    # Should have logged the block
    block_logs = [L for L in logs if "[DIVERGENCE_BLOCK]" in L]
    assert len(block_logs) >= 1, (
        "Expected at least 1 [DIVERGENCE_BLOCK] log line, got %d"
        % len(block_logs))
    print("  PASS: test_filenotfound_blocks_command")


# =====================================================================
# §F: No active_files context → graceful no-op
# =====================================================================

def test_no_active_files_context_graceful():
    """Without active_files context, we can't disprove a file
    choice.  The detector should be conservative and NOT flag
    (no information → no claim).  This avoids false positives
    in early-cycle situations where active_files isn't
    populated yet."""
    reasoning = (
        "I will use the MTZ file from the last refinement."
    )
    command = (
        "phenix.refine /path/model.pdb /path/some_data.mtz"
    )
    suspect_count, blocking_error = cross_check_reasoning_vs_command(
        reasoning=reasoning,
        command=command,
        # No active_files, no best_files
        active_files=None,
        best_files=None,
        file_exists_check=lambda p: True,
        log=None,
    )
    # Without context, we can't prove divergence — no flag.
    # (Note: the heuristic in _file_matches_category for
    # last_refine_mtz checks basename 'refine' first; here
    # the basename is "some_data.mtz" which has no 'refine'.
    # Without active_files to confirm divergence, the function
    # falls through to "assume match" → no flag.)
    assert suspect_count == 0, (
        "Without active_files context, should be conservative; "
        "got %d." % suspect_count)
    assert blocking_error is None
    print("  PASS: test_no_active_files_context_graceful")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    test_toms_case_detects_last_refinement_mismatch()
    test_reasoning_and_command_agree()
    test_no_categorical_phrase_no_detection()
    test_negation_cue_suppresses_detection()
    test_filenotfound_blocks_command()
    test_no_active_files_context_graceful()


if __name__ == "__main__":
    print("K_H15_ITEM_3: Reasoning/command divergence detection (v119.H15)")
    print("=" * 70)
    run_all_tests()
    print("=" * 70)
    print("K_H15_ITEM_3 complete.")
