"""
Fix Verification Tests
======================
Verify that v115.05 fixes are active in the installed code.
Each test checks a specific source pattern that distinguishes
the fixed code from the prior version.

Usage:
  libtbx.python tests/tst_fix_verification.py
"""
from __future__ import print_function
import sys, os

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_THIS = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_THIS)  # tests/ -> langchain/
if not os.path.isdir(os.path.join(_ROOT, "agent")):
    _PHENIX = os.environ.get("PHENIX", "")
    _ROOT = os.path.join(_PHENIX, "modules", "cctbx_project",
                         "libtbx", "langchain")

def _read(relpath):
    with open(os.path.join(_ROOT, relpath)) as f:
        return f.read()

def _read_abs(abspath):
    with open(abspath) as f:
        return f.read()


def run(verbose=False):
    """Run all fix verification tests.

    Returns True if all tests pass, False otherwise.
    """
    if not os.path.isdir(os.path.join(_ROOT, "agent")):
        print("ERROR: Cannot find agent/ directory.")
        print("  Tried: %s" % os.path.join(_ROOT, "agent"))
        print("  Set $PHENIX or run from the langchain/tests/ directory.")
        return False
    print("Using: %s" % _ROOT)

    passed = [0]
    failed = [0]

    def check(name, got, want_label):
        if got:
            print("  v115.05  %s" % name)
            passed[0] += 1
        else:
            print("  OLD CODE %s  <<<" % name)
            failed[0] += 1

    # ---------------------------------------------------------------------------
    # Test 1: workflow_engine.py — _is_at_target clashscore requires refine_count>=1
    #
    # OLD: clashscore < 10 with r_free=None → return True (even if refine_count=0)
    # NEW: clashscore < 10 with r_free=None → return True ONLY if refine_count >= 1
    # ---------------------------------------------------------------------------
    src_we = _read("agent/workflow_engine.py")

    # The new code has "refine_count >= 1" near the clashscore < 10 check
    idx_clash = src_we.find("clashscore < 10")
    assert idx_clash >= 0, "Cannot find 'clashscore < 10' in workflow_engine.py"
    block = src_we[idx_clash:idx_clash+200]
    check("clashscore guard requires refine_count >= 1",
          "refine_count >= 1" in block or "refine_count>=1" in block,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Test 2: workflow_engine.py — _is_at_target hopeless R-free requires autobuild_done
    #
    # OLD: r_free > 0.50 and refine_count >= 1 → return True
    # NEW: r_free > 0.50 and refine_count >= 1 → return True ONLY if autobuild_done
    # ---------------------------------------------------------------------------
    idx_hopeless = src_we.find("r_free > 0.50")
    assert idx_hopeless >= 0, "Cannot find 'r_free > 0.50' in workflow_engine.py"
    block2 = src_we[idx_hopeless:idx_hopeless+200]
    check("hopeless R-free requires autobuild_done",
          "autobuild_done" in block2,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Test 3: workflow_engine.py — Bug F routing requires autobuild_done
    #
    # OLD: r_free >= 0.45 + phaser_done + has_model_for_mr → obtain_model
    # NEW: same conditions BUT also requires autobuild_done
    # ---------------------------------------------------------------------------
    idx_bugf = src_we.find("Bug F fix")
    assert idx_bugf >= 0, "Cannot find 'Bug F fix' in workflow_engine.py"
    block3 = src_we[idx_bugf:idx_bugf+600]
    check("Bug F obtain_model requires autobuild_done",
          "autobuild_done" in block3,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Test 4: plan_generator.py — wants_polder context flag exists
    #
    # OLD: no wants_polder in context
    # NEW: wants_polder detected from advice or after_program
    # ---------------------------------------------------------------------------
    src_pg = _read("agent/plan_generator.py")
    check("plan_generator has wants_polder",
          "wants_polder" in src_pg,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Test 5: plan_generator.py — polder override fires despite "solve"
    #
    # OLD: "solve" in _mr_keywords blocks model_is_placed for polder
    # NEW: explicit polder override after the MR/placed check
    # ---------------------------------------------------------------------------
    # Find the polder override block
    idx_polder = src_pg.find("polder")
    has_override = False
    # Look for the pattern: if "polder" in advice → model_is_placed = True
    # This must appear AFTER the _mentions_mr/_mentions_placed check
    idx_mr_check = src_pg.find("_mentions_placed and not _mentions_mr")
    if idx_mr_check >= 0:
        after_check = src_pg[idx_mr_check:]
        # The override should set model_is_placed when polder is in advice
        has_override = ('if "polder" in advice' in after_check and
                        "model_is_placed" in after_check[after_check.find('if "polder" in advice'):
                                                         after_check.find('if "polder" in advice')+200])
    check("polder override fires despite 'solve' keyword",
          has_override,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Test 6: plan_templates.yaml — refine_placed_polder template exists
    # ---------------------------------------------------------------------------
    src_pt = _read("knowledge/plan_templates.yaml")
    check("refine_placed_polder template exists",
          "refine_placed_polder" in src_pt,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Test 7: plan_templates.yaml — early rebuild gate in mr_refine
    #
    # OLD: first gate is "r_free > 0.45 after 2 cycles"
    # NEW: first gate is "r_free > 0.50 after 1 cycles → try_rebuilding"
    # ---------------------------------------------------------------------------
    check("early rebuild gate (r_free > 0.50 after 1 cycles)",
          "r_free > 0.50 after 1 cycles" in src_pt,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Test 8: gate_evaluator.py — try_rebuilding advances to model_rebuilding
    #
    # OLD: try_rebuilding returns action="fallback" (weak hint, ignored)
    # NEW: try_rebuilding returns action="advance" when model_rebuilding exists
    # ---------------------------------------------------------------------------
    src_ge = _read("agent/gate_evaluator.py")
    idx_try = src_ge.find("try_rebuilding")
    assert idx_try >= 0, "Cannot find 'try_rebuilding' in gate_evaluator.py"
    block_try = src_ge[idx_try:idx_try+500]
    check("try_rebuilding advances to model_rebuilding stage",
          "model_rebuilding" in block_try and '"advance"' in block_try,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Test 9: workflow_engine.py — negligible anomalous guard removes autosol
    #
    # OLD: autosol stays in valid_programs even when measurability < 0.05
    # NEW: autosol removed when measurability < 0.05 and has_anomalous=False
    # ---------------------------------------------------------------------------
    idx_neg = src_we.find("Negligible-anomalous guard")
    check("negligible-anomalous guard removes autosol",
          idx_neg >= 0 and "0.05" in src_we[idx_neg:idx_neg+400],
          "v115.05")

    # ---------------------------------------------------------------------------
    # Tests 10-12: ai_agent.py fixes (programs/ai_agent.py)
    # These are in phenix/phenix/programs/ai_agent.py, not langchain.
    # ---------------------------------------------------------------------------
    _ai_agent_path = os.path.join(
        os.environ.get("PHENIX", ""),
        "modules", "phenix", "phenix", "programs", "ai_agent.py")
    if os.path.isfile(_ai_agent_path):
        src_ai = _read_abs(_ai_agent_path)

        # Test 10: unregistered explicit_program no longer raises Sorry
        idx_prog = src_ai.find("prog_info is None")
        if idx_prog >= 0:
            block_prog = src_ai[idx_prog:idx_prog+400]
            check("unregistered explicit_program downgrades to warning (no Sorry)",
                  "raise Sorry" not in block_prog,
                  "v115.05")
        else:
            check("unregistered explicit_program (not found)", False, "v115.05")

        # Test 11: preprocessing and needs-plan programs don't skip plan generation
        check("preprocessing/needs-plan programs generate full plan",
              "_preprocessing_programs" in src_ai and
              "_needs_plan_programs" in src_ai and
              "phenix.polder" in src_ai[
                  src_ai.find("_needs_plan_programs"):
                  src_ai.find("_needs_plan_programs")+200
              ],
              "v115.05")

        # Test 12: failed programs don't track output files
        check("skip output file tracking on FAILED programs",
              "skip_if_failed" in src_ai,
              "v115.05")

        # Test 14: metrics-based report selection (not INCOMPLETE when solved)
        check("report uses metrics to override INCOMPLETE status",
              "_metrics_good" in src_ai and "_best_rfree" in src_ai,
              "v115.05")
    else:
        print("  (skipped tests 10-12 — ai_agent.py not found at %s)" % _ai_agent_path)

    # ---------------------------------------------------------------------------
    # Test 13: graph_nodes.py — copies injection reads log_analysis
    #
    # OLD: BUILD only reads session_info["asu_copies"] (1-cycle delay)
    # NEW: BUILD also reads state["log_analysis"]["n_copies"] (same-cycle)
    # ---------------------------------------------------------------------------
    try:
        src_gn = _read("agent/graph_nodes.py")
        idx_copies = src_gn.find("PHASER COPIES INJECTION")
        if idx_copies >= 0:
            block_copies = src_gn[idx_copies:idx_copies+600]
            check("copies injection reads log_analysis (same-cycle)",
                  "log_analysis" in block_copies and "n_copies" in block_copies,
                  "v115.05")
        else:
            check("copies injection block (not found)", False, "v115.05")
    except Exception:
        print("  (skipped test 13 — graph_nodes.py not readable)")

    # ---------------------------------------------------------------------------
    # Test 15: workflow_state.py — no_ligand excluded from ligand classification
    # ---------------------------------------------------------------------------
    src_ws = _read("agent/workflow_state.py")
    check("no_ligand excluded from ligand file classification",
          "_anti_ligand_patterns" in src_ws and "no_ligand" in src_ws,
          "v115.05")

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    print()
    print("=" * 60)
    total = passed[0] + failed[0]
    if failed[0] == 0:
        print("  ALL %d TESTS CONFIRM v115.05 CODE IS INSTALLED" % total)
    else:
        print("  %d/%d tests show OLD CODE — v115.05 NOT fully installed" % (failed[0], total))
        print("  The fixes in the tarball are NOT active.")
    print("=" * 60)

    return failed[0] == 0

if __name__ == "__main__":
    success = run()
    sys.exit(0 if success else 1)
