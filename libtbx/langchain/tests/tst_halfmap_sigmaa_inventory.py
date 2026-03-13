#!/usr/bin/env python3
"""
Tests for half-map pair detection, sigma-A MTZ classification,
file inventory in expert prompt, and anomalous gate logic.

Run:
  python3 tests/tst_halfmap_sigmaa_inventory.py
"""
from __future__ import absolute_import, division, print_function
import os
import sys
import unittest

# ---------------------------------------------------------------------------
# Import helpers — work from both the dev tree and installed PHENIX
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
  sys.path.insert(0, _ROOT)


# =========================================================================
# Bug 1: Half-map pair detection
# =========================================================================

class TestHalfMapPairDetection(unittest.TestCase):
  """Verify that numbered CCP4 file pairs are promoted to half_map."""

  def _categorize(self, filenames):
    """Thin wrapper around _categorize_files."""
    from agent.workflow_state import _categorize_files
    # Provide full paths (the function uses os.path.basename)
    paths = [os.path.join("/data", fn) for fn in filenames]
    return _categorize_files(paths, files_local=False)

  def test_numbered_pair_detected(self):
    """Two _1/_2 ccp4 files with companion full map → half_map."""
    files = self._categorize([
      "7mjs_23883_H_1.ccp4",
      "7mjs_23883_H_2.ccp4",
      "7mjs_23883_H_full.ccp4",
      "7mjs_23883_H.fa",
    ])
    hm = [os.path.basename(f)
          for f in files.get("half_map", [])]
    self.assertEqual(sorted(hm),
      ["7mjs_23883_H_1.ccp4", "7mjs_23883_H_2.ccp4"],
      "Numbered pair should be in half_map")
    fm = [os.path.basename(f)
          for f in files.get("full_map", [])]
    self.assertEqual(fm, ["7mjs_23883_H_full.ccp4"],
      "Companion full map should remain in full_map")

  def test_explicit_half_still_works(self):
    """Files with 'half' in name still detected."""
    files = self._categorize([
      "half_map_1.ccp4",
      "half_map_2.ccp4",
      "seq.fa",
    ])
    hm = [os.path.basename(f)
          for f in files.get("half_map", [])]
    self.assertEqual(len(hm), 2)

  def test_three_numbered_maps_pair_promoted(self):
    """Three numbered maps: _1/_2 pair promoted, _3 stays
    (acts as companion full map)."""
    files = self._categorize([
      "map_1.ccp4",
      "map_2.ccp4",
      "map_3.ccp4",
      "seq.fa",
    ])
    hm = [os.path.basename(f)
          for f in files.get("half_map", [])]
    fm = [os.path.basename(f)
          for f in files.get("full_map", [])]
    # _1/_2 pair promoted (map_3 is the companion)
    self.assertEqual(sorted(hm),
      ["map_1.ccp4", "map_2.ccp4"])
    self.assertIn("map_3.ccp4", fm)

  def test_pair_with_full_map_promoted(self):
    """Pair + companion full_map → pair promoted, full
    stays."""
    files = self._categorize([
      "data_1.ccp4",
      "data_2.ccp4",
      "data_full.ccp4",
      "seq.fa",
    ])
    hm = [os.path.basename(f)
          for f in files.get("half_map", [])]
    fm = [os.path.basename(f)
          for f in files.get("full_map", [])]
    # _1/_2 pair should be promoted (companion full exists)
    self.assertEqual(sorted(hm),
      ["data_1.ccp4", "data_2.ccp4"])
    self.assertIn("data_full.ccp4", fm)

  def test_pair_without_companion_not_promoted(self):
    """Two _1/_2 maps without a companion full_map →
    stay in full_map (could be segmented maps)."""
    files = self._categorize([
      "resolve_cryo_em_1.ccp4",
      "resolve_cryo_em_2.ccp4",
      "seq.fa",
    ])
    hm = [os.path.basename(f)
          for f in files.get("half_map", [])]
    fm = [os.path.basename(f)
          for f in files.get("full_map", [])]
    self.assertEqual(hm, [],
      "Without companion full map, pair should NOT "
      "be promoted to half_map")
    self.assertEqual(len(fm), 2)

  def test_no_false_positive_on_non_map(self):
    """Non-map files with _1 suffix are not affected."""
    files = self._categorize([
      "model_1.pdb",
      "model_2.pdb",
      "data.mtz",
      "seq.fa",
    ])
    hm = files.get("half_map", [])
    self.assertEqual(hm, [],
      "PDB files should never be in half_map")

  def test_single_numbered_map_stays_full(self):
    """A single _1 map with no _2 partner stays in full_map."""
    files = self._categorize([
      "scan_1.ccp4",
      "seq.fa",
    ])
    fm = [os.path.basename(f)
          for f in files.get("full_map", [])]
    hm = [os.path.basename(f)
          for f in files.get("half_map", [])]
    self.assertIn("scan_1.ccp4", fm)
    self.assertEqual(hm, [])

  def test_explicit_half_blocks_pair_promotion(self):
    """When real half-maps (with 'half' in name) exist,
    numbered maps are NOT promoted even with companion."""
    files = self._categorize([
      "half_map_1.ccp4",
      "half_map_2.ccp4",
      "data_1.ccp4",
      "data_2.ccp4",
      "data_full.ccp4",
      "seq.fa",
    ])
    hm = [os.path.basename(f)
          for f in files.get("half_map", [])]
    fm = [os.path.basename(f)
          for f in files.get("full_map", [])]
    # Only the explicit half-maps should be in half_map
    self.assertEqual(sorted(hm),
      ["half_map_1.ccp4", "half_map_2.ccp4"])
    # The numbered pair + companion stay in full_map
    self.assertIn("data_1.ccp4", fm)
    self.assertIn("data_2.ccp4", fm)
    self.assertIn("data_full.ccp4", fm)


# =========================================================================
# Bug 3a: sigma-A MTZ classification
# =========================================================================

class TestSigmaAClassification(unittest.TestCase):
  """sigma-A MTZ files stay in data_mtz (they contain
  both Fobs and map coefficients).  We do NOT move them
  to map_coeffs_mtz because phenix.refine's data_mtz
  slot has exclude_categories: [map_coeffs_mtz] — dual
  membership would make refine unable to find its input
  data.  The file inventory annotates the filename for
  the expert instead."""

  def _categorize(self, filenames):
    from agent.workflow_state import _categorize_files
    paths = [os.path.join("/data", fn) for fn in filenames]
    return _categorize_files(paths, files_local=False)

  def test_sigmaa_stays_in_data_mtz(self):
    """a2u-sigmaa.mtz stays in data_mtz (needed by
    xtriage and refine)."""
    files = self._categorize([
      "a2u-sigmaa.mtz",
      "model.pdb",
      "seq.fa",
    ])
    dm = [os.path.basename(f)
          for f in files.get("data_mtz", [])]
    self.assertIn("a2u-sigmaa.mtz", dm,
      "sigmaa MTZ must remain in data_mtz")

  def test_plain_data_mtz_unaffected(self):
    """Regular data MTZ stays in data_mtz."""
    files = self._categorize([
      "data.mtz",
      "seq.fa",
    ])
    dm = [os.path.basename(f)
          for f in files.get("data_mtz", [])]
    self.assertIn("data.mtz", dm)


# =========================================================================
# Bug 3c: File inventory in expert prompt
# =========================================================================

class TestFileInventoryInPrompt(unittest.TestCase):
  """Verify the expert prompt includes categorized file info."""

  def test_inventory_rendered(self):
    """build_thinking_prompt includes AVAILABLE FILES."""
    from knowledge.thinking_prompts import (
      build_thinking_prompt)
    context = {
      "cycle_number": 2,
      "experiment_type": "xray",
      "workflow_state": "xray_analyzed",
      "valid_programs": ["phenix.refine"],
      "program_name": "phenix.xtriage",
      "metrics": {},
      "user_advice": "Solve the structure",
      "history_summary": "",
      "r_free_trend": [],
      "file_inventory": (
        "  Models: mup_mr_solution.pdb (model)\n"
        "  Reflection data: a2u-sigmaa.mtz "
        "(map coeffs mtz)\n"
        "  Sequences: a2u-globulin.seq (sequence)"
      ),
    }
    _sys, user = build_thinking_prompt(context)
    self.assertIn("AVAILABLE FILES", user,
      "Prompt should contain AVAILABLE FILES header")
    self.assertIn("mup_mr_solution.pdb", user,
      "Prompt should show the model filename")
    self.assertIn("map coeffs mtz", user,
      "Prompt should show the MTZ category")

  def test_no_inventory_when_empty(self):
    """No AVAILABLE FILES section when inventory is empty."""
    from knowledge.thinking_prompts import (
      build_thinking_prompt)
    context = {
      "cycle_number": 1,
      "experiment_type": "xray",
      "workflow_state": "xray_initial",
      "valid_programs": ["phenix.xtriage"],
      "program_name": "",
      "metrics": {},
      "user_advice": "",
      "history_summary": "",
      "r_free_trend": [],
      "file_inventory": "",
    }
    _sys, user = build_thinking_prompt(context)
    self.assertNotIn("AVAILABLE FILES", user)


class TestFileInventoryContext(unittest.TestCase):
  """Verify _build_thinking_context populates file_inventory."""

  def test_context_has_inventory(self):
    """With categorized_files in workflow_state, context
    has a non-empty file_inventory."""
    try:
      from agent.thinking_agent import (
        _build_thinking_context)
    except (ImportError, ModuleNotFoundError):
      self.skipTest(
        "thinking_agent requires strategy_memory "
        "(not in standalone tarball)")
    state = {
      "workflow_state": {
        "state": "xray_analyzed",
        "valid_programs": ["phenix.refine"],
        "categorized_files": {
          "model": ["/data/placed.pdb"],
          "sequence": ["/data/seq.fa"],
          "data_mtz": ["/data/obs.mtz"],
          "map_coeffs_mtz": [
            "/data/sigmaa.mtz"],
          "half_map": [],
          "full_map": [],
          "ligand_pdb": [],
          "ligand_cif": [],
          "refined": [],
          "phaser_output": [],
          "predicted": [],
          "rsr_output": [],
          "autobuild_output": [],
          "phased_data_mtz": [],
          "refine_map_coeffs": [],
          "denmod_map_coeffs": [],
        },
      },
      "session_info": {"experiment_type": "xray"},
      "history": [],
      "log_text": "",
      "log_analysis": {},
      "user_advice": "",
      "cycle_number": 1,
      "directives": {},
    }
    ctx = _build_thinking_context(state, "advanced")
    inv = ctx.get("file_inventory", "")
    self.assertIn("placed.pdb", inv,
      "Inventory should contain model filename")
    self.assertIn("sigmaa.mtz", inv,
      "Inventory should contain map_coeffs filename")
    self.assertIn("seq.fa", inv,
      "Inventory should contain sequence filename")


# =========================================================================
# Bug 2B: Anomalous gate guard — real gate_evaluator tests
# =========================================================================

class TestAnomalousGateGuard(unittest.TestCase):
  """Gate evaluator should stop before experimental_phasing
  when anomalous signal is negligible."""

  def _make_stage(self, stage_id, cycles_used=0):
    """Create a minimal stage-like object."""
    class FakeStage:
      def __init__(self, sid, cu):
        self.id = sid
        self.cycles_used = cu
        self.max_cycles = 1
        self.skip_if = None
        self.success_criteria = {}
        self.gate_conditions = []
        self.fallbacks = []
        self.programs = ["phenix.autosol"]
        self.start_cycle = 1
      def is_exhausted(self):
        return self.cycles_used >= self.max_cycles
    return FakeStage(stage_id, cycles_used)

  def _make_plan(self, stage):
    """Create a minimal plan-like object."""
    class FakePlan:
      def __init__(self, s):
        self._stage = s
      def is_complete(self):
        return False
      def current_stage(self):
        return self._stage
    return FakePlan(stage)

  def _try_import(self):
    try:
      from agent.gate_evaluator import (
        GateEvaluator, _get_metric_value)
      return GateEvaluator, _get_metric_value
    except ImportError:
      self.skipTest(
        "gate_evaluator not importable")

  def test_weak_signal_stops(self):
    """measurability=0.03, cycles_used=0 → gate returns
    stop."""
    GateEvaluator, _ = self._try_import()
    stage = self._make_stage(
      "experimental_phasing", cycles_used=0)
    plan = self._make_plan(stage)
    sm = {
      "anomalous_measurability": 0.03,
      "has_anomalous": False,
    }
    ev = GateEvaluator()
    result = ev.evaluate(plan, sm, None, 2)
    self.assertEqual(result.action, "stop",
      "Weak anomalous should trigger stop")
    self.assertIn("0.030", result.reason)

  def test_weak_signal_after_first_cycle(self):
    """measurability=0.03, cycles_used=1 (autosol just
    ran) → gate still stops.  This is the realistic
    case: gate evaluates after first autosol attempt."""
    GateEvaluator, _ = self._try_import()
    stage = self._make_stage(
      "experimental_phasing", cycles_used=1)
    plan = self._make_plan(stage)
    sm = {
      "anomalous_measurability": 0.03,
      "has_anomalous": False,
    }
    ev = GateEvaluator()
    result = ev.evaluate(plan, sm, None, 3)
    self.assertEqual(result.action, "stop",
      "Weak anomalous should stop after first "
      "autosol attempt")

  def test_strong_signal_continues(self):
    """measurability=0.15 → gate continues."""
    GateEvaluator, _ = self._try_import()
    stage = self._make_stage(
      "experimental_phasing", cycles_used=0)
    plan = self._make_plan(stage)
    sm = {
      "anomalous_measurability": 0.15,
      "has_anomalous": True,
    }
    ev = GateEvaluator()
    result = ev.evaluate(plan, sm, None, 2)
    self.assertNotEqual(result.action, "stop",
      "Strong anomalous should not stop")

  def test_has_anomalous_flag_overrides(self):
    """measurability=0.04 but has_anomalous=True → no stop."""
    GateEvaluator, _ = self._try_import()
    stage = self._make_stage(
      "experimental_phasing", cycles_used=0)
    plan = self._make_plan(stage)
    sm = {
      "anomalous_measurability": 0.04,
      "has_anomalous": True,
    }
    ev = GateEvaluator()
    result = ev.evaluate(plan, sm, None, 2)
    self.assertNotEqual(result.action, "stop",
      "has_anomalous=True should override weak "
      "measurability")

  def test_no_metrics_continues(self):
    """No anomalous data → gate continues."""
    GateEvaluator, _ = self._try_import()
    stage = self._make_stage(
      "experimental_phasing", cycles_used=0)
    plan = self._make_plan(stage)
    sm = {}  # no anomalous metrics
    ev = GateEvaluator()
    result = ev.evaluate(plan, sm, None, 2)
    self.assertNotEqual(result.action, "stop")

  def test_missing_has_anomalous_continues(self):
    """measurability known but has_anomalous unknown →
    fire guard (absent has_anomalous is not True, so
    weak measurability blocks autosol).
    v115.05: Changed from 'require explicit False' to
    'has_anomalous is not True' to fix AF_exoV_PredictAndBuild
    where has_anomalous was absent from context."""
    GateEvaluator, _ = self._try_import()
    stage = self._make_stage(
      "experimental_phasing", cycles_used=0)
    plan = self._make_plan(stage)
    # Only measurability set, has_anomalous absent
    sm = {"anomalous_measurability": 0.02}
    ev = GateEvaluator()
    result = ev.evaluate(plan, sm, None, 2)
    self.assertEqual(result.action, "stop",
      "Guard should fire when has_anomalous is absent "
      "and measurability < 0.05 (AF_exoV fix)")

  def test_other_stage_not_affected(self):
    """Non-experimental_phasing stage unaffected."""
    GateEvaluator, _ = self._try_import()
    stage = self._make_stage(
      "data_assessment", cycles_used=0)
    plan = self._make_plan(stage)
    sm = {
      "anomalous_measurability": 0.01,
      "has_anomalous": False,
    }
    ev = GateEvaluator()
    result = ev.evaluate(plan, sm, None, 1)
    self.assertNotEqual(result.action, "stop",
      "Guard should only fire for "
      "experimental_phasing")

  def test_already_running_not_stopped(self):
    """Stage with cycles_used > 1 → no stop (guard
    only fires on first evaluation)."""
    GateEvaluator, _ = self._try_import()
    stage = self._make_stage(
      "experimental_phasing", cycles_used=2)
    plan = self._make_plan(stage)
    sm = {
      "anomalous_measurability": 0.01,
      "has_anomalous": False,
    }
    ev = GateEvaluator()
    result = ev.evaluate(plan, sm, None, 4)
    self.assertNotEqual(result.action, "stop",
      "Guard should not fire after first evaluation")


# =========================================================================
# Bug 2B: structure_model.get_metric for anomalous
# =========================================================================

class TestStructureModelAnomalousMetric(unittest.TestCase):
  """get_metric should return anomalous_measurability."""

  def test_anomalous_measurability(self):
    try:
      from agent.structure_model import (
        StructureModel)
    except ImportError:
      self.skipTest("structure_model not importable")
    sm = StructureModel()
    sm.data_characteristics["anomalous"] = {
      "has_anomalous_data": False,
      "anomalous_measurability": 0.032,
    }
    val = sm.get_metric("anomalous_measurability")
    self.assertAlmostEqual(val, 0.032, places=3)
    val2 = sm.get_metric("has_anomalous")
    self.assertEqual(val2, False)


# =========================================================================
# Bug 3b: MR solution file categorization
# =========================================================================

class TestMRSolutionDetection(unittest.TestCase):
  """PDB files with solution in the name → phaser_output."""

  def _categorize(self, filenames):
    from agent.workflow_state import _categorize_files
    paths = [os.path.join("/data", fn)
             for fn in filenames]
    return _categorize_files(paths, files_local=False)

  def test_mr_solution_pdb(self):
    """mup_mr_solution.pdb → phaser_output."""
    files = self._categorize([
      "mup_mr_solution.pdb",
      "data.mtz",
      "seq.fa",
    ])
    po = [os.path.basename(f)
          for f in files.get("phaser_output", [])]
    self.assertIn("mup_mr_solution.pdb", po)

  def test_generic_pdb_not_phaser_output(self):
    """model.pdb → NOT phaser_output."""
    files = self._categorize([
      "model.pdb",
      "data.mtz",
      "seq.fa",
    ])
    po = [os.path.basename(f)
          for f in files.get("phaser_output", [])]
    self.assertNotIn("model.pdb", po)

  def test_nmr_solution_not_phaser_output(self):
    """nmr_solution.pdb → NOT phaser_output (no 'mr_solution')."""
    files = self._categorize([
      "nmr_solution.pdb",
      "data.mtz",
      "seq.fa",
    ])
    po = [os.path.basename(f)
          for f in files.get("phaser_output", [])]
    self.assertNotIn("nmr_solution.pdb", po,
      "nmr_solution should not match mr_solution")


# =========================================================================
# Bug 3b: Plan context placement detection
# =========================================================================

class TestPlanContextPlacement(unittest.TestCase):
  """Plan context should detect placed models from filenames."""

  def test_mr_solution_sets_placed(self):
    try:
      from agent.plan_generator import _build_context
    except ImportError:
      self.skipTest("plan_generator not importable")
    ctx = _build_context(
      available_files=[
        "/data/mup_mr_solution.pdb",
        "/data/a2u-sigmaa.mtz",
        "/data/seq.fa",
      ],
    )
    self.assertTrue(ctx["model_is_placed"],
      "MR solution + sigmaa MTZ should set "
      "model_is_placed")

  def test_generic_model_not_placed(self):
    try:
      from agent.plan_generator import _build_context
    except ImportError:
      self.skipTest("plan_generator not importable")
    ctx = _build_context(
      available_files=[
        "/data/model.pdb",
        "/data/data.mtz",
        "/data/seq.fa",
      ],
    )
    self.assertFalse(ctx["model_is_placed"],
      "Generic PDB + data MTZ should NOT set "
      "model_is_placed")

  def test_nmr_solution_not_placed(self):
    """nmr_solution.pdb should NOT trigger model_is_placed."""
    try:
      from agent.plan_generator import _build_context
    except ImportError:
      self.skipTest("plan_generator not importable")
    ctx = _build_context(
      available_files=[
        "/data/nmr_solution.pdb",
        "/data/data.mtz",
        "/data/seq.fa",
      ],
    )
    self.assertFalse(ctx["model_is_placed"],
      "nmr_solution should not match mr_solution "
      "pattern")


# =========================================================================
# Runner
# =========================================================================

def run():
  loader = unittest.TestLoader()
  suite = unittest.TestSuite()
  for cls in (
    TestHalfMapPairDetection,
    TestSigmaAClassification,
    TestFileInventoryInPrompt,
    TestFileInventoryContext,
    TestAnomalousGateGuard,
    TestStructureModelAnomalousMetric,
    TestMRSolutionDetection,
    TestPlanContextPlacement,
  ):
    suite.addTests(loader.loadTestsFromTestCase(cls))
  runner = unittest.TextTestRunner(verbosity=2)
  result = runner.run(suite)
  return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
  sys.exit(run())
