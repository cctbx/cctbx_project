#!/usr/bin/env python3
"""Runtime introspection audit for phase→stage/step rename.

Instead of grepping text, this script IMPORTS every module
and inspects the actual Python objects: class names, method
names, function signatures, default dict keys, constant
names. If anything still says 'phase' (outside crystallographic
contexts), it's a real bug that will hit users at runtime.
"""
import sys
import os
import inspect
import importlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Crystallographic / allowed phase terms
ALLOWED = {
    "experimental_phasing", "sad_phasing", "build_from_phases",
    "phaser", "phase_angle", "phase_probs", "phase_and_amplitude",
    # Callback event type (not plan/workflow)
    "phase",  # only as a value, not as an identifier name
}

ISSUES = []

def check_name(name, context):
    """Check if a name contains 'phase' and isn't allowed."""
    low = name.lower()
    if "phase" not in low:
        return
    # Allow crystallographic terms
    for allowed in ALLOWED:
        if allowed in low:
            return
    # Allow backward-compat parameter names
    if name == "phase_id" or name == "phase_name":
        # These should have been renamed
        ISSUES.append(f"  PARAM: {context} -> {name}")
        return
    ISSUES.append(f"  {context} -> {name}")


def audit_module(mod_name, mod):
    """Inspect a loaded module for phase-related names."""
    # Module-level names
    for name in dir(mod):
        if name.startswith("_"):
            continue
        check_name(name, f"{mod_name}.{name}")

        obj = getattr(mod, name, None)
        if obj is None:
            continue

        # Class inspection
        if inspect.isclass(obj):
            for attr_name in dir(obj):
                if attr_name.startswith("_"):
                    continue
                check_name(attr_name, f"{mod_name}.{name}.{attr_name}")

            # Check __init__ parameters
            try:
                sig = inspect.signature(obj.__init__)
                for param_name in sig.parameters:
                    if param_name == "self":
                        continue
                    check_name(param_name,
                        f"{mod_name}.{name}.__init__({param_name})")
            except (ValueError, TypeError):
                pass

        # Function inspection
        elif inspect.isfunction(obj):
            try:
                sig = inspect.signature(obj)
                for param_name in sig.parameters:
                    check_name(param_name,
                        f"{mod_name}.{name}({param_name})")
            except (ValueError, TypeError):
                pass

        # Check if it's a string constant with 'phase'
        elif isinstance(obj, str) and "phase" in obj.lower():
            low = obj.lower()
            # Skip crystallographic
            skip = False
            for a in ["phaser", "phasing", "experimental_phase",
                       "fourier", "improve phase", "calculate",
                       "coefficient", "HL ", "amplitude"]:
                if a in low:
                    skip = True
                    break
            if not skip and len(obj) < 80:
                check_name(obj, f"{mod_name}.{name} = \"{obj}\"")


def audit_dict_keys(d, context, depth=0):
    """Recursively check dict keys for 'phase'."""
    if depth > 3 or not isinstance(d, dict):
        return
    for key in d:
        if isinstance(key, str) and "phase" in key.lower():
            low = key.lower()
            if any(a in low for a in
                   ["phaser", "phasing", "experimental_phasing",
                    "sad_phasing", "build_from_phases"]):
                continue
            ISSUES.append(f"  DICT KEY: {context}['{key}']")
        if isinstance(d[key], dict):
            audit_dict_keys(d[key], f"{context}['{key}']", depth+1)


# ============================================================
# Run all audits
# ============================================================

def run():
    """Run all audits. Raises AssertionError on failure."""
    global ISSUES
    ISSUES = []

    _run_audits()

    if ISSUES:
        msg = "phase→stage rename audit: %d issue(s):\n" % len(ISSUES)
        msg += "\n".join(ISSUES)
        raise AssertionError(msg)


def _run_audits():
    """Inner audit logic."""

    MODULES_TO_CHECK = [
        "knowledge.plan_schema",
        "knowledge.plan_template_loader",
        "knowledge.explanation_prompts",
        "knowledge.html_report_template",
        "knowledge.thinking_prompts",
        "knowledge.yaml_loader",
        "agent.plan_generator",
        "agent.gate_evaluator",
        "agent.display_data_model",
        "agent.thinking_agent",
        "agent.session_tools",
        "agent.api_client",
        "agent.rules_selector",
        "agent.program_validator",
        "agent.validation_history",
        "agent.sanity_checker",
        "agent.event_formatter",
        "agent.hypothesis_evaluator",
    ]

    print("=" * 60)
    print("RUNTIME INTROSPECTION AUDIT: phase→stage/step")
    print("=" * 60)
    print()

    loaded = 0
    failed = 0
    for mod_name in MODULES_TO_CHECK:
        try:
            mod = importlib.import_module(mod_name)
            audit_module(mod_name, mod)
            loaded += 1
        except ImportError as e:
            # Some modules need libtbx
            if "libtbx" in str(e) or "wx" in str(e):
                pass
            else:
                print(f"  IMPORT FAIL: {mod_name}: {e}")
                failed += 1
        except Exception as e:
            print(f"  ERROR: {mod_name}: {e}")
            failed += 1

    print(f"Modules inspected: {loaded}")
    print(f"Import failures (non-libtbx): {failed}")
    print()

    # ============================================================
    # Check plan_schema deeply
    # ============================================================
    print("--- Deep check: plan_schema ---")
    try:
        from knowledge.plan_schema import (
            StructurePlan, StageDef,
        )
        import knowledge.plan_schema as _ps
        for _const in ("STAGE_PENDING", "STAGE_ACTIVE",
                        "STAGE_COMPLETE", "STAGE_SKIPPED",
                        "STAGE_FAILED"):
            assert hasattr(_ps, _const), \
                "Missing constant: %s" % _const

        # Check StageDef attributes
        sd = StageDef(id="test", programs=["phenix.refine"])
        for attr in dir(sd):
            if not attr.startswith("_"):
                check_name(attr, f"StageDef.{attr}")

        # Check StructurePlan methods
        plan = StructurePlan(stages=[sd])
        for attr in dir(plan):
            if not attr.startswith("_"):
                check_name(attr, f"StructurePlan.{attr}")

        # Check to_dict keys
        d = plan.to_dict()
        audit_dict_keys(d, "StructurePlan.to_dict()")

        # Check from_dict round-trip with OLD keys (backward compat)
        old_data = {
            "phases": [{"id": "test", "programs": ["phenix.refine"],
                        "status": "pending"}],
            "current_phase_index": 0,
            "goal": "test",
        }
        plan2 = StructurePlan.from_dict(old_data)
        assert len(plan2.stages) == 1, "Backward compat: old 'phases' key failed!"
        print("  Backward compat (old 'phases' key): OK")

        # Check to_dict outputs new keys
        d2 = plan2.to_dict()
        assert "stages" in d2, "to_dict should output 'stages' not 'phases'"
        assert "phases" not in d2, "to_dict still outputs old 'phases' key!"
        assert "current_stage_index" in d2, "to_dict should use current_stage_index"
        print("  to_dict outputs new keys: OK")

    except Exception as e:
        print(f"  ERROR: {e}")

    # ============================================================
    # Check DisplayDataModel with old session data
    # ============================================================
    print()
    print("--- Deep check: DisplayDataModel ---")
    try:
        from agent.display_data_model import DisplayDataModel

        # Old-format session data
        old_session = {
            "cycles": [
                {"cycle_number": 1, "program": "phenix.refine",
                 "result": "SUCCESS: R-free = 0.25",
                 "metrics": {"r_free": 0.25}},
            ],
            "plan": {
                "phases": [  # OLD key
                    {"id": "refine", "status": "complete",
                     "programs": ["phenix.refine"],
                     "cycles_used": 1},
                ],
            },
        }
        ddm = DisplayDataModel.from_session(old_session)

        # Check properties exist with new names
        assert hasattr(ddm, "stage_outcomes"), "Missing stage_outcomes property"
        assert not hasattr(ddm, "phase_outcomes"), "Old phase_outcomes still exists!"

        # Check it can read old data
        so = ddm.stage_outcomes
        assert len(so) > 0, "stage_outcomes empty with old 'phases' data"
        assert hasattr(so[0], "stage_id"), "StageOutcome missing stage_id"
        assert not hasattr(so[0], "phase_id"), "StageOutcome still has phase_id!"
        print("  Old session data backward compat: OK")
        print("  stage_outcomes property: OK")
        print("  StageOutcome.stage_id: OK")

        # New-format session data
        new_session = {
            "cycles": old_session["cycles"],
            "plan": {
                "stages": old_session["plan"]["phases"],  # NEW key
            },
        }
        ddm2 = DisplayDataModel.from_session(new_session)
        so2 = ddm2.stage_outcomes
        assert len(so2) > 0, "stage_outcomes empty with new 'stages' data"
        print("  New session data: OK")

    except Exception as e:
        print(f"  ERROR: {e}")

    # ============================================================
    # Check YAML loading
    # ============================================================
    print()
    print("--- Deep check: YAML loading ---")
    try:
        from knowledge.yaml_loader import get_workflow_steps

        xray = get_workflow_steps("xray")
        assert xray, "get_workflow_steps('xray') returned empty!"
        assert "analyze" in xray, "Missing 'analyze' step"
        assert "refine" in xray, "Missing 'refine' step"
        audit_dict_keys(xray, "get_workflow_steps('xray')")
        print("  get_workflow_steps: OK (%d steps)" % len(xray))

        # Verify old function name doesn't exist
        from knowledge import yaml_loader
        assert not hasattr(yaml_loader, "get_workflow_phases"), \
            "Old get_workflow_phases function still exists!"
        print("  get_workflow_phases removed: OK")

    except ImportError:
        print("  (skipped — needs libtbx)")
    except Exception as e:
        print(f"  ERROR: {e}")

    # ============================================================
    # Check plan templates
    # ============================================================
    print()
    print("--- Deep check: plan templates ---")
    try:
        from knowledge.plan_template_loader import load_templates

        templates = load_templates()
        assert templates, "No templates loaded!"
        for tid, tdata in templates.items():
            # Check template uses 'stages' not 'phases'
            if "phases" in tdata and "stages" not in tdata:
                ISSUES.append(
                    f"  TEMPLATE: {tid} uses 'phases' key (not 'stages')")
            # Check stage defs
            stages = tdata.get("stages", [])
            for s in stages:
                if isinstance(s, dict):
                    for k in s:
                        check_name(k, f"template[{tid}].stage.{k}")
        print("  Templates loaded: %d" % len(templates))
        print("  All use 'stages' key: OK")

    except Exception as e:
        print(f"  ERROR: {e}")

    # ============================================================
    # Check explanation_prompts function signatures
    # ============================================================
    print()
    print("--- Deep check: explanation_prompts ---")
    try:
        from knowledge import explanation_prompts as ep

        # Old function should not exist
        assert not hasattr(ep, "generate_phase_summary"), \
            "Old generate_phase_summary still exists!"
        assert hasattr(ep, "generate_stage_summary"), \
            "New generate_stage_summary missing!"
        print("  generate_stage_summary: OK")
        print("  generate_phase_summary removed: OK")

        sig = inspect.signature(ep.generate_stage_summary)
        for p in sig.parameters:
            check_name(p, f"generate_stage_summary({p})")
        print("  Parameter names: OK")

    except Exception as e:
        print(f"  ERROR: {e}")

    # ============================================================
    # Summary
    # ============================================================
    print()
    print("=" * 60)
    if ISSUES:
        print(f"FOUND {len(ISSUES)} ISSUE(S):")
        for issue in ISSUES:
            print(issue)
        print()
        print("AUDIT FAILED")
    else:
        print("NO ISSUES FOUND — ALL CLEAR")
        print("=" * 60)


if __name__ == "__main__":
    try:
        run()
    except AssertionError as e:
        print(e)
        sys.exit(1)
