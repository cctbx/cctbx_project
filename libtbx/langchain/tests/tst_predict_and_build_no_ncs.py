"""K_H14_2: predict_and_build no-NCS-input config fix (v119.H14.2).

Tom's 2026-05-26 1029B-sad ollama run (run after the H14.1 ship,
with the PDB file removed so predict_and_build became the workflow
choice) surfaced a programs.yaml configuration bug:

The directory contained leftover find_ncs.ncs_spec from a prior
workflow execution.  The agent's auto-discovery picked it up, looked
up the .ncs_spec extension in programs.yaml, found that
phenix.predict_and_build had it documented under inputs.optional
with flag `map_model.ncs_file=`, and emitted that PHIL line on the
predict_and_build command:

  phenix.predict_and_build ... map_model.ncs_file=.../find_ncs.ncs_spec ...

predict_and_build rejected the command with "Some PHIL parameters
are not recognized."  Cycle-level retry re-emitted the same bad
parameter and the workflow looped.

Per Tom's domain confirmation: predict_and_build does NOT accept
any NCS-file PHIL parameter.  It runs NCS detection internally if
needed.  The fix is purely configuration — remove the ncs_spec
entry from predict_and_build's inputs.optional in programs.yaml,
and remove {ncs_spec} from its command template.

This file pins the configuration:

  - predict_and_build's inputs.optional must NOT contain ncs_spec
  - predict_and_build's command template must NOT reference {ncs_spec}
  - map_symmetry's hint must NOT advertise predict_and_build as a
    downstream consumer of .ncs_spec output (otherwise users may be
    misled into expecting NCS-injection support that doesn't exist)
  - resolve_cryo_em's ncs_spec entry must remain intact (it DOES
    accept ncs_file= — only predict_and_build is being corrected)

Total: 5 tests.
"""
from __future__ import absolute_import, division, print_function

import os


def _load_programs_yaml():
    """Load the programs.yaml that ships with this directive_extractor."""
    try:
        import yaml
    except ImportError:
        return None, "PyYAML not available"

    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        # Sandbox layout (tests/ is sibling to knowledge/)
        os.path.join(os.path.dirname(here), "knowledge", "programs.yaml"),
        # PHENIX deployment layout
        "/net/cci-filer3/home/terwill/phenix_2026-02-20/modules/"
        "cctbx_project/libtbx/langchain/knowledge/programs.yaml",
    ]
    for path in candidates:
        if os.path.exists(path):
            try:
                with open(path) as f:
                    return yaml.safe_load(f), None
            except Exception as e:
                return None, "Error loading %s: %s" % (path, e)
    return None, "programs.yaml not found in any expected location"


# =====================================================================
# §A: predict_and_build no longer declares ncs_spec
# =====================================================================

def test_predict_and_build_has_no_ncs_spec_input():
    """predict_and_build's inputs.optional must NOT contain
    ncs_spec.  Pre-H14.2 it declared
    `flag: "map_model.ncs_file="` which predict_and_build rejects
    as an unknown PHIL parameter.
    """
    data, err = _load_programs_yaml()
    if data is None:
        print("  SKIP: %s" % err)
        return

    pb = data.get("phenix.predict_and_build", {})
    optional_inputs = pb.get("inputs", {}).get("optional", {})

    assert "ncs_spec" not in optional_inputs, (
        "predict_and_build.inputs.optional MUST NOT contain ncs_spec "
        "(H14.2 fix — predict_and_build does not accept an NCS file). "
        "Got optional keys: %r" % list(optional_inputs.keys()))
    print("  PASS: test_predict_and_build_has_no_ncs_spec_input")


def test_predict_and_build_command_has_no_ncs_spec_placeholder():
    """predict_and_build's command template must NOT reference
    {ncs_spec}.  If it did, the command builder would try to
    substitute an NCS file path into the command string."""
    data, err = _load_programs_yaml()
    if data is None:
        print("  SKIP: %s" % err)
        return

    pb = data.get("phenix.predict_and_build", {})
    cmd = pb.get("command", "")

    assert "{ncs_spec}" not in cmd, (
        "predict_and_build's command template MUST NOT include "
        "{ncs_spec} placeholder (H14.2 fix). Got: %r" % cmd)
    print("  PASS: test_predict_and_build_command_has_no_ncs_spec_placeholder")


def test_predict_and_build_command_template_intact():
    """The command template should still include the other expected
    placeholders — this is a regression guard that the H14.2 surgery
    didn't accidentally drop more than intended."""
    data, err = _load_programs_yaml()
    if data is None:
        print("  SKIP: %s" % err)
        return

    pb = data.get("phenix.predict_and_build", {})
    cmd = pb.get("command", "")

    # Must start with the program name
    assert cmd.startswith("phenix.predict_and_build"), (
        "Command template should start with 'phenix.predict_and_build'. "
        "Got: %r" % cmd)
    # Sequence is the only REQUIRED input — must be in the template
    assert "{sequence}" in cmd, (
        "Command template must include {sequence}. Got: %r" % cmd)
    # data_mtz, full_map, half_map are still expected optional inputs
    for placeholder in ("{data_mtz}", "{full_map}", "{half_map}"):
        assert placeholder in cmd, (
            "Command template must still include %s. Got: %r"
            % (placeholder, cmd))
    print("  PASS: test_predict_and_build_command_template_intact")


# =====================================================================
# §B: map_symmetry hint no longer advertises predict_and_build
# =====================================================================

def test_map_symmetry_hint_does_not_advertise_predict_and_build():
    """phenix.map_symmetry outputs .ncs_spec files.  Its hint that
    documents downstream consumers of those files must NOT include
    predict_and_build (since predict_and_build doesn't accept them)."""
    data, err = _load_programs_yaml()
    if data is None:
        print("  SKIP: %s" % err)
        return

    ms = data.get("phenix.map_symmetry", {})
    hints = ms.get("hints", [])

    # Find the hint about .ncs_spec downstream consumers
    consumer_hints = [
        h for h in hints
        if isinstance(h, str) and ".ncs_spec" in h
    ]
    assert consumer_hints, (
        "map_symmetry should have a hint mentioning .ncs_spec output. "
        "Got hints: %r" % hints)

    # The consumer hint must not list predict_and_build as a consumer
    # (use a word-boundary check to avoid false positives if the hint
    # explicitly mentions predict_and_build's NON-support)
    for h in consumer_hints:
        # Allow mentions of predict_and_build only if they're in the
        # "NOT predict_and_build" form
        lower = h.lower()
        # Find any "predict_and_build" mention
        if "predict_and_build" in lower:
            # It must be in the form that says it's NOT a consumer
            assert ("not predict_and_build" in lower
                    or "not.*predict_and_build" in lower
                    or "except predict_and_build" in lower
                    or "predict_and_build" + " — see" in h
                    or "NOT predict_and_build" in h), (
                "If map_symmetry hint mentions predict_and_build, it "
                "must explicitly say it's NOT a consumer (H14.2). "
                "Got hint: %r" % h)
    print("  PASS: test_map_symmetry_hint_does_not_advertise_predict_and_build")


# =====================================================================
# §C: resolve_cryo_em's ncs_spec entry is intact (sibling regression
# guard — the H14.2 surgery must not touch resolve_cryo_em)
# =====================================================================

def test_resolve_cryo_em_still_accepts_ncs_spec():
    """phenix.resolve_cryo_em DOES accept ncs_file=.  H14.2 only
    touches predict_and_build — resolve_cryo_em's ncs_spec input must
    remain intact."""
    data, err = _load_programs_yaml()
    if data is None:
        print("  SKIP: %s" % err)
        return

    rce = data.get("phenix.resolve_cryo_em", {})
    optional_inputs = rce.get("inputs", {}).get("optional", {})

    assert "ncs_spec" in optional_inputs, (
        "resolve_cryo_em must still accept ncs_spec input (only "
        "predict_and_build is being corrected in H14.2). Got: %r"
        % list(optional_inputs.keys()))

    ncs_entry = optional_inputs["ncs_spec"]
    assert ncs_entry.get("flag") == "ncs_file=", (
        "resolve_cryo_em's ncs_spec.flag should remain 'ncs_file=' "
        "(no map_model. prefix). Got: %r" % ncs_entry.get("flag"))

    # And the command template should still include {ncs_spec}
    cmd = rce.get("command", "")
    assert "{ncs_spec}" in cmd, (
        "resolve_cryo_em's command template must still include "
        "{ncs_spec}. Got: %r" % cmd)
    print("  PASS: test_resolve_cryo_em_still_accepts_ncs_spec")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: predict_and_build cleanup (3)
    test_predict_and_build_has_no_ncs_spec_input()
    test_predict_and_build_command_has_no_ncs_spec_placeholder()
    test_predict_and_build_command_template_intact()
    # §B: map_symmetry hint (1)
    test_map_symmetry_hint_does_not_advertise_predict_and_build()
    # §C: resolve_cryo_em sibling guard (1)
    test_resolve_cryo_em_still_accepts_ncs_spec()


if __name__ == "__main__":
    print("K_H14_2: predict_and_build no-NCS-input config fix (v119.H14.2)")
    print("=" * 70)
    run_all_tests()
    print("=" * 70)
    print("K_H14_2 complete.")
