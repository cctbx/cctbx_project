#!/usr/bin/env python
"""Tests for the P0 harness of the Phenix log extractor.

Run with: python tests/tst_log_structure_extractor.py

Corpus tests gate on PHENIX_LOG_CORPUS (corpus2/work, agent-stripped) and
PHENIX_LOG_CORPUS_AGENT (an agent-shape corpus, for region detection). They
SKIP when unset and FAIL LOUDLY when set but wrong -- a skip that reads as a
pass has already cost this project twice: tst_conformance reads
PHENIX_LOG_CORPUS while the prototype's own suite reads PHIL_CORPUS_DIR, so
running with the documented name silently skipped 4 of its 19 tests and still
printed a clean result.

Fixture-only tests encode the same assumptions as the code; that is how four
defects survived 50 passing tests in the consuming project. So every feature
here carries at least one CORPUS-LEVEL invariant and at least one NEGATIVE
CONTROL.
"""

from __future__ import absolute_import, division, print_function

import inspect
import os
import re
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

try:
    from tests.tst_utils import run_tests_with_fail_fast
except ImportError:
    # Same fallback the house tests carry, so this module still runs in a
    # partial checkout. Note tst_error_analyzer.py has no such fallback and
    # dies on import when tst_utils is absent.
    def run_tests_with_fail_fast():
        g = sys._getframe(1).f_globals
        tests = sorted(k for k in g if k.startswith("test_"))
        for name in tests:
            print("  Running %s..." % name)
            g[name]()
            print("  PASS: %s" % name)
        print("\nAll %d tests passed." % len(tests))

from log_extraction.log_structure_extractor import (  # noqa: E402
  AGENT_METRICS_MARKER, Candidate, CompletionRecord, ControlSkip, Cycle, Decision, Exclusion,
  Identification, Item, LabeledValue, LogStructure, Measurement, Phase,
  Section, Stage, Unparsed, SCREEN_RULES, SOURCE_AGENT, SOURCE_PROGRAM,
  SOURCE_UNKNOWN, UNCLAIMED, RULE_EXCLUDED, UNPARSED_STATUSES,
  GENERIC_ONLY, LABEL_DISTINCT_LIMIT, REFUSAL_BANNER, REFUSAL_LABEL_RUNAWAY,
  REFUSAL_TABLE_HEADER, REFUSAL_TOO_LONG, SECTION_MAX_TITLE,
  CYCLE_SENTINEL, REFUSAL_CYCLE_NO_METRICS, REFUSAL_SKIP_UNRECOGNISED,
  BANNER_BLACKLIST, COMPLETION_CLUSTER_GAP, SIGNAL_BANNER,
  SIGNAL_PHENIX_HEADER, TerminalFailure, find_preamble,
  find_terminal_failures,
  identify_program, REFUSAL_SETTING_NARRATION, REFUSAL_SETTING_NO_REASON,
  REFUSAL_SETTING_UP, find_completion_records,
  REFUSAL_STAGE_SUMMARY, find_decisions,
  find_cycles, find_labeled_values, find_phases_and_skips, find_regions,
  find_stage_tables,
  find_sections, main, scan, screen_line, split_lines, strip_agent_prefix)

# NO EXPORTS REQUIRED. corpus_paths finds the corpus shipped inside the
# package; the environment variables remain as an override for anyone who
# keeps it elsewhere. A suite that needs two exports before it checks
# anything real is a suite that gets run without them.
try:
    from log_extraction import corpus_paths
    CORPUS = corpus_paths.working_corpus()
    CORPUS_AGENT = corpus_paths.agent_corpus()
except ImportError:                                     # pragma: no cover
    CORPUS = os.environ.get("PHENIX_LOG_CORPUS")
    CORPUS_AGENT = os.environ.get("PHENIX_LOG_CORPUS_AGENT")
# Every Nth log. Default 1 = the whole corpus, which is what any reported
# number must be measured against. The NEGATIVE-CONTROL runner sets a stride,
# because it re-runs the whole suite once per control and 21 x 16 s is a
# runner nobody runs. A control only has to DISCRIMINATE, not to measure.
CORPUS_STRIDE = int(os.environ.get("PHENIX_LOG_CORPUS_STRIDE", "1"))


def _corpus_logs(root, stride=None):
  """Every .log under root, recursively. Never gates on the extension for
  behaviour -- only for collection -- and never keys anything on the name."""
  out = []
  for base, _, files in os.walk(root):
    for name in sorted(files):
      if name.endswith(".log") and not name.startswith("._"):
        out.append(os.path.join(base, name))
  if not out:
    raise AssertionError(
      "%s is set but contains no logs -- refusing to pass silently" % root)
  if stride is None:
    # THE STRIDE APPLIES ONLY TO THE LARGE WORKING CORPUS. FOUND IN P8
    # REVIEW: it was being applied to the 20-log wave-1 corpora too, so at
    # stride 40 those became one log and four tests failed for want of data
    # rather than for want of correctness. Sampling a 20-log reference set
    # is never what anyone meant.
    stride = CORPUS_STRIDE if os.path.abspath(root) == os.path.abspath(
      CORPUS or "") else 1
  return out[::stride]


def _read(path):
  handle = open(path, "rb")
  try:
    return handle.read().decode("utf-8", "replace")
  finally:
    handle.close()


_SCAN_CACHE = {}


def _scanned(path):
  """A cached scan, for the READ-ONLY corpus tests.

  Added when the suite reached 27 s and the negative-control runner -- which
  re-runs the whole suite once per control -- would have taken nine minutes.
  Tests that are ABOUT scanning (purity, the hint-independence invariant, the
  optimisation-equivalence check) call scan() directly and must not use this,
  or they would be asserting against a cached answer instead of a fresh one."""
  if path not in _SCAN_CACHE:
    _SCAN_CACHE[path] = scan(_read(path))
  return _SCAN_CACHE[path]


# ------------------------------------------------------------- contract shape


def test_contract_deviations_are_declared_not_silent():
    """Two v1 contract items differ in the implementation, and both are now
    DECLARED in requirements v2 rather than left as silent drift:

      Measurement.unit is permanently None in v1 (populating it is
      interpretation), and is RETAINED because the consumer built against
      the stub;

      Unparsed's v1 `kind` field is `screen_rule` here, because `kind` names
      the CHANNEL on every other item and is what items() sorts by.

    This test exists so that changing either without updating the spec
    fails."""
    assert Measurement("R free", 0.23, 1).unit is None
    item = Unparsed("x", 1, screen_rule="verb")
    assert item.screen_rule == "verb"
    assert item.kind == "unparsed"          # the CHANNEL, not the rule


def test_every_item_kind_accepts_a_line_and_defaults():
    items = [
        Section("Collecting inputs", 3, end=5),
        Phase("Running dock_and_rebuild", 7),
        Stage("1_bss", 9, metrics={"r_free": 0.488}, description="bulk"),
        Cycle(2, 11, metrics={"Built": 0}, sentinel=True),
        Decision("rebuild_in_place", "False", "maps_only=True", 13),
        ControlSkip("model-building", "'maps_only' is set", 15),
        Exclusion("9GSD.A", "not x-ray and not computational", 17),
        Measurement("Resolution", "2.50", 19, unit="A"),
        LabeledValue("Space group", "P 21 21 21", 21),
        CompletionRecord("Job complete", 23),
        Unparsed("weird", 25, screen_rule="verb"),
    ]
    for item in items:
        assert item.line >= 1, item
        assert item.source == SOURCE_PROGRAM, item
        assert item.section_id is None, item
        assert item.end >= item.line, item
        assert repr(item), item


def test_line_numbers_are_one_based_and_validated():
    # bool is an int subclass -- True would silently become line 1
    for bad in (0, -1, None, "3", 1.0, True):
        try:
            Phase("x", bad)
        except (ValueError, TypeError):
            continue
        raise AssertionError("accepted a bad line number: %r" % (bad,))
    assert Phase("x", 1).line == 1


def test_source_unknown_is_a_real_value():
    """Defaulting an ambiguous region to `program` would convert uncertainty
    into a positive misattribution -- the same error as reading an .eff's
    presence as evidence for every value in it."""
    assert Phase("x", 1, source=SOURCE_UNKNOWN).source == SOURCE_UNKNOWN
    assert Phase("x", 1, source=SOURCE_AGENT).source == SOURCE_AGENT
    try:
        Phase("x", 1, source="probably_program")
    except ValueError:
        return
    raise AssertionError("accepted an undefined source")


def test_section_id_none_is_normal_for_every_kind():
    """If items required an enclosing section, headerless programs would
    become second-class inputs and implicit sections would get invented."""
    for cls, args in ((Phase, ("x", 1)), (Decision, ("a", "b", "c", 1)),
                      (LabeledValue, ("k", "v", 1)),
                      (CompletionRecord, ("Job complete", 1))):
        item = cls(*args)
        assert item.section_id is None, item


def test_identification_abstains_by_default():
    ident = Identification()
    assert ident.is_unknown is True
    assert ident.candidates == []
    assert ident.name is None
    ident = Identification([Candidate("phenix.xtriage", 0.95,
                                      Candidate.SELF_IDENTIFICATION, 1)],
                           signals_fired=["banner"])
    assert ident.is_unknown is False
    assert ident.name == "phenix.xtriage"
    assert ident.candidates[0].evidence_kind == Candidate.SELF_IDENTIFICATION


def test_empty_and_degenerate_input():
    for text in ("", "   ", "\n\n\n", None):
        structure = scan(text)
        assert structure.is_flat is True, text
        assert structure.items() == []


def test_truncated_mid_line_and_no_trailing_newline():
    structure = scan("Processing files:\n" + "=" * 40 + "\n  Found phil, /x")
    assert structure is not None


def test_scan_is_pure():
    """No global state, no mutation of the input, same answer twice."""
    text = "Setting a=b as c\n" + "=" * 30 + "\nResolution: 2.5\n"
    original = text[:]
    first, second = scan(text), scan(text)
    assert text == original
    assert [repr(i) for i in first.items()] == [repr(i) for i in second.items()]


# ----------------------------------------------------------- the frozen screen


def test_screen_recognises_its_four_shapes():
    """REVISION (P0 review): one assertion here was a disjunction where the
    right-hand side was separately asserted on the next line, so it could not
    fail. Every assertion is now a single exact expectation."""
    assert screen_line("=" * 40) == "rule_line"
    assert screen_line("Skipping 9GSD.A - not x-ray") == "verb"
    assert screen_line("Setting n_cycle_build=1 as nbatch >1") == "verb"
    assert screen_line("            Running build_chains on segment 1") == "verb"
    assert screen_line("Resolution: 2.50") == "key_value"
    assert screen_line("  0.4064 0.4880 0.010 1.282 26.9") == "numeric_row"
    # a refine stage row is indented 6, and key_value requires indent <= 3,
    # so it lands in numeric_row. Recorded because it decides which rule P4's
    # stage parser has to reconcile its claims against.
    assert screen_line("      1_bss: 0.4064 0.4880 0.010 1.282") == "numeric_row"


def test_numeric_row_keeps_label_bearing_table_rows():
    """MEASURED IN REVIEW: a numeric-MAJORITY reading of the rule drops 11,363
    candidates including `| 4_0 (c) | 0 | 0.00 |`, a real table row with label
    cells. Over-selecting prose that carries four numbers is the cheaper
    error for a screen whose job is to over-select."""
    real = "  | 1      | 0.635           | 3.37            | 0.879              |"
    assert screen_line(real) == "numeric_row", screen_line(real)
    assert screen_line("  1  2  3") is None            # fewer than four
    # REVISION (twice): the first fixture was invented and carried only three
    # numeric cells, so it asserted the opposite of what it claimed; the
    # second was a real line that also had three. Both would have passed a
    # reading of the test and failed the data. This one is from
    # corpus2/work and is dropped by a majority rule (6 numerics, 13 tokens),
    # which is the case the test exists to pin down.


def test_strip_agent_prefix():
    """The wave-1 GUI-shape corpus still carries `LOG TEXT:` on line 1 of 15
    of its 20 logs, so tolerating the prefix is not optional."""
    assert strip_agent_prefix("LOG TEXT: Starting AutoBuild") == \
        "Starting AutoBuild"
    assert strip_agent_prefix("Starting AutoBuild") == "Starting AutoBuild"
    assert strip_agent_prefix("") == ""


def test_forms_and_is_flat_can_disagree_deliberately():
    """`forms` lists every non-empty channel INCLUDING unparsed; `is_flat`
    ignores unparsed. They therefore disagree on a log that yielded only
    unreadable structure, and that is intended, not a bug."""
    structure = scan("=" * 40 + "\n")
    assert structure.forms == ["unparsed"], structure.forms
    assert structure.is_flat is True


def test_cli_summary_mode_still_prints_one_line():
    """`--summary` is the machine-friendly form and must stay stable; the
    readable report is what someone running it on a log actually wants, so
    that is the default."""
    import tempfile
    handle, path = tempfile.mkstemp(suffix=".log")
    os.write(handle, b"Starting phenix.xtriage\nResolution: 2.5\n")
    os.close(handle)
    saved = sys.stdout
    try:
        import io
        sys.stdout = io.StringIO()
        code = main(["--summary", path])
        text = sys.stdout.getvalue()
    finally:
        sys.stdout = saved
        os.unlink(path)
    assert code == 0, code
    assert len(text.strip().splitlines()) == 1, text
    assert "forms=" in text and "unparsed=" in text, text


def test_cli_report_cites_line_numbers():
    """The report exists so a claim can be checked, so every finding it
    prints carries the log line it came from."""
    import tempfile
    handle, path = tempfile.mkstemp(suffix=".log")
    os.write(handle, b"Starting phenix.autobuild\n"
                     b"Skipping 9GSD.A - not x-ray\n"
                     b"Setting rebuild_in_place=False as maps_only=True\n")
    os.close(handle)
    saved = sys.stdout
    try:
        import io
        sys.stdout = io.StringIO()
        main([path])
        text = sys.stdout.getvalue()
    finally:
        sys.stdout = saved
        os.unlink(path)
    assert "phenix.autobuild" in text, text
    assert "EXCLUDED CANDIDATES" in text and "not x-ray" in text, text
    assert "rebuild_in_place = False   because maps_only=True" in text, text
    assert "NOT UNDERSTOOD" in text, text


def test_cli_reads_a_file_without_gating_on_the_extension():
    """A production defect traced to `if ext != ".log": return None`, after
    which the program name went unset and the pipeline printed a confident
    report about a program that never ran."""
    import tempfile
    handle, path = tempfile.mkstemp(suffix=".txt")
    os.write(handle, b"Starting phenix.xtriage\nResolution: 2.5\n")
    os.close(handle)
    # CAPTURE the CLI's own output. It used to print a genuine-looking
    # `/no/such/file: No such file or directory` into the middle of a passing
    # test log -- deliberate behaviour that reads exactly like a failure, and
    # a log that cries wolf trains a reader to skim past real errors.
    # Capturing also lets us assert on WHAT it printed, not just its code.
    saved_out, saved_err = sys.stdout, sys.stderr
    try:
        import io
        sys.stdout, sys.stderr = io.StringIO(), io.StringIO()
        ok_code = main([path])          # the readable report is the default
        reported = sys.stdout.getvalue()
        missing_code = main(["/no/such/file"])
        complaint = sys.stderr.getvalue()
    finally:
        sys.stdout, sys.stderr = saved_out, saved_err
        os.unlink(path)
    assert ok_code == 0, ok_code
    assert ".txt" in reported, reported
    assert "program:" in reported, reported
    assert missing_code == 1, missing_code
    assert "/no/such/file" in complaint, complaint


def test_screen_rejects_prose_and_bare_headings():
    assert screen_line("") is None
    assert screen_line("the refinement then proceeded normally") is None
    assert screen_line("Processing files:") is None      # no value after colon
    assert screen_line("== short ==") is None            # rule too short


def test_screen_does_not_swallow_timestamps_or_file_line_refs():
    """Raised in review. Neither reaches the screen, by construction: a bare
    timestamp starts with a digit and the screen needs a leading letter;
    file.py:124 has no colon-SPACE. Asserted rather than assumed."""
    assert screen_line("14:23:01") is None
    assert screen_line("  09:15") is None
    assert screen_line("/path/to/module.py:124") is None
    assert screen_line("module.py:124: warning") is None


def test_screen_keeps_log_level_lines():
    """The reviewer proposed filtering Note:/Warning:/Error: as console
    furniture. MEASURED: they are content, and only 3 of 30 failed runs carry
    ANY error keyword in what the program itself wrote, so this filter would
    delete the scarcest signal in the corpus."""
    assert screen_line("Sorry: number of groups of duplicate atom labels:  5")
    assert screen_line("WARNING: please remember the possibility of twin laws")
    assert screen_line("Warning: very small nonbonded interaction distances.")


# ---------------------------------------------------------------- the regions


def test_agent_header_region_is_found_and_marked():
    lines = ["WORKING DIRECTORY: /tmp/x",
             "COMMAND THAT WAS RUN: phenix.xtriage data.sca",
             "LOG TEXT: Starting phenix.xtriage",
             "Processing files:"]
    regions = find_regions(lines)
    assert [r.kind for r in regions] == ["agent_header", "agent_prefix"], regions
    header, prefix = regions
    assert (header.start, header.end) == (1, 2), header
    assert header.source == SOURCE_AGENT
    # REVISION (P0 review): the header first ran to line 3, swallowing the
    # program's own banner. `LOG TEXT: Starting phenix.xtriage` is nine
    # characters of agent in front of program content, and in 9 of the 20
    # agent-shape logs that content is the best identification signal there
    # is. The line is program-sourced; only the prefix is not.
    assert (prefix.start, prefix.end) == (3, 3), prefix
    assert prefix.source == SOURCE_PROGRAM, prefix
    assert not header.contains(4) and not prefix.contains(4)


def test_agent_footer_region_includes_its_rule():
    lines = ["Job complete", "*" * 50, "FINAL QUALITY METRICS REPORT:",
             "Resolution: 2.50", "*" * 50]
    regions = find_regions(lines)
    assert len(regions) == 1, regions
    assert (regions[0].start, regions[0].end) == (2, 5), regions
    assert regions[0].source == SOURCE_AGENT


def test_a_midfile_wrapper_block_is_unknown_not_agent():
    """ai_agent_*.log files QUOTE their children's wrapper blocks mid-file. We
    cannot tell the agent wrapping this log from this log quoting an agent, so
    we say unknown rather than assert agent."""
    lines = (["Running cycle 1"] + ["MTZ LABEL INFO:", "  FOBS 100%"]
             + ["still going"] * 40)
    regions = find_regions(lines)
    assert len(regions) == 1, regions
    assert regions[0].source == SOURCE_UNKNOWN, regions
    assert regions[0].start == 2, regions
    # and the tail-position case still reads as the agent's own footer
    tail = find_regions(["line"] * 40 + ["FINAL QUALITY METRICS REPORT:",
                                         "Resolution: 2.5"])
    assert tail[0].source == SOURCE_AGENT, tail


def test_a_log_with_no_wrapper_has_no_regions():
    assert find_regions(["Starting phenix.xtriage", "Processing files:"]) == []


def test_no_header_means_no_positional_arguments():
    """parse_command_line_files once fell through to parsing line 1 when no
    header existed, returning `libtbx.start_process` as a file the user
    supplied. No header, no header-derived anything."""
    structure = scan("Starting libtbx.start_process\nProcessing files:\n")
    assert structure.regions == []
    assert all(i.source == SOURCE_PROGRAM for i in structure.items())


# --------------------------------------------------- identification-independence


def test_identification_cannot_gate_extraction():
    """THE layering invariant. scan() never requires program identification;
    failure to identify cannot suppress otherwise extractable structure."""
    text = ("Starting phenix.refine\n" + "=" * 40 + "\nResolution: 2.5\n"
            "Skipping model-building as 'maps_only' is set\n")
    without = scan(text)
    for hint in ("phenix.refine", "phenix.autobuild", "not.a.program", ""):
        with_hint = scan(text, program_name=hint)
        assert [repr(i) for i in without.items()] == \
               [repr(i) for i in with_hint.items()], hint


# ------------------------------------------------------------ derived views


def test_is_flat_ignores_unparsed():
    """A log whose only output is `I saw things I could not read` is flat.
    Counting unparsed as structure would make every log non-flat and the
    reach metrics meaningless."""
    structure = scan("=" * 40 + "\n")
    assert structure.unparsed, "the screen should have caught the rule"
    assert structure.is_flat is True


def test_unparsed_counts_report_all_three_statuses():
    counts = scan("=" * 40 + "\n").unparsed_counts()
    assert set(counts) == set(UNPARSED_STATUSES), counts
    assert counts[UNCLAIMED] == 1, counts


def test_reach_reports_three_metrics():
    reach = scan("").reach()
    assert set(reach) == set(["structural", "basic", "diagnostic"]), reach
    assert reach == dict(structural=False, basic=False, diagnostic=False)


def test_items_are_returned_in_line_order():
    """When the program is unknown and section_id is None -- the normal state
    for a headerless log -- spatial proximity is the consumer's only remaining
    context, and it needs document order across kinds."""
    structure = LogStructure(n_lines=10)
    structure.phases.append(Phase("late", 9))
    structure.labeled_values.append(LabeledValue("k", "v", 2))
    structure.sections.append(Section("mid", 5))
    assert [i.line for i in structure.items()] == [2, 5, 9]


def test_metric_moves_and_exclusion_groups():
    structure = LogStructure()
    for name, value in (("0", 0.5371), ("1_bss", 0.4880), ("1_xyz", 0.4915)):
        structure.stages.append(
            Stage(name, len(structure.stages) + 1, metrics={"r_free": value}))
    moves = structure.metric_moves("r_free", 0.001)
    assert [m["stage"] for m in moves] == ["1_bss", "1_xyz"], moves
    assert abs(moves[1]["delta"] - 0.0035) < 1e-9, moves[1]

    for i, (item, reason) in enumerate([("A", "not x-ray"), ("B", "not x-ray"),
                                        (None, "too few residues")]):
        structure.exclusions.append(Exclusion(item, reason, i + 1))
    groups = structure.exclusion_groups()
    assert len(groups) == 2, groups
    assert groups[0]["count"] == 2, groups
    assert groups[1]["items"] == [None], groups


# --------------------------------------------------------------- P1 sections


def test_P1_form_a_title_inside_a_rule():
    for raw, title in (("=" * 12 + " Collecting inputs " + "=" * 12,
                        "Collecting inputs"),
                       ("-" * 10 + "Processing PDB file(s)" + "-" * 10,
                        "Processing PDB file(s)")):
        sections, claimed, refusals = find_sections([raw])
        assert [x[0] for x in sections] == [title], (raw, sections)
        assert claimed == {1: "sections"}, claimed
        assert refusals == [], refusals


def test_P1_form_b_title_over_a_rule():
    """The Program Template's own shape, which the prototype missed entirely
    -- most of why it found sections in only 15 of the 20 wave-1 logs."""
    sections, claimed, _ = find_sections(["Processing files:", "-" * 40,
                                          "  Found phil, /x.eff"])
    assert [x[0] for x in sections] == ["Processing files:"], sections
    assert sections[0][1] == 1 and sections[0][2] == 3, sections
    assert sorted(claimed) == [1, 2], claimed   # title AND its rule


def test_P1_a_bare_rule_is_not_a_section():
    """REVISION (P1): the first fixture here was `["text", "="*70]` -- taken
    from the wave-1 conformance suite -- and it FAILED, correctly. A text
    line above an `=` rule is form B, and in corpus2/work that shape is a
    real heading 1,818 times (`Starting job`, `OUTPUT FILES`). The fixture,
    not the parser, was wrong: it meant to test a rule with nothing above it.
    See HANDOFF 6 for the conformance-suite consequence."""
    sections, _, _ = find_sections(["", "=" * 70, ""])
    assert sections == [], sections
    assert scan("\n" + "=" * 70 + "\nmore\n").sections == []
    # a pure rule must not be read as an inline title made of dashes
    assert scan("-" * 40 + "\n").sections == []


def test_P1_a_blank_line_above_a_rule_is_not_a_refused_title():
    """DEFECT FOUND BY THE PAIRED-SHAPE INVARIANT, not by a fixture: a blank
    line above a rule was emitted as a rule_excluded record. A blank line is
    not a title we refused; it is not a title. The bad record showed up as a
    single extra item in one of the 20 wave-1 pairs."""
    structure = scan("\n" + "*" * 40 + "\n")
    blank_records = [u for u in structure.unparsed if not u.text.strip()]
    assert blank_records == [], [repr(u) for u in blank_records]
    assert structure.sections == []


def test_P1_a_banner_rule_does_not_underline_a_title():
    """`*` rules wrap a block above and below; they are delimiters, not
    underlines. MEASURED: accepting them would add 824 titles including
    `hydrogens.refine=riding` and `Unset PHENIX_OVERWRITE_ALL ...`."""
    structure = scan("hydrogens.refine=riding\n" + "*" * 50 + "\n")
    assert structure.sections == [], structure.sections
    refused = [u for u in structure.unparsed if u.excluded_by]
    assert REFUSAL_BANNER in [u.excluded_by for u in refused], refused
    # but a - or = underline does make a section
    assert scan("Outlier Rejection\n" + "-" * 50 + "\n").sections
    assert scan("Starting job\n" + "=" * 50 + "\n").sections


def test_P1_over_long_and_table_titles_are_refused_by_name():
    """Requirements 6 names an over-long section title as something to
    REPORT rather than accept. A refusal that does not say which rule
    refused it is a number that moved for unknown reasons."""
    long_title = "R" * (SECTION_MAX_TITLE + 1)
    structure = scan(long_title + "\n" + "-" * 40 + "\n")
    assert structure.sections == [], structure.sections
    refused = [u for u in structure.unparsed if u.status == RULE_EXCLUDED]
    assert [u.excluded_by for u in refused] == [REFUSAL_TOO_LONG], refused

    table = "| d_spacing | Expected rel. I | Data Z-score | Completeness |"
    structure = scan(table + "\n" + "-" * 40 + "\n")
    assert structure.sections == [], structure.sections
    refused = [u for u in structure.unparsed if u.status == RULE_EXCLUDED]
    assert [u.excluded_by for u in refused] == [REFUSAL_TABLE_HEADER], refused


def test_P1_sections_extend_to_the_next_one():
    text = ("First:\n" + "-" * 40 + "\nbody\nbody\n"
            "Second:\n" + "-" * 40 + "\ntail\n")
    sections = scan(text).sections
    assert [(s.title, s.line, s.end) for s in sections] == [
        ("First:", 1, 4), ("Second:", 5, 7)], sections
    assert [s.start for s in sections] == [1, 5]


def test_P1_items_inherit_the_section_they_sit_in():
    """section_id must still be able to be None -- that is the point of
    test_section_id_none_is_normal_for_every_kind -- but where a section
    does enclose an item, the item should say so."""
    text = ("Processing files:\n" + "-" * 40 + "\nResolution: 2.5\n")
    structure = scan(text)
    inside = [u for u in structure.unparsed if u.line == 3]
    assert inside and inside[0].section_id == 0, structure.unparsed


def test_P1_claimed_lines_leave_the_unclaimed_pile():
    """The whole point of the ledger: a line a parser understood must stop
    being reported as structure we could not read."""
    text = "Some Heading\n" + "=" * 40 + "\n"
    structure = scan(text)
    assert structure.sections, structure.sections
    counts = structure.unparsed_counts()
    assert counts[UNCLAIMED] == 0, [str(u) for u in structure.unparsed]


def test_P1_claiming_outside_the_frozen_screen_is_counted():
    """A parser claiming a line the frozen screen never proposed is how we
    learn the screen is too narrow. Freezing the screen without this number
    would be freezing it blind."""
    # `Processing files:` ends in a colon, so key_value does not take it and
    # no other screen rule matches -- yet the section parser claims it.
    assert screen_line("Processing files:") is None
    structure = scan("Processing files:\n" + "-" * 40 + "\n")
    assert structure.claimed_outside_screen == [1], \
        structure.claimed_outside_screen


def test_corpus_P1_sections_and_the_ledger_balance():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    logs = _corpus_logs(CORPUS)
    with_sections = total = outside = 0
    refusals = {}
    for path in logs:
        text = _read(path)
        lines = split_lines(text)
        structure = scan(text)
        if structure.sections:
            with_sections += 1
        total += len(structure.sections)
        outside += len(structure.claimed_outside_screen)
        previous_end = 0
        for section in structure.sections:
            assert 1 <= section.line <= section.end <= len(lines), \
                (path, section)
            assert section.line > previous_end, (path, section)
            previous_end = section.line
            assert len(section.title) <= SECTION_MAX_TITLE, (path, section)
        for item in structure.unparsed:
            if item.status == RULE_EXCLUDED:
                assert item.excluded_by, (path, item)
                refusals[item.excluded_by] = refusals.get(
                    item.excluded_by, 0) + 1
    assert with_sections > len(logs) // 2, with_sections
    print("      (%d/%d logs, %d sections, %d claims outside the screen,"
          " refusals %s)" % (with_sections, len(logs), total, outside,
                             refusals))


# ----------------------------------- P10 terminal failures and the preamble


def test_P11_compound_pairs_split_only_when_the_line_is_all_pairs():
    """`R:   0.42  Rfree:   0.48` is two numbers, and capturing it as label
    "R" value "0.42  Rfree: 0.48" buries the second. 9% of key:value lines
    across 554 logs carry two or more pairs.

    But a naive split turns `Space group: C 1 2 1 (No. 5)` into `C`, so the
    line is split ONLY when a sequence of `label: single-token` pairs
    accounts for it entirely -- a shape, not a judgement."""
    got = [(x.label, x.value) for x in scan("R:   0.42  Rfree:   0.48\n"
                                            ).labeled_values]
    assert got == [("R", "0.42"), ("Rfree", "0.48")], got
    for text, expected in [
            ("Space group: C 1 2 1 (No. 5)", "C 1 2 1 (No. 5)"),
            ("Map-model CC:  0.47 (overall) and 0.59 (local)",
             "0.47 (overall) and 0.59 (local)"),
            ("Working directory: /Users/t/Documents/x", "/Users/t/Documents/x")]:
        values = scan(text + "\n").labeled_values
        assert len(values) == 1, (text, values)
        assert values[0].value == expected, (text, values[0].value)


def test_P11_the_report_shows_the_numbers_not_a_count():
    """A REGRESSION I INTRODUCED: collapsing labeled values to one BACKGROUND
    line hid `Map-model CC: 0.47` and `R: 0.42` -- captured all along, and
    the two things Tom actually wanted from an AutoSol run."""
    import io
    from log_extraction.log_structure_extractor import report
    out = io.StringIO()
    report(scan("Map-model CC:  0.47 (overall)\nR:   0.42  Rfree:   0.48\n"),
           "x.log", out=out)
    text = out.getvalue()
    assert "VALUES THE RUN REPORTED" in text, text
    assert "Map-model CC" in text, text
    assert "Rfree" in text, text


def test_P11_the_values_are_ordered_most_recent_first():
    """Chronology, not a judgement about importance -- but ascending order
    buried the end-of-run results under 30 lines of setup."""
    import io
    from log_extraction.log_structure_extractor import report
    out = io.StringIO()
    report(scan("Early: 1\n" + "filler\n" * 5 + "Late: 2\n"), "x.log", out=out)
    text = out.getvalue()
    assert text.index("Late") < text.index("Early"), text


def test_docs_do_not_contradict_the_code_on_channel_count():
    """DOCUMENT DRIFT IS A DEFECT. Four documents quoted 11 channels, 30%
    diagnostic reach and superseded test counts after the changes that moved
    them -- drift caused by my own edits, and it took a reader to notice.

    A test cannot check every number, but it can check the ones that are
    structural facts about the code, so those at least cannot rot. The rest
    is handled by the record/current-state split declared in
    MEASUREMENTS.md: records are never edited, current-state documents are
    reconciled, and MEASUREMENTS.md is regenerated."""
    import log_extraction.log_structure_extractor as module
    root = os.path.dirname(os.path.dirname(os.path.abspath(module.__file__)))
    docs = os.path.join(root, "log_extraction", "docs")
    if not os.path.isdir(docs):
        SKIPPED.append("docs directory absent")
        print("      SKIP (docs directory absent)")
        return
    count = len(module.ITEM_KINDS)
    wrong = []
    for name in ("EXTRACTOR_ARCHITECTURE.md", "MEASUREMENTS.md",
                 "HANDOFF.md"):
        path = os.path.join(docs, name)
        if not os.path.exists(path):
            continue
        handle = open(path)
        try:
            text = handle.read()
        finally:
            handle.close()
        # only the CURRENT-STATE portion; records below the marker are dated
        for marker in ("# Parts VI-b onwards are RECORDS",
                       "## 7. The full defect record"):
            text = text.split(marker)[0]
        # BOTH PHRASINGS. The first version looked only for "N channels"
        # and the handoff writes "Channels: N" -- so the guard passed while
        # checking nothing in that file. Verified by deliberately writing a
        # wrong count into each phrasing and confirming this fails.
        for other in range(8, 20):
            if other == count:
                continue
            for phrasing in ("%d channels" % other, "Channels: %d" % other,
                             "%d output channels" % other):
                if phrasing.lower() in text.lower():
                    wrong.append((name, phrasing))
    assert not wrong, ("documents contradict ITEM_KINDS=%d: %s"
                       % (count, wrong))
    print("      (%d channels, agreed by code and documents)" % count)


def test_P12_options_are_consumed_by_position_not_by_value():
    """FOUND IN REVIEW: the first version dropped any path equal to the
    pattern, so `--grep Sorry Sorry.log` searched nothing and said so about
    no file."""
    import io
    import tempfile
    handle, path = tempfile.mkstemp(suffix=".log")
    os.write(handle, b"Skipping x - bad\n")
    os.close(handle)
    saved = sys.stdout
    try:
        sys.stdout = io.StringIO()
        code = main(["--grep", os.path.basename(path), path])
        text = sys.stdout.getvalue()
    finally:
        sys.stdout = saved
        os.unlink(path)
    assert code == 0, code
    assert "LINES MATCHING" in text, text


def test_P12_summary_refuses_flags_it_would_ignore():
    """Silently ignoring a flag is how a user comes to distrust output."""
    import io
    import tempfile
    handle, path = tempfile.mkstemp(suffix=".log")
    os.write(handle, b"x\n")
    os.close(handle)
    saved_out, saved_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = io.StringIO(), io.StringIO()
        code = main(["--summary", "--grep", "x", path])
        complaint = sys.stderr.getvalue()
    finally:
        sys.stdout, sys.stderr = saved_out, saved_err
        os.unlink(path)
    assert code == 2, code
    assert "--summary ignores" in complaint, complaint


def test_P12_grep_searches_the_RAW_lines_not_the_ledger():
    """A CORRECTION TO SOMETHING I SAID. I claimed unread prose "sits in the
    unparsed ledger". IT DOES NOT: `unparsed` holds only lines the frozen
    SCREEN proposed -- 737 of one AutoSol log's 3,877 -- and `Model is in
    /path` matches no screen rule, so it is in neither the channels nor the
    ledger. A --grep over `unparsed` would have missed exactly the lines it
    exists to find."""
    import io
    from log_extraction.log_structure_extractor import report
    text = ("Model is in /Users/t/x/y/model.pdb\n"
            "Space group: C 1 2 1\n")
    structure = scan(text)
    assert not [u for u in structure.unparsed if u.line == 1], structure.unparsed
    out = io.StringIO()
    report(structure, "x.log", out=out, pattern="Model is in",
           structure_text=text)
    shown = out.getvalue()
    assert "LINES MATCHING 'Model is in' (1)" in shown, shown
    assert "[not read]" in shown, shown


def test_P12_grep_says_what_read_each_hit():
    """The tag is the point: a reader can tell "the tool read this and put it
    somewhere" from "the tool never saw it"."""
    import io
    from log_extraction.log_structure_extractor import report
    text = "Space group: C 1 2 1\nModel is in /a/b/c/d.pdb\n"
    out = io.StringIO()
    report(scan(text), "x.log", out=out, pattern="e", structure_text=text)
    shown = out.getvalue()
    assert "[labeled_values]" in shown, shown
    assert "[not read]" in shown, shown


def test_P12_grep_without_the_text_says_so_rather_than_finding_nothing():
    """Silently reporting zero hits because the caller forgot an argument is
    the false-negative this project keeps meeting."""
    import io
    from log_extraction.log_structure_extractor import report
    out = io.StringIO()
    report(scan("Model is in /a/b/c/d.pdb\n"), "x.log", out=out,
           pattern="Model")
    assert "--grep needs the log text" in out.getvalue(), out.getvalue()


def test_P11_label_filter_prefers_an_exact_match():
    """`--label R` matched every label containing the letter r -- 196 of
    them, which is not a series, it is the whole log again."""
    import io
    from log_extraction.log_structure_extractor import report
    out = io.StringIO()
    report(scan("R: 0.42\nResolution: 2.5\nR: 0.44\n"), "x.log", out=out,
           only_label="R")
    text = out.getvalue()
    assert "EVERY VALUE FOR 'R' (2)" in text, text
    assert "Resolution" not in text.split("EVERY VALUE")[1], text


def test_P10_an_abort_marker_is_reported():
    """THE motivating case. An AutoSol run ended `STOPWIZARD IN
    AutoSol_run_1_/STOPWIZARD` on line 3877 and the report showed only a
    CHILD completion from line 1994, 1,883 lines earlier."""
    structure = scan("work\nSTOPWIZARD IN AutoSol_run_1_/STOPWIZARD\n")
    got = [(f.failure_kind, f.line) for f in structure.terminal_failures]
    assert got == [(TerminalFailure.ABORT, 2)], got


def test_P10_a_sorry_block_keeps_its_remedy():
    """The remedy is the most useful text in a failed log and it spans blank
    lines, so the block ends at TWO blanks or at a line that resumes the run
    -- not at the first blank."""
    text = ("Traceback (most recent call last):\n"
            "  File \"x.py\", line 3\n"
            "    raise Sorry(msg)\n"
            "Sorry: No array of R-free flags found.\n"
            "run this command again with\n"
            "\n"
            "xray_data.r_free_flags.generate=True\n"
            "\n"
            "Please try again.\n")
    failures = scan(text).terminal_failures
    kinds = [f.failure_kind for f in failures]
    assert kinds == [TerminalFailure.TRACEBACK, TerminalFailure.SORRY], kinds
    sorry = failures[1]
    assert sorry.text.startswith("No array of R-free flags"), sorry.text
    assert "generate=True" in sorry.remedy, sorry.remedy
    assert "Please try again." in sorry.remedy, sorry.remedy


def test_P10_a_capped_block_says_it_was_capped():
    """FOUND IN REVIEW: one Sorry block ran to the 40-line cap and stopped
    there BY LUCK -- the real end was a box-drawing border one line later.
    A cap that silently decides a boundary is the silent-truncation class
    this project keeps meeting, so the boundary rule was widened and
    cap-truncation is now flagged."""
    from log_extraction.log_structure_extractor import SORRY_MAX_LINES
    text = "Sorry: broke\n" + "".join("  detail %d\n" % i
                                      for i in range(SORRY_MAX_LINES + 10))
    failure = scan(text).terminal_failures[0]
    assert failure.truncated is True, (failure.line, failure.end)
    ordinary = scan("Sorry: broke\n  detail\n\n\nnext\n").terminal_failures[0]
    assert ordinary.truncated is False, ordinary


def test_P10_a_box_border_ends_the_block():
    text = "Sorry: broke\n  detail\n\u2514\u2500\u2500\u2500\nmore output\n"
    failure = scan(text).terminal_failures[0]
    assert failure.end == 2, (failure.line, failure.end)


def test_P10_a_terminal_failure_is_not_a_verdict():
    """Non-goal 2.8 stands. A log may carry BOTH a completion record and a
    traceback; both are reported and the reader sees the conflict."""
    text = ("Job complete\n\n\nTraceback (most recent call last):\n"
            "  File \"x.py\", line 3\n")
    structure = scan(text)
    assert structure.completion_records, structure.completion_records
    assert structure.terminal_failures, structure.terminal_failures


def test_P10_terminal_failures_are_diagnostic_reach():
    assert scan("STOPWIZARD here\n").reach()["diagnostic"] is True


def test_P10_the_phenix_preamble_is_a_region_and_yields_no_sections():
    """`Texas Engineering Experiment Station` was a section: a copyright
    continuation line that happens to sit above a rule."""
    text = ("               PHENIX autosol  Sun Jul 19 10:20:17 2026\n"
            "Phenix developers include:\n"
            "  Paul Adams, Pavel Afonine\n"
            "Phenix home page:\n"
            "  http://www.phenix-online.org/\n"
            + "-" * 79 + "\n"
            "  - Texas Agricultural Experiment Station &\n"
            "    Texas Engineering Experiment Station\n"
            + "-" * 79 + "\n"
            "Adams, P.D. et al. (2010) Acta Cryst. D66, 213-221.\n"
            "\n"
            "Real content:\n" + "-" * 40 + "\n")
    structure = scan(text)
    regions = [r for r in structure.regions if r.kind == "phenix_preamble"]
    assert regions, structure.regions
    titles = [x.title for x in structure.sections]
    assert "Texas Engineering Experiment Station" not in titles, titles
    assert "Real content:" in titles, titles


def test_P10_the_report_leads_with_how_it_ended():
    """The three things a reader wants were 3rd, 4th and 6th, under 59
    sections and 342 labeled values."""
    import io
    from log_extraction.log_structure_extractor import report
    out = io.StringIO()
    report(scan("work\nSTOPWIZARD IN run/STOPWIZARD\n"), "x.log", out=out)
    text = out.getvalue()
    assert "HOW IT ENDED" in text, text
    assert text.index("HOW IT ENDED") < text.index("NOT UNDERSTOOD"), text
    # REVISION (P11): this asserted a `BACKGROUND` line that collapsed
    # sections AND labeled values into one count. Collapsing the values was
    # a regression -- they are where the numbers are -- so the report now
    # lists them grouped by label and only the SECTIONS are summarised.
    quiet = io.StringIO()
    report(scan("A heading\n" + "-" * 40 + "\nK: v\n"), "y.log", out=quiet)
    text = quiet.getvalue()
    assert "SECTIONS  1" in text, text
    assert "VALUES THE RUN REPORTED" in text, text


def test_P10_long_paths_are_elided_in_the_report():
    import io
    from log_extraction.log_structure_extractor import report
    out = io.StringIO()
    report(scan("Running refinement with /Users/t/Documents/x/y/z/model.pdb\n"),
           "z.log", out=out)
    text = out.getvalue()
    assert ".../model.pdb" in text, text
    assert "/Users/t/Documents/x/y/z" not in text, text


def test_P10_unclaimed_shapes_are_shown_not_just_counted():
    """`unclaimed=175` is a number nobody can act on."""
    import io
    from log_extraction.log_structure_extractor import report
    out = io.StringIO()
    report(scan("".join("SCALE%d  1.0  2.0  3.0  4.0\n" % i
                        for i in range(4))), "w.log", out=out)
    text = out.getvalue()
    assert "commonest unclaimed shapes" in text, text
    assert "SCALE#" in text, text


# -------------------------------------------------------- P8 measurements


def test_P8_the_agent_block_is_extracted_and_marked():
    """Requirements 4.8.3: extracting the agent's summary is useful;
    presenting it as the program's own measurement is not."""
    # a block at the TAIL is the agent's own wrapper -> source=agent
    tail = scan("program output\n" * 20 + "*" * 40 + "\n"
                "FINAL QUALITY METRICS REPORT:\n" + "-" * 40 + "\n"
                "Resolution: 1.74\nSpace Group: I 4 (No. 79)\n" + "*" * 40)
    got = [(m.name, m.value, m.source) for m in tail.measurements]
    assert got == [("Resolution", "1.74", SOURCE_AGENT),
                   ("Space Group", "I 4 (No. 79)", SOURCE_AGENT)], got

    # the same block QUOTED MID-FILE cannot be attributed -> source=unknown.
    # REVISION (P8 review): this test first asserted `agent` for both, which
    # is the hardcoding the review removed. 76 measurements in corpus2/work
    # come from the agent's own logs quoting their children mid-file.
    middle = scan("FINAL QUALITY METRICS REPORT:\n" + "-" * 40 + "\n"
                  "Resolution: 1.74\n" + "more program output\n" * 40)
    assert [m.source for m in middle.measurements] == [SOURCE_UNKNOWN], \
        [(m.name, m.source) for m in middle.measurements]


def test_P8_the_opening_rule_does_not_close_the_block():
    """DEFECT FOUND IN P8. The block opens with the marker, then a `-----`
    rule, then the metrics, and closes with `*****`. Treating any rule as the
    close ended the block immediately and the channel emitted NOTHING -- on
    every log, silently."""
    structure = scan("FINAL QUALITY METRICS REPORT:\n" + "-" * 40 + "\n"
                     "Resolution: 1.74\n")
    assert len(structure.measurements) == 1, structure.measurements


def test_P8_attached_numbers_are_re_emitted_flat():
    """A number reaches this channel because it arrived inside a structure
    already parsed. The context says which one."""
    structure = scan(" stage r-work r-free\n  1_bss: 0.4064 0.4880\n"
                     "Cycle 2  R/Rfree=0.20/0.23  Built=125\n")
    by_context = {}
    for m in structure.measurements:
        by_context.setdefault(m.context, []).append(m.name)
    assert by_context["stage:1_bss"] == ["r_free", "r_work"], by_context
    assert "cycle:2" in by_context, by_context
    assert all(m.source == SOURCE_PROGRAM for m in structure.measurements)


def test_P8_no_vocabulary_is_carried():
    """Non-goal 2.9: v1 does not grow per-program metric patterns. The module
    must contain no list of metric names to look for."""
    import log_extraction.log_structure_extractor as module
    source = inspect.getsource(module)
    # comments may DISCUSS the metrics; code must not name them. Strip
    # comment text first -- the first version flagged the module's own
    # explanation of why it has no vocabulary.
    code = "\n".join(line.split("#")[0] for line in source.split("\n"))
    for token in ("clashscore", "molprobity", "ramachandran", "rotamer",
                  "bayes", "plddt"):
        assert token not in code.lower(), token


def test_corpus_P8_attached_measurements_sit_on_their_parent_line():
    """An attached measurement claims no literal evidence of its own (its
    name is derived), so the span verifier skips it. Its provenance is the
    parent stage or cycle, which IS verified -- but only if it actually cites
    the same line. 9,433 measurements would otherwise be unchecked."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    checked = 0
    # BOTH corpora, for the same reason as the summary-row test above
    roots = [r for r in (CORPUS, CORPUS_AGENT) if r]
    for path in [p for r in roots for p in _corpus_logs(r)]:
        structure = _scanned(path)
        parents = {}
        for stage in structure.stages:
            parents.setdefault("stage:" + stage.name, set()).add(stage.line)
        for cycle in structure.cycles:
            parents.setdefault("cycle:%d" % cycle.number, set()).add(cycle.line)
        for m in structure.measurements:
            if not m.context or m.context == AGENT_METRICS_MARKER:
                continue
            assert m.line in parents.get(m.context, set()), (path, m)
            checked += 1
    assert checked > 0, "no attached measurement was seen -- checked nothing"
    print("      (%d attached measurements, all on their parent line)"
          % checked)


def test_corpus_P8_no_line_is_both_a_measurement_and_a_labeled_value():
    """Plan 2.9 requires it and nothing tested it until now."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    for path in _corpus_logs(CORPUS):
        structure = _scanned(path)
        agent_lines = set(m.line for m in structure.measurements
                          if m.source == SOURCE_AGENT)
        label_lines = set(x.line for x in structure.labeled_values)
        assert not (agent_lines & label_lines), (path,
                                                 agent_lines & label_lines)



def test_P10_the_phenix_program_header_identifies_the_program():
    """The second self-identification signal, and the shape comes from the
    Phenix SOURCE rather than from guessing at logs: ~60 programs print
    `60*"*" + 10*" " + "PHENIX <name>"`. The name is whatever precedes the
    weekday, so multi-word, CamelCase and dotted forms all survive.

    Added because an AutoSol run Tom fed the tool went UNIDENTIFIED while
    line 3 read `PHENIX autosol` -- something any casual reader gets right."""
    cases = [("               PHENIX autosol  Sun Jul 19 10:20:17 2026",
              "phenix.autosol"),
             ("          PHENIX predict_and_build  Wed Jun 17 10:06:16 2026",
              "phenix.predict_and_build"),
             ("          PHENIX Plan SAD experiment  Wed Jun 17 10:06:16 2026",
              "phenix.plan_sad_experiment"),
             ("          PHENIX AutoSol  Wed Aug 12 15:11:04 2026",
              "phenix.autosol"),
             ("          PHENIX phenix.ai_agent Wed Aug 12 15:06:26 2026",
              "phenix.ai_agent")]
    for text, expected in cases:
        ident = scan(text + "\n").identification
        assert ident.name == expected, (text, ident.name)
        assert SIGNAL_PHENIX_HEADER in ident.signals_fired, ident
    print("      (%d header forms)" % len(cases))


def test_P10_the_header_needs_a_date_after_the_name():
    """Without the weekday anchor a multi-word name has no end, and
    `PHENIX components are copyrighted by:` -- which is in the preamble of
    every Phenix log -- would become a program name."""
    for text in ("PHENIX components are copyrighted by:",
                 "PHENIX home page:",
                 "Enter phenix.acknowledgments for details."):
        assert scan(text + "\n").identification.is_unknown, text


def test_P7_a_positional_banner_identifies_the_program():
    ident = scan("Starting phenix.xtriage\non Wed by terwill\n").identification
    assert ident.name == "phenix.xtriage", ident
    assert ident.candidates[0].evidence_kind == Candidate.SELF_IDENTIFICATION
    assert ident.candidates[0].line == 1, ident
    assert ident.signals_fired == [SIGNAL_BANNER], ident


def test_P7_the_agent_prefix_does_not_hide_the_banner():
    """15 of the 20 wave-1 GUI-shape logs begin with a bare `LOG TEXT:`, and
    in 9 the payload is the banner. Marking that line as the agent's would
    have buried the single best identification signal there is."""
    ident = scan("LOG TEXT: Starting AutoBuild\n").identification
    assert ident.name == "phenix.autobuild", ident


def test_P9_a_decorated_banner_is_still_a_banner():
    """FOUND WHEN THE HOLDOUT WAS OPENED. 18 of the 87 held-out logs write
    `**Starting phenix.molprobity **` -- a form absent from the working
    corpus -- and the extractor abstained on every one. Abstention, not a
    wrong name, so the hard gate held; but it is the largest coverage miss
    the holdout exposed. The 32% holdout coverage measured BEFORE this fix
    is the number of record; anything after it is in-sample."""
    assert scan("**Starting phenix.molprobity **\n").identification.name == \
        "phenix.molprobity"
    assert scan("Starting phenix.xtriage\n").identification.name == \
        "phenix.xtriage"
    # and the guards still hold
    for text in ("Starting CC of ligand as input to map: 0.714\n",
                 "Starting AutoBuild with the command:\n",
                 "Starting phenix\n"):
        assert scan(text).identification.is_unknown, text


def test_P7_degenerate_banners_abstain_but_report_the_signal():
    """`Starting phenix` with no name (15 logs) and `Starting
    libtbx.start_process` (11 logs) identify nothing. Abstaining while
    REPORTING that the signal fired is the difference between a fingerprint
    that fails loudly and one that rots quietly."""
    for text in ("Starting phenix\non Wed\n", "Starting libtbx.start_process\n"):
        ident = scan(text).identification
        assert ident.is_unknown, (text, ident)
        assert ident.candidates == [], ident
        assert ident.signals_fired == [SIGNAL_BANNER], ident


def test_P7_a_banner_with_trailing_text_abstains():
    """LATENT DEFECT FOUND IN P7 REVIEW, before it bit. The first pattern
    captured everything to end of line. The corpus proves the unsafe shape
    exists: `Starting AutoBuild with the command:` at line ~55 in four logs,
    and ligandfit's `Starting CC of ligand as input to map: 0.714` -- a
    metric wearing a verb. Whole-line capture gives a junk name; a
    first-token rule gives `phenix.cc`, which is WRONG and fails the hard
    gate. A bare single token abstains on both."""
    for text in ("Starting CC of ligand as input to map: 0.714\n",
                 "Starting AutoBuild with the command:\n",
                 "Starting job\n"):
        ident = scan(text).identification
        assert ident.is_unknown, (text, ident.name)
    assert scan("Starting AutoBuild\n").identification.name == \
        "phenix.autobuild"


def test_P7_identification_is_computed_last():
    """Layer B is computed after every layer-A parser and used by none of
    them. The ORDER is the enforcement -- a parser cannot consult a program
    name that does not exist yet. REVISION: the first version sat mid-way
    through the parsers with a comment claiming it was last."""
    import log_extraction.log_structure_extractor as module
    body = inspect.getsource(module.scan)
    assert body.rindex("structure.identification = identify_program") > \
        body.rindex("find_labeled_values("), \
        "identification must be computed after the last layer-A parser"


def test_P7_no_signal_means_no_signal_reported():
    ident = scan("Processing files:\n").identification
    assert ident.is_unknown and ident.signals_fired == [], ident


def test_P7_a_mention_is_not_a_self_identification():
    """predict_and_build says `phenix.refine` 38 times and never names
    itself. Any frequency rule identifies every wizard as its child."""
    text = ("Running refinement\n" + "phenix.refine step\n" * 40)
    assert scan(text).identification.is_unknown


def test_P7_the_filename_is_never_used():
    """ai_analysis derived the program from the filename, and when derivation
    failed it carried on and printed a confident report about a program that
    never ran. scan() cannot see a filename at all -- asserted, because the
    API shape is what enforces it."""
    import inspect
    names = list(inspect.signature(scan).parameters)
    assert names == ["text", "program_name"], names


def _program_matches(truth, got):
    """Filename-derived truth is coarser than the banner (`phaser_mr` vs
    `phaser`), so a prefix match either way counts. Extracted into a helper
    because the D9 meta-test correctly refused a three-way disjunction inside
    an assert -- and the fix is to name the condition, not to weaken the
    guard."""
    return truth == got or truth.startswith(got) or got.startswith(truth)



def test_corpus_P7_wave1_ground_truth():
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    truth_path = os.path.join(CORPUS_AGENT, "program_truth.json")
    if not os.path.exists(truth_path):
        SKIPPED.append("no program_truth.json")
        print("      SKIP (no program_truth.json)")
        return
    import json
    handle = open(truth_path)
    try:
        truth = json.load(handle)
    finally:
        handle.close()
    named = wrong = 0
    for path in _corpus_logs(CORPUS_AGENT):
        ident = _scanned(path).identification
        if ident.is_unknown:
            continue
        expected = truth[os.path.basename(path)]
        if ident.name == expected:
            named += 1
        else:
            wrong += 1
    assert wrong == 0, wrong
    print("      (%d of %d agent-shape logs named, 0 wrong; the coverage"
          " target is a full-corpus claim, see validate.py)"
          % (named, len(_corpus_logs(CORPUS_AGENT))))


# ----------------------------------------------------- P6 completion records


def test_P6_a_three_line_ending_is_one_record():
    """The Program Template closes with `Job complete` / `usr+sys time:` /
    `wall clock time:` on consecutive lines. Counting three would treble
    every completion count in the corpus."""
    structure = scan("Job complete\nusr+sys time: 5 seconds\n"
                     "wall clock time: 6 seconds\n")
    assert len(structure.completion_records) == 1, \
        structure.completion_records
    assert structure.completion_records[0].line == 1


def test_P6_finished_with_names_its_child():
    """The only positive evidence of attribution the corpus offers -- found
    by reading the shapes, not assumed. 113 such lines."""
    records = scan("Finished with predict_chain\nmore output\n"
                   ).completion_records
    assert [(r.applies_to) for r in records] == ["predict_chain"], records


def test_P6_repeated_markers_of_one_kind_stay_separate():
    """CORRECTION FOUND IN P6 REVIEW. The planning note said phaser writes
    `EXIT STATUS: SUCCESS` per module and that one log carried 59 of them.
    Re-measured: `EXIT STATUS` appears EXACTLY ONCE in each of the 24 logs
    carrying it; the per-module marker is `Finished: <timestamp>` -- 59 in
    one phaser log, 34 and 29 in others. The conclusion (a list, not a
    scalar) was right; the evidence cited for it was wrong, and it was
    repeated in the plan and three handoffs before being checked."""
    text = "".join("Finished: Wed Aug 12 15:16:1%d 2026\nwork\n" % i
                   for i in range(5))
    records = scan(text).completion_records
    assert len(records) == 5, records
    assert all(r.applies_to == "unknown" for r in records), records


def test_P6_a_child_ending_does_not_absorb_the_runs_own():
    """FOUND BY A MISSED PREDICTION (P6.3). In 26 logs a child's `Finished
    with process_predicted_model` sat within six lines of the program's own
    `Job complete` / `usr+sys` / `wall clock` ending. Every marker kind
    differed, so they clustered -- and the run's own ending vanished into a
    record attributed to the child. `Finished with X` is self-contained."""
    structure = scan("Finished with process_predicted_model\n\n"
                     + "=" * 40 + "\nJob complete\n"
                     "usr+sys time: 2.43 seconds\n"
                     "wall clock time: 2.69 seconds\n")
    got = [(r.text, r.applies_to) for r in structure.completion_records]
    assert len(got) == 2, got
    assert got[0][1] == "process_predicted_model", got
    assert got[1] == ("Job complete", "top_level"), got


def test_P6_top_level_requires_positive_evidence():
    """A child completing is not the run completing. `top_level` only when
    the record is the last event with nothing but blank lines after it."""
    last = scan("Job complete\n\n\n").completion_records
    assert last[0].applies_to == "top_level", last

    followed = scan("Job complete\nstill working\nJob complete\nmore\n"
                    ).completion_records
    assert [r.applies_to for r in followed] == ["unknown", "unknown"], followed


def test_P6_an_empty_list_is_not_a_verdict():
    """Non-goal 2.8. Of 30 failed runs, 27 carry no error keyword at all in
    what the program itself wrote, so absence of a completion record is not
    evidence of failure -- and absence of error is not evidence of success.
    The API must expose no way to read it as either."""
    structure = scan("some output with no ending\n")
    assert structure.completion_records == []

    # REVISION (P10): this banned any attribute containing "failure", and
    # fired on the new `terminal_failures` CHANNEL. Loosening a non-goal
    # guard is exactly how a non-goal dies, so the distinction is written
    # down rather than assumed: what 2.8 forbids is a VERDICT about the run,
    # not a channel of things the program said. The ban is now on
    # verdict-shaped names, and every observation channel must be a LIST --
    # a verdict would be a bool or a string.
    banned = ("success", "succeeded", "failed", "is_failure", "outcome",
              "verdict", "status", "passed", "ok")
    for name in dir(structure):
        if name.startswith("_"):
            continue
        lowered = name.lower()
        for word in banned:
            assert word not in lowered, (name, word)
    assert isinstance(structure.terminal_failures, list)
    assert isinstance(structure.completion_records, list)


def test_P6_completion_is_basic_reach_not_diagnostic():
    reach = scan("Job complete\n").reach()
    assert reach["basic"] is True, reach
    assert reach["diagnostic"] is False, reach


def test_corpus_P6_completion_records_and_their_attribution():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    tally = {}
    # existence claims below ("some record names a child") need the whole
    # corpus; the per-log assertions are sampling-safe either way
    for path in _corpus_logs(CORPUS, stride=1):
        structure = _scanned(path)
        lines = split_lines(_read(path))
        top = [r for r in structure.completion_records
               if r.applies_to == CompletionRecord.TOP_LEVEL]
        assert len(top) <= 1, (path, top)
        for record in structure.completion_records:
            assert 1 <= record.line <= len(lines), (path, record)
            assert record.text.strip(), record
            tally[record.applies_to] = tally.get(record.applies_to, 0) + 1
    assert CompletionRecord.TOP_LEVEL in tally, tally
    assert CompletionRecord.UNKNOWN in tally, tally
    named = [k for k in tally if k not in (CompletionRecord.TOP_LEVEL,
                                           CompletionRecord.UNKNOWN)]
    assert named, "no record was attributed to a named child"
    print("      (top_level=%d, unknown=%d, %d named children)"
          % (tally[CompletionRecord.TOP_LEVEL],
             tally[CompletionRecord.UNKNOWN], len(named)))



def test_P5_setting_value_and_reason_are_separate():
    decisions, claims, _ = find_decisions(
        ["Setting rebuild_in_place=False as maps_only=True",
         "Setting n_cycle_build=1 as nbatch >1 (nbatch =3)"])
    assert [(d[0], d[1], d[2]) for d in decisions] == [
        ("rebuild_in_place", "False", "maps_only=True"),
        ("n_cycle_build", "1", "nbatch >1 (nbatch =3)")], decisions
    assert sorted(claims) == [1, 2], claims


def test_P5_a_setting_with_no_reason_is_not_a_decision():
    """Requirements 4.4: a decision is the program changing its own
    configuration AND STATING THE REASON -- announcing a branch. A value with
    no reason is not a branch. MEASURED: 1,075 of the 1,283 `Setting` lines
    in corpus2/work carry no reason; they are reported, not pattern-fitted."""
    structure = scan("Setting output.overwrite to True\n"
                     "Setting macro_cycles to  3\n")
    assert structure.decisions == [], structure.decisions
    refused = [u for u in structure.unparsed
               if u.excluded_by == REFUSAL_SETTING_NO_REASON]
    assert len(refused) == 2, [repr(u) for u in structure.unparsed]


def test_P5_setting_up_is_a_gerund_not_a_decision():
    structure = scan("Setting up prediction keywords...\n")
    assert structure.decisions == [], structure.decisions


def test_P5_refusals_say_accurately_why():
    """FOUND IN P5 REVIEW. Every refused `Setting` line carried
    `setting_states_no_reason`, and for 490 of the 1,075 that is the wrong
    diagnosis -- a reader auditing the ledger would hunt for a missing reason
    on lines that are not assignments at all. Same failure as D20 and D31."""
    cases = [("Setting output.overwrite to True", REFUSAL_SETTING_NO_REASON),
             ("Setting macro_cycles=3", REFUSAL_SETTING_NO_REASON),
             ("Setting up prediction keywords...", REFUSAL_SETTING_UP),
             ("Setting default value of  True  for  find_ncs",
              REFUSAL_SETTING_NARRATION),
             ("Setting parameters", REFUSAL_SETTING_NARRATION)]
    for text, expected in cases:
        refused = [u.excluded_by for u in scan(text + "\n").unparsed
                   if u.excluded_by]
        assert refused == [expected], (text, refused)


def test_corpus_P5_every_refused_setting_carries_one_of_the_three():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    counts = {}
    # "all three reasons occur" is an existence claim -> whole corpus
    for path in _corpus_logs(CORPUS, stride=1):
        for item in _scanned(path).unparsed:
            if item.excluded_by in (REFUSAL_SETTING_NO_REASON,
                                    REFUSAL_SETTING_UP,
                                    REFUSAL_SETTING_NARRATION):
                counts[item.excluded_by] = counts.get(item.excluded_by, 0) + 1
    assert len(counts) == 3, counts
    print("      (%s)" % ", ".join("%s=%d" % kv for kv in sorted(counts.items())))


def test_P5_decisions_are_diagnostic():
    assert scan("Setting a=b as c\n").reach()["diagnostic"] is True


def test_corpus_P5_the_failed_build_states_its_branch():
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    path = os.path.join(CORPUS_AGENT, "autobuild_4_bromodomain.log")
    if not os.path.exists(path):
        SKIPPED.append("not in this corpus")
        print("      SKIP (not in this corpus)")
        return
    structure = scan(_read(path))
    settings = dict((d.setting, (d.value, d.reason)) for d in structure.decisions)
    assert "rebuild_in_place" in settings, settings
    assert settings["rebuild_in_place"] == ("False", "maps_only=True"), settings
    print("      (%d decisions)" % len(structure.decisions))


def test_corpus_P5_every_decision_states_a_reason():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    total = 0
    for path in _corpus_logs(CORPUS):
        for decision in _scanned(path).decisions:
            assert decision.reason and decision.reason.strip(), (path, decision)
            assert decision.setting and decision.value, (path, decision)
            total += 1
    # a floor, not a count. The COUNT is a full-corpus claim and lives in
    # validate.py; here we only need the invariant to have had work to do.
    assert total > 0, total
    print("      (%d decisions, all with a stated reason)" % total)


# ------------------------------------------------ P4 stage tables and cycles


def test_P4_rows_are_anchored_to_their_header():
    """A global row pattern matches 40 lines in refine_5_beta_blip where only
    29 belong to the table. Rows are a contiguous run after a header."""
    lines = [" stage r-work r-free bonds",
             "      0    : 0.4747 0.5371 0.010",
             "      1_bss: 0.4064 0.4880 0.010",
             "----------",
             "  other  : 0.1 0.2 0.3"]
    stages, claims, _ = find_stage_tables(lines)
    assert [x[0] for x in stages] == ["0", "1_bss"], stages


def test_P4_the_end_row_is_a_summary_not_a_stage():
    """THE trap. The table's last row is `end` with r_free 0.4942, while the
    last real stage carries 0.4935 -- and 0.4942 is exactly the number the
    consuming project identified as HIDING the finding. A parser that keeps
    `end` reports the summary and loses the trajectory."""
    lines = [" stage r-work r-free",
             "      0    : 0.4747 0.5371",
             "    3_adp: 0.3966 0.4935",
             "      end: 0.3966 0.4942"]
    stages, _, _ = find_stage_tables(lines)
    assert [x[0] for x in stages] == ["0", "3_adp"], stages
    assert stages[-1][2]["r_free"] == 0.4935, stages[-1]


def test_P4_the_end_row_is_reported_with_its_own_reason():
    """FOUND IN P4 REVIEW. The summary row was claimed and emitted nowhere.
    It surfaced anyway only because it sits above a rule and the SECTION
    parser refused it as `section_title_is_a_numeric_row` -- an accident, and
    a misleading diagnosis. 77 such rows across the two corpora were relying
    on that coincidence."""
    structure = scan(" stage r-work r-free\n"
                     "  3_adp: 0.3966 0.4935\n"
                     "    end: 0.3966 0.4942\n"
                     "next line, no rule below\n")
    assert [x.name for x in structure.stages] == ["3_adp"], structure.stages
    refused = [u for u in structure.unparsed
               if u.excluded_by == REFUSAL_STAGE_SUMMARY]
    assert len(refused) == 1, [repr(u) for u in structure.unparsed]
    assert refused[0].line == 3, refused[0]


def test_corpus_P4_every_end_row_is_accounted_for():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    import log_extraction.log_structure_extractor as module
    total = 0
    # BOTH corpora. On the smoke set the stage tables live in the
    # agent-shape logs, so reading CORPUS alone made this report
    # "0 summary rows, all named correctly" -- a test that verified nothing
    # and said PASS. Same class as the vacuous probe in P4.
    roots = [r for r in (CORPUS, CORPUS_AGENT) if r]
    for path in [p for r in roots for p in _corpus_logs(r)]:
        text = _read(path)
        lines = split_lines(text)
        structure = _scanned(path)
        reported = dict((u.line, u.excluded_by) for u in structure.unparsed)
        for index, raw in enumerate(lines):
            row = module._STAGE_ROW_RE.match(raw)
            if not row or row.group(1) != module.STAGE_SUMMARY_ROW:
                continue
            total += 1
            assert reported.get(index + 1) == REFUSAL_STAGE_SUMMARY, \
                (path, index + 1, reported.get(index + 1))
    assert total > 0, "no stage summary row was seen -- this checked nothing"
    print("      (%d summary rows, all named correctly)" % total)


def test_P4_the_legend_names_the_stage_and_absence_is_not_invented():
    lines = ["1_bss: bulk solvent correction and/or scaling",
             " stage r-work r-free",
             "  1_bss: 0.40 0.48",
             "  9_zzz: 0.41 0.49"]
    stages, _, _ = find_stage_tables(lines)
    named = dict((x[0], x[3]) for x in stages)
    assert named["1_bss"] == "bulk solvent correction and/or scaling", named
    assert named["9_zzz"] is None, named


def test_P4_r_over_rfree_is_one_key_and_two_values():
    cycles, _, _ = find_cycles(
        ["Cycle 2  R/Rfree=0.20/0.23  Built=125  Placed=125 Resolution=2.5 A"])
    metrics = cycles[0][2]
    assert metrics["r_work"] == 0.20 and metrics["r_free"] == 0.23, metrics
    assert metrics["Built"] == 125 and metrics["Placed"] == 125, metrics


def test_P4_the_sentinel_is_flagged_and_never_quoted():
    """999.90 means `no usable result`. Passing it through as an R-factor is
    worse than omitting it -- a consumer would quote a nonsense number."""
    cycles, _, _ = find_cycles(
        ["Cycle 2  R/Rfree=999.90/999.90  Built=0  Placed=0 Resolution=2.1 A"])
    number, line, metrics, sentinel = cycles[0]
    assert sentinel is True, cycles[0]
    assert "r_work" not in metrics and "r_free" not in metrics, metrics
    assert metrics["Built"] == 0, metrics
    assert CYCLE_SENTINEL not in metrics.values(), metrics


def test_P4_a_bare_counter_is_not_a_cycle_record():
    structure = scan("Cycle 3 of morphing chain trace\n")
    assert structure.cycles == [], structure.cycles
    refused = [u for u in structure.unparsed
               if u.excluded_by == REFUSAL_CYCLE_NO_METRICS]
    assert len(refused) == 1, structure.unparsed


def test_corpus_P4_the_refine_trajectory():
    """The finding that made the case for findings-first evidence: bulk
    solvent gained 0.049 and coordinate and ADP refinement gave it back."""
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    path = os.path.join(CORPUS_AGENT, "refine_5_beta_blip.log")
    if not os.path.exists(path):
        SKIPPED.append("not in this corpus")
        print("      SKIP (not in this corpus)")
        return
    structure = scan(_read(path))
    assert len(structure.stages) == 28, len(structure.stages)
    assert structure.stages[0].metrics["r_free"] == 0.5371
    assert structure.stages[-1].metrics["r_free"] == 0.4935
    worse = [m["stage"] for m in structure.metric_moves("r_free", 0.001)
             if m["delta"] > 0]
    assert worse == ["1_xyzrec", "1_realsrl2", "1_adp"], worse
    assert [x for x in structure.stages if x.description], "no legend applied"


def test_corpus_P4_the_failed_build_is_visible():
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    path = os.path.join(CORPUS_AGENT, "autobuild_4_bromodomain.log")
    if not os.path.exists(path):
        SKIPPED.append("not in this corpus")
        print("      SKIP (not in this corpus)")
        return
    structure = scan(_read(path))
    sentinels = [c for c in structure.cycles if c.sentinel]
    assert sentinels, structure.cycles
    assert sentinels[-1].metrics.get("Built") == 0, sentinels[-1]
    assert any("model-building" in x.what for x in structure.control_skips), \
        structure.control_skips


def test_corpus_P4_no_stage_metric_is_the_sentinel():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    for path in _corpus_logs(CORPUS):
        for cycle in _scanned(path).cycles:
            assert CYCLE_SENTINEL not in cycle.metrics.values(), (path, cycle)


# ------------------------------------------- P3 phases, skips and exclusions


def test_P3_running_is_a_phase_and_starting_is_not():
    """An early draft of the requirements had this backwards. autobuild
    writes `Starting job`, `Starting mtz file`, `Starting model list`;
    ligandfit writes `Starting CC of ligand as input to map: 0.714`, which is
    a metric wearing a verb."""
    structure = scan("Running dock_and_rebuild on total of 1 chains\n"
                     "Starting job\n"
                     "Starting CC of ligand as input to map: 0.714\n")
    assert [x.name for x in structure.phases] == [
        "dock_and_rebuild on total of 1 chains"], structure.phases


def test_P3_a_control_skip_is_not_an_item_exclusion():
    """Merging the two either invents a phase that never ran or loses the
    user complaint the consuming project exists to answer."""
    structure = scan(
        "Skipping model-building as 'maps_only' is set\n"
        "Skipping 9GSD.A - not x-ray and not computational\n")
    assert [(x.what, x.reason) for x in structure.control_skips] == [
        ("model-building", "'maps_only' is set")], structure.control_skips
    assert [(x.item, x.reason) for x in structure.exclusions] == [
        ("9GSD.A", "not x-ray and not computational")], structure.exclusions


def test_P3_an_item_name_is_optional():
    """`Skipping - not protein or too few residues.` names no candidate.
    Requiring a name silently dropped 3 of find_reference's 18."""
    structure = scan("Skipping - not protein or too few residues.\n")
    assert len(structure.exclusions) == 1, structure.exclusions
    assert structure.exclusions[0].item is None, structure.exclusions[0]
    assert structure.exclusions[0].reason == "not protein or too few residues."


def test_P3_a_control_skip_with_no_subject_is_not_invented():
    """FOUND BY A MISSED PREDICTION (P3.3, band zero). An exclusion may be
    anonymous -- a rejected candidate with no name is still a rejection. A
    control skip may NOT: `Skipping as CC_mask is < 1/2 min_cc` does not say
    which phase did not run, and emitting what=None would invent a phase that
    never existed. One such line in corpus2/work."""
    structure = scan("Skipping as CC_mask is < 1/2 min_cc\n")
    assert structure.control_skips == [], structure.control_skips
    assert structure.exclusions == [], structure.exclusions
    refused = [u for u in structure.unparsed
               if u.excluded_by == REFUSAL_SKIP_UNRECOGNISED]
    assert len(refused) == 1, structure.unparsed
    # but an anonymous EXCLUSION is still emitted
    assert len(scan("Skipping - not protein.\n").exclusions) == 1


def test_P3_an_unrecognised_skip_is_reported_not_dropped():
    """MEASURED: 233 of the 381 `Skipping` lines in corpus2/work match
    neither form -- `Skipping remainder region 1 (already written out)`.
    Requirements 6 asks that they be reported. They are NOT fitted with a
    third pattern: a control-versus-item distinction guessed from a corpus I
    have already read is the kind of wrong that stays invisible."""
    structure = scan("Skipping ...already tested this set\n")
    assert structure.control_skips == [] and structure.exclusions == []
    refused = [u for u in structure.unparsed
               if u.excluded_by == REFUSAL_SKIP_UNRECOGNISED]
    assert len(refused) == 1, structure.unparsed


def test_P3_a_line_can_be_both_a_phase_and_a_section():
    """Requirements 3: nothing is deduplicated across kinds. SPECIFIC parsers
    ignore other parsers' claims; the GENERIC channel respects them. A
    heading that is also a phase is two true facts; a heading also captured
    as a nameless fact is one fact counted twice."""
    structure = scan("Running refine_ca to refine and make full model\n"
                     + "-" * 40 + "\n")
    assert len(structure.phases) == 1, structure.phases
    assert len(structure.sections) == 1, structure.sections
    assert structure.phases[0].line == structure.sections[0].line
    assert not [x for x in structure.labeled_values
                if x.line == structure.phases[0].line]


def test_P3_skips_and_exclusions_are_diagnostic_phases_are_not():
    """Announcing that something began is not saying what happened."""
    assert scan("Running something\n").reach()["diagnostic"] is False
    assert scan("Skipping x - bad\n").reach()["diagnostic"] is True
    assert scan("Skipping x as y is set\n").reach()["diagnostic"] is True


def test_corpus_P3_the_complaint_case_is_exclusions_not_skips():
    """THE case the consuming project exists to answer: chain A's search
    returned 9 hits and 8 were skipped, including the user's own input
    structure. 18 exclusions in 2 reason groups, 0 control skips. No band --
    an off-by-one here is a defect."""
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    path = os.path.join(CORPUS_AGENT, "find_reference_14.log")
    if not os.path.exists(path):
        SKIPPED.append("find_reference_14.log not in this corpus")
        print("      SKIP (find_reference_14.log not in this corpus)")
        return
    structure = scan(_read(path))
    assert len(structure.exclusions) == 18, len(structure.exclusions)
    assert structure.control_skips == [], structure.control_skips
    groups = structure.exclusion_groups()
    assert [g["count"] for g in groups] == [15, 3], groups
    assert "not x-ray" in groups[0]["reason"], groups[0]
    assert any(x.item is None for x in structure.exclusions), \
        "the unnamed exclusions were dropped"


def test_corpus_P3_no_phase_comes_from_a_starting_line():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    for path in _corpus_logs(CORPUS):
        lines = split_lines(_read(path))
        for phase in _scanned(path).phases:
            raw = lines[phase.line - 1].strip()
            assert not raw.startswith("Starting"), (path, phase)


# --------------------------------------------------------- P2 labeled values


def test_SCREEN_v2_single_character_labels_are_admitted():
    """SCREEN VERSION 2, signed off. `{1,48}` after the leading letter meant
    2-49 characters, so `R:` could not match -- and 115 lines across 33 logs
    carrying R and R-free were invisible to the screen ENTIRELY, not merely
    unclaimed. Known limitation: the line holds two facts and is captured as
    one; splitting compound lines is interpretation."""
    assert screen_line("R:   0.45  Rfree:   0.46") == "key_value"
    values = scan("R:   0.45  Rfree:   0.46\n").labeled_values
    # REVISION (P11): this used to assert the whole rest of the line as one
    # value, and recorded that as a known limitation. The compound-pair
    # split fixed it -- the line is now two values, which is what it says.
    assert [(x.label, x.value) for x in values] == [
        ("R", "0.45"), ("Rfree", "0.46")], values
    # and the upper bound did not move: a 49-character label still matches
    assert screen_line("R" * 49 + ": x") == "key_value"


def test_corpus_SCREEN_v2_the_r_factor_lines_are_now_seen():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    seen = 0
    # existence test, no sampling -- and it must look in BOTH corpora, since
    # the R/Rfree lines in the smoke set are in agent-shape logs
    for root in [r for r in (CORPUS, CORPUS_AGENT) if r]:
        for path in _corpus_logs(root, stride=1):
            for item in _scanned(path).labeled_values:
                if item.label == "R" and "Rfree" in item.value:
                    seen += 1
    assert seen > 0, "no R/Rfree line was captured"
    print("      (%d R/Rfree lines captured)" % seen)


def test_P2_a_key_value_line_becomes_a_labeled_value():
    structure = scan("Space group: P 21 21 21\nResolution: 2.10 A\n")
    got = [(x.label, x.value, x.line) for x in structure.labeled_values]
    assert got == [("Space group", "P 21 21 21", 1),
                   ("Resolution", "2.10 A", 2)], got


def test_P2_values_are_not_interpreted():
    """No unit stripping, no number parsing, no name normalisation. That
    mapping is the treadmill this tool exists to escape."""
    values = dict((x.label, x.value) for x in scan(
        "R Free: 0.2321\nr_free: 0.2321\nResolution: 2.10 A\n"
    ).labeled_values)
    assert values["R Free"] == "0.2321", values
    assert values["r_free"] == "0.2321", values      # NOT merged with R Free
    assert values["Resolution"] == "2.10 A", values  # unit not stripped


def test_P2_identical_label_value_pairs_collapse():
    """The same fact stated twice is one fact. repeat_count holds the number
    of occurrences and line/end hold the first and last, so nothing about
    where it happened is lost."""
    structure = scan("CC: 0.75\nfiller\nCC: 0.75\nCC: 0.80\n")
    got = [(x.label, x.value, x.line, x.end, x.repeat_count)
           for x in structure.labeled_values]
    assert got == [("CC", "0.75", 1, 3, 2), ("CC", "0.80", 4, 4, 1)], got


def test_P2_a_runaway_label_collapses_and_says_so():
    """MEASURED: `Target left` occurs 1,999 times in one log. A runaway
    collapses to one record and the lines it would have produced are
    REPORTED as rule_excluded, never silently dropped."""
    text = "".join("Target left: %d\n" % i
                   for i in range(LABEL_DISTINCT_LIMIT + 5))
    structure = scan(text)
    assert len(structure.labeled_values) == 1, structure.labeled_values
    record = structure.labeled_values[0]
    assert record.repeat_count == LABEL_DISTINCT_LIMIT + 5, record
    assert record.line == 1 and record.end == LABEL_DISTINCT_LIMIT + 5, record
    refused = [u for u in structure.unparsed
               if u.excluded_by == REFUSAL_LABEL_RUNAWAY]
    assert len(refused) == LABEL_DISTINCT_LIMIT + 4, len(refused)


def test_P2_a_series_just_below_the_limit_is_kept_whole():
    """The limit is a guard against runaway, not a filter. A per-cycle series
    -- `CC` with 26 distinct values in one corpus log -- must survive."""
    text = "".join("CC: 0.%02d\n" % i for i in range(LABEL_DISTINCT_LIMIT))
    structure = scan(text)
    assert len(structure.labeled_values) == LABEL_DISTINCT_LIMIT, \
        len(structure.labeled_values)


def test_P2_a_section_title_is_not_also_a_labeled_value():
    """Lines claimed by another parser are skipped. Without this a heading
    would be counted twice and the reach metrics would double-count it.

    REVISION (P2 controls): the first fixture used `Processing files:` as the
    heading -- which ends in a colon and so is not a key_value candidate at
    all. The test therefore passed with the claim check DISABLED, and the
    control reported no effect. The heading here is deliberately one that IS
    a key_value candidate."""
    structure = scan("Space group: P 21 21 21\n" + "-" * 40 + "\n"
                     "Found: /x.eff\n")
    assert [x.title for x in structure.sections] == ["Space group: P 21 21 21"]
    assert [x.label for x in structure.labeled_values] == ["Found"], \
        structure.labeled_values


def test_P2_generic_only_is_reported_separately():
    """THE ANTI-GAMING SPLIT. A generic capture channel that claimed every
    Key: value line would drive `unclaimed` to near zero by construction --
    the metric defeated by the feature meant to help the long tail. A line
    recorded but not understood is `generic_only`, not `unclaimed`."""
    structure = scan("Space group: P 21 21 21\n")
    counts = structure.unparsed_counts()
    assert counts[GENERIC_ONLY] == 1, counts
    assert counts[UNCLAIMED] == 0, counts
    assert structure.labeled_values, structure.labeled_values


def test_P2_labeled_values_lift_basic_reach_not_diagnostic():
    """Plan 6: three reach metrics, always together. A labeled value is
    basic information; it is not a stage, cycle or decision."""
    reach = scan("Space group: P 21 21 21\n").reach()
    assert reach["basic"] is True, reach
    assert reach["diagnostic"] is False, reach


def test_corpus_P2_every_labeled_value_is_on_its_cited_line():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    total = collapsed = 0
    for path in _corpus_logs(CORPUS):
        text = _read(path)
        lines = split_lines(text)
        for item in _scanned(path).labeled_values:
            assert 1 <= item.line <= item.end <= len(lines), (path, item)
            assert item.label in lines[item.line - 1], (path, item)
            assert item.value in lines[item.line - 1], (path, item)
            assert item.repeat_count >= 1, item
            total += 1
            if item.repeat_count > 1:
                collapsed += 1
    assert total > 0, total
    print("      (%d labeled values, %d collapsed groups)"
          % (total, collapsed))


def test_corpus_P2_no_line_is_both_claimed_and_unclaimed():
    """The ledger must balance: a line is claimed, generically claimed, or
    unclaimed -- never two of those."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    for path in _corpus_logs(CORPUS):
        structure = _scanned(path)
        seen = {}
        for item in structure.unparsed:
            assert item.line not in seen, (path, item.line, seen[item.line],
                                           item.status)
            seen[item.line] = item.status


# ------------------------------------------- regressions: one test per defect
#
# Every defect found in a review or by a control gets a test here, named for
# what it protects, so the same mistake cannot come back quietly. A defect
# without a regression test is a defect we have merely moved.


def test_D5_bool_is_not_a_line_number():
    """bool is a subclass of int, so Phase("x", True) silently became line 1."""
    try:
        Phase("x", True)
    except ValueError:
        return
    raise AssertionError("True was accepted as a line number")


def test_D6_module_is_free_of_python3_only_syntax():
    """The module must run under libtbx.python. The docstring's own example
    used open(path, errors="replace"), which is Python 3 only -- in a file
    whose header promises Python 2 compatibility."""
    path = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "log_extraction", "log_structure_extractor.py")
    handle = open(path)
    try:
        source = handle.read()
    finally:
        handle.close()
    banned = [("f-string", re.compile(r"""[^\w]f["']""")),
              ("dataclass", re.compile(r"\bdataclass\b")),
              # REVISION (P5 review): `(?:=` in a non-capturing group
              # contains `:=` and tripped this. A meta-test with a false
              # positive is worse than none -- same class as D30.
              ("walrus", re.compile(r"(?<!\?):=")),
              ("open(..., errors=)", re.compile(r"open\([^)]*errors\s*=")),
              ("py3 print-file kwarg in docstring",
               re.compile(r"print\([^)]*, file=")),
              ("nonlocal", re.compile(r"\bnonlocal\b"))]
    offenders = [name for name, pattern in banned if pattern.search(source)]
    assert not offenders, offenders
    assert "from __future__ import" in source


def test_D7_section_exposes_both_start_and_line():
    """Requirements 3 names a section's first line `start`; every other item
    calls it `line`. The frozen conformance suite reads `start`."""
    section = Section("Collecting inputs", 5, end=9)
    assert section.start == section.line == 5
    assert section.end == 9


def test_D9_no_test_in_this_file_asserts_a_disjunction():
    """A test assertion of the form `assert A or B` can pass for the wrong
    reason. One such assertion sat in this file and could not fail.

    Meta-tests are usually a bad idea; this one earns its place because the
    failure it guards against is invisible to a green suite by definition."""
    path = os.path.abspath(__file__)
    handle = open(path)
    try:
        body = handle.read()
    finally:
        handle.close()
    # REVISION (P3): the first version searched the raw line, so it fired
    # on `assert x.reason == "not protein or too few residues."` -- the word
    # `or` inside a STRING LITERAL. A meta-test with a false positive is
    # worse than none, because the tempting fix is to weaken the assertion
    # it flags. String content is removed before the check.
    offenders = []
    for number, line in enumerate(body.split("\n"), 1):
        stripped = line.strip()
        if not stripped.startswith("assert "):
            continue
        code = re.sub(r"'[^']*'", "", re.sub(r'"[^"]*"', "", stripped))
        code = code.split("#")[0]
        if re.search(r"\bor\b", code) and "offenders" not in code:
            offenders.append((number, stripped[:60]))
    assert not offenders, offenders


def test_D11_screen_fixtures_are_real_corpus_lines():
    """Two invented fixtures asserted the opposite of what they claimed --
    a hand-written pipe-table row given three numeric cells instead of four.
    Fixtures that claim to represent the corpus are now checked against it."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    # Both come from the SMOKE corpus, which is what ships. The previous
    # pair came from the full corpus and vanished when it moved out of the
    # tree -- exactly the fixture drift this test exists to catch, caught by
    # this test.
    fixtures = [
        "  | 1      | 0.635           | 3.37            | 0.879              |",
        "            Running build_chains on segment 1",
    ]
    seen = dict((f, False) for f in fixtures)
    # EXISTENCE tests must see the WHOLE corpus. FOUND IN P8: the control
    # runner samples with a stride, and at stride 12 this test and the
    # R-factor one failed in the BASELINE -- so every control number that run
    # was measured against a red baseline. A test that asserts "this line
    # exists somewhere" cannot be sampled.
    roots = [r for r in (CORPUS, CORPUS_AGENT) if r]
    for root in roots:
        for path in _corpus_logs(root, stride=1):
            text = _read(path)
            for fixture in fixtures:
                if not seen[fixture] and fixture in text:
                    seen[fixture] = True
            if all(seen.values()):
                break
    missing = [f for f, found in seen.items() if not found]
    assert not missing, ("invented fixtures: %r" % missing)
    print("      (%d fixtures found in corpus2/work)" % len(fixtures))


def test_D12_every_negative_control_still_patches_something():
    """Control C8 silently reported 'patch target not found' after an
    unrelated edit moved the line it patched. Controls need the same
    maintenance as tests, and nothing was checking them."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, root)
    from log_extraction import run_controls
    handle = open(os.path.join(root, "log_extraction", "log_structure_extractor.py"))
    try:
        source = handle.read()
    finally:
        handle.close()
    # EXACTLY ONCE, not merely present. FOUND IN REVIEW: two controls had
    # patterns occurring 5 and 2 times, and `replace(..., 1)` silently
    # patched the FIRST site -- C34 was disabling section refusals while
    # claiming to disable decision refusals, and had been giving a false
    # result since P5. This is the third instance of the same failure
    # (D28 was the first), so the check is now uniqueness, not existence.
    stale = [name for name, _, (pattern, _) in run_controls.CONTROLS
             if pattern not in source]
    assert not stale, ("patch target absent: %s" % stale)
    ambiguous = [(name, source.count(pattern))
                 for name, _, (pattern, _) in run_controls.CONTROLS
                 if source.count(pattern) != 1]
    assert not ambiguous, ("patch target is not unique: %s" % ambiguous)
    print("      (%d controls, all targets present)"
          % len(run_controls.CONTROLS))


def test_D13_line_numbers_match_grep_not_splitlines():
    """str.splitlines() breaks on bare CR, VT, FF, NEL and U+2028 as well as
    newline. corpus2/work/ok/fem_407.log carries 100 bare CRs: splitlines()
    gives 539 lines where grep -n gives 439, so every cited line after the
    first CR pointed a human at the wrong place."""
    for breaker in ("\r", "\x0b", "\x0c", "\x85", "\u2028"):
        text = "alpha%sbeta\ngamma\n" % breaker
        assert len(split_lines(text)) == 2, (breaker, split_lines(text))
        assert scan(text).n_lines == 2, breaker
    assert split_lines("a\r\nb\r\n") == ["a", "b"]      # CRLF
    assert split_lines("a\nb") == ["a", "b"]            # no trailing newline
    assert split_lines("") == []


def test_D14_the_agent_prefix_needs_no_header_above_it():
    """15 of the 20 wave-1 GUI-shape logs begin with a bare `LOG TEXT:` and no
    header at all. The first version only emitted a prefix region when a
    header preceded it, so the commonest case in the shipped GUI corpus went
    unmarked."""
    regions = scan("LOG TEXT: Starting AutoBuild\nProcessing files:\n").regions
    assert [r.kind for r in regions] == ["agent_prefix"], regions
    assert regions[0].source == SOURCE_PROGRAM, regions
    assert (regions[0].start, regions[0].end) == (1, 1), regions


def test_D13_corpus_no_log_disagrees_with_newline_counting():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    disagree = []
    for path in _corpus_logs(CORPUS):
        text = _read(path)
        if len(text.splitlines()) != scan(text).n_lines:
            disagree.append((os.path.basename(path), len(text.splitlines()),
                             scan(text).n_lines))
    # This is expected to be NON-empty: it names the logs where the naive
    # reading would have been wrong. It must be reported, not asserted away.
    print("      (%d logs where splitlines() would have been wrong: %s)"
          % (len(disagree), disagree or "none"))


# ------------------------------------------------------------ corpus invariants


def test_corpus_every_log_scans_without_raising():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    logs = _corpus_logs(CORPUS)
    for path in logs:
        scan(_read(path))
    print("      (%d logs)" % len(logs))


def test_corpus_every_item_cites_a_line_that_exists_and_matches():
    """SPAN VERIFICATION. Hand-authored spans in the consuming project cited
    one log's line numbers against another log and validated cleanly; only a
    mechanical verifier caught it."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    checked = 0
    for path in _corpus_logs(CORPUS):
        text = _read(path)
        # MUST use the extractor's own splitting. The first version used
        # str.splitlines() here and the verifier then validated spans against
        # a DIFFERENT line map than the one the items were numbered in --
        # a verifier checking the wrong ruler. It failed on fem_407.log,
        # the one corpus log carrying bare CR characters.
        lines = split_lines(text)
        for item in _scanned(path).items():
            assert 1 <= item.line <= len(lines), (path, item)
            assert item.line <= item.end <= len(lines), (path, item)
            evidence = item.evidence()
            if evidence:
                assert evidence in lines[item.line - 1], (path, item)
            checked += 1
    assert checked > 0, checked
    print("      (%d items verified)" % checked)


def test_corpus_identification_never_gates_extraction():
    """The layering invariant, at corpus scale rather than on one fixture."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    logs = _corpus_logs(CORPUS)
    for path in logs:
        text = _read(path)
        a = [repr(i) for i in scan(text).items()]
        b = [repr(i) for i in scan(text, program_name="phenix.refine").items()]
        assert a == b, path
    print("      (%d logs)" % len(logs))


def test_corpus_is_gui_shape_with_no_wrapper_left():
    """The working corpus is the test input and must carry no evidence of its
    own outcome. 486 lines of agent footer leaked into 25 of 29 error logs as
    first delivered, because the stripping script removed only the lines that
    BEGAN with a wrapper key and not the blocks beneath them.

    REVISION: the first version asserted NO regions at all and failed on two
    logs. Those two are the agent's OWN logs, which quote their children's
    wrapper blocks mid-file -- content, not residue. The invariant that
    actually catches leakage is about the ENDS of the file: no agent header,
    and no footer region in the last quarter.
    """
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    logs = _corpus_logs(CORPUS)
    dirty = []
    for path in logs:
        for region in _scanned(path).regions:
            if region.kind == "agent_header" or region.source == SOURCE_AGENT:  # noqa
                dirty.append((os.path.basename(path), region))
    assert dirty == [], dirty[:5]
    print("      (%d logs, 0 with wrapper residue at either end)"
          % len(logs))


def test_corpus_agent_shape_regions_are_found():
    """Region detection has nothing to detect on the stripped corpus, so it is
    tested against an agent-shape one. Without this the feature would be
    'verified' by a corpus that cannot exercise it."""
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    logs = _corpus_logs(CORPUS_AGENT)
    # counted BY KIND. The first version counted "not a header" as a footer,
    # so when agent_prefix regions were added the footer count jumped from 14
    # to 29 and the test happily passed while reporting a wrong number.
    kinds = {}
    for path in logs:
        for region in scan(_read(path)).regions:
            kinds[region.kind] = kinds.get(region.kind, 0) + 1
    assert kinds.get("agent_header", 0) >= 1, kinds
    assert kinds.get("agent_footer", 0) >= 1, kinds
    print("      (%d logs: %s)" % (
        len(logs), ", ".join("%s=%d" % kv for kv in sorted(kinds.items()))))


def test_agent_region_marks_the_items_inside_it():
    """FOUND BY A NEGATIVE CONTROL, not by review: control C4 forced every
    item to source=program and the whole suite still passed, which means the
    marking was implemented and never tested. find_regions was covered; the
    propagation of a region's source onto the items inside it was not."""
    text = ("WORKING DIRECTORY: /tmp/x\n"
            "COMMAND THAT WAS RUN: phenix.xtriage data.sca opt=1\n"
            "LOG TEXT: Starting phenix.xtriage\n"
            "  Found phil, /tmp/x.eff\n"
            + "*" * 50 + "\n"
            "FINAL QUALITY METRICS REPORT:\n"
            "Space Group: I 4 (No. 79)\n")
    structure = scan(text)
    by_line = dict((item.line, item.source) for item in structure.unparsed)
    assert by_line.get(2) == SOURCE_AGENT, by_line   # the command line
    assert by_line.get(4) == SOURCE_PROGRAM, by_line  # the program's own line
    # line 7 sits inside the footer; since P8 it is an agent MEASUREMENT
    # rather than an unparsed record, so check the item stream instead
    inside_footer = [i for i in structure.items() if i.line == 7]
    assert inside_footer, structure.items()
    assert all(i.source == SOURCE_AGENT for i in inside_footer), inside_footer
    # REVISION (P1): this fixture used `Resolution: 2.50` on line 4, which
    # sits above the footer's `****` rule and so became a section title once
    # P1 landed. Replaced with a line that is not adjacent to a rule.


def test_corpus_agent_shape_items_inside_wrappers_are_marked():
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    agent = program = 0
    for path in _corpus_logs(CORPUS_AGENT):  # noqa
        structure = scan(_read(path))
        for item in structure.items():
            if item.source == SOURCE_AGENT:
                agent += 1
            elif item.source == SOURCE_PROGRAM:
                program += 1
    assert agent > 0, "no item was attributed to the agent in an agent corpus"
    assert program > agent, (program, agent)
    print("      (%d agent-sourced items, %d program-sourced)"
          % (agent, program))


def test_corpus_paired_shapes_agree_once_agent_items_are_removed():
    """Requirements 4.8: an extractor that only works on agent-shape logs is
    the wrong tool. The wave-1 conformance suite tests this as a raw count
    equality, which P1 breaks -- correctly, because the agent's own wrapper
    now yields one section per wrapped log. MEASURED on all 14 wrapped pairs:
    the delta is exactly -1 and in every case the extra section is
    source='agent'. Excluding agent-sourced items, as 4.8 itself demands,
    the shapes agree exactly."""
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    gui_root = corpus_paths.agent_gui_corpus() or (
        CORPUS_AGENT.rstrip("/") + "_gui")
    if not os.path.isdir(gui_root):
        print("      SKIP (no paired GUI-shape corpus beside %s)" % CORPUS_AGENT)
        return
    pairs = differed = 0
    for path in _corpus_logs(CORPUS_AGENT):
        twin = os.path.join(gui_root, os.path.basename(path))
        if not os.path.exists(twin):
            continue
        pairs += 1
        a = [(i.kind, i.evidence()) for i in scan(_read(path)).items()
             if i.source != SOURCE_AGENT]
        b = [(i.kind, i.evidence()) for i in scan(_read(twin)).items()
             if i.source != SOURCE_AGENT]
        # SEMANTIC IDENTITY, not equal counts. Two logs can each yield 17
        # items and disagree about all 17. Raised in review; the weaker
        # assertion had been described as semantic in the architecture
        # document and in requirements 4.8, and was not.
        if a != b:
            differed += 1
            first = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y),
                         min(len(a), len(b)))
            raise AssertionError(
                "%s: agent and GUI shapes disagree at item %d: %r vs %r"
                % (os.path.basename(path), first,
                   a[first:first + 1], b[first:first + 1]))
    assert pairs >= 5, pairs
    print("      (%d pairs, semantically identical)" % pairs)


def test_P1_a_section_does_not_cross_a_provenance_boundary():
    """FOUND IN P1 REVIEW. A section runs to the next heading or to EOF, so
    the last program-written section of a wrapped log swallowed the agent's
    footer. molprobity's `Summary` section then "contained" the agent's
    FINAL QUALITY METRICS REPORT, and a consumer grouping by section would
    attribute the agent's numbers to the program."""
    text = ("Summary\n" + "-" * 40 + "\nClashscore 3.56\n"
            + "*" * 50 + "\nFINAL QUALITY METRICS REPORT:\n"
            "Resolution: 2.50\n")
    structure = scan(text)
    assert structure.sections, structure.sections
    section = structure.sections[0]
    footer = [r for r in structure.regions if r.kind == "agent_footer"]
    assert footer, structure.regions
    assert section.end < footer[0].start, (section, footer[0])


def test_no_corpus_test_reports_that_it_checked_nothing():
    """A test that verifies nothing and prints PASS is the failure this
    project keeps meeting -- the vacuous probe in P4, the flat-log audit
    whose sample definition drifted in P9, and two smoke tests that reported
    "0 summary rows" and "0 attached measurements" on Tom's first run.

    This walks the corpus tests' own reported counts and fails if any of the
    ones that MUST have material found none. It is cheap insurance against
    the same mistake in a tier whose whole point is to be small."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    roots = [r for r in (CORPUS, CORPUS_AGENT) if r]
    paths = [p for r in roots for p in _corpus_logs(r)]
    material = dict(stages=0, cycles=0, decisions=0, exclusions=0,
                    control_skips=0, completion_records=0, measurements=0,
                    labeled_values=0, sections=0, phases=0)
    for path in paths:
        structure = _scanned(path)
        for kind in material:
            material[kind] += len(getattr(structure, kind))
    empty = sorted(k for k, v in material.items() if v == 0)
    assert not empty, (
        "the corpus available to the smoke tier contains NO %s, so every"
        " invariant about them passes vacuously" % ", ".join(empty))
    print("      (every channel has material: %s)"
          % ", ".join("%s=%d" % kv for kv in sorted(material.items())))


def test_corpus_P1_no_section_spans_a_provenance_boundary():
    if not CORPUS_AGENT:
        SKIPPED.append("PHENIX_LOG_CORPUS_AGENT unset")
        print("      SKIP (PHENIX_LOG_CORPUS_AGENT unset)")
        return
    crossings = []
    for path in _corpus_logs(CORPUS_AGENT):
        structure = _scanned(path)
        for section in structure.sections:
            for region in structure.regions:
                if region.source == section.source:
                    continue
                if region.contains(section.end) != region.contains(section.line):
                    crossings.append((os.path.basename(path), section.title))
    assert crossings == [], crossings[:5]
    print("      (0 crossings)")


def test_corpus_optimised_paths_agree_with_naive_reference():
    """An optimisation that changes behaviour is a defect wearing a
    performance argument. Both P1 optimisations -- the early-exit numeric
    scan and the cursor-based section lookup -- are checked against naive
    reference implementations over the whole corpus, not argued about."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    import log_extraction.log_structure_extractor as module

    def naive_numeric_row(text, minimum=4):
        tokens = text.split()
        if len(tokens) < minimum:
            return False
        numeric = sum(1 for t in tokens
                      if module._NUMBER_RE.match(t.rstrip(":,")))
        return numeric >= minimum

    lines_checked = sections_checked = 0
    for path in _corpus_logs(CORPUS):
        text = _read(path)
        lines = split_lines(text)
        for raw in lines:
            assert module._is_numeric_row(raw) == naive_numeric_row(raw), raw
            lines_checked += 1
        structure = _scanned(path)
        bounds = [(x.line, x.end, x.section_id) for x in structure.sections]
        for item in structure.unparsed:
            naive = None
            for start, end, ident in bounds:
                if start <= item.line <= end:
                    naive = ident
                    break
            assert item.section_id == naive, (path, item, naive)
            sections_checked += 1
    print("      (%d lines, %d section lookups)"
          % (lines_checked, sections_checked))


def test_corpus_screen_rules_all_fire():
    """A screen rule that never fires is either dead code or a wrong pattern.
    Either way it should not sit in a frozen screen unnoticed."""
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    fired = {}
    for path in _corpus_logs(CORPUS):
        for item in _scanned(path).unparsed:
            fired[item.screen_rule] = fired.get(item.screen_rule, 0) + 1
    missing = [r for r in SCREEN_RULES if r not in fired]
    assert not missing, ("screen rules that never fired: %s" % missing)
    print("      (%s)" % ", ".join(
        "%s=%d" % (r, fired[r]) for r in SCREEN_RULES))


def test_order_is_preserved_within_each_kind():
    """Requirements 3: order is preserved -- stage and cycle sequence IS the
    finding. Never tested until the review noticed nothing asserted it."""
    structure = scan("Aa: 1\nBb: 2\nCc: 3\n")
    assert [x.line for x in structure.labeled_values] == [1, 2, 3], \
        structure.labeled_values
    text = ("One:\n" + "-" * 40 + "\nbody\nTwo:\n" + "-" * 40 + "\n"
            "Three:\n" + "-" * 40 + "\n")
    assert [x.title for x in scan(text).sections] == ["One:", "Two:",
                                                     "Three:"]



def test_corpus_purity_at_scale():
    if not CORPUS:
        SKIPPED.append("PHENIX_LOG_CORPUS unset")
        print("      SKIP (PHENIX_LOG_CORPUS unset)")
        return
    logs = _corpus_logs(CORPUS)
    for path in logs:          # the smoke set is small enough to take whole
        text = _read(path)
        assert [repr(i) for i in scan(text).items()] == \
               [repr(i) for i in scan(text).items()], path
    print("      (%d logs)" % len(logs))


SKIPPED = []


def run_all_tests():
    """House protocol: cctbx-style FAIL-FAST.

    The first assertion failure raises an uncaught exception with a full
    traceback rather than collecting failures, which is what
    tests/run_all_tests.py expects: it calls this, catches whatever comes
    out, and marks the module FAILED.

    Corpus tests SKIP when PHENIX_LOG_CORPUS is unset and FAIL LOUDLY when it
    is set but wrong -- a skip that reads as a pass has cost this project
    twice.
    """
    del SKIPPED[:]
    run_tests_with_fail_fast()
    _report_skips()
    _report_scope()


def _report_scope():
    """Tier reporting. NOT a skip and NOT a failure -- explicit scope.

    The in-tree suite is Tier 1 (fixtures) + Tier 2 (a 1.25 MB smoke corpus).
    Tier 3, the full-population claims, needs the 21 MB validation corpus
    which is deliberately not in cctbx_project. Saying so every run is how it
    stays hard to forget; a suite that silently omits a third of its scope is
    the false-green this project has been bitten by three times."""
    print("")
    print("  Tier 1 (fixtures) + Tier 2 (smoke corpus, %s): RUN" %
          ("present" if CORPUS else "ABSENT"))
    print("  Tier 3 (full 21 MB corpus, population claims): NOT REQUESTED")
    print("      phenix.python langchain/log_extraction/validate.py")


def _report_skips():
    """A SKIP THAT READS AS A PASS IS THE FAILURE THIS PROJECT KEEPS MAKING.

    Without a corpus this module prints "All 130 tests passed" and the house
    runner shows a green tick -- while 34 corpus-level invariants, which are
    most of the value, did not run. Twice already a suite has reported clean
    while silently skipping (tst_conformance vs the prototype's own suite
    reading different env vars; the ligandfit suite counting skips as passes).
    So say so, loudly, every time."""
    if not SKIPPED:
        return
    print("")
    print("!" * 68)
    print("!! %d of the tests above did NOT RUN -- they were SKIPPED." % len(SKIPPED))
    print("!! These are the CORPUS-LEVEL INVARIANTS, and they are most of")
    print("!! the value: fixture-only tests encode the same assumptions as")
    print("!! the code. A green result here does NOT mean the extractor was")
    print("!! checked against real logs.")
    print("!!")
    print("!! The corpus normally ships INSIDE the package and is found")
    print("!! automatically. It is missing or in the wrong place:")
    print("!!")
    for line in corpus_paths.describe().splitlines():
        print("!! %s" % line)
    print("!!")
    print("!! Expected at: %s" % corpus_paths.DEFAULT_ROOT)
    print("!" * 68)


def run_all_tests_counting():
    """Non-fail-fast mode, for the NEGATIVE CONTROL runner only.

    A control result is the DIFFERENCE a disabled feature makes, so the
    magnitude matters: "C1 disabled -> 24 passed, 10 failed" says more than
    "failed". Fail-fast would stop at the first one and throw that away.

    Not the house protocol, and not what run_all_tests.py calls. Reached with
    `tst_log_structure_extractor.py --count`.
    """
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = errors = 0
    for test in tests:
        try:
            test()
        except AssertionError as error:
            failed += 1
            print("  FAILED: %s" % (error,))
        except Exception as error:                          # noqa: BLE001
            failed += 1
            errors += 1
            print("  ERROR: %s: %s" % (type(error).__name__, error))
    _report_skips()
    note = ""
    if errors:
        note = " (%d were errors, not assertion failures)" % errors
    print("\n%d passed, %d failed, %d total%s"
          % (len(tests) - failed, failed, len(tests), note))
    return failed


if __name__ == "__main__":
    if "--count" in sys.argv:
        sys.exit(1 if run_all_tests_counting() else 0)
    run_all_tests()
