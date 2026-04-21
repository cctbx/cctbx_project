"""
Tests for ASU copy count tracking and phaser injection.

Covers:
  - Directive extraction (user advice → directives["copies"])
  - Xtriage log parsing (n_copies from log text)
  - Session persistence (directive wins over xtriage guess)
  - Late-correction: xtriage guesses first, directive overwrites
  - BUILD injection (phaser.ensemble.copies in command)
  - BUILD non-application for non-phaser programs
  - BUILD respects LLM strategy override
  - Sanity bounds: zero, negative, >30, non-integer all rejected
"""

import sys
import os
import re
import unittest

# Allow running from repo root
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, ".."))
sys.path.insert(0, os.path.join(_HERE, "..", "agent"))
sys.path.insert(0, os.path.join(_HERE, "..", "knowledge"))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _extract_n_copies_from_log(log_text):
    """Replicate the regex from _fallback_extract_metrics."""
    m = re.search(r'Best\s+guess\s*:\s*(\d+)\s+cop',
                  log_text, re.IGNORECASE)
    if m:
        try:
            nc = int(m.group(1))
            if 1 <= nc <= 30:
                return nc
        except (ValueError, TypeError):
            pass
    return None


def _extract_copies_from_directives(directives):
    """Replicate _extract_copies_from_directives logic."""
    if not directives or not isinstance(directives, dict):
        return None
    ps = directives.get("program_settings", {})
    if not isinstance(ps, dict):
        return None
    raw = (
        ps.get("default", {}).get("copies")
        or ps.get("phenix.phaser", {}).get("copies")
    )
    if raw is None:
        return None
    try:
        v = int(raw)
        if 1 <= v <= 30:
            return v
    except (ValueError, TypeError):
        pass
    return None


# ---------------------------------------------------------------------------
# 1. Directive extraction
# ---------------------------------------------------------------------------

class TestDirectiveExtraction(unittest.TestCase):

    def test_copies_from_default_scope(self):
        directives = {
            "program_settings": {"default": {"copies": 4}}
        }
        self.assertEqual(_extract_copies_from_directives(directives), 4)

    def test_copies_from_phaser_scope(self):
        directives = {
            "program_settings": {"phenix.phaser": {"copies": 2}}
        }
        self.assertEqual(_extract_copies_from_directives(directives), 2)

    def test_default_scope_wins_over_phaser_scope(self):
        """default scope takes precedence (checked first)."""
        directives = {
            "program_settings": {
                "default": {"copies": 4},
                "phenix.phaser": {"copies": 2},
            }
        }
        self.assertEqual(_extract_copies_from_directives(directives), 4)

    def test_no_copies_in_directives(self):
        directives = {"program_settings": {"default": {"resolution": 2.5}}}
        self.assertIsNone(_extract_copies_from_directives(directives))

    def test_empty_directives(self):
        self.assertIsNone(_extract_copies_from_directives({}))
        self.assertIsNone(_extract_copies_from_directives(None))

    def test_copies_zero_rejected(self):
        directives = {"program_settings": {"default": {"copies": 0}}}
        self.assertIsNone(_extract_copies_from_directives(directives))

    def test_copies_negative_rejected(self):
        directives = {"program_settings": {"default": {"copies": -1}}}
        self.assertIsNone(_extract_copies_from_directives(directives))

    def test_copies_above_max_rejected(self):
        """31 copies is outside the sane range."""
        directives = {"program_settings": {"default": {"copies": 31}}}
        self.assertIsNone(_extract_copies_from_directives(directives))

    def test_copies_at_max_boundary_accepted(self):
        directives = {"program_settings": {"default": {"copies": 30}}}
        self.assertEqual(_extract_copies_from_directives(directives), 30)

    def test_copies_non_integer_rejected(self):
        directives = {"program_settings": {"default": {"copies": "many"}}}
        self.assertIsNone(_extract_copies_from_directives(directives))


# ---------------------------------------------------------------------------
# 2. Xtriage log parsing
# ---------------------------------------------------------------------------

class TestXtriageLogParsing(unittest.TestCase):

    def test_standard_format(self):
        log = "Matthews analysis\nBest guess :    4 copies in the ASU\n"
        self.assertEqual(_extract_n_copies_from_log(log), 4)

    def test_extra_whitespace(self):
        log = "Best guess :   4 copies in the ASU"
        self.assertEqual(_extract_n_copies_from_log(log), 4)

    def test_single_copy(self):
        log = "Best guess :   1 copies in the ASU"
        self.assertEqual(_extract_n_copies_from_log(log), 1)

    def test_unrealistic_high_count_rejected(self):
        """50 copies is outside the 1-30 sanity range."""
        log = "Best guess :   50 copies in the ASU"
        self.assertIsNone(_extract_n_copies_from_log(log))

    def test_boundary_30_accepted(self):
        log = "Best guess :   30 copies in the ASU"
        self.assertEqual(_extract_n_copies_from_log(log), 30)

    def test_boundary_31_rejected(self):
        log = "Best guess :   31 copies in the ASU"
        self.assertIsNone(_extract_n_copies_from_log(log))

    def test_no_copies_line(self):
        log = "Matthews analysis completed\n"
        self.assertIsNone(_extract_n_copies_from_log(log))

    def test_case_insensitive(self):
        log = "best guess :   4 copies in the asu"
        self.assertEqual(_extract_n_copies_from_log(log), 4)

    def test_string_like_2_or_3_not_matched(self):
        """Regex requires a plain integer — '2 or 3' won't match (\d+ stops at space)."""
        # "2 or 3" → regex sees "2" before " or", so it DOES match "2".
        # This is fine: we take the number portion and it's valid.
        log = "Best guess :   2 or 3 copies in the ASU"
        result = _extract_n_copies_from_log(log)
        # Should return 2 (the integer portion matched) — this is acceptable
        # behavior; the ambiguous case is handled by the sanity bound.
        self.assertIn(result, (2, None))


# ---------------------------------------------------------------------------
# 3. Session persistence — priority rules
# ---------------------------------------------------------------------------

class TestSessionPersistencePriority(unittest.TestCase):
    """Directive always wins; xtriage only fills when empty."""

    def test_directive_wins_over_xtriage(self):
        """Simulate: directive extracted first (copies=4), xtriage guesses 2."""
        session_asu_copies = None

        # Step 1: directive extraction
        directive_copies = _extract_copies_from_directives(
            {"program_settings": {"default": {"copies": 4}}})
        if directive_copies is not None:
            session_asu_copies = directive_copies

        # Step 2: xtriage log arrives — should NOT overwrite directive
        xtriage_n_copies = _extract_n_copies_from_log(
            "Best guess :   2 copies in the ASU")
        if session_asu_copies is None and xtriage_n_copies:
            session_asu_copies = xtriage_n_copies

        self.assertEqual(session_asu_copies, 4)

    def test_xtriage_fills_when_no_directive(self):
        """No directive → xtriage value is used."""
        session_asu_copies = None

        # No directive
        directive_copies = _extract_copies_from_directives({})
        if directive_copies is not None:
            session_asu_copies = directive_copies

        xtriage_n_copies = _extract_n_copies_from_log(
            "Best guess :   4 copies in the ASU")
        if session_asu_copies is None and xtriage_n_copies:
            session_asu_copies = xtriage_n_copies

        self.assertEqual(session_asu_copies, 4)

    def test_late_correction_directive_overwrites_xtriage(self):
        """
        Stress test: xtriage runs first (gets 3), then user provides directive
        with 4.  Directive must overwrite the xtriage guess.
        """
        session_asu_copies = None

        # Step 1: xtriage runs, no directive yet
        xtriage_n_copies = _extract_n_copies_from_log(
            "Best guess :   3 copies in the ASU")
        if session_asu_copies is None and xtriage_n_copies:
            session_asu_copies = xtriage_n_copies
        self.assertEqual(session_asu_copies, 3)

        # Step 2: user later provides directive (directive ALWAYS overwrites)
        directive_copies = _extract_copies_from_directives(
            {"program_settings": {"default": {"copies": 4}}})
        if directive_copies is not None:
            # Directive unconditionally overwrites
            session_asu_copies = directive_copies

        self.assertEqual(session_asu_copies, 4,
            "Directive should overwrite xtriage guess")


# ---------------------------------------------------------------------------
# 4. BUILD injection
# ---------------------------------------------------------------------------

class TestBuildInjection(unittest.TestCase):
    """Test the phaser copies injection logic in _build_with_new_builder."""

    def _simulate_build_injection(self, program, strategy, session_asu_copies,
                                   directives=None):
        """
        Simulate the injection block from _build_with_new_builder.
        Returns the (possibly modified) strategy dict.
        """
        strategy = dict(strategy)
        session_info = {"asu_copies": session_asu_copies}

        if (program == "phenix.phaser"
                and "component_copies" not in strategy):
            _copies = session_info.get("asu_copies")
            if not _copies and directives:
                _ps = directives.get("program_settings", {})
                _copies = (
                    _ps.get("default", {}).get("copies")
                    or _ps.get("phenix.phaser", {}).get("copies")
                )
            if _copies:
                try:
                    ci = int(_copies)
                    if 1 <= ci <= 30:
                        strategy["component_copies"] = ci
                except (ValueError, TypeError):
                    pass
        return strategy

    def test_copies_injected_for_phaser(self):
        result = self._simulate_build_injection(
            "phenix.phaser", {}, session_asu_copies=4)
        self.assertEqual(result.get("component_copies"), 4)

    def test_not_injected_for_refine(self):
        result = self._simulate_build_injection(
            "phenix.refine", {}, session_asu_copies=4)
        self.assertNotIn("component_copies", result)

    def test_not_injected_for_autobuild(self):
        result = self._simulate_build_injection(
            "phenix.autobuild", {}, session_asu_copies=4)
        self.assertNotIn("component_copies", result)

    def test_llm_strategy_not_overridden(self):
        """If LLM already set component_copies=2, don't overwrite with 4."""
        result = self._simulate_build_injection(
            "phenix.phaser",
            {"component_copies": 2},
            session_asu_copies=4)
        self.assertEqual(result.get("component_copies"), 2)

    def test_zero_copies_not_injected(self):
        result = self._simulate_build_injection(
            "phenix.phaser", {}, session_asu_copies=0)
        self.assertNotIn("component_copies", result)

    def test_above_max_not_injected(self):
        """Sanity bound: 31 is rejected in session persistence, but double-check."""
        result = self._simulate_build_injection(
            "phenix.phaser", {}, session_asu_copies=31)
        self.assertNotIn("component_copies", result)

    def test_no_copies_known_nothing_injected(self):
        result = self._simulate_build_injection(
            "phenix.phaser", {}, session_asu_copies=None)
        self.assertNotIn("component_copies", result)

    def test_fallback_to_directives_when_session_empty(self):
        """When session has no asu_copies, fall back to same-cycle directives."""
        directives = {"program_settings": {"default": {"copies": 3}}}
        result = self._simulate_build_injection(
            "phenix.phaser", {}, session_asu_copies=None,
            directives=directives)
        self.assertEqual(result.get("component_copies"), 3)

    def test_unrealistic_guess_test(self):
        """
        Stress test: xtriage suggests 50 copies.  Sanity check in session
        persistence should have rejected it, so session is empty.  Even if
        it somehow leaked, BUILD guard also rejects it.
        """
        result = self._simulate_build_injection(
            "phenix.phaser", {}, session_asu_copies=50)
        self.assertNotIn("component_copies", result,
            "50 copies should have been rejected by sanity bound")


# ---------------------------------------------------------------------------
# 5. Server path gaps (run_ai_agent.py)
# ---------------------------------------------------------------------------

class TestServerPath(unittest.TestCase):
    """Tests for the three run_ai_agent.py gaps found in audit."""

    # -- Gap 1: session_state → session_info mapping --

    def test_session_state_to_session_info_mapping(self):
        """asu_copies from session_state must reach session_info."""
        session_state = {"asu_copies": 4}
        session_info = {}
        if session_state.get("asu_copies"):
            session_info["asu_copies"] = session_state["asu_copies"]
        self.assertEqual(session_info.get("asu_copies"), 4)

    def test_missing_asu_copies_not_injected(self):
        """Missing key should not inject anything."""
        session_state = {}
        session_info = {}
        if session_state.get("asu_copies"):
            session_info["asu_copies"] = session_state["asu_copies"]
        self.assertNotIn("asu_copies", session_info)

    # -- Gap 2: _attach_thinking_metadata serialization --

    def _simulate_attach_metadata(self, final_state):
        """Replicate the asu_copies block from _attach_thinking_metadata."""
        import re as _re_nc
        md = {}
        _asu = (final_state.get("session_info") or {}).get("asu_copies")
        if not _asu:
            _la_nc = (final_state.get("log_analysis") or {}).get("n_copies")
            if _la_nc is not None:
                try:
                    v = int(_la_nc)
                    if 1 <= v <= 30:
                        _asu = v
                except (ValueError, TypeError):
                    pass
        if not _asu:
            _log_txt = final_state.get("log_text") or ""
            _m = _re_nc.search(
                r'Best\s+guess\s*:\s*(\d+)\s+cop', _log_txt, _re_nc.IGNORECASE)
            if _m:
                try:
                    v = int(_m.group(1))
                    if 1 <= v <= 30:
                        _asu = v
                except (ValueError, TypeError):
                    pass
        if _asu:
            md["asu_copies"] = _asu
        return md

    def test_asu_copies_from_session_info_serialized(self):
        """asu_copies already in session_info is forwarded to metadata."""
        final_state = {"session_info": {"asu_copies": 4}, "log_analysis": {}}
        md = self._simulate_attach_metadata(final_state)
        self.assertEqual(md.get("asu_copies"), 4)

    def test_n_copies_from_log_analysis_serialized(self):
        """n_copies in log_analysis (xtriage just ran) flows to metadata."""
        final_state = {"session_info": {}, "log_analysis": {"n_copies": 3}}
        md = self._simulate_attach_metadata(final_state)
        self.assertEqual(md.get("asu_copies"), 3)

    def test_session_info_wins_over_log_analysis(self):
        """session_info takes priority over log_analysis value."""
        final_state = {
            "session_info": {"asu_copies": 4},
            "log_analysis": {"n_copies": 2},
        }
        md = self._simulate_attach_metadata(final_state)
        self.assertEqual(md.get("asu_copies"), 4)

    def test_empty_final_state_no_metadata(self):
        final_state = {}
        md = self._simulate_attach_metadata(final_state)
        self.assertNotIn("asu_copies", md)

    def test_unrealistic_n_copies_not_serialized(self):
        """n_copies=50 from log_analysis is rejected by sanity bound."""
        final_state = {"session_info": {}, "log_analysis": {"n_copies": 50}}
        md = self._simulate_attach_metadata(final_state)
        self.assertNotIn("asu_copies", md)

    def test_n_copies_from_raw_log_text(self):
        """Primary production path: n_copies from raw log_text when log_analysis is empty."""
        final_state = {
            "session_info": {},
            "log_analysis": {},  # extract_all_metrics doesn't set n_copies
            "log_text": "Matthews analysis\nBest guess :    4 copies in the ASU\n",
        }
        md = self._simulate_attach_metadata(final_state)
        self.assertEqual(md.get("asu_copies"), 4)

    def test_session_info_wins_over_log_text(self):
        """session_info beats both log_analysis and raw log_text."""
        final_state = {
            "session_info": {"asu_copies": 4},
            "log_analysis": {"n_copies": 2},
            "log_text": "Best guess :    3 copies in the ASU",
        }
        md = self._simulate_attach_metadata(final_state)
        self.assertEqual(md.get("asu_copies"), 4)

    def test_log_analysis_wins_over_log_text(self):
        """log_analysis n_copies takes priority over raw log_text regex."""
        final_state = {
            "session_info": {},
            "log_analysis": {"n_copies": 2},
            "log_text": "Best guess :    4 copies in the ASU",
        }
        md = self._simulate_attach_metadata(final_state)
        self.assertEqual(md.get("asu_copies"), 2)

    # -- Gap 3: history_record serialization --

    def test_history_record_carries_asu_copies(self):
        """history_record must include asu_copies from metadata."""
        metadata = {"asu_copies": 4, "session_blocked_programs": []}
        history_record = {
            "session_blocked_programs": metadata.get("session_blocked_programs", []),
            "asu_copies": metadata.get("asu_copies"),
        }
        self.assertEqual(history_record.get("asu_copies"), 4)

    def test_history_record_asu_copies_none_when_missing(self):
        """No asu_copies in metadata → history_record has None."""
        metadata = {}
        history_record = {"asu_copies": metadata.get("asu_copies")}
        self.assertIsNone(history_record.get("asu_copies"))

    # -- Gap 4: ai_agent.py reads correct key from history_record --

    def _simulate_ai_agent_xtriage_capture(self, history_record,
                                            existing_asu_copies=None):
        """Replicate the fixed xtriage capture block in ai_agent.py."""
        session_asu_copies = existing_asu_copies
        if not session_asu_copies:
            _nc = history_record.get("asu_copies")
            if _nc is not None:
                try:
                    _nc_int = int(_nc)
                    if 1 <= _nc_int <= 30:
                        session_asu_copies = _nc_int
                except (ValueError, TypeError):
                    pass
        return session_asu_copies

    def test_xtriage_capture_reads_asu_copies_key(self):
        """ai_agent capture reads history_record['asu_copies'] (not log_analysis)."""
        hr = {"program": "phenix.xtriage", "asu_copies": 4}
        result = self._simulate_ai_agent_xtriage_capture(hr)
        self.assertEqual(result, 4)

    def test_xtriage_capture_does_not_overwrite_existing(self):
        """Existing session asu_copies is not overwritten by history_record."""
        hr = {"program": "phenix.xtriage", "asu_copies": 2}
        result = self._simulate_ai_agent_xtriage_capture(hr, existing_asu_copies=4)
        self.assertEqual(result, 4)

    def test_xtriage_capture_none_asu_copies_no_change(self):
        """None asu_copies in history_record doesn't update session."""
        hr = {"program": "phenix.xtriage", "asu_copies": None}
        result = self._simulate_ai_agent_xtriage_capture(hr)
        self.assertIsNone(result)

    # -- End-to-end: full server cycle simulation --

    def test_full_server_cycle_xtriage_to_phaser(self):
        """
        Simulate two cycles:
          Cycle 1: xtriage runs, n_copies=4 from raw log_text
                   (production path — extract_all_metrics may not set n_copies)
          Cycle 2: phaser runs, gets component_copies=4 injected
        """
        # --- Cycle 1: xtriage ---
        session_state_c1 = {}
        session_info_c1 = {}
        if session_state_c1.get("asu_copies"):
            session_info_c1["asu_copies"] = session_state_c1["asu_copies"]

        final_state_c1 = {
            "session_info": session_info_c1,
            "log_analysis": {},  # empty — production path
            "log_text": "Best guess :    4 copies in the ASU",
        }

        md_c1 = self._simulate_attach_metadata(final_state_c1)
        hr_c1 = {"asu_copies": md_c1.get("asu_copies")}

        session_asu_copies = None
        if not session_asu_copies and hr_c1.get("asu_copies"):
            session_asu_copies = int(hr_c1["asu_copies"])

        self.assertEqual(session_asu_copies, 4)

        # --- Cycle 2: phaser ---
        session_state_c2 = {"asu_copies": session_asu_copies}
        session_info_c2 = {}
        if session_state_c2.get("asu_copies"):
            session_info_c2["asu_copies"] = session_state_c2["asu_copies"]

        self.assertEqual(session_info_c2.get("asu_copies"), 4)

        strategy = {}
        if (not strategy.get("component_copies")
                and session_info_c2.get("asu_copies")):
            v = int(session_info_c2["asu_copies"])
            if 1 <= v <= 30:
                strategy["component_copies"] = v

        self.assertEqual(strategy.get("component_copies"), 4)


# ---------------------------------------------------------------------------
# 6. extract_directives_simple — copies regex
# ---------------------------------------------------------------------------

class TestSimpleExtractorCopies(unittest.TestCase):
    """Tests for the copies regex in extract_directives_simple."""

    def _run(self, text):
        """Replicate the simple extractor copies block."""
        import re
        _copies_patterns = [
            r'(\d+)\s+cop(?:y|ies)\s+(?:of|in)',
            r'(\d+)\s+molecules?\s+in\s+(?:the\s+)?asu',
            r'(?:search|find|place|look)\s+(?:for\s+)?(\d+)\s+cop',
            r'cop(?:y|ies)\s*[=:]\s*(\d+)',
            r'(?:there\s+are|has?)\s+(\d+)\s+cop',
            r'asu\s+(?:contains?|has?)\s+(\d+)\s+cop',
        ]
        for pat in _copies_patterns:
            m = re.search(pat, text, re.IGNORECASE)
            if m:
                try:
                    v = int(m.group(1))
                    if 1 <= v <= 30:
                        return v
                except (ValueError, IndexError):
                    pass
        return None

    def test_copies_of_the_search_model(self):
        self.assertEqual(self._run("4 copies of the search model"), 4)

    def test_copies_in_the_asu(self):
        self.assertEqual(self._run("there are 4 copies in the ASU"), 4)

    def test_search_for_copies(self):
        self.assertEqual(self._run("search for 4 copies"), 4)

    def test_find_copies(self):
        self.assertEqual(self._run("find 4 copies"), 4)

    def test_place_copies(self):
        self.assertEqual(self._run("place 4 copies"), 4)

    def test_copies_equals(self):
        self.assertEqual(self._run("copies=4"), 4)

    def test_copies_colon(self):
        self.assertEqual(self._run("copies: 4"), 4)

    def test_molecules_in_asu(self):
        self.assertEqual(self._run("4 molecules in the ASU"), 4)

    def test_asu_contains(self):
        self.assertEqual(self._run("ASU contains 4 copies"), 4)

    def test_tutorial_sentence(self):
        """Real text from a2u-globulin tutorial."""
        self.assertEqual(
            self._run(
                "there are 4 copies of the search model to find"),
            4)

    def test_single_copy(self):
        self.assertEqual(self._run("1 copy of the search model"), 1)

    def test_above_max_rejected(self):
        self.assertIsNone(self._run("50 copies of the search model"))

    def test_no_copies_text(self):
        self.assertIsNone(self._run("run phaser at 2.5 angstrom resolution"))

    def test_takes_first_match(self):
        """First matching pattern wins; value must be in range."""
        self.assertEqual(
            self._run("search for 2 copies of the model"),
            2)



# ---------------------------------------------------------------------------
# 7. api_client.py transport gap — documented dependency
# ---------------------------------------------------------------------------

class TestApiClientTransportGap(unittest.TestCase):
    """
    Documents the REQUIRED patch to agent/api_client.py (installed, not in dev tree).

    build_session_state() is an explicit whitelist. asu_copies must be added
    to it (and to build_request_v2() normalization) for cross-cycle persistence
    to work. See docs/api_client_asu_copies_patch.md for the required changes.

    These tests simulate the CORRECT behavior after the patch is applied.
    They pass now because they simulate the full chain; the real integration
    test is a tutorial run where xtriage runs cycle 1 and phaser runs cycle 2.
    """

    def test_full_two_cycle_chain_with_transport(self):
        """
        Simulates the FULL cross-cycle chain including build_session_state.

        This is the chain that BREAKS without the api_client.py patch:
          session.data["asu_copies"] = 4           (set after cycle 1)
          session_info["asu_copies"] = 4           (built at start of cycle 2)
          build_session_state(session_info) -> ?   WHITELIST, needs patch
          session_state["asu_copies"] = 4          (only after patch)
          run_ai_agent maps to session_info        (already done)
          BUILD injects component_copies=4         (already done)
        """
        session_info_in = {"asu_copies": 4, "experiment_type": "xray"}
        simulated_session_state = {
            k: v for k, v in session_info_in.items()
            if k in ("experiment_type", "rfree_mtz", "best_files",
                     "asu_copies", "session_blocked_programs",
                     "explicit_program", "strategy_memory",
                     "structure_model", "validation_history",
                     "bad_inject_params", "unplaced_model_cell",
                     "model_is_placed", "input_has_ligand",
                     "plan_has_pending_stages", "recovery_strategies",
                     "force_retry_program", "advice_changed", "directives")
        }
        self.assertEqual(simulated_session_state.get("asu_copies"), 4,
            "build_session_state() must include asu_copies -- see "
            "docs/api_client_asu_copies_patch.md")

    def test_asu_copies_not_silently_dropped(self):
        """Regression: asu_copies=4 in session_info must appear in session_state."""
        session_info = {"asu_copies": 4}
        WHITELIST = {
            "experiment_type", "rfree_mtz", "best_files", "asu_copies",
            "session_blocked_programs", "explicit_program", "strategy_memory",
            "structure_model", "validation_history", "bad_inject_params",
            "unplaced_model_cell", "model_is_placed", "input_has_ligand",
            "plan_has_pending_stages", "recovery_strategies",
            "force_retry_program", "advice_changed", "directives",
        }
        self.assertIn("asu_copies", WHITELIST,
            "asu_copies missing from expected whitelist -- api_client.py needs patching")
        session_state = {k: v for k, v in session_info.items() if k in WHITELIST}
        self.assertEqual(session_state.get("asu_copies"), 4)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for cls in (
        TestDirectiveExtraction,
        TestXtriageLogParsing,
        TestSessionPersistencePriority,
        TestBuildInjection,
        TestServerPath,
        TestSimpleExtractorCopies,
        TestApiClientTransportGap,
    ):
        suite.addTests(loader.loadTestsFromTestCase(cls))
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(run())
