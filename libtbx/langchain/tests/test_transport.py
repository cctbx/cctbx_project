"""
Unit tests for the transport module.

Run with: python -m pytest tests/test_transport.py -v
"""

from __future__ import absolute_import, division, print_function

import unittest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.transport import (
    remove_markers,
    truncate_quoted_strings,
    replace_tabs,
    remove_control_chars,
    truncate_string,
    sanitize_string,
    sanitize_for_transport,
    sanitize_dict_recursive,
    encode_for_rest,
    decode_from_rest,
    get_transport_config,
    get_field_settings,
    clear_config_cache,
    sanitize_field,
    sanitize_request,
    sanitize_response,
    prepare_request_for_transport,
    process_request_from_transport,
    prepare_response_for_transport,
    process_response_from_transport,
    verify_roundtrip,
)


class TestRemoveMarkers(unittest.TestCase):
    """Tests for remove_markers function."""
    
    def test_removes_zztazz(self):
        """ZZTAZZ markers are removed."""
        text = "Hello ZZTAZZ world"
        result = remove_markers(text)
        self.assertEqual(result, "Hello  world")
    
    def test_removes_multiple_markers(self):
        """Multiple different markers are removed."""
        text = "ZZCRZZ line ZZTAZZ tab ZZLBZZ brace"
        result = remove_markers(text)
        self.assertEqual(result, " line  tab  brace")
    
    def test_removes_all_marker_types(self):
        """All ZZxxZZ patterns are removed."""
        markers = ["ZZCRZZ", "ZZLBZZ", "ZZRBZZ", "ZZSCZZ", "ZZHAZZ", 
                   "ZZTAZZ", "ZZDQZZ", "ZZSQZZ", "ZZSRZZ", "ZZEXZZ", "ZZDSZZ"]
        for marker in markers:
            text = f"before {marker} after"
            result = remove_markers(text)
            self.assertEqual(result, "before  after", f"Failed for {marker}")
    
    def test_preserves_non_markers(self):
        """Non-marker text is preserved."""
        text = "Hello world, no markers here!"
        result = remove_markers(text)
        self.assertEqual(result, text)
    
    def test_handles_empty_string(self):
        """Empty string returns empty string."""
        self.assertEqual(remove_markers(""), "")
    
    def test_handles_none(self):
        """None returns empty string."""
        self.assertEqual(remove_markers(None), "")
    
    def test_handles_non_string(self):
        """Non-string is converted to string."""
        self.assertEqual(remove_markers(123), "123")


class TestTruncateQuotedStrings(unittest.TestCase):
    """Tests for truncate_quoted_strings function."""
    
    def test_truncates_long_quoted_string(self):
        """Long quoted strings are truncated."""
        content = "a" * 1000
        text = f"data='{content}'"
        result = truncate_quoted_strings(text, max_len=100)
        self.assertIn("...[truncated]", result)
        self.assertLess(len(result), len(text))
    
    def test_preserves_short_quoted_string(self):
        """Short quoted strings are unchanged."""
        text = "data='short'"
        result = truncate_quoted_strings(text, max_len=100)
        self.assertEqual(result, text)
    
    def test_preserves_exact_max_len(self):
        """String at exact max_len is unchanged."""
        content = "a" * 100
        text = f"'{content}'"
        result = truncate_quoted_strings(text, max_len=100)
        self.assertEqual(result, text)
    
    def test_truncates_multiple_quoted_strings(self):
        """Multiple quoted strings are each truncated."""
        long_content = "x" * 500
        text = f"a='{long_content}' b='{long_content}'"
        result = truncate_quoted_strings(text, max_len=50)
        # Both should be truncated
        self.assertEqual(result.count("...[truncated]"), 2)
    
    def test_handles_nested_context(self):
        """Quoted strings with context are handled."""
        text = "pdb70_text='101\t7mjs_H\t1.000\t132' other='short'"
        result = truncate_quoted_strings(text, max_len=10)
        self.assertIn("pdb70_text=", result)
        self.assertIn("other='short'", result)  # Short one unchanged
    
    def test_handles_empty_quoted_string(self):
        """Empty quoted strings are unchanged."""
        text = "data=''"
        result = truncate_quoted_strings(text, max_len=100)
        self.assertEqual(result, text)
    
    def test_handles_no_quotes(self):
        """Text without quotes is unchanged."""
        text = "no quotes here"
        result = truncate_quoted_strings(text, max_len=100)
        self.assertEqual(result, text)


class TestReplaceTabs(unittest.TestCase):
    """Tests for replace_tabs function."""
    
    def test_replaces_single_tab(self):
        """Single tab is replaced with space."""
        self.assertEqual(replace_tabs("hello\tworld"), "hello world")
    
    def test_replaces_multiple_tabs(self):
        """Multiple tabs are replaced."""
        self.assertEqual(replace_tabs("a\tb\tc"), "a b c")
    
    def test_preserves_non_tab_text(self):
        """Text without tabs is unchanged."""
        text = "no tabs here"
        self.assertEqual(replace_tabs(text), text)
    
    def test_handles_empty_string(self):
        """Empty string returns empty string."""
        self.assertEqual(replace_tabs(""), "")


class TestRemoveControlChars(unittest.TestCase):
    """Tests for remove_control_chars function."""
    
    def test_removes_null(self):
        """Null character is removed."""
        self.assertEqual(remove_control_chars("hello\x00world"), "helloworld")
    
    def test_removes_bell(self):
        """Bell character is removed."""
        self.assertEqual(remove_control_chars("hello\x07world"), "helloworld")
    
    def test_preserves_newline(self):
        """Newline is preserved."""
        self.assertEqual(remove_control_chars("hello\nworld"), "hello\nworld")
    
    def test_preserves_normal_text(self):
        """Normal text is unchanged."""
        text = "Hello, World! 123"
        self.assertEqual(remove_control_chars(text), text)
    
    def test_removes_delete_char(self):
        """DEL character (0x7f) is removed."""
        self.assertEqual(remove_control_chars("hello\x7fworld"), "helloworld")


class TestTruncateString(unittest.TestCase):
    """Tests for truncate_string function."""
    
    def test_truncates_long_string(self):
        """Long string is truncated."""
        text = "a" * 1000
        result = truncate_string(text, max_len=100)
        self.assertLessEqual(len(result), 100)
        self.assertIn("[truncated]", result)
    
    def test_preserves_short_string(self):
        """Short string is unchanged."""
        text = "short"
        result = truncate_string(text, max_len=100)
        self.assertEqual(result, text)
    
    def test_preserves_exact_max_len(self):
        """String at exact max_len is unchanged."""
        text = "a" * 100
        result = truncate_string(text, max_len=100)
        self.assertEqual(result, text)
    
    def test_none_max_len_no_truncation(self):
        """None max_len means no truncation."""
        text = "a" * 10000
        result = truncate_string(text, max_len=None)
        self.assertEqual(result, text)
    
    def test_no_marker_option(self):
        """Can truncate without marker."""
        text = "a" * 100
        result = truncate_string(text, max_len=50, add_marker=False)
        self.assertEqual(result, "a" * 50)
        self.assertNotIn("[truncated]", result)


class TestSanitizeString(unittest.TestCase):
    """Tests for sanitize_string function."""
    
    def test_removes_markers_and_tabs(self):
        """Combines marker removal and tab replacement."""
        text = "Hello ZZTAZZ world\twith tabs"
        result = sanitize_string(text)
        self.assertNotIn("ZZTAZZ", result)
        self.assertNotIn("\t", result)
    
    def test_removes_control_chars(self):
        """Control characters are removed."""
        text = "Hello\x00\x07world"
        result = sanitize_string(text)
        self.assertEqual(result, "Helloworld")
    
    def test_applies_max_len(self):
        """Max length is applied."""
        text = "a" * 1000
        result = sanitize_string(text, max_len=100)
        self.assertLessEqual(len(result), 100)
    
    def test_full_sanitization_pipeline(self):
        """Full pipeline works correctly."""
        text = "data ZZTAZZ\twith\x00control ZZCRZZ chars"
        result = sanitize_string(text, max_len=50)
        self.assertNotIn("ZZ", result)
        self.assertNotIn("\t", result)
        self.assertNotIn("\x00", result)
        self.assertLessEqual(len(result), 50)


class TestSanitizeForTransport(unittest.TestCase):
    """Tests for sanitize_for_transport function."""
    
    def test_basic_sanitization(self):
        """Basic sanitization without quote truncation."""
        text = "Hello ZZTAZZ\tworld"
        result = sanitize_for_transport(text)
        self.assertNotIn("ZZTAZZ", result)
        self.assertNotIn("\t", result)
    
    def test_with_quote_truncation(self):
        """Quote truncation is applied when enabled."""
        long_content = "x" * 1000
        text = f"data='{long_content}'"
        result = sanitize_for_transport(text, truncate_quotes=True, quote_max_len=100)
        self.assertIn("...[truncated]", result)
        self.assertLess(len(result), len(text))
    
    def test_without_quote_truncation(self):
        """Quotes not truncated when disabled."""
        long_content = "x" * 1000
        text = f"data='{long_content}'"
        result = sanitize_for_transport(text, truncate_quotes=False)
        self.assertNotIn("...[truncated]", result)
    
    def test_markers_removed_before_truncation(self):
        """ZZxxZZ markers are removed before quote truncation."""
        # This ensures we don't split a marker in half
        text = "data='ZZTAZZ ZZTAZZ ZZTAZZ' other"
        result = sanitize_for_transport(text, truncate_quotes=True, quote_max_len=10)
        # Should not have partial markers like "ZZT" or "AZZ"
        self.assertNotIn("ZZT", result)
        self.assertNotIn("AZZ", result)
    
    def test_order_of_operations(self):
        """Correct order: remove markers, truncate quotes, sanitize."""
        # Create text with markers inside quotes
        markers = "ZZTAZZ " * 100  # Would be ~600 chars
        text = f"data='{markers}' end"
        result = sanitize_for_transport(text, truncate_quotes=True, quote_max_len=50)
        # Markers should be removed first, then truncation
        self.assertNotIn("ZZTAZZ", result)
        self.assertIn("data=", result)
        self.assertIn("end", result)


class TestSanitizeDictRecursive(unittest.TestCase):
    """Tests for sanitize_dict_recursive function."""
    
    def test_sanitizes_nested_dict(self):
        """Nested dict strings are sanitized."""
        obj = {"a": {"b": "Hello\tZZTAZZ"}}
        result = sanitize_dict_recursive(obj)
        self.assertEqual(result["a"]["b"], "Hello ")
    
    def test_sanitizes_list_items(self):
        """List items are sanitized."""
        obj = ["Hello\tworld", "foo\x00bar"]
        result = sanitize_dict_recursive(obj)
        self.assertEqual(result[0], "Hello world")
        self.assertEqual(result[1], "foobar")
    
    def test_preserves_non_strings(self):
        """Non-string values are preserved."""
        obj = {"num": 123, "flag": True, "none": None}
        result = sanitize_dict_recursive(obj)
        self.assertEqual(result["num"], 123)
        self.assertEqual(result["flag"], True)
        self.assertIsNone(result["none"])
    
    def test_applies_max_len(self):
        """Max length is applied to strings."""
        obj = {"long": "a" * 1000}
        result = sanitize_dict_recursive(obj, max_len_per_string=100)
        self.assertLessEqual(len(result["long"]), 100)
    
    def test_max_depth_protection(self):
        """Maximum depth prevents infinite recursion."""
        # Create deeply nested structure
        obj = {"a": {"a": {"a": {"a": {"a": {"a": {"a": {"a": {"a": {"a": {"a": "deep"}}}}}}}}}}}
        # Should not raise, and should handle gracefully
        result = sanitize_dict_recursive(obj, max_depth=5)
        self.assertIsNotNone(result)


class TestEncodeDecodeRoundtrip(unittest.TestCase):
    """Tests for encode/decode functions."""
    
    def test_roundtrip_simple(self):
        """Simple text survives roundtrip."""
        text = "Hello world"
        encoded = encode_for_rest(text)
        decoded = decode_from_rest(encoded)
        self.assertEqual(decoded, text)
    
    def test_roundtrip_with_special_chars(self):
        """Text with special chars survives roundtrip."""
        text = '{"key": "value"}\nnewline'
        encoded = encode_for_rest(text)
        decoded = decode_from_rest(encoded)
        self.assertEqual(decoded, text)
    
    def test_roundtrip_with_tabs(self):
        """Text with tabs survives roundtrip."""
        text = "col1\tcol2\tcol3"
        encoded = encode_for_rest(text)
        decoded = decode_from_rest(encoded)
        self.assertEqual(decoded, text)
    
    def test_sanitize_then_roundtrip(self):
        """Sanitized text survives roundtrip."""
        text = "data ZZTAZZ\twith markers"
        sanitized = sanitize_for_transport(text)
        encoded = encode_for_rest(sanitized)
        decoded = decode_from_rest(encoded)
        self.assertEqual(decoded, sanitized)
        # And the result should be clean
        self.assertNotIn("ZZTAZZ", decoded)
        self.assertNotIn("\t", decoded)


class TestEdgeCases(unittest.TestCase):
    """Tests for edge cases and real-world scenarios."""
    
    def test_pdb70_text_scenario(self):
        """Real-world pdb70_text scenario is handled."""
        # Simulates the actual problematic data
        pdb70_data = "101\t7mjs_H\t1.000\t132\t" * 100  # Tab-separated data
        text = f"pdb70_text='{pdb70_data}'"
        
        result = sanitize_for_transport(text, truncate_quotes=True, quote_max_len=100)
        
        # Should be much shorter
        self.assertLess(len(result), len(text))
        # Should have truncation marker
        self.assertIn("[truncated]", result)
        # Should not have tabs (they were in the quoted string that got truncated)
        # The truncation happens after marker removal but tabs in remaining text become spaces
        self.assertNotIn("\t", result)
    
    def test_log_with_embedded_markers(self):
        """Log content with embedded markers from previous runs."""
        log = """
        Running program...
        Previous output: data=ZZDQZZ101ZZTAZZ7mjsZZDQZZ
        More output ZZCRZZ here
        """
        result = sanitize_for_transport(log)
        self.assertNotIn("ZZDQZZ", result)
        self.assertNotIn("ZZTAZZ", result)
        self.assertNotIn("ZZCRZZ", result)
    
    def test_large_log_performance(self):
        """Large log is handled efficiently."""
        # Create a large log with many quoted strings
        log_parts = [f"metric_{i}='{'x' * 1000}'" for i in range(100)]
        log = "\n".join(log_parts)
        
        # Should complete without hanging
        result = sanitize_for_transport(log, truncate_quotes=True, quote_max_len=100)
        
        # Should be much smaller
        self.assertLess(len(result), len(log) / 5)


class TestConfigLoading(unittest.TestCase):
    """Tests for YAML config loading."""
    
    def setUp(self):
        """Clear config cache before each test."""
        clear_config_cache()
    
    def tearDown(self):
        """Clear config cache after each test."""
        clear_config_cache()
    
    def test_get_default_config(self):
        """Default config is returned when YAML not available."""
        config = get_transport_config()
        self.assertIsNotNone(config)
        self.assertIn('fields', config)
        self.assertIn('request', config['fields'])
    
    def test_config_has_log_content_settings(self):
        """Config has log_content field settings."""
        config = get_transport_config()
        settings = get_field_settings(config, 'request', 'log_content')
        self.assertTrue(settings.get('sanitize', False))
        self.assertIsNotNone(settings.get('max_length'))
    
    def test_config_has_history_settings(self):
        """Config has history item field settings."""
        config = get_transport_config()
        settings = get_field_settings(config, 'request', 'history.result')
        self.assertTrue(settings.get('sanitize', False))


class TestSanitizeField(unittest.TestCase):
    """Tests for config-driven field sanitization."""
    
    def test_sanitize_log_content(self):
        """Log content is sanitized with quote truncation."""
        text = "data='x" * 100 + "' end ZZTAZZ"
        result = sanitize_field(text, 'log_content', 'request')
        self.assertNotIn("ZZTAZZ", result)
        # Should have truncation if quotes enabled in config
    
    def test_sanitize_user_advice(self):
        """User advice is sanitized."""
        text = "Hello\tworld ZZTAZZ"
        result = sanitize_field(text, 'user_advice', 'request')
        self.assertNotIn("ZZTAZZ", result)
        self.assertNotIn("\t", result)
    
    def test_unknown_field_passes_through(self):
        """Unknown fields pass through unchanged."""
        text = "unchanged"
        result = sanitize_field(text, 'unknown_field', 'request')
        self.assertEqual(result, text)


class TestSanitizeRequest(unittest.TestCase):
    """Tests for full request sanitization."""
    
    def test_sanitize_request_log_content(self):
        """Request log_content is sanitized."""
        request = {
            'api_version': '2.0',
            'log_content': 'Hello\tZZTAZZ world',
            'files': ['/path/to/file.mtz'],
        }
        result = sanitize_request(request)
        self.assertNotIn("ZZTAZZ", result['log_content'])
        self.assertNotIn("\t", result['log_content'])
        # Files should be unchanged
        self.assertEqual(result['files'], request['files'])
    
    def test_sanitize_request_history(self):
        """Request history items are sanitized."""
        request = {
            'api_version': '2.0',
            'history': [
                {
                    'program': 'phenix.refine',
                    'result': 'Output\tZZTAZZ here',
                }
            ],
        }
        result = sanitize_request(request)
        self.assertNotIn("ZZTAZZ", result['history'][0]['result'])
    
    def test_sanitize_request_preserves_structure(self):
        """Request structure is preserved."""
        request = {
            'api_version': '2.0',
            'cycle_number': 1,
            'files': ['/a.mtz', '/b.pdb'],
            'session_state': {'resolution': 2.0},
        }
        result = sanitize_request(request)
        self.assertEqual(result['api_version'], '2.0')
        self.assertEqual(result['cycle_number'], 1)
        self.assertEqual(result['files'], ['/a.mtz', '/b.pdb'])
        self.assertEqual(result['session_state']['resolution'], 2.0)


class TestSanitizeResponse(unittest.TestCase):
    """Tests for full response sanitization."""
    
    def test_sanitize_response_decision(self):
        """Response decision is sanitized."""
        response = {
            'api_version': '2.0',
            'decision': {
                'program': 'phenix.refine',
                'reasoning': 'Based on\tZZTAZZ analysis',
            },
        }
        result = sanitize_response(response)
        self.assertNotIn("ZZTAZZ", result['decision']['reasoning'])
    
    def test_sanitize_response_debug_log(self):
        """Response debug_log is sanitized."""
        response = {
            'api_version': '2.0',
            'debug_log': ['Entry 1\tZZTAZZ', 'Entry 2'],
        }
        result = sanitize_response(response)
        self.assertNotIn("ZZTAZZ", result['debug_log'][0])


class TestUnifiedTransport(unittest.TestCase):
    """Tests for unified transport functions."""
    
    def test_prepare_request_basic(self):
        """Basic request preparation works."""
        request = {
            'api_version': '2.0',
            'files': ['/path/to/file.mtz'],
            'log_content': 'Hello world',
        }
        encoded, original = prepare_request_for_transport(request, do_encode=True)
        self.assertIsInstance(encoded, str)
        self.assertIsInstance(original, str)
        # Encoded should have markers
        self.assertIn('ZZLBZZ', encoded)  # Left brace
    
    def test_prepare_request_no_encode(self):
        """Request preparation without encoding."""
        request = {
            'api_version': '2.0',
            'files': [],
        }
        result, original = prepare_request_for_transport(request, do_encode=False)
        self.assertEqual(result, original)
        # Should be valid JSON
        import json
        parsed = json.loads(result)
        self.assertEqual(parsed['api_version'], '2.0')
    
    def test_process_request_basic(self):
        """Basic request processing works."""
        request = {
            'api_version': '2.0',
            'files': ['/path/to/file.mtz'],
            'cycle_number': 1,
        }
        encoded, _ = prepare_request_for_transport(request, do_encode=True)
        decoded = process_request_from_transport(encoded, was_encoded=True)
        self.assertNotIn('error', decoded)
        self.assertEqual(decoded['api_version'], '2.0')
        self.assertEqual(decoded['cycle_number'], 1)
    
    def test_process_request_not_encoded(self):
        """Request processing without encoding."""
        import json
        request = {'api_version': '2.0', 'files': []}
        json_str = json.dumps(request)
        decoded = process_request_from_transport(json_str, was_encoded=False)
        self.assertEqual(decoded['api_version'], '2.0')
    
    def test_process_request_invalid_json(self):
        """Invalid JSON returns error."""
        result = process_request_from_transport("not json", was_encoded=False)
        self.assertIn('error', result)
    
    def test_prepare_response_basic(self):
        """Basic response preparation works."""
        response = {
            'api_version': '2.0',
            'decision': {'program': 'phenix.refine'},
        }
        encoded, original = prepare_response_for_transport(response, do_encode=True)
        self.assertIsInstance(encoded, str)
        self.assertIn('ZZLBZZ', encoded)
    
    def test_process_response_basic(self):
        """Basic response processing works."""
        response = {
            'api_version': '2.0',
            'decision': {'program': 'phenix.refine'},
        }
        encoded, _ = prepare_response_for_transport(response, do_encode=True)
        decoded = process_response_from_transport(encoded, was_encoded=True)
        self.assertNotIn('error', decoded)
        self.assertEqual(decoded['decision']['program'], 'phenix.refine')
    
    def test_full_roundtrip(self):
        """Full request roundtrip preserves data."""
        request = {
            'api_version': '2.0',
            'files': ['/a.mtz', '/b.pdb'],
            'cycle_number': 3,
            'session_state': {
                'resolution': 2.5,
                'experiment_type': 'xray',
            },
            'history': [
                {'program': 'phenix.xtriage', 'result': 'OK'},
            ],
            'log_content': 'Some log output',
        }
        
        # Prepare
        encoded, _ = prepare_request_for_transport(request, do_encode=True)
        
        # Process
        decoded = process_request_from_transport(encoded, was_encoded=True)
        
        # Verify
        self.assertNotIn('error', decoded)
        self.assertEqual(decoded['api_version'], '2.0')
        self.assertEqual(len(decoded['files']), 2)
        self.assertEqual(decoded['cycle_number'], 3)
        self.assertEqual(decoded['session_state']['experiment_type'], 'xray')
        self.assertEqual(len(decoded['history']), 1)
    
    def test_roundtrip_with_problematic_content(self):
        """Roundtrip handles content that previously caused issues."""
        request = {
            'api_version': '2.0',
            'files': [],
            'log_content': "pdb70_text='101\t7mjs_H\t1.000' ZZTAZZ marker",
            'session_state': {'experiment_type': 'cryoem'},
        }
        
        # Prepare
        encoded, _ = prepare_request_for_transport(request, do_encode=True)
        
        # Process
        decoded = process_request_from_transport(encoded, was_encoded=True)
        
        # Should not fail
        self.assertNotIn('error', decoded)
        # Markers should be removed
        self.assertNotIn('ZZTAZZ', decoded.get('log_content', ''))
        # experiment_type preserved
        self.assertEqual(decoded['session_state']['experiment_type'], 'cryoem')


class TestVerifyRoundtrip(unittest.TestCase):
    """Tests for verify_roundtrip helper."""
    
    def test_verify_roundtrip_success(self):
        """Successful roundtrip verification."""
        request = {
            'api_version': '2.0',
            'files': ['/a.mtz'],
            'session_state': {'experiment_type': 'xray'},
            'history': [],
        }
        success, message, details = verify_roundtrip(request)
        self.assertTrue(success)
        self.assertIn('OK', message)
    
    def test_verify_roundtrip_with_content(self):
        """Roundtrip verification with substantial content."""
        request = {
            'api_version': '2.0',
            'files': [f'/file_{i}.mtz' for i in range(10)],
            'session_state': {'experiment_type': 'cryoem', 'resolution': 2.5},
            'history': [{'program': f'prog_{i}'} for i in range(5)],
            'log_content': 'x' * 1000,
        }
        success, message, details = verify_roundtrip(request)
        self.assertTrue(success, f"Roundtrip failed: {message}")


def run_all_tests():
    """Run all transport tests (for run_all_tests.py compatibility)."""
    # Create a test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestRemoveMarkers))
    suite.addTests(loader.loadTestsFromTestCase(TestTruncateQuotedStrings))
    suite.addTests(loader.loadTestsFromTestCase(TestReplaceTabs))
    suite.addTests(loader.loadTestsFromTestCase(TestRemoveControlChars))
    suite.addTests(loader.loadTestsFromTestCase(TestTruncateString))
    suite.addTests(loader.loadTestsFromTestCase(TestSanitizeString))
    suite.addTests(loader.loadTestsFromTestCase(TestSanitizeForTransport))
    suite.addTests(loader.loadTestsFromTestCase(TestSanitizeDictRecursive))
    suite.addTests(loader.loadTestsFromTestCase(TestEncodeDecodeRoundtrip))
    suite.addTests(loader.loadTestsFromTestCase(TestEdgeCases))
    suite.addTests(loader.loadTestsFromTestCase(TestConfigLoading))
    suite.addTests(loader.loadTestsFromTestCase(TestSanitizeField))
    suite.addTests(loader.loadTestsFromTestCase(TestSanitizeRequest))
    suite.addTests(loader.loadTestsFromTestCase(TestSanitizeResponse))
    suite.addTests(loader.loadTestsFromTestCase(TestUnifiedTransport))
    suite.addTests(loader.loadTestsFromTestCase(TestVerifyRoundtrip))
    
    # Run with verbosity
    runner = unittest.TextTestRunner(verbosity=1)
    result = runner.run(suite)
    
    # Raise exception if any tests failed (for run_all_tests.py)
    if not result.wasSuccessful():
        raise Exception(f"{len(result.failures)} test(s) failed, {len(result.errors)} error(s)")


if __name__ == '__main__':
    unittest.main()
