"""
Transport module for PHENIX AI Agent client-server communication.

This module provides consistent encoding/decoding for both LocalAgent
and RemoteAgent, ensuring identical behavior regardless of transport.

Sanitization Order (important!):
1. Remove ZZxxZZ markers first (before truncation to avoid splitting markers)
2. Truncate long quoted strings
3. Replace tabs and control characters
4. Apply max length truncation

Usage:
    from agent.transport import sanitize_string, sanitize_for_transport

    # Sanitize a single string
    clean = sanitize_string(text, max_len=500)

    # Sanitize with quoted string truncation (for log content)
    clean = sanitize_for_transport(log_text, max_len=50000, truncate_quotes=True)

    # Get config-driven field settings
    config = get_transport_config()
    settings = get_field_settings(config, 'request', 'log_content')
"""

from __future__ import absolute_import, division, print_function

import os
import re

# =============================================================================
# CONFIGURATION LOADING
# =============================================================================

# Cache for loaded config
_config_cache = None

def get_transport_config():
    """
    Load transport configuration from YAML.

    Returns cached config if already loaded.

    Returns:
        dict: Transport configuration, or default config if YAML not found
    """
    global _config_cache

    if _config_cache is not None:
        return _config_cache

    # Try to load YAML config
    config = _load_yaml_config()

    if config is None:
        # Use default config
        config = _get_default_config()

    _config_cache = config
    return config


def _load_yaml_config():
    """
    Attempt to load transport.yaml configuration.

    Returns:
        dict or None: Parsed config, or None if not found/failed
    """
    try:
        import yaml
    except ImportError:
        return None

    # Find the YAML file
    yaml_paths = [
        os.path.join(os.path.dirname(os.path.dirname(__file__)), 'knowledge', 'transport.yaml'),
        os.path.join(os.path.dirname(__file__), '..', 'knowledge', 'transport.yaml'),
    ]

    # Also try PHENIX paths
    try:
        import libtbx.load_env
        phenix_path = os.path.join(
            libtbx.env.find_in_repositories("phenix"),
            'phenix', 'phenix_ai', 'knowledge', 'transport.yaml'
        )
        yaml_paths.append(phenix_path)
    except Exception:
        pass

    for yaml_path in yaml_paths:
        if os.path.exists(yaml_path):
            try:
                with open(yaml_path, 'r') as f:
                    return yaml.safe_load(f)
            except Exception as e:
                import sys
                print(f"Warning: Failed to load transport.yaml: {e}", file=sys.stderr)
                return None

    return None


def _get_default_config():
    """
    Return default configuration when YAML is not available.

    Returns:
        dict: Default transport configuration
    """
    return {
        'marker_pattern': r'ZZ[A-Z]{2}ZZ',
        'replacements': {
            '\t': ' ',
            '\r': '',
        },
        'quoted_string_truncation': {
            'enabled': True,
            'quote_chars': ["'"],
            'default_max_length': 500,
            'truncation_marker': '...[truncated]',
        },
        'fields': {
            'request': {
                'log_content': {
                    'max_length': 50000,
                    'sanitize': True,
                    'truncate_quotes': True,
                    'quote_max_length': 500,
                },
                'user_advice': {
                    'max_length': 10000,
                    'sanitize': True,
                    'truncate_quotes': False,
                },
                'history': {
                    'type': 'list',
                    'item_fields': {
                        'program': {'max_length': 100, 'sanitize': True},
                        'command': {'max_length': 2000, 'sanitize': True},
                        'result': {
                            'max_length': 1000,
                            'sanitize': True,
                            'truncate_quotes': True,
                            'quote_max_length': 200,
                        },
                        'metrics': {
                            'type': 'dict',
                            'sanitize': 'recursive',
                            'max_length_per_value': 500,
                        },
                    },
                },
            },
            'response': {
                'decision': {
                    'type': 'dict',
                    'fields': {
                        'reasoning': {'max_length': 5000, 'sanitize': True},
                        'command': {'max_length': 5000, 'sanitize': True},
                    },
                },
                'debug_log': {
                    'type': 'list',
                    'item_max_length': 1000,
                    'sanitize': True,
                },
            },
        },
        'performance': {
            'max_recursion_depth': 10,
        },
    }


def get_field_settings(config, direction, field_path):
    """
    Get sanitization settings for a specific field.

    Args:
        config: Transport configuration dict
        direction: 'request' or 'response'
        field_path: Dot-separated path like 'log_content' or 'history.result'

    Returns:
        dict: Field settings, or empty dict if not found
    """
    fields = config.get('fields', {}).get(direction, {})

    parts = field_path.split('.')
    current = fields

    for part in parts:
        if isinstance(current, dict):
            if part in current:
                current = current[part]
            elif 'item_fields' in current and part in current['item_fields']:
                current = current['item_fields'][part]
            elif 'fields' in current and part in current['fields']:
                current = current['fields'][part]
            else:
                return {}
        else:
            return {}

    return current if isinstance(current, dict) else {}


def clear_config_cache():
    """Clear the cached configuration (useful for testing)."""
    global _config_cache
    _config_cache = None

# =============================================================================
# CONSTANTS
# =============================================================================

# Pattern to match ZZxxZZ markers (where xx is any two uppercase letters)
# These are REST encoding markers that must be removed before transport
ZZXXZZ_PATTERN = re.compile(r'ZZ[A-Z]{2}ZZ')

# Control characters to remove (except newline \n which is \x0a)
# Includes: \x00-\x08, \x0b, \x0c, \x0e-\x1f, \x7f
CONTROL_CHAR_PATTERN = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]')

# Pattern to match single-quoted strings
# Captures content between single quotes (non-greedy for nested quotes)
SINGLE_QUOTED_PATTERN = re.compile(r"'([^']*)'")


# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def remove_markers(text):
    """
    Remove ZZxxZZ REST encoding markers from text.

    These markers appear in logs from previous server runs and would be
    decoded back to control characters, breaking JSON parsing.

    Args:
        text: String that may contain ZZxxZZ markers

    Returns:
        String with all ZZxxZZ markers removed

    Example:
        >>> remove_markers("Hello ZZTAZZ world ZZCRZZ")
        "Hello  world "
    """
    if not isinstance(text, str):
        return str(text) if text else ""
    return ZZXXZZ_PATTERN.sub('', text)


def truncate_quoted_strings(text, max_len=500, quote_char="'"):
    """
    Truncate long quoted strings in text.

    This is useful for log content that contains large data dumps
    like HHblits pdb70_text output. These are enclosed in quotes
    and can be safely truncated without losing important information.

    Args:
        text: String containing quoted substrings
        max_len: Maximum length for quoted content (default 500)
        quote_char: Quote character to look for (default single quote)

    Returns:
        String with long quoted sections truncated

    Example:
        >>> truncate_quoted_strings("data='abcd...1000 chars...'", max_len=10)
        "data='abcd...[truncated]...'"
    """
    if not isinstance(text, str):
        return str(text) if text else ""

    if quote_char == "'":
        pattern = SINGLE_QUOTED_PATTERN
    else:
        # For other quote chars, build pattern dynamically
        pattern = re.compile(r"{0}([^{0}]*){0}".format(re.escape(quote_char)))

    def truncate_match(match):
        content = match.group(1)
        if len(content) <= max_len:
            return match.group(0)  # Return unchanged

        # Keep beginning with truncation marker
        # We keep the start because it usually has the variable/field name context
        truncated = content[:max_len] + "...[truncated]"
        return quote_char + truncated + quote_char

    return pattern.sub(truncate_match, text)


def replace_tabs(text):
    """
    Replace tab characters with spaces.

    Tabs cause issues with REST transport encoding/decoding.

    Args:
        text: String that may contain tabs

    Returns:
        String with tabs replaced by spaces
    """
    if not isinstance(text, str):
        return str(text) if text else ""
    return text.replace('\t', ' ')


def remove_control_chars(text):
    """
    Remove control characters except newline.

    Control characters (ASCII 0-31 except newline, and 127) can break
    JSON encoding and REST transport.

    Args:
        text: String that may contain control characters

    Returns:
        String with control characters removed (newlines preserved)
    """
    if not isinstance(text, str):
        return str(text) if text else ""
    return CONTROL_CHAR_PATTERN.sub('', text)


def truncate_string(text, max_len, add_marker=True):
    """
    Truncate string to maximum length.

    Args:
        text: String to truncate
        max_len: Maximum length (None for no limit)
        add_marker: Whether to add "... [truncated]" marker

    Returns:
        Truncated string
    """
    if not isinstance(text, str):
        return str(text) if text else ""

    if max_len is None or len(text) <= max_len:
        return text

    if add_marker:
        marker = "... [truncated]"
        return text[:max_len - len(marker)] + marker
    else:
        return text[:max_len]


# =============================================================================
# HIGH-LEVEL FUNCTIONS
# =============================================================================

def sanitize_string(text, max_len=None):
    """
    Apply standard sanitization to a string.

    This applies the basic sanitization steps:
    1. Remove ZZxxZZ markers
    2. Replace tabs with spaces
    3. Remove control characters
    4. Truncate to max_len (if specified)

    Does NOT truncate quoted strings - use sanitize_for_transport() for that.

    Args:
        text: String to sanitize
        max_len: Maximum length (None for no limit)

    Returns:
        Sanitized string
    """
    if not isinstance(text, str):
        return str(text) if text else ""

    # Step 1: Remove ZZxxZZ markers
    text = remove_markers(text)

    # Step 2: Replace tabs
    text = replace_tabs(text)

    # Step 3: Remove control characters
    text = remove_control_chars(text)

    # Step 4: Truncate if needed
    if max_len is not None:
        text = truncate_string(text, max_len)

    return text


def sanitize_for_transport(text, max_len=None, truncate_quotes=False, quote_max_len=500):
    """
    Full sanitization for transport, including quoted string truncation.

    This is the main entry point for sanitizing text that will be sent
    through the REST transport layer.

    Sanitization order:
    1. Remove ZZxxZZ markers (before truncation to avoid splitting)
    2. Truncate long quoted strings (if enabled)
    3. Replace tabs with spaces
    4. Remove control characters
    5. Truncate to max_len

    Args:
        text: String to sanitize
        max_len: Maximum total length (None for no limit)
        truncate_quotes: Whether to truncate long quoted strings
        quote_max_len: Max length for quoted string content (default 500)

    Returns:
        Sanitized string ready for transport
    """
    if not isinstance(text, str):
        return str(text) if text else ""

    # Step 1: Remove ZZxxZZ markers FIRST (before any truncation)
    text = remove_markers(text)

    # Step 2: Truncate long quoted strings (optional)
    if truncate_quotes:
        text = truncate_quoted_strings(text, max_len=quote_max_len)

    # Step 3: Replace tabs
    text = replace_tabs(text)

    # Step 4: Remove control characters
    text = remove_control_chars(text)

    # Step 5: Truncate total length if needed
    if max_len is not None:
        text = truncate_string(text, max_len)

    return text


def sanitize_dict_recursive(obj, max_len_per_string=500, depth=0, max_depth=10):
    """
    Recursively sanitize all strings in a dict/list structure.

    Used for sanitizing metrics and other nested data structures.

    Args:
        obj: Dict, list, or scalar value
        max_len_per_string: Max length for each string value
        depth: Current recursion depth
        max_depth: Maximum recursion depth to prevent infinite loops

    Returns:
        Sanitized structure with same shape as input
    """
    if depth > max_depth:
        return obj

    if isinstance(obj, dict):
        return {k: sanitize_dict_recursive(v, max_len_per_string, depth + 1, max_depth)
                for k, v in obj.items()}
    elif isinstance(obj, list):
        return [sanitize_dict_recursive(item, max_len_per_string, depth + 1, max_depth)
                for item in obj]
    elif isinstance(obj, str):
        return sanitize_string(obj, max_len=max_len_per_string)
    else:
        return obj


# =============================================================================
# CONFIG-DRIVEN SANITIZATION
# =============================================================================

def sanitize_field(value, field_path, direction='request', config=None):
    """
    Sanitize a field value using YAML configuration.

    This looks up the field settings from the config and applies
    the appropriate sanitization.

    Args:
        value: The value to sanitize
        field_path: Dot-separated path like 'log_content' or 'history.result'
        direction: 'request' or 'response'
        config: Transport config (loads default if None)

    Returns:
        Sanitized value
    """
    if config is None:
        config = get_transport_config()

    settings = get_field_settings(config, direction, field_path)

    if not settings:
        # No config for this field - return unchanged
        return value

    # Check if sanitization is enabled for this field
    if not settings.get('sanitize', False):
        return value

    # Handle different types
    field_type = settings.get('type', 'string')

    if field_type == 'list':
        if isinstance(value, list):
            item_max = settings.get('item_max_length')
            return [sanitize_string(item, max_len=item_max) if isinstance(item, str) else item
                    for item in value]
        return value

    elif field_type == 'dict' or settings.get('sanitize') == 'recursive':
        if isinstance(value, dict):
            max_per_value = settings.get('max_length_per_value', 500)
            max_depth = config.get('performance', {}).get('max_recursion_depth', 10)
            return sanitize_dict_recursive(value, max_per_value, max_depth=max_depth)
        return value

    elif isinstance(value, str):
        # String field
        max_len = settings.get('max_length')
        truncate_quotes = settings.get('truncate_quotes', False)
        quote_max_len = settings.get('quote_max_length', 500)

        if truncate_quotes:
            return sanitize_for_transport(
                value,
                max_len=max_len,
                truncate_quotes=True,
                quote_max_len=quote_max_len
            )
        else:
            return sanitize_string(value, max_len=max_len)

    return value


def sanitize_request(request, config=None):
    """
    Sanitize a full request dict using YAML configuration.

    This applies field-specific sanitization rules from the config.

    Args:
        request: Request dict to sanitize
        config: Transport config (loads default if None)

    Returns:
        Sanitized request dict
    """
    if config is None:
        config = get_transport_config()

    if not isinstance(request, dict):
        return request

    result = {}

    for key, value in request.items():
        if key == 'log_content':
            result[key] = sanitize_field(value, 'log_content', 'request', config)

        elif key == 'user_advice':
            result[key] = sanitize_field(value, 'user_advice', 'request', config)

        elif key == 'history' and isinstance(value, list):
            # Sanitize each history item
            sanitized_history = []
            for item in value:
                if isinstance(item, dict):
                    sanitized_item = {}
                    for item_key, item_value in item.items():
                        sanitized_item[item_key] = sanitize_field(
                            item_value, f'history.{item_key}', 'request', config
                        )
                    sanitized_history.append(sanitized_item)
                else:
                    sanitized_history.append(item)
            result[key] = sanitized_history

        else:
            # Pass through unchanged
            result[key] = value

    return result


def sanitize_response(response, config=None):
    """
    Sanitize a full response dict using YAML configuration.

    This applies field-specific sanitization rules from the config.

    Args:
        response: Response dict to sanitize
        config: Transport config (loads default if None)

    Returns:
        Sanitized response dict
    """
    if config is None:
        config = get_transport_config()

    if not isinstance(response, dict):
        return response

    result = {}

    for key, value in response.items():
        if key == 'decision' and isinstance(value, dict):
            sanitized_decision = {}
            for dec_key, dec_value in value.items():
                sanitized_decision[dec_key] = sanitize_field(
                    dec_value, f'decision.{dec_key}', 'response', config
                )
            result[key] = sanitized_decision

        elif key == 'debug_log':
            result[key] = sanitize_field(value, 'debug_log', 'response', config)

        elif key == 'history_record' and isinstance(value, dict):
            # Use request history settings
            sanitized_record = {}
            for rec_key, rec_value in value.items():
                sanitized_record[rec_key] = sanitize_field(
                    rec_value, f'history.{rec_key}', 'request', config
                )
            result[key] = sanitized_record

        else:
            # Pass through unchanged
            result[key] = value

    return result


# =============================================================================
# ENCODING/DECODING (Phase 3 - placeholders for now)
# =============================================================================

def encode_for_rest(text):
    """
    Encode text for REST transport using text_as_simple_string.

    This is a wrapper that will be used by both LocalAgent and RemoteAgent
    to ensure consistent encoding.

    Args:
        text: String to encode (typically JSON)

    Returns:
        Encoded string safe for REST transport
    """
    # Import here to avoid circular imports and allow standalone testing
    try:
        from phenix.rest import text_as_simple_string
        return text_as_simple_string(text)
    except ImportError:
        # Fallback for testing without PHENIX
        # This mimics the basic behavior
        replacements = [
            ('\n', 'ZZCRZZ'),
            ('{', 'ZZLBZZ'),
            ('}', 'ZZRBZZ'),
            (';', 'ZZSCZZ'),
            ('#', 'ZZHAZZ'),
            ('\t', 'ZZTAZZ'),
            ('"', 'ZZDQZZ'),
            ("'", 'ZZSQZZ'),
            ('`', 'ZZSRZZ'),
            ('!', 'ZZEXZZ'),
            ('$', 'ZZDSZZ'),
        ]
        result = text
        for char, marker in replacements:
            result = result.replace(char, marker)
        return result


def decode_from_rest(encoded):
    """
    Decode text from REST transport using simple_string_as_text.

    This is a wrapper that will be used by both LocalAgent and RemoteAgent
    to ensure consistent decoding.

    Args:
        encoded: Encoded string from REST transport

    Returns:
        Decoded original string
    """
    # Import here to avoid circular imports and allow standalone testing
    try:
        from phenix.rest import simple_string_as_text
        return simple_string_as_text(encoded)
    except ImportError:
        # Fallback for testing without PHENIX
        replacements = [
            ('ZZCRZZ', '\n'),
            ('ZZLBZZ', '{'),
            ('ZZRBZZ', '}'),
            ('ZZSCZZ', ';'),
            ('ZZHAZZ', '#'),
            ('ZZTAZZ', '\t'),
            ('ZZDQZZ', '"'),
            ('ZZSQZZ', "'"),
            ('ZZSRZZ', '`'),
            ('ZZEXZZ', '!'),
            ('ZZDSZZ', '$'),
        ]
        result = encoded
        for marker, char in replacements:
            result = result.replace(marker, char)
        return result


# =============================================================================
# UNIFIED REQUEST/RESPONSE PROCESSING
# =============================================================================

def prepare_request_for_transport(request, do_encode=True):
    """
    Prepare a request dict for transport.

    This is the unified entry point for both LocalAgent and RemoteAgent.
    It applies sanitization and encoding in the correct order.

    Steps:
    1. Sanitize the request (using YAML config)
    2. Serialize to JSON
    3. Encode for REST transport (if do_encode=True)

    Args:
        request: Request dict to prepare
        do_encode: Whether to apply REST encoding (default True)

    Returns:
        tuple: (prepared_string, original_json)
            - prepared_string: Ready for transport (encoded if do_encode)
            - original_json: The JSON before encoding (for debugging)
    """
    import json

    # Step 1: Sanitize the request
    sanitized = sanitize_request(request)

    # Step 2: Serialize to JSON
    json_str = json.dumps(sanitized, indent=None, separators=(',', ':'))

    # Final safety: replace any remaining tab escape sequences
    json_str = json_str.replace('\\t', ' ')

    # Step 3: Encode if requested
    if do_encode:
        encoded = encode_for_rest(json_str)
        return encoded, json_str
    else:
        return json_str, json_str


def process_request_from_transport(encoded_or_json, was_encoded=True):
    """
    Process a received request from transport.

    This is the unified entry point for receiving requests.
    It decodes and parses the request.

    Steps:
    1. Decode from REST encoding (if was_encoded=True)
    2. Parse JSON to dict

    Args:
        encoded_or_json: The received string
        was_encoded: Whether REST encoding was applied (default True)

    Returns:
        dict: Parsed request, or dict with 'error' key on failure
    """
    import json

    try:
        # Step 1: Decode if needed
        if was_encoded:
            json_str = decode_from_rest(encoded_or_json)
        else:
            json_str = encoded_or_json

        # Step 2: Parse JSON
        request = json.loads(json_str)
        return request

    except json.JSONDecodeError as e:
        return {'error': f'JSON decode error: {e}'}
    except Exception as e:
        return {'error': f'Request processing error: {e}'}


def prepare_response_for_transport(response, do_encode=True):
    """
    Prepare a response dict for transport back to client.

    This applies sanitization and encoding.

    Args:
        response: Response dict to prepare
        do_encode: Whether to apply REST encoding (default True)

    Returns:
        tuple: (prepared_string, original_json)
    """
    import json

    # Step 1: Sanitize the response
    sanitized = sanitize_response(response)

    # Step 2: Serialize to JSON
    json_str = json.dumps(sanitized, indent=None, separators=(',', ':'))

    # Final safety: replace any remaining tab escape sequences
    json_str = json_str.replace('\\t', ' ')

    # Step 3: Encode if requested
    if do_encode:
        encoded = encode_for_rest(json_str)
        return encoded, json_str
    else:
        return json_str, json_str


def process_response_from_transport(encoded_or_json, was_encoded=True):
    """
    Process a received response from transport.

    Args:
        encoded_or_json: The received string
        was_encoded: Whether REST encoding was applied (default True)

    Returns:
        dict: Parsed response, or dict with 'error' key on failure
    """
    import json

    try:
        # Step 1: Decode if needed
        if was_encoded:
            json_str = decode_from_rest(encoded_or_json)
        else:
            json_str = encoded_or_json

        # Step 2: Parse JSON
        response = json.loads(json_str)
        return response

    except json.JSONDecodeError as e:
        return {'error': f'JSON decode error: {e}'}
    except Exception as e:
        return {'error': f'Response processing error: {e}'}


def verify_roundtrip(request, logger=None):
    """
    Verify that a request survives the encode/decode roundtrip.

    This is useful for debugging transport issues.

    Args:
        request: Request dict to verify
        logger: Optional logger for output

    Returns:
        tuple: (success, message, details)
    """
    import json

    try:
        # Prepare for transport
        encoded, original_json = prepare_request_for_transport(request, do_encode=True)

        # Process back
        decoded_request = process_request_from_transport(encoded, was_encoded=True)

        if 'error' in decoded_request:
            return False, f"Decode failed: {decoded_request['error']}", {
                'original_length': len(original_json),
                'encoded_length': len(encoded),
            }

        # Compare key fields
        original = json.loads(original_json)

        checks = []

        # Check files count
        orig_files = len(original.get('files', []))
        dec_files = len(decoded_request.get('files', []))
        if orig_files != dec_files:
            checks.append(f"files count: {orig_files} -> {dec_files}")

        # Check session_state
        orig_exp = original.get('session_state', {}).get('experiment_type')
        dec_exp = decoded_request.get('session_state', {}).get('experiment_type')
        if orig_exp != dec_exp:
            checks.append(f"experiment_type: {orig_exp} -> {dec_exp}")

        # Check history count
        orig_hist = len(original.get('history', []))
        dec_hist = len(decoded_request.get('history', []))
        if orig_hist != dec_hist:
            checks.append(f"history count: {orig_hist} -> {dec_hist}")

        if checks:
            return False, "Roundtrip data mismatch", {
                'differences': checks,
                'original_length': len(original_json),
                'encoded_length': len(encoded),
            }

        return True, "Roundtrip OK", {
            'original_length': len(original_json),
            'encoded_length': len(encoded),
        }

    except Exception as e:
        return False, f"Roundtrip error: {e}", {}
