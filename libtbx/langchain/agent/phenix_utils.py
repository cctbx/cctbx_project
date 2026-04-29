"""
Utility functions previously imported from phenix modules.

This module provides standalone implementations of functions that were
previously imported from phenix.rest and phenix.command_line.doc,
allowing the langchain agent to work without phenix installed.

Functions:
- text_as_simple_string: Encode text for REST transport
- simple_string_as_text: Decode text from REST transport
- load_url: Open a URL in the system browser
"""

from __future__ import absolute_import, division, print_function

import os
import sys


# =============================================================================
# REST ENCODING/DECODING
# =============================================================================

# Special marker strings for REST transport encoding
SPECIAL_STRING_LIST = [
    "ZZCRZZ",  # newline
    "ZZLBZZ",  # left brace {
    "ZZRBZZ",  # right brace }
    "ZZSCZZ",  # semicolon ;
    "ZZHAZZ",  # hash #
    "ZZTAZZ",  # tab
    "ZZDQZZ",  # double quote "
    "ZZSQZZ",  # single quote '
    "ZZSRZZ",  # backtick `
    "ZZEXZZ",  # exclamation !
    "ZZDSZZ",  # dollar sign $
]

TEXT_TO_REPLACE_LIST = [
    "\n", "{", "}", ";", "#",
    '\t', '"', "'", '`', '!', '$'
]


def allowed_chars():
    """Return all allowed characters for simple_string."""
    text = 'abcdefghijklmnopqrstuvwxyz0123456789'
    text += '~@%^&*()_-–[]+=|?<>,./:; \\'
    text += "".join(['"', "'", '`', '!', '#', '$', '{', '}', '\t', '\n'])
    return text


def unique_chars(text):
    """Return unique chars in text. Apply lower case first."""
    text = text.lower()

    if len(text) <= 1:
        return text

    first_char = text[0]
    t0 = first_char
    for i in range(len(text)):  # max tries
        text = text[1:].replace(t0, "") + t0
        if len(text) <= 1:
            return text
        t0 = text[0]
        if t0 == first_char:  # went all the way through
            break
    chars = [c for c in text]
    chars.sort()
    text = "".join(chars)
    return text


def non_allowed_chars(text):
    """Return any non-allowed characters in text.

    Allowed are: A-Z, a-z, 0-9, ~@%^&*()_-[]+=|?<>,./: and space
    Special are: ' " ` ! # $ { } \\t \\n
    All others are non-allowed (i.e., backslash, backslash-a)
    """
    text = unique_chars(text)
    for c in allowed_chars():
        text = text.replace(c, "")
    return text


def remove_non_allowed_chars(text):
    """Remove non-allowed characters from text."""
    if not text:
        return text
    non_allowed = non_allowed_chars(text)
    for c in non_allowed:
        text = text.replace(c, "")
    return text


def as_text(s):
    """Returns s as a string whether it is utf or not."""
    if hasattr(s, 'decode'):
        return s.decode(encoding='utf-8')
    else:
        return s


def text_as_simple_string(text):
    """Encode text for REST transport by replacing special characters.

    Replaces newline characters, braces, and other special chars with
    marker strings (ZZCRZZ, ZZLBZZ, etc.) and removes non-allowed chars.

    Args:
        text (str): The input string.

    Returns:
        str: The encoded string safe for REST transport.
    """
    if not text:
        return text
    working_text = remove_non_allowed_chars(text)
    for text_to_replace, special_string in zip(
            TEXT_TO_REPLACE_LIST, SPECIAL_STRING_LIST):
        working_text = working_text.replace(text_to_replace, special_string)
    return as_text(working_text)


def simple_string_as_text(simple_string):
    """Decode text from REST transport by restoring special characters.

    Replaces marker strings (ZZCRZZ, ZZLBZZ, etc.) back to their
    original characters (newlines, braces, etc.).

    Args:
        simple_string (str): The encoded string from REST transport.

    Returns:
        str: The decoded original string.
    """
    if not simple_string:
        return simple_string
    working_text = simple_string
    for text_to_replace, special_string in zip(
            TEXT_TO_REPLACE_LIST, SPECIAL_STRING_LIST):
        working_text = working_text.replace(special_string, text_to_replace)
    return as_text(working_text)


# =============================================================================
# BROWSER/URL UTILITIES
# =============================================================================

# Module-level cache for browser path
_browser_path = None
_browser_name = None


def load_url(url, reset_ld_library_path=True):
    """Open a URL in the system's default browser.

    Args:
        url (str): The URL to open (typically a file:// URL).
        reset_ld_library_path (bool): On Linux, whether to clear
            LD_LIBRARY_PATH before launching browser (avoids library
            conflicts with PHENIX). Defaults to True.

    Raises:
        RuntimeError: If no suitable browser can be found on Linux.
    """
    global _browser_path, _browser_name

    # Try to import easy_run from libtbx
    try:
        from libtbx import easy_run
    except ImportError:
        # Fallback: use subprocess directly
        import subprocess

        class _EasyRunFallback:
            @staticmethod
            def call(cmd):
                if isinstance(cmd, list):
                    subprocess.Popen(cmd)
                else:
                    subprocess.Popen(cmd, shell=True)

        easy_run = _EasyRunFallback()

    if sys.platform == 'win32':
        easy_run.call(["explorer", "file://%s" % url])

    elif sys.platform == "darwin":
        easy_run.call("open 'file://%s'" % url)

    else:  # linux, etc.
        if _browser_path is None:
            browser_names = [
                "firefox",
                "mozilla",
                "netscape",
            ]

            for browser_name in browser_names:
                for path in os.environ.get("PATH", "").split(os.pathsep):
                    browser_candidate = os.path.join(path, browser_name)
                    if os.path.isfile(browser_candidate):
                        _browser_path = browser_candidate
                        _browser_name = browser_name
                        break
                if _browser_path is not None:
                    break

            if _browser_path is None:
                print("Could not find a browser in PATH")
                print("\nBrowser choices:")
                for browser_name in browser_names:
                    print("\t%s" % browser_name)
                print("\nPath choices:")
                for path in os.environ.get("PATH", "").split(os.pathsep):
                    print("\t%s" % path)
                raise RuntimeError("Couldn't find a valid browser.")

        print("  Using browser %s" % _browser_path)
        if reset_ld_library_path:
            easy_run.call(
                "env LD_LIBRARY_PATH='' LD_PRELOAD='' %s '%s' &" % (_browser_path, url))
        else:
            easy_run.call("%s '%s' &" % (_browser_path, url))
