"""
Test utilities for cctbx-style test execution.

This module provides:
1. Assert helper functions (assert_equal, assert_true, etc.)
2. A fail-fast test runner that crashes on first failure with full traceback

This matches the cctbx testing convention where tests are plain functions
using simple assert statements, and failures crash immediately.

Usage in test files:

    from tests.test_utils import (
        assert_equal, assert_true, assert_false, assert_none,
        assert_not_none, assert_in, assert_raises,
        run_tests_with_fail_fast
    )

    def test_something():
        assert_equal(1, 1)
        assert_true(True)

    def test_another():
        assert_in("a", ["a", "b", "c"])

    if __name__ == "__main__":
        run_tests_with_fail_fast()
"""

from __future__ import absolute_import, division, print_function

import unittest
import inspect


# =============================================================================
# ASSERT HELPERS
# =============================================================================

def assert_equal(actual, expected, message=""):
    """Assert two values are equal."""
    if actual != expected:
        msg = f"expected {expected!r}, got {actual!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_not_equal(actual, expected, message=""):
    """Assert two values are not equal."""
    if actual == expected:
        msg = f"expected not equal to {expected!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_true(value, message=""):
    """Assert value is True."""
    if not value:
        msg = f"expected True, got {value!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_false(value, message=""):
    """Assert value is False."""
    if value:
        msg = f"expected False, got {value!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_none(value, message=""):
    """Assert value is None."""
    if value is not None:
        msg = f"expected None, got {value!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_not_none(value, message=""):
    """Assert value is not None."""
    if value is None:
        msg = "expected not None"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_in(item, container, message=""):
    """Assert item is in container."""
    if item not in container:
        msg = f"{item!r} not found in {container!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_not_in(item, container, message=""):
    """Assert item is not in container."""
    if item in container:
        msg = f"{item!r} unexpectedly found in {container!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_is(actual, expected, message=""):
    """Assert actual is expected (identity check)."""
    if actual is not expected:
        msg = f"{actual!r} is not {expected!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_is_not(actual, expected, message=""):
    """Assert actual is not expected (identity check)."""
    if actual is expected:
        msg = f"{actual!r} is {expected!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_is_instance(obj, cls, message=""):
    """Assert obj is an instance of cls."""
    if not isinstance(obj, cls):
        msg = f"{obj!r} is not an instance of {cls}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_greater(a, b, message=""):
    """Assert a > b."""
    if not (a > b):
        msg = f"{a!r} is not greater than {b!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_greater_equal(a, b, message=""):
    """Assert a >= b."""
    if not (a >= b):
        msg = f"{a!r} is not greater than or equal to {b!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_less(a, b, message=""):
    """Assert a < b."""
    if not (a < b):
        msg = f"{a!r} is not less than {b!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_less_equal(a, b, message=""):
    """Assert a <= b."""
    if not (a <= b):
        msg = f"{a!r} is not less than or equal to {b!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_almost_equal(actual, expected, places=7, message=""):
    """Assert two floats are equal to given decimal places."""
    if round(abs(actual - expected), places) != 0:
        msg = f"{actual!r} != {expected!r} within {places} places"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_raises(exception_type, callable_func, *args, **kwargs):
    """Assert that callable raises the specified exception type."""
    try:
        callable_func(*args, **kwargs)
    except exception_type:
        return  # Expected exception was raised
    except Exception as e:
        raise AssertionError(
            f"Expected {exception_type.__name__}, got {type(e).__name__}: {e}"
        )
    else:
        raise AssertionError(f"Expected {exception_type.__name__} but no exception was raised")


def assert_startswith(string, prefix, message=""):
    """Assert string starts with prefix."""
    if not string.startswith(prefix):
        msg = f"{string!r} does not start with {prefix!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_endswith(string, suffix, message=""):
    """Assert string ends with suffix."""
    if not string.endswith(suffix):
        msg = f"{string!r} does not end with {suffix!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


def assert_contains(string, substring, message=""):
    """Assert string contains substring."""
    if substring not in string:
        msg = f"{substring!r} not found in {string!r}"
        if message:
            msg = f"{message}: {msg}"
        raise AssertionError(msg)


# =============================================================================
# TEST RUNNER
# =============================================================================


def run_tests_with_fail_fast(test_classes=None, verbose=True):
    """
    Run tests with fail-fast behavior.

    Supports both unittest.TestCase classes and plain test_* functions.
    Unlike the standard unittest runner which catches assertion failures
    and continues, this runner lets the first failure raise an uncaught
    exception with full traceback.

    Args:
        test_classes: List of TestCase classes to run. If None, discovers
                     all TestCase subclasses AND test_* functions in the
                     calling module.
        verbose: If True, print test names as they run.

    Raises:
        AssertionError: On first test failure (with full traceback)
        Exception: On first test error (with full traceback)
    """
    # Get the calling module
    frame = inspect.currentframe()
    try:
        caller_module = inspect.getmodule(frame.f_back)
    finally:
        del frame

    total_tests = 0

    # If no classes specified, discover from calling module
    if test_classes is None:
        # Find all TestCase subclasses in that module
        test_classes = []
        for name, obj in inspect.getmembers(caller_module):
            if (inspect.isclass(obj) and
                issubclass(obj, unittest.TestCase) and
                obj is not unittest.TestCase):
                test_classes.append(obj)

    # Sort classes by name for consistent ordering
    test_classes = sorted(test_classes, key=lambda c: c.__name__)

    # Run TestCase classes
    for test_class in test_classes:
        # Get test methods (methods starting with 'test')
        test_methods = [m for m in dir(test_class)
                       if m.startswith('test') and
                       callable(getattr(test_class, m))]
        test_methods.sort()

        # Run setUpClass if present
        if hasattr(test_class, 'setUpClass') and test_methods:
            try:
                test_class.setUpClass()
            except Exception:
                if verbose:
                    print(f"  {test_class.__name__}.setUpClass ... ERROR")
                raise

        for method_name in test_methods:
            total_tests += 1

            # Create instance
            instance = test_class(method_name)

            if verbose:
                # Print test name
                print(f"  {test_class.__name__}.{method_name} ... ", end="", flush=True)

            try:
                instance.setUp()
            except Exception:
                if verbose:
                    print("ERROR (setUp)")
                raise

            # Run the test - let exceptions propagate!
            try:
                getattr(instance, method_name)()
                if verbose:
                    print("ok")
            except AssertionError:
                if verbose:
                    print("FAIL")
                raise
            except Exception:
                if verbose:
                    print("ERROR")
                raise
            finally:
                # Run tearDown
                try:
                    instance.tearDown()
                except Exception:
                    pass  # Don't mask test failures

        # Run tearDownClass after all test methods
        if hasattr(test_class, 'tearDownClass') and test_methods:
            try:
                test_class.tearDownClass()
            except Exception:
                pass

    # Also discover and run plain test_* functions from the caller's module
    test_functions = []
    for name, obj in inspect.getmembers(caller_module):
        if (name.startswith('test_') and
            inspect.isfunction(obj)):
            # Check if function belongs to this module (not imported)
            # For __main__, function's __module__ is also '__main__'
            if hasattr(obj, '__module__'):
                if obj.__module__ == caller_module.__name__:
                    test_functions.append((name, obj))
            else:
                # No __module__ attribute, assume it's local
                test_functions.append((name, obj))

    # Sort by name for consistent ordering
    test_functions.sort(key=lambda x: x[0])

    for func_name, func in test_functions:
        total_tests += 1

        if verbose:
            print(f"  {func_name} ... ", end="", flush=True)

        try:
            func()
            if verbose:
                print("ok")
        except AssertionError:
            if verbose:
                print("FAIL")
            raise
        except Exception:
            if verbose:
                print("ERROR")
            raise

    if verbose:
        print(f"\n{total_tests} tests passed.")


def run_tests_from_module(module, verbose=True):
    """
    Run all tests from a module with fail-fast behavior.

    Args:
        module: The module containing TestCase classes and/or test functions
        verbose: If True, print test names as they run
    """
    test_classes = []
    for name, obj in inspect.getmembers(module):
        if (inspect.isclass(obj) and
            issubclass(obj, unittest.TestCase) and
            obj is not unittest.TestCase):
            test_classes.append(obj)

    # Note: run_tests_with_fail_fast will also discover plain functions
    run_tests_with_fail_fast(test_classes, verbose=verbose)
