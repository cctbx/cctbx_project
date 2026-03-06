"""Minimal test utilities for the agent test suite."""

import traceback


def assert_equal(a, b, msg=""):
    if a != b:
        raise AssertionError("%s: %r != %r" % (msg or "assert_equal", a, b))


def assert_not_equal(a, b, msg=""):
    if a == b:
        raise AssertionError("%s: %r == %r (expected different)" % (msg or "assert_not_equal", a, b))


def assert_true(val, msg=""):
    if not val:
        raise AssertionError(msg or "Expected True, got %r" % val)


def assert_false(val, msg=""):
    if val:
        raise AssertionError(msg or "Expected False, got %r" % val)


def assert_in(item, container, msg=""):
    if item not in container:
        raise AssertionError(msg or "%r not in %r" % (item, container))


def assert_not_in(item, container, msg=""):
    if item in container:
        raise AssertionError(msg or "%r unexpectedly found in %r" % (item, container))


def assert_none(val, msg=""):
    if val is not None:
        raise AssertionError(msg or "Expected None, got %r" % (val,))


def assert_not_none(val, msg=""):
    if val is None:
        raise AssertionError(msg or "Expected non-None value")


def assert_greater(a, b, msg=""):
    if not (a > b):
        raise AssertionError("%s: %r is not greater than %r" % (msg or "assert_greater", a, b))


def assert_greater_equal(a, b, msg=""):
    if not (a >= b):
        raise AssertionError("%s: %r is not >= %r" % (msg or "assert_greater_equal", a, b))


def assert_less(a, b, msg=""):
    if not (a < b):
        raise AssertionError("%s: %r is not less than %r" % (msg or "assert_less", a, b))


def assert_less_equal(a, b, msg=""):
    if not (a <= b):
        raise AssertionError("%s: %r is not <= %r" % (msg or "assert_less_equal", a, b))


def assert_is_instance(obj, cls, msg=""):
    if not isinstance(obj, cls):
        raise AssertionError(msg or "Expected instance of %s, got %s" % (cls.__name__, type(obj).__name__))


def assert_raises(exc_type, fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except exc_type:
        return
    except Exception as e:
        raise AssertionError("Expected %s, got %s: %s" % (exc_type.__name__, type(e).__name__, e))
    raise AssertionError("Expected %s but no exception raised" % exc_type.__name__)


def run_tests_with_fail_fast(test_functions=None, verbose=True):
    """Run test functions, stopping on first failure.

    If *test_functions* is None, auto-discover all ``test_*`` functions in the
    caller's global namespace.
    """
    if test_functions is None:
        import inspect
        frame = inspect.stack()[1]
        caller_globals = frame[0].f_globals
        test_functions = [
            v for k, v in sorted(caller_globals.items())
            if k.startswith("test_") and callable(v)
        ]

    passed = 0
    failed = 0
    errors = []
    for fn in test_functions:
        name = fn.__name__
        try:
            fn()
            passed += 1
            if verbose:
                print("  PASS: %s" % name)
        except Exception as e:
            failed += 1
            errors.append((name, e))
            if verbose:
                print("  FAIL: %s — %s" % (name, e))
                traceback.print_exc()
            break  # fail fast

    print("\n%d passed, %d failed" % (passed, failed))
    if errors:
        first_name, first_exc = errors[0]
        summary = "First failure: %s — %s" % (
            first_name, first_exc)
        print(summary)
        # Raise so callers (including run_all_tests.py) detect
        # the failure.  The full traceback was already printed
        # above at the point of failure.
        raise AssertionError(summary)
    return True
