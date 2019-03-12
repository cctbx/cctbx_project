#!/usr/bin/env python
# coding: utf-8

"""
Test new features of the find_unused_imports_crude tool.

At the moment this only tests new features to verify that it's
not breaking *the existing behaviour*. Any unwanted existing
behaviour is not reinforced with this test.
"""

from __future__ import absolute_import, division, print_function
from libtbx.command_line.find_unused_imports_crude import inspect

test_cases_to_catch = {
    "from x import y": {"y"},
    "import x": {"x"},
    "import a,b,c\na, b": {"c"},
    "import x # noqa: W334": {"x"},
    "import x, y, z # noqa: W334": {"x", "y", "z"},
}

test_cases_to_ignore = [
    "from x import y # noqa",
    "import x\nx",
    "import x # noqa",
    "import x # noqa: F401",
    "import x, y, z # noqa: F401",
    "import x, y, z # noqa - all of these are essential",
]


def test_find_unused_imports_crude():
    # Check the cases where we expect to catch
    for case, results in test_cases_to_catch.items():
        assert set(inspect(case.splitlines())) == results, case

    # Test cases where we expect to pass
    for case in test_cases_to_ignore:
        assert not inspect(case.splitlines())


if __name__ == "__main__":
    test_find_unused_imports_crude()
    print("OK")
